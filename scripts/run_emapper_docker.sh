#!/usr/bin/env bash
set -euo pipefail

# eggNOG-mapper via Docker helper
# - Downloads databases to a local dir
# - Runs emapper.py on a protein FASTA
# - Registra el comando Docker usado en un archivo para reproducibilidad

SCRIPT_NAME=$(basename "$0")

usage() {
  cat <<USAGE
Usage: $SCRIPT_NAME <download|run|both> [options]

Actions:
  download                Download eggNOG databases to --data_dir
  run                     Run emapper.py on --input (proteins FASTA)
  both                    Download (if needed) and run

Options:
  --image IMG             Docker image (default: quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_1)
  --data_dir DIR          Database dir on host (default: ./eggnog_data)
  --input FILE            Input proteins FASTA (default: data/genome/proteins.faa)
  --out_dir DIR           Output dir on host (default: reports)
  --out_prefix NAME       Output file prefix (default: emapper)
  --threads N             CPU threads (default: 8)
  --mode {diamond|mmseqs} Search backend (default: diamond)
  --tax_scope TAXID       Limit to taxonomy (e.g., 33090 Viridiplantae). Omit to search all.
  --resume                Resume a previous run if outputs exist
  --override              Overwrite existing outputs
  --extra ARGS            Extra args passed verbatim to emapper.py (quoted string)
  -h, --help              Show this help

Examples:
  $SCRIPT_NAME download --data_dir ~/eggnog_data
  $SCRIPT_NAME run --input data/genome/proteins.faa --data_dir ~/eggnog_data --threads 16
  $SCRIPT_NAME both --data_dir ~/eggnog_data --threads 12 --tax_scope 33090
Nota: el script guarda el comando ejecutado en <out_dir>/<out_prefix>.docker_cmd.txt
USAGE
}

ACTION=${1-}
if [[ -z "$ACTION" || "$ACTION" == "-h" || "$ACTION" == "--help" ]]; then
  usage; exit 0
fi
shift

IMG="${EGGNOG_DOCKER_IMAGE:-quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_1}"
DOCKER_NET="${DOCKER_NETWORK_MODE:-host}"
DATA_DIR="${EGGNOG_DATA_DIR:-./eggnog_data}"
INPUT="data/genome/proteins.faa"
OUT_DIR="reports"
OUT_PREFIX="emapper"
THREADS=8
MODE="diamond"
TAX_SCOPE=""
EXTRA_ARGS=""
RESUME=0
OVERRIDE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --image) IMG=$2; shift 2;;
    --data_dir) DATA_DIR=$2; shift 2;;
    --input) INPUT=$2; shift 2;;
    --out_dir) OUT_DIR=$2; shift 2;;
    --out_prefix) OUT_PREFIX=$2; shift 2;;
    --threads) THREADS=$2; shift 2;;
    --mode) MODE=$2; shift 2;;
    --tax_scope) TAX_SCOPE=$2; shift 2;;
    --resume) RESUME=1; shift 1;;
    --override) OVERRIDE=1; shift 1;;
    --extra) EXTRA_ARGS=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 2;;
  esac
done

# Resolve absolute paths
WORK_DIR=$(pwd)
DATA_DIR_ABS=$(mkdir -p "$DATA_DIR" && cd "$DATA_DIR" && pwd)
OUT_DIR_ABS=$(mkdir -p "$OUT_DIR" && cd "$OUT_DIR" && pwd)
INPUT_ABS=$(cd "$(dirname "$INPUT")" && pwd)/"$(basename "$INPUT")"

if [[ "$ACTION" != "download" ]]; then
  if [[ ! -f "$INPUT_ABS" ]]; then
    echo "ERROR: input not found: $INPUT_ABS" >&2
    exit 2
  fi
fi

docker_pull_if_needed() {
  if ! docker image inspect "$IMG" >/dev/null 2>&1; then
    echo "Pulling image: $IMG" >&2
    docker pull "$IMG"
  fi
}

run_download() {
  echo "[download] Using data dir: $DATA_DIR_ABS" >&2
  DL_OPTS="-y"
  # If user plans to use mmseqs, prefetch that DB as well
  if [[ "$MODE" == "mmseqs" ]]; then
    DL_OPTS="$DL_OPTS -M"
  fi
  docker run --rm --network "$DOCKER_NET" \
    -v "$DATA_DIR_ABS":/db \
    "$IMG" \
    download_eggnog_data.py $DL_OPTS --data_dir /db
}

run_emapper() {
  echo "[run] Input: $INPUT_ABS" >&2
  echo "[run] Output dir: $OUT_DIR_ABS (prefix: $OUT_PREFIX)" >&2
  echo "[run] Data dir: $DATA_DIR_ABS" >&2
  BACKEND_OPT="-m $MODE"
  TAX_OPT=""
  if [[ -n "$TAX_SCOPE" ]]; then
    TAX_OPT="--tax_scope $TAX_SCOPE"
  fi
  # Build container-visible paths and mounts
  MOUNTS="-v \"$WORK_DIR\":/work -w /work -v \"$DATA_DIR_ABS\":/db"

  # Input path inside container
  INPUT_CONT=""
  if [[ "$INPUT_ABS" == "$WORK_DIR"/* ]]; then
    INPUT_REL=${INPUT_ABS#"$WORK_DIR"/}
    INPUT_CONT="/work/$INPUT_REL"
  else
    INPUT_DIR_HOST=$(dirname "$INPUT_ABS")
    INPUT_BASE=$(basename "$INPUT_ABS")
    MOUNTS="$MOUNTS -v \"$INPUT_DIR_HOST\":/in"
    INPUT_CONT="/in/$INPUT_BASE"
  fi

  # Output dir inside container
  OUT_DIR_CONT=""
  if [[ "$OUT_DIR_ABS" == "$WORK_DIR"/* ]]; then
    OUT_REL=${OUT_DIR_ABS#"$WORK_DIR"/}
    OUT_DIR_CONT="/work/$OUT_REL"
  else
    MOUNTS="$MOUNTS -v \"$OUT_DIR_ABS\":/out"
    OUT_DIR_CONT="/out"
  fi

  echo "[run] Container input path: $INPUT_CONT" >&2
  echo "[run] Container output dir: $OUT_DIR_CONT" >&2

  # Handle resume/override flags (override takes precedence if both are set)
  RESUME_OPT=""
  OVERRIDE_OPT=""
  if [[ "$RESUME" -eq 1 && "$OVERRIDE" -eq 1 ]]; then
    echo "[run] Both --resume and --override set; using --override" >&2
    OVERRIDE_OPT="--override"
  elif [[ "$RESUME" -eq 1 ]]; then
    RESUME_OPT="--resume"
  elif [[ "$OVERRIDE" -eq 1 ]]; then
    OVERRIDE_OPT="--override"
  fi

  CMD="docker run --rm --network \"$DOCKER_NET\" \
    $MOUNTS \
    \"$IMG\" \
    emapper.py -i \"$INPUT_CONT\" --itype proteins $BACKEND_OPT \
      --data_dir /db --cpu \"$THREADS\" -o \"$OUT_PREFIX\" --output_dir \"$OUT_DIR_CONT\" \
      $TAX_OPT $RESUME_OPT $OVERRIDE_OPT ${EXTRA_ARGS}"

  # Log exact command for reproducibility
  echo "$CMD" | sed 's/\\\\/\\n\\\\/g' > "$OUT_DIR_ABS/$OUT_PREFIX.docker_cmd.txt"
  echo "[run] Logged docker command to $OUT_DIR_ABS/$OUT_PREFIX.docker_cmd.txt" >&2

  eval "$CMD"
}

docker_pull_if_needed

case "$ACTION" in
  download)
    run_download ;;
  run)
    run_emapper ;;
  both)
    # Download only if core DB files are missing
    DMND_UNDERSCORE="$DATA_DIR_ABS/eggnog_proteins.dmnd"
    DMND_DOT="$DATA_DIR_ABS/eggnog.proteins.dmnd"
    if [[ ! -s "$DATA_DIR_ABS/eggnog.db" || ( ! -s "$DMND_UNDERSCORE" && ! -s "$DMND_DOT" ) ]]; then
      run_download
    else
      echo "[both] DB appears present in $DATA_DIR_ABS; skipping download." >&2
    fi
    run_emapper ;;
  *)
    echo "Unknown action: $ACTION" >&2; usage; exit 2 ;;
esac
