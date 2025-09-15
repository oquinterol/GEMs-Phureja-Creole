#!/usr/bin/env bash
set -euo pipefail

# Uso: rxn2cpd.sh RXN_LIST OUT_CMP_LIST
# RXN_LIST: archivo con IDs tipo R00036 (uno por línea)
# OUT_CMP_LIST: salida con CIDs tipo C00001 (únicos, ordenados)

if [[ $# -lt 2 ]]; then
  echo "Uso: $0 RXN_LIST OUT_CMP_LIST" >&2
  exit 1
fi

RXN_LIST="$1"
OUT_LIST="$2"
mkdir -p "$(dirname "$OUT_LIST")"

# Prepara lotes rn:Rxxxx + rn:Ryyyy (10 por request para no pasarse de URL)
mk_batches() {
  awk 'NF{gsub(/^rn:/,""); print "rn:"$0}' "$RXN_LIST" \
  | paste -d'+' - - - - - - - - - - 2>/dev/null || true
}

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

: > "$tmp"
mk_batches | while IFS= read -r chunk; do
  [[ -z "$chunk" ]] && continue
  curl -fsS "https://rest.kegg.jp/link/cpd/${chunk}" || {
    echo "ERROR: fallo link cpd para: ${chunk}" >&2
    exit 2
  }
  sleep 0.5
done >> "$tmp"

# Extrae segunda columna (CIDs), normaliza y deduplica
awk -F'\t' 'NF>=2{print $2}' "$tmp" \
| sed -E 's/^cpd://; s/^(C[0-9]{5}).*/\1/' \
| awk 'length($0)==6' \
| sort -u > "$OUT_LIST"

# Sanidad
if [[ ! -s "$OUT_LIST" ]]; then
  echo "ERROR: OUT_CMP_LIST quedó vacío. ¿RXN_LIST tiene RIDs válidos?" >&2
  exit 3
fi
