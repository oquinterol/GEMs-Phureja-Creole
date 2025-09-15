#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Uso: $0 KO_LIST OUT_DIR" >&2; exit 1
fi
KO_LIST="$1"
OUT_DIR="$2"
mkdir -p "$OUT_DIR"

VALID_TMP="$(mktemp)"
INVALID_TMP="${OUT_DIR}/ko_links.invalid"
: > "$INVALID_TMP"

# Normaliza: trim, quita 'ko:' si viene, valida K + 5 dígitos
while IFS= read -r ko; do
  ko_trim="$(echo "$ko" | tr -d '\r' | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//; s/^ko://I')"
  [[ -z "$ko_trim" ]] && continue
  if [[ "$ko_trim" =~ ^K[0-9]{5}$ ]]; then
    echo "$ko_trim" >> "$VALID_TMP"
  else
    echo "$ko" >> "$INVALID_TMP"
  fi
done < "$KO_LIST"

if [[ ! -s "$VALID_TMP" ]]; then
  echo "ERROR: No hay KOs válidos en ${KO_LIST}. Revisa ${INVALID_TMP}" >&2
  exit 1
fi

KO2RXN="${OUT_DIR}/ko2rxn.tsv"
KO2MOD="${OUT_DIR}/ko2mod.tsv"
KO2MAP="${OUT_DIR}/ko2map.tsv"
: > "$KO2RXN"; : > "$KO2MOD"; : > "$KO2MAP"
ERR_LOG="${OUT_DIR}/ko_links.errors"; : > "$ERR_LOG"

# Batches de 10: ko:Kxxxxx+ko:Kyyyy...
sed 's/^/ko:/' "$VALID_TMP" | xargs -n10 | sed 's/ /+/g' | while IFS= read -r chunk; do
  [[ -z "$chunk" ]] && continue
  curl -fsS "https://rest.kegg.jp/link/reaction/${chunk}" >> "$KO2RXN" || echo "[link/reaction] 400: $chunk" >> "$ERR_LOG"
  curl -fsS "https://rest.kegg.jp/link/module/${chunk}"  >> "$KO2MOD" || echo "[link/module] 400: $chunk"  >> "$ERR_LOG"
  curl -fsS "https://rest.kegg.jp/link/pathway/${chunk}" >> "$KO2MAP" || echo "[link/pathway] 400: $chunk" >> "$ERR_LOG"
  sleep 0.5
done

for f in "$KO2RXN" "$KO2MOD" "$KO2MAP"; do sort -u -o "$f" "$f"; done

if [[ -s "$INVALID_TMP" ]]; then echo "Aviso: se ignoraron KOs inválidos → $INVALID_TMP" >&2; fi
if [[ -s "$ERR_LOG" ]]; then echo "Algunos batches fallaron → $ERR_LOG" >&2; fi

rm -f "$VALID_TMP"
