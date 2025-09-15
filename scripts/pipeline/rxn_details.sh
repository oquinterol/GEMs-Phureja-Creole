#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Uso: $0 RXN_LIST OUT_FLAT" >&2; exit 1
fi
RXN_LIST="$1"
OUT_FLAT="$2"
mkdir -p "$(dirname "$OUT_FLAT")"

mk_batches() {
  awk '{print "rn:"$0}' "$RXN_LIST" | paste -d'+' - - - - - - - - - - 2>/dev/null || true
}

: > "$OUT_FLAT"
mk_batches | while IFS= read -r chunk; do
  [[ -z "$chunk" ]] && continue
  curl -fsS "https://rest.kegg.jp/get/${chunk}" >> "$OUT_FLAT" || true
  printf "\n" >> "$OUT_FLAT"
  sleep 0.5
done
