#!/usr/bin/env bash
# Scarica un file processato del 4DN Data Portal (es. un .mcool) in data/,
# leggendo URL pubblico e MD5 dall'API 4DN. Il download riprende se interrotto.
# Uso: scripts/fetch_4dn.sh 4DNFIXVAKX9Q
set -euo pipefail

acc="${1:?uso: $0 <accession 4DN, es. 4DNFIXVAKX9Q>}"
[[ "$acc" =~ ^4DNFI[A-Z0-9]{7}$ ]] || { echo "accession non valida: $acc" >&2; exit 1; }
dest="$(cd "$(dirname "$0")/.." && pwd)/data"

read -r url md5 fname size < <(
  curl -fsSL -H 'Accept: application/json' "https://data.4dnucleome.org/files-processed/${acc}/" |
  python3 -c '
import json, sys
f = json.load(sys.stdin)
if f.get("status") != "released" or not f.get("open_data_url"):
    acc, status = f.get("accession"), f.get("status")
    sys.exit(f"{acc}: non pubblico (status={status})")
print(f["open_data_url"], f["md5sum"], f["display_title"], f["file_size"])
')

out="${dest}/${fname}"
if [ -f "$out" ] && echo "${md5}  ${out}" | md5sum -c --status; then
  echo "${fname}: già presente e verificato"; exit 0
fi

free=$(df --output=avail -B1 "$dest" | tail -1)
(( free > size )) || { echo "spazio insufficiente: servono $((size / 10**9)) GB" >&2; exit 1; }

echo "scarico ${fname} ($((size / 10**9)) GB) da ${url}"
# S3 a volte si blocca a metà: sotto 10 kB/s per 60 s si interrompe e si
# riprende dal punto raggiunto (-C -), finché il file non è completo.
for attempt in $(seq 50); do
  have=$(stat -c %s "${out}.part" 2>/dev/null || echo 0)
  (( have >= size )) && break
  (( attempt > 1 )) && echo "ripresa n. $((attempt - 1)) da $((have / 10**6)) MB"
  curl -fL -C - --speed-limit 10000 --speed-time 60 --retry 5 --retry-all-errors \
    "$url" -o "${out}.part" || sleep 10
done
(( $(stat -c %s "${out}.part") == size )) || { echo "download incompleto" >&2; exit 1; }
echo "verifico MD5..."
echo "${md5}  ${out}.part" | md5sum -c --quiet
mv "${out}.part" "$out"
echo "pronto: ${out}"
