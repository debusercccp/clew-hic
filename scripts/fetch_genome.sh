#!/usr/bin/env bash
# Scarica un genoma UCSC (fasta + chrom.sizes) in data/genomes/<asm>/,
# verifica l'MD5, decomprime e crea l'indice .fai (serve a bioframe/pyfaidx).
# Uso: scripts/fetch_genome.sh hg38   (oppure mm10)
set -euo pipefail

asm="${1:?uso: $0 <hg38|mm10>}"
base="https://hgdownload.soe.ucsc.edu/goldenPath/${asm}/bigZips"
dest="$(cd "$(dirname "$0")/.." && pwd)/data/genomes/${asm}"
mkdir -p "$dest"
cd "$dest"

curl -fsSL "${base}/md5sum.txt" -o md5sum.txt
for f in "${asm}.chrom.sizes" "${asm}.fa.gz"; do
  if [ -f "$f" ] && grep " ${f}$" md5sum.txt | md5sum -c --status; then
    echo "${f}: già presente e verificato"
    continue
  fi
  echo "scarico ${f}..."
  curl -fL -C - --retry 5 "${base}/${f}" -o "$f"
  grep " ${f}$" md5sum.txt | md5sum -c
done

if [ ! -f "${asm}.fa" ]; then
  echo "decomprimo ${asm}.fa.gz..."
  gzip -dc "${asm}.fa.gz" > "${asm}.fa.tmp" && mv "${asm}.fa.tmp" "${asm}.fa"
fi

if [ ! -f "${asm}.fa.fai" ]; then
  echo "indicizzo ${asm}.fa..."
  docker run --rm -u "$(id -u):$(id -g)" -v "$dest:/g" clew-env:latest \
    python -c "import pyfaidx; pyfaidx.Faidx('/g/${asm}.fa')"
fi
echo "pronto: ${dest}/${asm}.fa"
