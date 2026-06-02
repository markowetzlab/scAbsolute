#!/bin/bash
# One-time script to lift mm9 FSU RepliChip BigWig files over to mm10 coordinates.
#
# Prerequisites (install once):
#   conda install -c bioconda ucsc-liftover ucsc-bigwigtobedgraph ucsc-bedgraphtobigwig
# OR download UCSC binaries directly:
#   https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/
#
# Run from: scAbsolute/data/replicationTiming/
# Output:   repliChip_data_mm10/ (used by replicationInfo.R)
# After running this script, re-run:
#   Rscript replicationInfo.R    -> regenerates replicationInfo_mouse.RDS
#   Rscript replicationInfo_per_bin.R -> regenerates replicationTiming_per_binSize_mouse.RDS

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MM9_DIR="${SCRIPT_DIR}/repliChip_data_mm9"
MM10_DIR="${SCRIPT_DIR}/repliChip_data_mm10"
BASE_URL="http://hgdownload.cse.ucsc.edu/goldenpath/mm9/encodeDCC/wgEncodeFsuRepliChip"
CHAIN="${SCRIPT_DIR}/mm9ToMm10.over.chain"
CHROM_SIZES="${SCRIPT_DIR}/mm10.chrom.sizes"

mkdir -p "$MM9_DIR" "$MM10_DIR"

echo "=== Downloading liftover chain and chrom sizes ==="
if [ ! -f "$CHAIN" ]; then
  curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz" \
    | gunzip > "$CHAIN"
fi
if [ ! -f "$CHROM_SIZES" ]; then
  curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes" \
    -o "$CHROM_SIZES"
fi

files=(
  "wgEncodeFsuRepliChipCh12FWaveSignalRep1"
  "wgEncodeFsuRepliChipCh12FWaveSignalRep2"
  "wgEncodeFsuRepliChipEpisc5MWaveSignalRep1"
  "wgEncodeFsuRepliChipEpisc5MWaveSignalRep2"
  "wgEncodeFsuRepliChipEpisc7FWaveSignalRep1"
  "wgEncodeFsuRepliChipEpisc7FWaveSignalRep2"
  "wgEncodeFsuRepliChipEs46cMDifff6dWaveSignalRep1"
  "wgEncodeFsuRepliChipEs46cMWaveSignalRep1"
  "wgEncodeFsuRepliChipEsd3MDiffe3dWaveSignalRep1"
  "wgEncodeFsuRepliChipEsd3MDiffe3dWaveSignalRep2"
  "wgEncodeFsuRepliChipEsd3MDiffe6dWaveSignalRep1"
  "wgEncodeFsuRepliChipEsd3MDiffe6dWaveSignalRep2"
  "wgEncodeFsuRepliChipEsd3MDiffe9dWaveSignalRep1"
  "wgEncodeFsuRepliChipEsd3MDiffe9dWaveSignalRep2"
  "wgEncodeFsuRepliChipEsd3MDiffg3dWaveSignalRep1"
  "wgEncodeFsuRepliChipEsd3MDiffg3dWaveSignalRep2"
  "wgEncodeFsuRepliChipEsd3MWaveSignalRep1"
  "wgEncodeFsuRepliChipEsd3MWaveSignalRep2"
  "wgEncodeFsuRepliChipEsem5sUDiffhsoxmWaveSignalRep1"
  "wgEncodeFsuRepliChipEsem5sUDiffhsoxmWaveSignalRep2"
  "wgEncodeFsuRepliChipEsem5sUDiffhsoxpWaveSignalRep1"
  "wgEncodeFsuRepliChipEsem5sUDiffhsoxpWaveSignalRep2"
  "wgEncodeFsuRepliChipEstt2MDifff9dWaveSignalRep1"
  "wgEncodeFsuRepliChipEstt2MDifff9dWaveSignalRep2"
  "wgEncodeFsuRepliChipEstt2MWaveSignalRep1"
  "wgEncodeFsuRepliChipEstt2MWaveSignalRep2"
  "wgEncodeFsuRepliChipJ185aUWaveSignalRep1"
  "wgEncodeFsuRepliChipJ185aUWaveSignalRep2"
  "wgEncodeFsuRepliChipL1210FWaveSignalRep1"
  "wgEncodeFsuRepliChipL1210FWaveSignalRep2"
  "wgEncodeFsuRepliChipMefMWaveSignalRep1"
  "wgEncodeFsuRepliChipMelMWaveSignalRep1"
  "wgEncodeFsuRepliChipMelMWaveSignalRep2"
)

for f in "${files[@]}"; do
  out="${MM10_DIR}/${f}.bigWig"
  if [ -f "$out" ]; then
    echo "Skipping $f (already lifted)"
    continue
  fi

  echo "=== Processing $f ==="

  mm9_bw="${MM9_DIR}/${f}.bigWig"
  if [ ! -f "$mm9_bw" ]; then
    echo "  Downloading mm9 BigWig..."
    curl -s "${BASE_URL}/${f}.bigWig" -o "$mm9_bw"
  fi

  echo "  Converting BigWig -> bedGraph..."
  bigWigToBedGraph "$mm9_bw" "${MM9_DIR}/${f}.bedGraph"

  echo "  Lifting over mm9 -> mm10..."
  liftOver \
    "${MM9_DIR}/${f}.bedGraph" \
    "$CHAIN" \
    "${MM10_DIR}/${f}_unsorted.bedGraph" \
    "${MM10_DIR}/${f}_unmapped.txt"

  echo "  Sorting and converting bedGraph -> BigWig..."
  sort -k1,1 -k2,2n "${MM10_DIR}/${f}_unsorted.bedGraph" \
    > "${MM10_DIR}/${f}_sorted.bedGraph"
  bedGraphToBigWig \
    "${MM10_DIR}/${f}_sorted.bedGraph" \
    "$CHROM_SIZES" \
    "$out"

  rm -f "${MM10_DIR}/${f}_unsorted.bedGraph" \
        "${MM10_DIR}/${f}_sorted.bedGraph"

  echo "  Done: $out"
done

echo ""
echo "=== Liftover complete ==="
echo "Next steps:"
echo "  Rscript replicationInfo.R"
echo "  Rscript replicationInfo_per_bin.R"
