#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 || $# -gt 7 ]]; then
  echo "Usage: $0 samples.tsv genes.bed variant_gene_pairs.tsv output_dir phaser_dir [threads=4] [blacklist.bed]" >&2
  exit 1
fi

samples=$1
features=$2
pairs=$3
out=$4
phaser_dir=$5
threads=${6:-4}
blacklist=${7:-}

phaser="$phaser_dir/phaser/phaser.py"
gene_ae="$phaser_dir/phaser_gene_ae/phaser_gene_ae.py"
expr_matrix="$phaser_dir/phaser_pop/phaser_expr_matrix.py"
cis_var="$phaser_dir/phaser_pop/phaser_cis_var.py"

for cmd in python3 bcftools bgzip tabix; do
  command -v "$cmd" >/dev/null || { echo "Missing command: $cmd" >&2; exit 1; }
done
for file in "$samples" "$features" "$pairs" "$phaser" "$gene_ae" "$expr_matrix" "$cis_var"; do
  [[ -f "$file" ]] || { echo "Missing file: $file" >&2; exit 1; }
done
[[ -z "$blacklist" || -f "$blacklist" ]] || { echo "Missing blacklist: $blacklist" >&2; exit 1; }

mkdir -p "$out/sample" "$out/tissue"
tail -n +2 "$samples" | while IFS=$'\t' read -r sample tissue bam vcf vcf_sample; do
  [[ -n "$sample" && -n "$tissue" && -n "$bam" && -n "$vcf" && -n "$vcf_sample" ]] || continue
  [[ -f "$bam" && ( -f "$bam.bai" || -f "${bam%.bam}.bai" ) ]] || { echo "BAM/index missing: $bam" >&2; exit 1; }
  [[ -f "$vcf" && ( -f "$vcf.tbi" || -f "$vcf.csi" ) ]] || { echo "VCF/index missing: $vcf" >&2; exit 1; }

  sample_dir="$out/sample/$sample"
  prefix="$sample_dir/$sample"
  mkdir -p "$sample_dir"
  blacklist_arg=()
  [[ -z "$blacklist" ]] || blacklist_arg=(--blacklist "$blacklist")

  python3 "$phaser" \
    --vcf "$vcf" --bam "$bam" --sample "$vcf_sample" --o "$prefix" \
    --paired_end 1 --mapq 255 --baseq 10 --gw_phase_vcf 1 \
    --pass_only 0 --unphased_vars 1 --threads "$threads" \
    "${blacklist_arg[@]}"

  python3 "$gene_ae" \
    --haplotypic_counts "$prefix.haplotypic_counts.txt" \
    --features "$features" --o "$prefix.gene_ae.txt"
done

cut -f2 "$samples" | tail -n +2 | sort -u | while read -r tissue; do
  [[ -n "$tissue" ]] || continue
  tissue_dir="$out/tissue/$tissue"
  mkdir -p "$tissue_dir/gene_ae"
  map="$tissue_dir/sample_map.tsv"
  vcf_list="$tissue_dir/phaser_vcfs.list"
  printf 'vcf_sample\tbed_sample\n' > "$map"
  : > "$vcf_list"

  awk -F '\t' -v t="$tissue" 'NR>1 && $2==t {print $1"\t"$3"\t"$5}' "$samples" |
  while IFS=$'\t' read -r sample bam vcf_sample; do
    ln -sf "$out/sample/$sample/$sample.gene_ae.txt" "$tissue_dir/gene_ae/$sample.txt"
    printf '%s\t%s\n' "$vcf_sample" "$bam" >> "$map"
    printf '%s\n' "$out/sample/$sample/$sample.vcf.gz" >> "$vcf_list"
  done

  merged_vcf="$tissue_dir/$tissue.phaser_merged.vcf.gz"
  if [[ $(wc -l < "$vcf_list") -eq 1 ]]; then
    bcftools view -Oz -o "$merged_vcf" "$(<"$vcf_list")"
  else
    bcftools merge -l "$vcf_list" -Oz -o "$merged_vcf"
  fi
  tabix -f -p vcf "$merged_vcf"

  matrix_prefix="$tissue_dir/$tissue.haplotype_expression"
  (cd "$tissue_dir" && python3 "$expr_matrix" \
    --gene_ae_dir gene_ae --features "$features" --o "$matrix_prefix" --t "$threads")

  python3 "$cis_var" \
    --bed "$matrix_prefix.gw_phased.bed.gz" --vcf "$merged_vcf" \
    --pairs "$pairs" --map "$map" --o "$tissue_dir/$tissue.ASE_aFC.tsv" \
    --min_cov 8 --bs 100 --t "$threads"
done

echo "ASE analysis completed: $out"
