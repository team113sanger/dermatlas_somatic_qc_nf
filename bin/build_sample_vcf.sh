#!/usr/bin/env bash

# Build a single SigProfiler-ready VCF for one sample by pulling PASS variants
# from the sample's caveman VCF (SBS/DBS targets) and pindel VCF (ID targets),
# then concatenating. Logs a WARNING if variant counts between the MAF slice
# and the resulting VCF disagree, but never fails.

set -euo pipefail

usage() {
    cat >&2 <<EOF
Usage: $0 -s <sample_id> -c <caveman.vcf.gz> -p <pindel.vcf.gz> \\
          -x <sbs-dbs_targets.tsv> -i <id_targets.tsv> [-m <keep.maf>]

  -s  Sample ID (Tumor_Sample_Barcode)
  -c  Caveman snpflagged VCF (bgzipped, tabix-indexed)
  -p  Pindel snpflagged VCF (bgzipped, tabix-indexed)
  -x  SBS/DBS targets TSV from maf2targets.R
  -i  ID targets TSV from maf2targets.R
  -m  Optional keep MAF for the subcohort (enables warn-only sanity check)
EOF
    exit 1
}

sample=""
caveman=""
pindel=""
sbs_targets=""
id_targets=""
maf=""

while getopts ":s:c:p:x:i:m:h" flag; do
    case "${flag}" in
        s) sample=$OPTARG ;;
        c) caveman=$OPTARG ;;
        p) pindel=$OPTARG ;;
        x) sbs_targets=$OPTARG ;;
        i) id_targets=$OPTARG ;;
        m) maf=$OPTARG ;;
        h|*) usage ;;
    esac
done

if [[ -z $sample || -z $caveman || -z $pindel || -z $sbs_targets || -z $id_targets ]]; then
    usage
fi

sbs_tmp="${sample}.sbs.vcf.gz"
id_tmp="${sample}.id.vcf.gz"
out_vcf="${sample}.vcf"

if [[ -s $sbs_targets ]]; then
    bcftools view -f PASS --targets-file "$sbs_targets" -O z -o "$sbs_tmp" "$caveman"
else
    bcftools view -h -O z -o "$sbs_tmp" "$caveman"
fi

if [[ -s $id_targets ]]; then
    bcftools view -f PASS --targets-file "$id_targets" -O z -o "$id_tmp" "$pindel"
else
    bcftools view -h -O z -o "$id_tmp" "$pindel"
fi

tabix -p vcf "$sbs_tmp"
tabix -p vcf "$id_tmp"

bcftools concat -a "$sbs_tmp" "$id_tmp" > "$out_vcf"
rm -f "$sbs_tmp" "$sbs_tmp".tbi "$id_tmp" "$id_tmp".tbi

if [[ -n $maf && -s $maf ]]; then
    chr_col=$(awk -v RS='\t' '/Chromosome/{print NR; exit}' "$maf")
    end_col=$(awk -v RS='\t' '/End_Position/{print NR; exit}' "$maf")
    ref_col=$(awk -v RS='\t' '/Reference_Allele/{print NR; exit}' "$maf")
    alt_col=$(awk -v RS='\t' '/Tumor_Seq_Allele2/{print NR; exit}' "$maf")
    pos_col=$(awk -v RS='\t' '/POS_VCF/{print NR; exit}' "$maf")
    bc_col=$(awk -v RS='\t' '/Tumor_Sample_Barcode/{print NR; exit}' "$maf")

    expected=$(awk -v s="$sample" -v bc="$bc_col" -v c1="$chr_col" -v c2="$end_col" \
                   -v c3="$ref_col" -v c4="$alt_col" -v c5="$pos_col" \
                   'BEGIN{FS=OFS="\t"} $bc==s && $0 ~ /\y(SNP|DNP|INS|DEL)\y/ {print $c1,$c2,$c3,$c4,$c5}' \
                   "$maf" | sort -u | wc -l)
    observed=$(grep -vc '^#' "$out_vcf" || true)

    if [[ $expected -ne $observed ]]; then
        unique_sites=$(grep -v '^#' "$out_vcf" | cut -f1,2 | sort -u | wc -l)
        echo "WARNING: ${sample} - MAF vs VCF variant counts differ (MAF=${expected}, VCF=${observed}, unique_sites=${unique_sites})" >&2
    else
        echo "OK: ${sample} - ${expected} variants match between MAF and VCF" >&2
    fi
fi
