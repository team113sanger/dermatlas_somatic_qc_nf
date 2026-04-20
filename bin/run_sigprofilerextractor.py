#!/usr/bin/env python
"""Run SigProfilerExtractor on a directory of per-sample VCFs.

Positional args:
  1) vcf_dir      - directory containing *.vcf files (one per sample)
  2) out_dir      - output directory for SigProfiler results
  3) genome_build - reference/opportunity genome (e.g. GRCh38, GRCh37)
  4) seed_file    - optional Seeds.txt for reproducible re-runs ("random" if omitted)
"""

import os
import sys

from SigProfilerExtractor import sigpro as sig


def main() -> int:
    if len(sys.argv) < 4:
        sys.stderr.write(
            "Usage: run_sigprofilerextractor.py <vcf_dir> <out_dir> <genome_build> [seed_file]\n"
        )
        return 1

    vcf_dir = sys.argv[1]
    out_dir = sys.argv[2]
    genome_build = sys.argv[3]
    seed_file = sys.argv[4] if len(sys.argv) >= 5 else "random"

    if not os.path.isdir(vcf_dir):
        sys.exit(f"VCF directory does not exist: {vcf_dir}")
    os.makedirs(out_dir, exist_ok=True)

    sig.sigProfilerExtractor(
        "vcf",
        out_dir,
        vcf_dir,
        reference_genome=genome_build,
        opportunity_genome=genome_build,
        context_type="96,DINUC,ID",
        exome=True,
        minimum_signatures=1,
        maximum_signatures=10,
        nmf_replicates=500,
        cpu=12,
        seeds=seed_file,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
