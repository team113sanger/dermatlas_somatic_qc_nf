#!/usr/bin/env python3
"""Stage pre-converted PNG plots as MultiQC custom content.

Plot dirs arrive already namespaced as ``{subcohort}__{plotdir}`` from the
upstream CONVERT_PLOTS_TO_PNG process, so every section gets a unique id.
"""

import argparse
import glob
import os
import shutil


PLOT_SECTIONS = {
    "AF_vs_depth_snv_samples": {
        "id": "af_vs_depth_snv",
        "section_name": "AF vs Depth (SNV)",
        "description": "Allele frequency versus read depth for SNVs across samples",
    },
    "AF_vs_depth_indel_samples": {
        "id": "af_vs_depth_indel",
        "section_name": "AF vs Depth (Indel)",
        "description": "Allele frequency versus read depth for indels across samples",
    },
    "mutation_types_barplot_samples": {
        "id": "mutation_types",
        "section_name": "Mutation Types",
        "description": "Distribution of mutation types across samples",
    },
    "mutation_types_barplot_proportion_samples": {
        "id": "mutation_types_proportion",
        "section_name": "Mutation Types (Proportion)",
        "description": "Proportional distribution of mutation types across samples",
    },
    "gene_tileplot": {
        "id": "gene_tileplot",
        "section_name": "Gene Tileplot",
        "description": "Gene-level mutation landscape across samples",
    },
    "AF_vs_depth_recurrent_genes": {
        "id": "af_vs_depth_recurrent_genes",
        "section_name": "AF vs Depth (Recurrent Genes)",
        "description": "Allele frequency versus read depth for recurrently mutated genes",
    },
    "AF_vs_depth_recurrent_sites": {
        "id": "af_vs_depth_recurrent_sites",
        "section_name": "AF vs Depth (Recurrent Sites)",
        "description": "Allele frequency versus read depth for recurrently mutated sites",
    },
}

FILTER_LABELS = {
    "plots_keep_vaf_size_filt_matched": ("keep", "Keep (all PASS)"),
    "plots_keepPA_vaf_size_filt_matched": ("keepPA", "KeepPA (protein-altering)"),
}


def parse_dir_name(dir_name):
    """Split ``{subcohort}__{plotdir}`` into its parts."""
    if "__" not in dir_name:
        return None, dir_name
    sub, _, plotdir = dir_name.partition("__")
    return sub, plotdir


def stage_png(png_path, outdir, section_config, subcohort, plotdir):
    filter_key, filter_label = FILTER_LABELS.get(plotdir, (plotdir, plotdir))
    basename = os.path.splitext(os.path.basename(png_path))[0]
    unique = f"{subcohort}_{filter_key}_{basename}"
    dest_png = os.path.join(outdir, f"{unique}_mqc.png")
    yaml_path = os.path.join(outdir, f"{unique}_mqc.yaml")

    shutil.copyfile(png_path, dest_png)
    with open(yaml_path, "w") as fh:
        fh.write(f"id: '{subcohort}_{filter_key}_{section_config['id']}'\n")
        fh.write(f"parent_id: '{subcohort}'\n")
        fh.write(f"parent_name: '{subcohort}'\n")
        fh.write(f"section_name: '{section_config['section_name']} — {filter_label}'\n")
        fh.write(f"description: '{section_config['description']} ({filter_label})'\n")
        fh.write("plot_type: 'image'\n")


def main():
    parser = argparse.ArgumentParser(description="Stage PNG plots for MultiQC custom content")
    parser.add_argument("--plot-dirs", nargs="+", required=True,
                        help="Directories named {subcohort}__{plotdir} containing PNGs")
    parser.add_argument("--outdir", required=True, help="Output dir for MultiQC inputs")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for plot_dir in args.plot_dirs:
        if not os.path.isdir(plot_dir):
            continue
        subcohort, plotdir = parse_dir_name(os.path.basename(os.path.normpath(plot_dir)))
        if subcohort is None:
            continue
        for png_file in glob.glob(os.path.join(plot_dir, "*.png")):
            basename = os.path.splitext(os.path.basename(png_file))[0]
            if basename in PLOT_SECTIONS:
                stage_png(png_file, args.outdir, PLOT_SECTIONS[basename], subcohort, plotdir)


if __name__ == "__main__":
    main()
