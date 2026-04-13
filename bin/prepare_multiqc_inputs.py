#!/usr/bin/env python3
"""Prepare pipeline outputs for MultiQC custom content.

Reformats TSV files with YAML front-matter headers expected by MultiQC's
custom content module, and converts PDF plots to PNG for embedding.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys


# Mapping from PDF base name to MultiQC section config
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


def write_yaml_header(fh, section_id, section_name, description, plot_type, pconfig=None):
    """Write MultiQC YAML front-matter comment block."""
    fh.write(f"# id: '{section_id}'\n")
    fh.write(f"# section_name: '{section_name}'\n")
    fh.write(f"# description: '{description}'\n")
    fh.write(f"# plot_type: '{plot_type}'\n")
    if pconfig:
        fh.write("# pconfig:\n")
        for key, val in pconfig.items():
            fh.write(f"#     {key}: '{val}'\n")


def prepare_tmb(tmb_files, outdir):
    """Combine TMB TSV files and add MultiQC bargraph header."""
    outpath = os.path.join(outdir, "tmb_mqc.tsv")
    rows = []
    for tmb_file in tmb_files:
        with open(tmb_file) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    rows.append(line)

    if not rows:
        return

    with open(outpath, "w") as fh:
        write_yaml_header(
            fh,
            section_id="tmb",
            section_name="Tumour Mutation Burden",
            description="Mutations per megabase for each sample",
            plot_type="bargraph",
            pconfig={
                "id": "tmb_bargraph",
                "title": "Tumour Mutation Burden",
                "ylab": "Mutations per Mb",
            },
        )
        fh.write("Sample\tMutations_per_Mb\n")
        for row in rows:
            fh.write(row + "\n")


def prepare_table_tsv(input_path, outdir, section_id, section_name, description):
    """Copy a TSV file with MultiQC table header added."""
    basename = os.path.splitext(os.path.basename(input_path))[0]
    outpath = os.path.join(outdir, f"{basename}_mqc.tsv")
    with open(input_path) as fin, open(outpath, "w") as fout:
        write_yaml_header(fout, section_id, section_name, description, "table")
        shutil.copyfileobj(fin, fout)


def convert_pdf_to_png(pdf_path, outdir, section_config):
    """Convert a PDF plot to PNG and create a companion YAML descriptor."""
    basename = os.path.splitext(os.path.basename(pdf_path))[0]
    png_path = os.path.join(outdir, f"{basename}_mqc.png")
    yaml_path = os.path.join(outdir, f"{basename}_mqc.yaml")

    # Try pdftoppm (poppler-utils) first, fall back to ghostscript
    converted = False
    if shutil.which("pdftoppm"):
        result = subprocess.run(
            ["pdftoppm", "-png", "-r", "150", "-singlefile", pdf_path, png_path.replace(".png", "")],
            capture_output=True,
        )
        if result.returncode == 0:
            converted = True
    if not converted and shutil.which("gs"):
        result = subprocess.run(
            [
                "gs", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m",
                "-r150", f"-sOutputFile={png_path}", pdf_path,
            ],
            capture_output=True,
        )
        if result.returncode == 0:
            converted = True

    if not converted:
        print(f"WARNING: Could not convert {pdf_path} to PNG (no pdftoppm or gs available)", file=sys.stderr)
        return

    # Write companion YAML descriptor for MultiQC custom content image
    with open(yaml_path, "w") as fh:
        fh.write(f"id: '{section_config['id']}'\n")
        fh.write(f"section_name: '{section_config['section_name']}'\n")
        fh.write(f"description: '{section_config['description']}'\n")
        fh.write("plot_type: 'image'\n")


def find_files(paths, pattern):
    """Find files matching a glob pattern across multiple paths."""
    results = []
    for p in paths:
        if os.path.isdir(p):
            results.extend(glob.glob(os.path.join(p, pattern)))
        elif os.path.isfile(p) and glob.fnmatch.fnmatch(os.path.basename(p), pattern):
            results.append(p)
    return sorted(results)


def main():
    parser = argparse.ArgumentParser(description="Prepare MultiQC custom content inputs")
    parser.add_argument("--qc-tsvs", nargs="+", required=True, help="QC TSV files from QC_VARIANTS")
    parser.add_argument("--tmb-tsvs", nargs="+", required=True, help="TMB TSV files from CALCULATE_SAMPLE_TMB")
    parser.add_argument("--plot-dirs", nargs="+", required=True, help="Plot directories from QC_VARIANTS")
    parser.add_argument("--outdir", required=True, help="Output directory for MultiQC inputs")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Prepare TMB bargraph
    tmb_files = [f for f in args.tmb_tsvs if os.path.isfile(f)]
    prepare_tmb(tmb_files, args.outdir)

    # 2. Prepare recurrent genes table
    for tsv in args.qc_tsvs:
        basename = os.path.basename(tsv)
        if "top_recurrently_mutated_genes" in basename:
            prepare_table_tsv(
                tsv, args.outdir,
                section_id="recurrent_genes",
                section_name="Recurrently Mutated Genes",
                description="Top recurrently mutated genes across the cohort",
            )
        elif "top_recurrently_mutated_sites" in basename:
            prepare_table_tsv(
                tsv, args.outdir,
                section_id="recurrent_sites",
                section_name="Recurrently Mutated Sites",
                description="Top recurrently mutated sites across the cohort",
            )

    # 3. Convert PDF plots to PNG
    for plot_dir in args.plot_dirs:
        if not os.path.isdir(plot_dir):
            continue
        for pdf_file in glob.glob(os.path.join(plot_dir, "*.pdf")):
            basename = os.path.splitext(os.path.basename(pdf_file))[0]
            if basename in PLOT_SECTIONS:
                convert_pdf_to_png(pdf_file, args.outdir, PLOT_SECTIONS[basename])


if __name__ == "__main__":
    main()
