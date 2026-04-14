#!/usr/bin/env python3
"""Stage plots + TMB for a single MultiQC report grouped by subcohort.

Writes a dynamic `dynamic_mqc_config.yaml` declaring every image as
`custom_data` with matching `sp` patterns, plus one TMB bargraph per
subcohort. Section ordering is controlled via `report_section_order`.
"""

import argparse
import glob
import os
import shutil


PLOT_SECTIONS = {
    "AF_vs_depth_snv_samples": ("af_vs_depth_snv", "AF vs Depth (SNV)",
                                "Allele frequency versus read depth for SNVs"),
    "AF_vs_depth_indel_samples": ("af_vs_depth_indel", "AF vs Depth (Indel)",
                                  "Allele frequency versus read depth for indels"),
    "mutation_types_barplot_samples": ("mutation_types", "Mutation Types",
                                       "Distribution of mutation types across samples"),
    "mutation_types_barplot_proportion_samples": ("mutation_types_prop", "Mutation Types (Proportion)",
                                                  "Proportional distribution of mutation types"),
    "gene_tileplot": ("gene_tileplot", "Gene Tileplot",
                      "Gene-level mutation landscape across samples"),
    "AF_vs_depth_recurrent_genes": ("af_vs_depth_rec_genes", "AF vs Depth (Recurrent Genes)",
                                    "AF vs depth for recurrently mutated genes"),
    "AF_vs_depth_recurrent_sites": ("af_vs_depth_rec_sites", "AF vs Depth (Recurrent Sites)",
                                    "AF vs depth for recurrently mutated sites"),
}

PLOT_ORDER = [v[0] for v in PLOT_SECTIONS.values()]

FILTER_LABELS = {
    "plots_keep_vaf_size_filt_matched": ("keep", "Keep (all PASS)"),
    "plots_keepPA_vaf_size_filt_matched": ("keepPA", "KeepPA (protein-altering)"),
}
FILTER_ORDER = ["keep", "keepPA"]


def parse_dir_name(dir_name):
    if "__" not in dir_name:
        return None, dir_name
    sub, _, plotdir = dir_name.partition("__")
    return sub, plotdir


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot-dirs", nargs="+", required=True)
    parser.add_argument("--tmb-files", nargs="*", default=[],
                        help="Per-subcohort TMB TSVs named <subcohort>.tmb.tsv")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    custom_data = {}
    sp = {}
    section_order = []
    subcohorts_seen = []

    # TMB bargraphs: one per (subcohort, filter). Filenames look like "sub__keepPA.tmb.tsv"
    for tmb in sorted(args.tmb_files):
        stem = os.path.basename(tmb).replace(".tmb.tsv", "")
        if "__" not in stem:
            continue
        sub, _, filter_key = stem.partition("__")
        filter_label = dict(FILTER_LABELS.values()).get(filter_key, filter_key)
        if sub not in subcohorts_seen:
            subcohorts_seen.append(sub)
        section_id = f"{sub}_{filter_key}_tmb"
        dest = os.path.join(args.outdir, f"{section_id}_mqc.tsv")
        with open(tmb) as fin, open(dest, "w") as fout:
            fout.write(f"# id: '{section_id}'\n")
            fout.write(f"# parent_id: '{sub}'\n")
            fout.write(f"# parent_name: 'Subcohort: {sub}'\n")
            fout.write(f"# section_name: 'Tumour Mutation Burden — {filter_label}'\n")
            fout.write("# description: 'Mutations per megabase per sample'\n")
            fout.write("# plot_type: 'bargraph'\n")
            fout.write("# pconfig:\n")
            fout.write(f"#     id: '{section_id}_plot'\n")
            fout.write(f"#     title: 'Tumour Mutation Burden ({filter_label})'\n")
            fout.write("#     ylab: 'Mutations per Mb'\n")
            fout.write("Sample\tMutations_per_Mb\n")
            for line in fin:
                line = line.strip()
                if line:
                    fout.write(line + "\n")
        section_order.append(section_id)

    # Images: copy + record custom_data/sp entry per PNG
    for plot_dir in sorted(args.plot_dirs):
        if not os.path.isdir(plot_dir):
            continue
        sub, plotdir = parse_dir_name(os.path.basename(os.path.normpath(plot_dir)))
        if sub is None:
            continue
        if sub not in subcohorts_seen:
            subcohorts_seen.append(sub)
        filter_key, filter_label = FILTER_LABELS.get(plotdir, (plotdir, plotdir))

        for png in sorted(glob.glob(os.path.join(plot_dir, "*.png"))):
            stem = os.path.splitext(os.path.basename(png))[0]
            if stem not in PLOT_SECTIONS:
                continue
            plot_id, plot_name, plot_desc = PLOT_SECTIONS[stem]
            section_id = f"{sub}_{filter_key}_{plot_id}"
            dest = os.path.join(args.outdir, f"{section_id}.png")
            shutil.copyfile(png, dest)
            custom_data[section_id] = {
                "parent_id": sub,
                "parent_name": f"Subcohort: {sub}",
                "section_name": f"{plot_name} — {filter_label}",
                "description": plot_desc,
                "plot_type": "image",
            }
            sp[section_id] = {"fn": f"{section_id}.png"}
            section_order.append(section_id)

    # Order sections: all of subcohort A (TMB → keep plots → keepPA plots), then B, ...
    def rank(section_id):
        for s_idx, sub in enumerate(subcohorts_seen):
            for f_idx, f in enumerate(FILTER_ORDER):
                if section_id == f"{sub}_{f}_tmb":
                    return (s_idx, f_idx, -1)
                for p_idx, pid in enumerate(PLOT_ORDER):
                    if section_id == f"{sub}_{f}_{pid}":
                        return (s_idx, f_idx, p_idx)
        return (99, 99, 99)

    section_order.sort(key=rank)

    config_path = os.path.join(args.outdir, "dynamic_mqc_config.yaml")
    with open(config_path, "w") as fh:
        fh.write("custom_data:\n")
        for sid, cfg in custom_data.items():
            fh.write(f"  {sid}:\n")
            for k, v in cfg.items():
                fh.write(f"    {k}: '{v}'\n")
        fh.write("sp:\n")
        for sid, pat in sp.items():
            fh.write(f"  {sid}:\n")
            fh.write(f"    fn: '{pat['fn']}'\n")
        fh.write("report_section_order:\n")
        # MultiQC reads this as a map of {section: {order: N}}; higher = earlier
        for i, sid in enumerate(section_order):
            fh.write(f"  {sid}:\n")
            fh.write(f"    order: {len(section_order) - i}\n")


if __name__ == "__main__":
    main()
