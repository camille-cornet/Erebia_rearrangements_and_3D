#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A2_hicExplorer
OUTDIR=$MAINDIR/HiC/A4b_callAB_cool_strength
cd $OUTDIR

### Using cooltools to calculate PC1 score, then A and B compartment
# Also calculate compartment strength as AA or BB interactions over AB or BA interactions (need expected cis values first)

### Python script for the compartment calling at 40 kb

    cat > "call_compartments.py" <<EOF
#!/usr/bin/env python
import os
import sys
import argparse
import cooler
import cooltools
import cooltools.api.expected
import cooltools.api.eigdecomp
import cooltools.api.saddle
from cooltools.api.saddle import saddle_strength
import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import bioframe

def run_pipeline(sample, outdir, merian, nproc=4):

    cool_path     = f"{outdir}/{sample}/{sample}_mapq30_corr_40kb_clean.cool"
    expected_path = f"{outdir}/{sample}/{sample}_expected_40kb.tsv"
    eigs_prefix   = f"{outdir}/{sample}/{sample}_AB_40kb"
    saddle_prefix = f"{outdir}/{sample}/{sample}_saddle_40kb"

    # 1 — Open cooler
    clr = cooler.Cooler(cool_path)

    # 2 — Create dummy weight column of ones (because the matrix was corrected with HiCExplorer)
    with h5py.File(cool_path, "r+") as f:
        n_bins = f['bins/chrom'].shape[0]
        f['bins'].create_dataset('weight', data=np.ones(n_bins, dtype=np.float64))
    clr = cooler.Cooler(cool_path)

    # 3 — Compute expected values
    expected = cooltools.api.expected.expected_cis(clr, nproc=nproc, clr_weight_name='weight')
    expected.to_csv(expected_path, sep="\t", index=False)

    # 4 — Eigenvectors (PC1 score and compartments)
    # Load gene density and align to cooler bins
    bins = clr.bins()[:][['chrom', 'start', 'end']]
    gene_density = pd.read_csv(
        f"{outdir}/{sample}/{sample}_gene_density_40kb.bedgraph",
        sep="\t", header=None,
        names=["chrom", "start", "end", "gene_density"]
    )
    phasing_track = bins.merge(gene_density, on=["chrom", "start", "end"], how="left")

    # Calling compartment
    _, eigvec_df = cooltools.api.eigdecomp.eigs_cis(clr, phasing_track=phasing_track, n_eigs=3, clr_weight_name='weight')
    eigvec_df.to_csv(f"{eigs_prefix}.cis.vecs.tsv", sep="\t", index=False)

    # 5 — Compartment strength (saddle-based)
    # Defining function to be later applied to both genome-wide, per chromosome and per Merian
    e1_track = eigvec_df[["chrom", "start", "end", "E1"]]
    N_BINS, QRANGE, EXTENT = 50, (0.02, 0.98), 5  # parameters to choose the most A and most B corners of the saddle plot to compute strength
    mids = (e1_track["start"].values + e1_track["end"].values) // 2
    chrom_order = {c: i for i, c in enumerate(clr.chromnames)}     # necessary for Merian matching

    # Function to assign bins to chromosomes or Merians
    def subset_track(view_df):
        keep = np.zeros(len(e1_track), dtype=bool)
        for _, r in view_df.iterrows():
            keep |= ((e1_track["chrom"].values == r["chrom"]) &
                        (mids >= r["start"]) & (mids < r["end"]))
        return e1_track[keep]

    # Function to run cooltools.saddle and saddle_strength on the regions matched with subset_track
    def unit_strength(view_df, expected_df, label):
        view_df = view_df.reset_index(drop=True)
        track = subset_track(view_df)
        S, C = cooltools.api.saddle.saddle(
            clr, expected_df, track, contact_type="cis",
            n_bins=N_BINS, qrange=QRANGE, view_df=view_df,
            clr_weight_name="weight", expected_value_col="count.avg",
            verbose=False,
        )
        return saddle_strength(S, C)[EXTENT]

    # Call compartment strength function per chromosome and genome-wide
    full_view = pd.DataFrame({"chrom": clr.chromnames, "start": 0,
                                "end": [clr.chromsizes[c] for c in clr.chromnames],
                                "name": clr.chromnames})

    chrom_rows = [{"sample": sample, "chrom": "genome_wide", "res_kb": "40kb",
                    "strength": unit_strength(full_view, expected, "genome_wide")}]

    for chrom in clr.chromnames:
        view = full_view[full_view["chrom"] == chrom]
        exp_sub = expected[expected["region1"] == chrom]
        chrom_rows.append({"sample": sample, "chrom": chrom, "res_kb": "40kb",
                            "strength": unit_strength(view, exp_sub, f"chrom {chrom}")})

    pd.DataFrame(chrom_rows).to_csv(f"{saddle_prefix}.strength.tsv", sep="\t", index=False)

    # Call compartment strength function per Merian
    merian_df = pd.read_csv(merian, sep="\t", header=None,
                            names=["chrom", "start", "end", "merian"])
    merian_df["end"] = [min(e, clr.chromsizes[c])
                        for c, e in zip(merian_df["chrom"], merian_df["end"])]
    merian_df = merian_df[merian_df["end"] > merian_df["start"]] # just to make sure coordinates match
    merian_df["_ord"] = merian_df["chrom"].map(chrom_order) 
    merian_df = merian_df.sort_values(["_ord", "start"]).drop(columns="_ord").reset_index(drop=True) # sort by chromosome order
    merian_df["name"] = merian_df["merian"] + "_" + merian_df.groupby("merian").cumcount().astype(str)  # needed in case of Merian fission

    merian_view = merian_df[["chrom", "start", "end", "name"]]
    expected_merian = cooltools.api.expected.expected_cis(
        clr, view_df=merian_view, nproc=nproc, clr_weight_name="weight")

    merian_rows = []
    for m, sub in merian_df.groupby("merian"):
        view = sub[["chrom", "start", "end", "name"]]
        exp_sub = expected_merian[expected_merian["region1"].isin(sub["name"])]
        merian_rows.append({"sample": sample, "merian": m, "res_kb": "40kb",
                            "strength": unit_strength(view, exp_sub, f"Merian {m}")})

    pd.DataFrame(merian_rows).to_csv(f"{saddle_prefix}.merian_strength.tsv", sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cooltools compartment calling pipeline")
    parser.add_argument("--sample",  required=True, help="Species code")
    parser.add_argument("--merian", required=True, help="Path to Merian correspondance bed (chrom start end merian)")
    parser.add_argument("--outdir",  required=True, help="Output directory")
    parser.add_argument("--nproc",   type=int, default=4, help="Number of processors (default: 4)")
    args = parser.parse_args()

    run_pipeline(args.sample, args.outdir, args.merian, args.nproc)

EOF


compartments() {
    SAMPLE=$1
    mkdir $OUTDIR/${SAMPLE}
    SIZES=$MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes

    source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
    conda activate cooler_env

    # Need to clean the matrix because of HiC explorer bin sizes
    HICCORR=$INDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_10kb.cool
    cooler dump --join $HICCORR | cooler load -f bg2 --count-as-float $SIZES:10000 - $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_10kb_clean.cool
    HICCORR=$INDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_40kb.cool
    cooler dump --join $HICCORR | cooler load -f bg2 --count-as-float $SIZES:40000 - $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_40kb_clean.cool

    # Make gene bedgraph, for AB phasing and plotting
    HICCLEAN=$OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_10kb_clean.cool
    cooler dump --table bins $HICCLEAN \
    | awk '{print $1"\t"$2"\t"$3}' > $OUTDIR/${SAMPLE}/${SAMPLE}_windows_10kb.bed
    HICCLEAN=$OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_40kb_clean.cool
    cooler dump --table bins $HICCLEAN \
    | awk '{print $1"\t"$2"\t"$3}' > $OUTDIR/${SAMPLE}/${SAMPLE}_windows_40kb.bed
    bedtools coverage -a $OUTDIR/${SAMPLE}/${SAMPLE}_windows_10kb.bed -b $OUTDIR/genes/${SAMPLE}.bed -counts \
    | awk '{print $1, $2, $3, $4}' OFS='\t' > $OUTDIR/${SAMPLE}/${SAMPLE}_gene_density_10kb.bedgraph
    bedtools coverage -a $OUTDIR/${SAMPLE}/${SAMPLE}_windows_40kb.bed -b $OUTDIR/genes/${SAMPLE}.bed -counts \
    | awk '{print $1, $2, $3, $4}' OFS='\t' > $OUTDIR/${SAMPLE}/${SAMPLE}_gene_density_40kb.bedgraph

    # Merian correspondence file to directly compute compartment strength per Merian
    MERIAN=$MAINDIR/HiC/D6_general_corr/merians/synteny_Merian_joined_${SAMPLE}.bed

    # Run cooltools to call compartments inside python script
    python call_compartments.py \
        --sample $SAMPLE \
        --merian $MERIAN \
        --outdir $OUTDIR \
        --nproc 4
    conda deactivate

    # Aggregate for plotting in HiGlass
    conda activate clodius_env
    clodius aggregate bedfile \
        --output-file $OUTDIR/genes/${SAMPLE}.multires.bed \
        --chromsizes-filename $MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes \
        --no-header $OUTDIR/genes/${SAMPLE}.bed
    clodius aggregate bedgraph          \
        --output-file $OUTDIR/${SAMPLE}/${SAMPLE}_gene_density_10kb.hitile \
        --chromosome-col 1              \
        --from-pos-col 2                \
        --to-pos-col 3                  \
        --value-col 4                   \
        --chromsizes-filename $MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes \
        --nan-value NA                  \
        --no-header $OUTDIR/${SAMPLE}/${SAMPLE}_gene_density_10kb.bedgraph
    awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}' \
        $OUTDIR/${SAMPLE}/${SAMPLE}_AB_40kb.cis.vecs.tsv \
        | awk 'NF==4' > $OUTDIR/${SAMPLE}/${SAMPLE}_40kb_AB.bedgraph
    clodius aggregate bedgraph          \
        --output-file $OUTDIR/${SAMPLE}/${SAMPLE}b_40kb_AB.hitile \
        --chromosome-col 1              \
        --from-pos-col 2                \
        --to-pos-col 3                  \
        --value-col 4                   \
        --chromsizes-filename $MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes \
        --nan-value NA                  \
        --no-header $OUTDIR/${SAMPLE}/${SAMPLE}_40kb_AB.bedgraph

    conda deactivate

    # Get separate beds for A and B
    awk '$4 > 0 {print $1, $2, $3}' OFS="\t" $OUTDIR/${SAMPLE}/${SAMPLE}_40kb_AB.bedgraph > $OUTDIR/${SAMPLE}/${SAMPLE}_40kb_A.bed
    awk '$4 < 0 {print $1, $2, $3}' OFS="\t" $OUTDIR/${SAMPLE}/${SAMPLE}_40kb_AB.bedgraph > $OUTDIR/${SAMPLE}/${SAMPLE}_40kb_B.bed

}

### Parallelize over individuals
export INDIR OUTDIR
export -f compartments
parallel --colsep '\t' 'compartments {1}' :::: $MAINDIR/HiC/A2_hicExplorer/species_list.txt


### Get AB compartment percentages per species and per chromosome
while read SAMPLE; do
    BEDGRAPH=$OUTDIR/${SAMPLE}/${SAMPLE}_40kb_AB.bedgraph
    FAI=../../assemblies/chromosomes_only/$SAMPLE.fasta.fai

    awk -v sample="$SAMPLE" \
        -v perchrom="$OUTDIR/AB_per_chrom.tsv" \
        -v perspecies="$OUTDIR/AB_per_species.tsv" '
        NR==FNR { chromsize[$1] = $2; genome += $2; next }
        {
            chrom = $1
            size  = $3 - $2
            seen[chrom] = 1
            if ($4 > 0) { sumA[chrom] += size; gA += size }
            if ($4 < 0) { sumB[chrom] += size; gB += size }
        }
        END {
            for (chrom in seen) {
                total = chromsize[chrom]
                if (total == 0) continue
                a = sumA[chrom]; b = sumB[chrom]
                printf "%s\t%s\t%.4f\t%.4f\t%.4f\n", sample, chrom,
                    a/total*100, b/total*100, (total - a - b)/total*100 >> perchrom
            }
            printf "%s\t%.4f\t%.4f\t%.4f\n", sample,
                gA/genome*100, gB/genome*100,
                (genome - gA - gB)/genome*100 >> perspecies
        }
    ' $FAI $BEDGRAPH

done < $MAINDIR/HiC/A2_hicExplorer/species_list.txt

### Get compartment strength (saddle-based) per species, per chromosome, and per Merian
echo -e "species\tchrom\tstrength" > $OUTDIR/compartment_strength_per_chrom.tsv
echo -e "species\tstrength" > $OUTDIR/compartment_strength_per_species.tsv
echo -e "species\tmerian\tstrength" > $OUTDIR/compartment_strength_per_merian.tsv

while read SAMPLE; do
    f="$OUTDIR/${SAMPLE}/${SAMPLE}_saddle_40kb.strength.tsv"
    tail -n +2 "$f" | awk -F'\t' -v OFS='\t' '$2!="genome_wide"{print $1,$2,$4}' >> $OUTDIR/compartment_strength_per_chrom.tsv
    tail -n +2 "$f" | awk -F'\t' -v OFS='\t' '$2=="genome_wide"{print $1,$4}' >> $OUTDIR/compartment_strength_per_species.tsv

    m="$OUTDIR/${SAMPLE}/${SAMPLE}_saddle_40kb.merian_strength.tsv"
    tail -n +2 "$m" | awk -F'\t' -v OFS='\t' '{print $1,$2,$4}' >> $OUTDIR/compartment_strength_per_merian.tsv
done < $MAINDIR/HiC/A2_hicExplorer/species_list.txt

### SOFTWARE VERSIONS
# cooler v0.10.4
# cooltools v0.7.1
# clodius v0.14.3
