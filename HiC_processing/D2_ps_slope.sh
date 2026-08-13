#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A2_hicExplorer
OUTDIR=$MAINDIR/HiC/D2_ps_slope
cd $OUTDIR

PSslope() {
    SP=$1
    mkdir $OUTDIR/${SP}
    cd $OUTDIR/${SP}
    
    source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
    conda activate cooler_env

    # Need to clean the matrix because of HiC explorer bin sizes
    HICCORR=$INDIR/${SP}/${SP}_mapq30_norm_corr_100kb.cool
    cat $MAINDIR/assemblies/chromosomes_only/${SP}.fasta.fai | cut -f 1,2 > $OUTDIR/${SP}/${SP}_chromsizes.txt
    SIZES=$OUTDIR/${SP}/${SP}_chromsizes.txt
    cooler dump --join $HICCORR | cooler load -f bg2 --count-as-float $SIZES:100000 - $OUTDIR/${SP}/${SP}_mapq30_norm_corr_100kb_clean.cool
    HICCORR=$OUTDIR/${SP}/${SP}_mapq30_norm_corr_100kb_clean.cool

    conda deactivate

    cat > "run_${SP}.py" <<EOF
import numpy as np
import pandas as pd
import cooler
import cooltools

clr = cooler.Cooler("$HICCORR")

# Remove sex chrom
def get_sex_chroms(path, species):
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == species:
                return {p for p in parts[1:] if p != "NONE"}
    return set()

sex_chroms = get_sex_chroms("$OUTDIR/sexchromlist.txt", "$SP")

chromsizes = pd.Series(
    dict(clr.chromsizes),  
).drop(labels=sex_chroms, errors='ignore')

view_df = pd.DataFrame({
    'chrom': chromsizes.index,
    'start': 0,
    'end': chromsizes.values,
    'name': chromsizes.index
})

# Smoothed curve per chromosome and per species (aggregated over chromosomes)
cvd = cooltools.expected_cis(
    clr=clr,
    view_df=view_df,
    smooth=True,
    aggregate_smoothed=True,
    smooth_sigma=0.1,
    nproc=1,
    clr_weight_name=None # because corrected in HiCexplorer
)
# Save
output_file = f"$OUTDIR/${SP}/${SP}_mapq30_expected_cis_100kb_perchrom_smooth.tsv"
cvd.to_csv(output_file, sep="\t", index=False)

# Take a single value for each genomic separation
cvd.loc[cvd['dist'] < 2, 'count.avg.smoothed.agg'] = np.nan
cvd.loc[cvd['dist'] < 2, 'count.avg.smoothed'] = np.nan
cvd_merged = cvd.drop_duplicates(subset=['dist'])[['dist_bp', 'count.avg.smoothed.agg']]

# Calculate slope as the derivative of the curve in log-log space (remove neg values)
cvd_plot = cvd_merged.copy()
cvd_plot = cvd_plot[cvd_plot['count.avg.smoothed.agg'] > 0]
der = np.gradient(
    np.log(cvd_plot['count.avg.smoothed.agg']),
    np.log(cvd_plot['dist_bp']))

# Save one dataframe with smoothed curve and slope per species
ps_slope_df = pd.DataFrame({
    'dist_bp': cvd_plot['dist_bp'],
    'contact_prob_smoothed': cvd_plot['count.avg.smoothed.agg'],
    'slope_loglog': der
})
output_file_slope = f"$OUTDIR/${SP}/${SP}_mapq30_expected_cis_100kb_smooth_slope.tsv"
ps_slope_df.to_csv(output_file_slope, sep="\t", index=False)

### Slope per chromosome (smoothed but not aggregated)
chrom_slopes = []
for chrom, grp in cvd.groupby('region1'):
    grp_plot = grp[['dist_bp', 'count.avg.smoothed']].dropna()
    grp_plot = grp_plot[grp_plot['count.avg.smoothed'] > 0].copy()
    if len(grp_plot) < 3:   # need enough points for gradient
        continue
    der = np.gradient(
        np.log(grp_plot['count.avg.smoothed']),
        np.log(grp_plot['dist_bp']))
    grp_plot['slope_loglog'] = der
    grp_plot['chrom'] = chrom
    chrom_slopes.append(grp_plot)

ps_slope_perchrom_df = pd.concat(chrom_slopes)[['chrom', 'dist_bp', 'count.avg.smoothed', 'slope_loglog']]
ps_slope_perchrom_df = ps_slope_perchrom_df.rename(columns={'count.avg.smoothed': 'contact_prob_smoothed'})

output_file_slope_perchrom = f"$OUTDIR/${SP}/${SP}_mapq30_expected_cis_100kb_smooth_slope_perchrom.tsv"
ps_slope_perchrom_df.to_csv(output_file_slope_perchrom, sep="\t", index=False)

EOF

    python3 run_${SP}.py
}

### Parallelize over individuals
export INDIR OUTDIR SP HICCORR
export -f PSslope
parallel --colsep '\t' 'PSslope {1}' :::: $MAINDIR/HiC/A2_hicExplorer/species_list.txt

### SOFTWARE VERSIONS
# cooler v0.10.4
# cooltools v0.7.1
