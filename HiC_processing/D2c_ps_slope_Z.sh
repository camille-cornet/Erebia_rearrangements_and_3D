#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A2_hicExplorer
OUTDIR=$MAINDIR/HiC/D2_ps_slope_Z
cd $OUTDIR

# P(s) slope calculation per chromosome but include the Z this time
### Limiting the distance to consider to 0 - 10 Mb 

PSslope() {
    SP=$1
    mkdir $OUTDIR/${SP}
    cd $OUTDIR/${SP}

    SIZES=$MAINDIR/HiC/D2_ps_slope/${SP}/${SP}_chromsizes.txt
    HICCORR=$MAINDIR/HiC/D2_ps_slope/${SP}/${SP}_mapq30_norm_corr_100kb_clean.cool

    cat > "run_${SP}.py" <<EOF
import numpy as np
import pandas as pd
import cooler
import cooltools

clr = cooler.Cooler("$HICCORR")

DIST_MAX_BP = 10_000_000

# All-chrom view (includes sex chroms)
chromsizes = pd.Series(dict(clr.chromsizes))
view_df = pd.DataFrame({
    "chrom": chromsizes.index,
    "start": 0,
    "end": chromsizes.values,
    "name": chromsizes.index,
})

cvd = cooltools.expected_cis(
    clr=clr, view_df=view_df, smooth=True, aggregate_smoothed=False,
    smooth_sigma=0.1, nproc=1, clr_weight_name=None,
)
cvd.loc[cvd["dist"] < 2, "count.avg.smoothed"] = np.nan

# Per-chrom log-log slope
chrom_slopes = []
for chrom, grp in cvd.groupby("region1"):
    grp_plot = grp[["dist_bp", "count.avg.smoothed"]].dropna()
    grp_plot = grp_plot[
        (grp_plot["count.avg.smoothed"] > 0) & (grp_plot["dist_bp"] <= DIST_MAX_BP)
    ]
    if len(grp_plot) < 3:
        continue
    slope = np.gradient(
        np.log(grp_plot["count.avg.smoothed"]), np.log(grp_plot["dist_bp"])
    )
    chrom_slopes.append(pd.DataFrame({"chrom": chrom, "slope_loglog": slope}))

# Average slope per chrom
avg_slope_all = (
    pd.concat(chrom_slopes)
    .groupby("chrom")["slope_loglog"]
    .mean()
    .reset_index()
    .rename(columns={"slope_loglog": "avg_slope_loglog"})
)
avg_slope_all.insert(0, "species", "$SP")

output_file = f"$OUTDIR/$SP/${SP}_mapq30_avg_slope_loglog_perchrom_allchroms.tsv"
avg_slope_all.to_csv(output_file, sep="\t", index=False)

EOF

    python3 run_${SP}.py
}

### Parallelize over species
export INDIR OUTDIR SP HICCORR
export -f PSslope
parallel --colsep '\t' 'PSslope {1}' :::: $MAINDIR/HiC/A2_hicExplorer/species_list.txt

### Combine avg slope per chrom across all species (all chroms)
header=1
while read SP; do
    f="$OUTDIR/${SP}/${SP}_mapq30_avg_slope_loglog_perchrom_allchroms.tsv"
    if [ "$header" -eq 1 ]; then
        cat "$f"
        header=0
    else
        tail -n +2 "$f"
    fi
done < $MAINDIR/HiC/A2_hicExplorer/species_list.txt \
> $OUTDIR/avg_slope_loglog_perchrom_allchroms_allspecies.tsv

### SOFTWARE VERSIONS
# cooler v0.10.4
# cooltools v0.7.1