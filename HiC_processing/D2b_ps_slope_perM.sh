#!/bin/bash

INDIR=$MAINDIR/HiC/D2_ps_slope
OUTDIR=$MAINDIR/HiC/D2b_ps_slope_perM
MERIANDIR=$MAINDIR/HiC/D6_general_corr/merians
cd $OUTDIR

PSslope() {
    SP=$1
    mkdir -p $OUTDIR/${SP}
    cd $OUTDIR/${SP}

    HICCORR=$INDIR/${SP}/${SP}_mapq30_norm_corr_100kb_clean.cool
    MERIAN_BED=$MERIANDIR/synteny_Merian_joined_${SP}.bed

    cat > "run_${SP}.py" <<EOF
import numpy as np
import pandas as pd
import cooler
import cooltools

clr = cooler.Cooler("$HICCORR")

# Load Merian BED
merian_bed = pd.read_csv("$MERIAN_BED", sep="\t", header=None,
                         names=["chrom", "start", "end", "merian"])

# P(s) slope per Merian — run expected_cis separately per Merian
merian_slopes = []

for merian_name, merian_grp in merian_bed.groupby("merian"):

    # Build view for this Merian only (its intervals are non-overlapping by definition)
    view = merian_grp[["chrom", "start", "end"]].copy().reset_index(drop=True)
    view["name"] = [f"{merian_name}_part{i+1}" if len(view) > 1 else merian_name
                    for i in range(len(view))]

    try:
        cvd = cooltools.expected_cis(
            clr=clr,
            view_df=view,
            smooth=True,
            aggregate_smoothed=False,
            smooth_sigma=0.1,
            nproc=1,
            clr_weight_name=None
        )
    except Exception as e:
        print(f"Skipping {merian_name}: {e}")
        continue

    cvd.loc[cvd['dist'] < 2, 'count.avg.smoothed'] = np.nan

    # Average smoothed contact probability across scaffold parts at each dist_bp
    agg = (cvd.groupby("dist_bp", as_index=False)["count.avg.smoothed"]
              .mean()
              .rename(columns={"count.avg.smoothed": "contact_prob_smoothed"}))
    agg = agg.dropna().query("contact_prob_smoothed > 0").copy()

    if len(agg) < 3:
        print(f"Skipping {merian_name}: too few points after filtering")
        continue

    der = np.gradient(np.log(agg["contact_prob_smoothed"]),
                      np.log(agg["dist_bp"]))
    agg["slope_loglog"] = der
    agg["merian"] = merian_name
    merian_slopes.append(agg)

ps_slope_perMerian_df = pd.concat(merian_slopes)[
    ["merian", "dist_bp", "contact_prob_smoothed", "slope_loglog"]
]
ps_slope_perMerian_df.to_csv(
    f"$OUTDIR/${SP}/${SP}_mapq30_expected_cis_100kb_smooth_slope_perMerian.tsv",
    sep="\t", index=False)
print(f"Done. {ps_slope_perMerian_df['merian'].nunique()} Merians written.")

# Global aggregated slope (across all Merians)
# Collect all per-Merian smoothed curves and average them
global_agg = (ps_slope_perMerian_df
              .groupby("dist_bp", as_index=False)["contact_prob_smoothed"]
              .mean())
global_agg = global_agg.query("contact_prob_smoothed > 0").copy()
der_global = np.gradient(np.log(global_agg["contact_prob_smoothed"]),
                         np.log(global_agg["dist_bp"]))
global_agg["slope_loglog"] = der_global
global_agg.to_csv(
    f"$OUTDIR/${SP}/${SP}_mapq30_expected_cis_100kb_smooth_slope.tsv",
    sep="\t", index=False)

EOF

    python3 run_${SP}.py
}

export INDIR OUTDIR MERIANDIR SP HICCORR
export -f PSslope
parallel --colsep '\t' 'PSslope {1}' :::: $MAINDIR/HiC/A2_hicExplorer/species_list.txt

# Restrict the range of distances from 0 to 10 Mb 
for f in */*_perMerian.tsv; do
    awk -F'\t' 'NR==1 || ($2 <= 10000000)' "$f" \
        > "${f%.tsv}_0_10Mb.tsv"
done

# Create output files with averages per Merian
awk -F'\t' '
NR==1 { next }
FNR==1 { next }
{
    # Extract species from path
    split(FILENAME, path, "/")
    sp = path[1]
    key = sp "\t" $1
    sum[key] += $4
    count[key]++
}
END {
    print "species\tMerian\tavg_slope_loglog_all"
    for (key in sum)
        print key "\t" sum[key]/count[key]
}
' */*_perMerian.tsv > all_species_avg_slope_perMerian.tsv

awk -F'\t' '
NR==1 { next }
FNR==1 { next }
{
    # Extract species from path
    split(FILENAME, path, "/")
    sp = path[1]
    key = sp "\t" $1
    sum[key] += $4
    count[key]++
}
END {
    print "species\tMerian\tavg_slope_loglog_0_10"
    for (key in sum)
        print key "\t" sum[key]/count[key]
}
' */*_perMerian_0_10Mb.tsv > all_species_avg_slope_perMerian_0_10Mb.tsv

### SOFTWARE VERSIONS
# cooltools v0.7.1
