#!/bin/bash

# Submit with: bsub -q long -n 6 -M60000 -J B1b_AB_regioneR_perbp -o $MAINDIR/HiC/B1b_AB_regioneR_perbp/B1b_AB_regioneR_perbp.out -e $MAINDIR/HiC/B1b_AB_regioneR_perbp/B1b_AB_regioneR_perbp.err -R "select[mem>60000] rusage[mem=60000] span[hosts=1]" bash B1b_AB_regioneR_perbp.sh

### Set directories
OUTDIR=$MAINDIR/HiC/D6_general_corr
cd $OUTDIR

### Convert A% from chromosomes to Merians
echo -e "species\tMerian\ttot_bp_A" > Aperc_per_Merian.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    comp=comp_beds/${SP}_40kb_A_norm_nosex.bed
    merian=merians/synteny_Merian_joined_${SP}.bed
    # Intersect: for each A-comp bin, get the Merian it overlaps
    # -wa prints the original A bin, -wb prints the Merian region
    bedtools intersect -a "$comp" -b "$merian" -wa -wb \
    | awk -v sp="$SP" '{
        merian = $7          # 4th column of the merian bed = M label
        bp     = $3 - $2     # size of the A-comp bin
        tot[merian] += bp
    }
    END {
        for (m in tot) print sp "\t" m "\t" tot[m]
    }'
done >> Aperc_per_Merian.tsv

### Convert B% from chromosomes to Merians
echo -e "species\tMerian\ttot_bp_B" > Bperc_per_Merian.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    comp=comp_beds/${SP}_40kb_B_norm_nosex.bed
    merian=merians/synteny_Merian_joined_${SP}.bed
    # Intersect: for each A-comp bin, get the Merian it overlaps
    # -wa prints the original A bin, -wb prints the Merian region
    bedtools intersect -a "$comp" -b "$merian" -wa -wb \
    | awk -v sp="$SP" '{
        merian = $7          # 4th column of the merian bed = M label
        bp     = $3 - $2     # size of the B-comp bin
        tot[merian] += bp
    }
    END {
        for (m in tot) print sp "\t" m "\t" tot[m]
    }'
done >> Bperc_per_Merian.tsv

### Convert TAD size from chromosomes to Merians
echo -e "species\tMerian\tavg_TAD_length" > TADs_length_per_Merian.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    tads=../A3b_callTADs_norm/${SP}/${SP}_10kb_TADs_norm_domains.bed
    merian=merians/synteny_Merian_joined_${SP}.bed
    # -f 1.0: TAD must be fully contained within a Merian
    bedtools intersect -a "$tads" -b "$merian" -wa -wb -f 1.0 \
    | awk 'BEGIN{OFS="\t"} {
        tad   = $1":"$2"-"$3
        m     = $13            # Merian label (4th col of merian bed, offset by 9 tads cols)
        len   = $3 - $2
        if (tad in tad_merian) {
            tad_merian[tad] = "MULTI"
        } else {
            tad_merian[tad] = m
            tad_len[tad]    = len
        }
    }
    END {
        for (tad in tad_merian) {
            if (tad_merian[tad] == "MULTI") continue
            m = tad_merian[tad]
            sum[m]   += tad_len[tad]
            count[m] += 1
        }
        for (m in sum) print sp, m, sum[m]/count[m]
    }' sp="$SP" OFS="\t"
done >> TADs_length_per_Merian.tsv

### Convert PC1 score from chromosomes to Merians
# If a window overlaps several Merians, skip it
echo -e "species\tMerian\tavg_PC1" > PC1_per_Merian.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    bg=comp_bedgraphs/${SP}_40kb_AB_norm_nosex.bedgraph
    merian=merians/synteny_Merian_joined_${SP}.bed
    # -wa -wb: get the bedgraph window + the Merian it overlaps
    # -f 1.0: the window must be 100% within the Merian (avoids boundary windows)
    bedtools intersect -a "$bg" -b "$merian" -wa -wb -f 1.0 \
    | awk 'BEGIN{OFS="\t"} {
        win  = $1":"$2"-"$3   # unique ID for this window
        m    = $8             # Merian label (4th col of merian bed)
        score = $4
        # track which merians each window overlaps
        if (win in win_merian) {
            win_merian[win] = "MULTI"   # flag if seen with >1 merian
        } else {
            win_merian[win] = m
            win_score[win]  = score
        }
    }
    END {
        for (win in win_merian) {
            if (win_merian[win] == "MULTI") continue
            m = win_merian[win]
            sum[m]   += win_score[win]
            count[m] += 1
        }
        for (m in sum) print sp, m, sum[m]/count[m]
    }' sp="$SP" OFS="\t"
done >> PC1_per_Merian.tsv

### Convert Insu score from chromosomes to Merians
# If a window overlaps several Merians, skip it
echo -e "species\tMerian\tavg_insulation" > TADs_insulation_per_Merian.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    bg=tads_bedgraphs/${SP}_10kb_TADs_norm_score_nosex.bedgraph
    merian=merians/synteny_Merian_joined_${SP}.bed
    # -wa -wb: get the bedgraph window + the Merian it overlaps
    # -f 1.0: the window must be 100% within the Merian (avoids boundary windows)
    bedtools intersect -a "$bg" -b "$merian" -wa -wb -f 1.0 \
    | awk 'BEGIN{OFS="\t"} {
        win  = $1":"$2"-"$3   # unique ID for this window
        m    = $8             # Merian label (4th col of merian bed)
        score = $4
        # track which merians each window overlaps
        if (win in win_merian) {
            win_merian[win] = "MULTI"   # flag if seen with >1 merian
        } else {
            win_merian[win] = m
            win_score[win]  = score
        }
    }
    END {
        for (win in win_merian) {
            if (win_merian[win] == "MULTI") continue
            m = win_merian[win]
            sum[m]   += win_score[win]
            count[m] += 1
        }
        for (m in sum) print sp, m, sum[m]/count[m]
    }' sp="$SP" OFS="\t"

done >> TADs_insulation_per_Merian.tsv

### Conver ratio inter/intra from chromosomes to Merian
echo -e "species\tMerian\tavg_ratio_inter_intra" > ratio_inter_intra_per_Merian.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    ratio=../D3_inter_intra/${SP}/${SP}_inter_intra_perbin.tsv
    merian=merians/synteny_Merian_joined_${SP}.bed
    # Skip header, intersect bins with Merian regions
    tail -n +2 "$ratio" \
    | bedtools intersect -a stdin -b "$merian" -wa -wb -f 1.0 \
    | awk 'BEGIN{OFS="\t"} {
        bin   = $1":"$2"-"$3
        m     = $10             # Merian label (4th col of merian bed, offset by 6 ratio cols)
        score = $6              # ratio_inter_intra
        if (score == "NA") next
        if (bin in bin_merian) {
            bin_merian[bin] = "MULTI"
        } else {
            bin_merian[bin] = m
            bin_score[bin]  = score
        }
    }
    END {
        for (bin in bin_merian) {
            if (bin_merian[bin] == "MULTI") continue
            m = bin_merian[bin]
            sum[m]   += bin_score[bin]
            count[m] += 1
        }
        for (m in sum) print sp, m, sum[m]/count[m]
    }' sp="$SP" OFS="\t"
done >> ratio_inter_intra_per_Merian.tsv

### VNE per Merian
echo -e "species\tMerian\tavg_S" > S_per_Merian.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    bed=../D5_VNE/${SP}/PYTHON/${SP}_S_per_40kb.bed
    merian=merians/synteny_Merian_joined_${SP}.bed
    bedtools intersect -a "$bed" -b "$merian" -wa -wb \
    | awk 'BEGIN{OFS="\t"} {
        bin   = $1":"$2"-"$3
        m     = $8              # Merian label (4th col of merian bed, offset by 4 bed cols)
        score = $4
        if (bin in bin_merian) {
            bin_merian[bin] = "MULTI"
        } else {
            bin_merian[bin] = m
            bin_score[bin]  = score
        }
    }
    END {
        for (bin in bin_merian) {
            if (bin_merian[bin] == "MULTI") continue
            m = bin_merian[bin]
            sum[m]   += bin_score[bin]
            count[m] += 1
        }
        for (m in sum) print sp, m, sum[m]/count[m]
    }' sp="$SP" OFS="\t"
done >> S_per_Merian.tsv

### Convert R1 interchromosomal contacts from windows to Merians
echo -e "species\tMerian\tavg_R1_interchrom" > R1_interchrom_per_Merian.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    bg=../D1b_R1_contacts_perM/${SP}/${SP}_R1_interchrom_avg_norm_100kb.bedgraph
    merian=merians/synteny_Merian_joined_${SP}.bed
    bedtools intersect -a "$bg" -b "$merian" -wa -wb -f 1.0 \
    | awk 'BEGIN{OFS="\t"} {
        win   = $1":"$2"-"$3
        m     = $8
        score = $4
        if (win in win_merian) {
            win_merian[win] = "MULTI"
        } else {
            win_merian[win] = m
            win_score[win]  = score
        }
    }
    END {
        for (win in win_merian) {
            if (win_merian[win] == "MULTI") continue
            m = win_merian[win]
            sum[m]   += win_score[win]
            count[m] += 1
        }
        for (m in sum) print sp, m, sum[m]/count[m]
    }' sp="$SP" OFS="\t"
done >> R1_interchrom_per_Merian.tsv

### Convert R1 interchromosomal contacts from windows to Merians (no rRNA clusters)
echo -e "species\tMerian\tavg_R1_interchrom_norRNA" > R1_interchrom_per_Merian_norRNA.tsv
for SP in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f1); do
    bg=../D1b_R1_contacts_perM/${SP}/${SP}_LINE_R1_interchrom_avg_no_rRNA_clusters_norm_100kb.bedgraph
    merian=merians/synteny_Merian_joined_${SP}.bed
    bedtools intersect -a "$bg" -b "$merian" -wa -wb -f 1.0 \
    | awk 'BEGIN{OFS="\t"} {
        win   = $1":"$2"-"$3
        m     = $8
        score = $4
        if (win in win_merian) {
            win_merian[win] = "MULTI"
        } else {
            win_merian[win] = m
            win_score[win]  = score
        }
    }
    END {
        for (win in win_merian) {
            if (win_merian[win] == "MULTI") continue
            m = win_merian[win]
            sum[m]   += win_score[win]
            count[m] += 1
        }
        for (m in sum) print sp, m, sum[m]/count[m]
    }' sp="$SP" OFS="\t"
done >> R1_interchrom_per_Merian_norRNA.tsv


### To correct for extant chromosome size and not Merian size,
# Get the size of the chromosome each Merian lives on 
# If fissions, take the average
for f in merians/synteny_Merian_joined_*.bed; do
    SP=$(basename $f | sed 's/synteny_Merian_joined_//;s/\.bed//')
    awk -v sp="$SP" '
    {
        if ($3 > chrom_max[$1]) chrom_max[$1] = $3   # max end coord per chrom
        merian_chrom[$4][$1] = 1                      # which chroms each merian is on
    }
    END {
        for (m in merian_chrom) {
            sum = 0
            n = 0
            for (chrom in merian_chrom[m]) {
                sum += chrom_max[chrom]
                n++
            }
            print sp, m, sum/n    # average chrom size across chroms for this merian
        }
    }' OFS="\t" $f
done > merian_chrom_sizes.tsv
