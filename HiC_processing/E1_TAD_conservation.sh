#!/bin/bash

# Of the gene pairs that colocalize in a TAD in species1, which proportion also colocalizes in a TAD in species2?
# And of TADs with more than one gene in species1, which % have the exact same genes colocalize in one TAD in species2? 

INDIR=$MAINDIR/HiC/E2_AB_conservation
OUTDIR=$MAINDIR/HiC/E1_TAD_conservation
ORTHODIR=$MAINDIR/Downstream/D1_orthofinder/results/Results_Erebia
GENEDIR=$MAINDIR/Downstream/D3_genespace/bed
TADDIR=$MAINDIR/HiC/A3_callTADs
cd $OUTDIR

# Function to assign genes to TADs
get_tad() {
  bedtools intersect -a $1 -b $2 -wa -wb \
  | awk 'OFS="\t" {
      overlap = (($3 < $7 ? $3 : $7) - ($2 > $5 ? $2 : $5))
      if (overlap > best[$4]) { best[$4] = overlap; tad[$4] = $8 }
    }
    END { for (g in tad) print g"\t"tad[g] }'
}

##### Species pair Eligea and Eottomana
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/C0056.bed $TADDIR/C0055/C0055_10kb_TADs_domains.bed > $OUTDIR/C0055_gene_tad.txt
get_tad $GENEDIR/Eligea.bed $TADDIR/Eligea/Eligea_10kb_TADs_domains.bed > $OUTDIR/Eligea_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $4!="" && $20!="" && $4!~/,/ && $20!~/,/ {print $20"\t"$4}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N5.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/C0055_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/Eligea_gene_tad.txt) \
  > $OUTDIR/Eott_Elig_gene_tad.tsv
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (C0055): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (Eligea): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (Eligea): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (C0055): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (C0055): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (Eligea): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Eott_Elig_gene_tad.tsv > $OUTDIR/Eott_Elig_TAD_summary.txt

##### Species pair Egraucasica and Eottomana
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/C0056.bed $TADDIR/C0055/C0055_10kb_TADs_domains.bed > $OUTDIR/C0055_gene_tad.txt
get_tad $GENEDIR/C0111.bed $TADDIR/C0100/C0100_10kb_TADs_domains.bed > $OUTDIR/C0100_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $4!="" && $9!="" && $4!~/,/ && $9!~/,/ {print $9"\t"$4}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N9.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/C0055_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/C0100_gene_tad.txt) \
  > $OUTDIR/Eott_Egra_gene_tad.tsv
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (C0055): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (C0100): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (C0100): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (C0055): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (C0055): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (C0100): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Eott_Egra_gene_tad.tsv > $OUTDIR/Eott_Egra_TAD_summary.txt

##### Species pair Erondoui and Ecassioides
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/X3529.bed $TADDIR/X3531/X3531_10kb_TADs_domains.bed > $OUTDIR/X3531_gene_tad.txt
get_tad $GENEDIR/Erondoui.bed $TADDIR/Erondoui/Erondoui_10kb_TADs_domains.bed > $OUTDIR/Erondoui_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $41!="" && $31!="" && $41!~/,/ && $31!~/,/ {print $31"\t"$41}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N23.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/X3531_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/Erondoui_gene_tad.txt) \
  > $OUTDIR/Ecas_Eron_gene_tad.tsv
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (X3531): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (Erondoui): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (Erondoui): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (X3531): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (X3531): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (Erondoui): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Ecas_Eron_gene_tad.tsv > $OUTDIR/Ecas_Eron_TAD_summary.txt

##### Species pair Enivalis and Ecalcaria
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/C0079.bed $TADDIR/C0080/C0080_10kb_TADs_domains.bed > $OUTDIR/C0080_gene_tad.txt
get_tad $GENEDIR/X3336.bed $TADDIR/X3258/X3258_10kb_TADs_domains.bed > $OUTDIR/X3336_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $5!="" && $38!="" && $5!~/,/ && $38!~/,/ {print $38"\t"$5}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N24.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/C0080_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/X3336_gene_tad.txt) \
  > $OUTDIR/Ecal_Eniv_gene_tad.tsv
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (C0080): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (X3336): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (X3336): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (C0080): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (C0080): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (X3336): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Ecal_Eniv_gene_tad.tsv > $OUTDIR/Ecal_Eniv_TAD_summary.txt

##### Species pair Egorge and Epluto
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/C0096.bed $TADDIR/C0001/C0001_10kb_TADs_domains.bed > $OUTDIR/C0001_gene_tad.txt
get_tad $GENEDIR/Egorge.bed $TADDIR/Egorge/Egorge_10kb_TADs_domains.bed > $OUTDIR/Egorge_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $8!="" && $19!="" && $8!~/,/ && $19!~/,/ {print $19"\t"$8}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N26.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/C0001_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/Egorge_gene_tad.txt) \
  > $OUTDIR/Eplu_Egor_gene_tad.tsv
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (C0001): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (Egorge): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (Egorge): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (C0001): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (C0001): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (Egorge): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Eplu_Egor_gene_tad.tsv > $OUTDIR/Eplu_Egor_TAD_summary.txt

##### Species pair Eepiphron and Epharte
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/Epharte.bed $TADDIR/Epharte/Epharte_10kb_TADs_domains.bed > $OUTDIR/Epharte_gene_tad.txt
get_tad $GENEDIR/Eepiphron.bed $TADDIR/Eepiphron/Eepiphron_10kb_TADs_domains.bed > $OUTDIR/Eepiphron_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $29!="" && $16!="" && $29!~/,/ && $16!~/,/ {print $16"\t"$29}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N35.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/Epharte_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/Eepiphron_gene_tad.txt) \
  > $OUTDIR/Epha_Eepi_gene_tad.tsv
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (Epharte): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (Eepiphron): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (Eepiphron): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (Epharte): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (Epharte): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (Eepiphron): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Epha_Eepi_gene_tad.tsv > $OUTDIR/Epha_Eepi_TAD_summary.txt


##### Species pair Emedusa and Etriaria
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/Emedusa.bed $TADDIR/Emedusa/Emedusa_10kb_TADs_domains.bed > $OUTDIR/Emedusa_gene_tad.txt
get_tad $GENEDIR/Etriaria.bed $TADDIR/Etriaria/Etriaria_10kb_TADs_domains.bed > $OUTDIR/Etriaria_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $34!="" && $22!="" && $34!~/,/ && $22!~/,/ {print $22"\t"$34}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N31.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/Emedusa_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/Etriaria_gene_tad.txt) \
  > $OUTDIR/Emed_Etri_gene_tad.tsv
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (Emedusa): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (Etriaria): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (Etriaria): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (Emedusa): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (Emedusa): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (Etriaria): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Emed_Etri_gene_tad.tsv > $OUTDIR/Emed_Etri_TAD_summary.txt


##### Species pair Estirius and Emontana
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/Estirius.bed $TADDIR/Estirius/Estirius_10kb_TADs_domains.bed > $OUTDIR/Estirius_gene_tad.txt
get_tad $GENEDIR/Emontana.bed $TADDIR/Emontana/Emontana_10kb_TADs_domains.bed > $OUTDIR/Emontana_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $32!="" && $27!="" && $32!~/,/ && $27!~/,/ {print $27"\t"$32}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N37.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/Estirius_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/Emontana_gene_tad.txt) \
  > $OUTDIR/Esti_Emon_gene_tad.tsv
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (Estirius): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (Emontana): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (Emontana): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (Estirius): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (Estirius): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (Emontana): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Esti_Emon_gene_tad.tsv > $OUTDIR/Esti_Emon_TAD_summary.txt

##### Species pair Etyndarus and Enivalis
# Step 1: Assign genes to TADs - take TAD with most overlap (same logic as compartments)
get_tad $GENEDIR/X3738.bed $TADDIR/X3737/X3737_10kb_TADs_domains.bed > $OUTDIR/X3738_gene_tad.txt
get_tad $GENEDIR/X3336.bed $TADDIR/X3258/X3258_10kb_TADs_domains.bed > $OUTDIR/X3336_gene_tad.txt
# Step2: Add TAD info for both species (same join logic as compartments)
join -t $'\t' -1 2 -2 1 \
  <(awk -F'\t' 'NR>1 && $42!="" && $38!="" && $42!~/,/ && $38!~/,/ {print $38"\t"$42}' \
      $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N29.tsv | sort -t$'\t' -k2,2) \
  <(sort -k1,1 $OUTDIR/X3738_gene_tad.txt) \
  | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
  | sort -k1,1 \
  | join -t $'\t' -1 1 -2 1 - \
      <(sort -k1,1 $OUTDIR/X3336_gene_tad.txt) \
  > $OUTDIR/Etyn_Eniv_gene_tad.tsv
# cols: X3336_gene | X3738_gene | X3738_tad | X3336_tad
# Step 3: count colocalizing pairs
awk -F'\t' '
{
  combo[$3"\t"$4]++
  tad_sp1[$3]++
  tad_sp2[$4]++
  sp1_to_sp2[$3][$4]++
  sp2_to_sp1[$4][$3]++
}
END {
  pairs_sp1_also_sp2 = 0
  pairs_sp2_also_sp1 = 0

  for (c in combo) {
    k = combo[c]
    pairs_sp1_also_sp2 += k*(k-1)/2
    pairs_sp2_also_sp1 += k*(k-1)/2
  }

  pairs_sp1 = 0
  for (tad in tad_sp1) {
    k = tad_sp1[tad]
    pairs_sp1 += k*(k-1)/2
  }

  pairs_sp2 = 0
  for (tad in tad_sp2) {
    k = tad_sp2[tad]
    pairs_sp2 += k*(k-1)/2
  }

  tads_sp1_total = 0
  tads_sp1_conserved = 0
  for (tad1 in sp1_to_sp2) {
    if (tad_sp1[tad1] < 2) continue   # skip TADs with only 1 gene
    tads_sp1_total++
    n_sp2_tads = 0
    for (tad2 in sp1_to_sp2[tad1]) n_sp2_tads++
    if (n_sp2_tads == 1) tads_sp1_conserved++
  }

  tads_sp2_total = 0
  tads_sp2_conserved = 0
  for (tad2 in sp2_to_sp1) {
    if (tad_sp2[tad2] < 2) continue   # skip TADs with only 1 gene
    tads_sp2_total++
    n_sp1_tads = 0
    for (tad1 in sp2_to_sp1[tad2]) n_sp1_tads++
    if (n_sp1_tads == 1) tads_sp2_conserved++
  }

  printf "Gene pairs co-localising in same TAD in sp1 (X3738): %d\n", pairs_sp1
  printf "Of those, also co-localising in sp2 (X3336): %d (%.1f%%)\n\n", pairs_sp1_also_sp2, (pairs_sp1 > 0 ? pairs_sp1_also_sp2/pairs_sp1*100 : 0)
  printf "Gene pairs co-localising in same TAD in sp2 (X3336): %d\n", pairs_sp2
  printf "Of those, also co-localising in sp1 (X3738): %d (%.1f%%)\n\n", pairs_sp2_also_sp1, (pairs_sp2 > 0 ? pairs_sp2_also_sp1/pairs_sp2*100 : 0)
  printf "TADs with 2+ orthologues in sp1 (X3738): %d\n", tads_sp1_total
  printf "Of those, orthologues co-localise in same sp2 TAD: %d (%.1f%%)\n\n", tads_sp1_conserved, (tads_sp1_total > 0 ? tads_sp1_conserved/tads_sp1_total*100 : 0)
  printf "TADs with 2+ orthologues in sp2 (X3336): %d\n", tads_sp2_total
  printf "Of those, orthologues co-localise in same sp1 TAD: %d (%.1f%%)\n", tads_sp2_conserved, (tads_sp2_total > 0 ? tads_sp2_conserved/tads_sp2_total*100 : 0)
}' $OUTDIR/Etyn_Eniv_gene_tad.tsv > $OUTDIR/Etyn_Eniv_TAD_summary.txt
