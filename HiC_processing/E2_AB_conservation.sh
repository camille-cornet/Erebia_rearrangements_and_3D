#!/bin/bash

### Use Orthofinder and check if the same genes are in A or B in different species

OUTDIR=$MAINDIR/HiC/E2_AB_conservation
ORTHODIR=$MAINDIR/Downstream/D1_orthofinder/results/Results_Erebia
GENEDIR=$MAINDIR/Downstream/D3_genespace/bed
COMPDIR=$MAINDIR/HiC/A4_callAB_cool
cd $OUTDIR

# Compartment lookup function: gene bed + bedgraph -> gene\tcomp
get_comp() {
  bedtools intersect -a $1 -b $2 -wa -wb \
  | awk 'OFS="\t" {
      overlap = (($3 < $7 ? $3 : $7) - ($2 > $5 ? $2 : $5))
      if (overlap > best[$4]) { best[$4] = overlap; comp[$4] = ($8 > 0) ? "A" : "B" }
    }
    END { for (g in comp) print g"\t"comp[g] }'
}

##### Species pair Eligea and Eottomana (N5 node and 20 and 4 col in Orthofinder)
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "Eligea\tC0056\tC0056_comp\tEligea_comp" > $OUTDIR/Eott_Elig_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $4!="" && $20!="" && $4!~/,/ && $20!~/,/ {print $20"\t"$4}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N5.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/C0056.bed $COMPDIR/C0055/C0055_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/Eligea.bed $COMPDIR/Eligea/Eligea_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Eott_Elig_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  Elig=$4; Eott=$3
  if (Elig=="A") Elig_A++
  if (Elig=="B") Elig_B++
  if (Eott=="A") Eott_A++
  if (Eott=="B") Eott_B++
  if (Elig=="A" && Eott=="A") both_A++
  if (Elig=="B" && Eott=="B") both_B++
  if (Elig=="A" && Eott=="B") Elig_A_Eott_B++
  if (Elig=="B" && Eott=="A") Elig_B_Eott_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Eligea A: %d (%.1f%% of total)\n", Elig_A, Elig_A/total*100
  printf "Eligea B: %d (%.1f%% of total)\n", Elig_B, Elig_B/total*100
  printf "Eottomana A: %d (%.1f%% of total)\n", Eott_A, Eott_A/total*100
  printf "Eottomana B: %d (%.1f%% of total)\n\n", Eott_B, Eott_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Eligea A / Eottomana B: %d (%.1f%% of total)\n", Elig_A_Eott_B, Elig_A_Eott_B/total*100
  printf "Eligea B / Eottomana A: %d (%.1f%% of total)\n\n", Elig_B_Eott_A, Elig_B_Eott_A/total*100
  printf "Of Eottomana A: %.1f%% also A in Eligea, %.1f%% B in Eligea\n", both_A/Eott_A*100, Elig_B_Eott_A/Eott_A*100
  printf "Of Eottomana B: %.1f%% also B in Eligea, %.1f%% A in Eligea\n", both_B/Eott_B*100, Elig_A_Eott_B/Eott_B*100
  printf "Of Eligea A: %.1f%% also A in Eottomana, %.1f%% B in Eottomana\n", both_A/Elig_A*100, Elig_A_Eott_B/Elig_A*100
  printf "Of Eligea B: %.1f%% also B in Eottomana, %.1f%% A in Eottomana\n", both_B/Elig_B*100, Elig_B_Eott_A/Elig_B*100
}' $OUTDIR/Eott_Elig_single_copy_AB.tsv > $OUTDIR/Eott_Elig_AB_summary.txt

##### Species pair Egraucasica and Eottomana 
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "C0111\tC0056\tC0056_comp\tC0111_comp" > $OUTDIR/Eott_Egra_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $4!="" && $9!="" && $4!~/,/ && $9!~/,/ {print $9"\t"$4}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N9.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/C0056.bed $COMPDIR/C0055/C0055_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/C0111.bed $COMPDIR/C0100/C0100_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Eott_Egra_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  Egra=$4; Eott=$3
  if (Egra=="A") Egra_A++
  if (Egra=="B") Egra_B++
  if (Eott=="A") Eott_A++
  if (Eott=="B") Eott_B++
  if (Egra=="A" && Eott=="A") both_A++
  if (Egra=="B" && Eott=="B") both_B++
  if (Egra=="A" && Eott=="B") Egra_A_Eott_B++
  if (Egra=="B" && Eott=="A") Egra_B_Eott_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Egraucasica A: %d (%.1f%% of total)\n", Egra_A, Egra_A/total*100
  printf "Egraucasica B: %d (%.1f%% of total)\n", Egra_B, Egra_B/total*100
  printf "Eottomana A: %d (%.1f%% of total)\n", Eott_A, Eott_A/total*100
  printf "Eottomana B: %d (%.1f%% of total)\n\n", Eott_B, Eott_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Egraucasica A / Eottomana B: %d (%.1f%% of total)\n", Egra_A_Eott_B, Egra_A_Eott_B/total*100
  printf "Egraucasica B / Eottomana A: %d (%.1f%% of total)\n\n", Egra_B_Eott_A, Egra_B_Eott_A/total*100
  printf "Of Eottomana A: %.1f%% also A in Egraucasica, %.1f%% B in Egraucasica\n", both_A/Eott_A*100, Egra_B_Eott_A/Eott_A*100
  printf "Of Eottomana B: %.1f%% also B in Egraucasica, %.1f%% A in Egraucasica\n", both_B/Eott_B*100, Egra_A_Eott_B/Eott_B*100
  printf "Of Egraucasica A: %.1f%% also A in Eottomana, %.1f%% B in Eottomana\n", both_A/Egra_A*100, Egra_A_Eott_B/Egra_A*100
  printf "Of Egraucasica B: %.1f%% also B in Eottomana, %.1f%% A in Eottomana\n", both_B/Egra_B*100, Egra_B_Eott_A/Egra_B*100
}' $OUTDIR/Eott_Egra_single_copy_AB.tsv > $OUTDIR/Eott_Egra_AB_summary.txt

##### Species pair Erondoui and Ecassioides (N23 node and 31 and 41 col in Orthofinder)
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "Erondoui\tX3529\tX3529_comp\tErondoui_comp" > $OUTDIR/Ecas_Eron_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $41!="" && $31!="" && $41!~/,/ && $31!~/,/ {print $31"\t"$41}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N23.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/X3529.bed $COMPDIR/X3531/X3531_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/Erondoui.bed $COMPDIR/Erondoui/Erondoui_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Ecas_Eron_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  eron=$4; ecas=$3
  if (eron=="A") eron_A++
  if (eron=="B") eron_B++
  if (ecas=="A") ecas_A++
  if (ecas=="B") ecas_B++
  if (eron=="A" && ecas=="A") both_A++
  if (eron=="B" && ecas=="B") both_B++
  if (eron=="A" && ecas=="B") eron_A_ecas_B++
  if (eron=="B" && ecas=="A") eron_B_ecas_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Erondoui A: %d (%.1f%% of total)\n", eron_A, eron_A/total*100
  printf "Erondoui B: %d (%.1f%% of total)\n", eron_B, eron_B/total*100
  printf "Ecassioides A: %d (%.1f%% of total)\n", ecas_A, ecas_A/total*100
  printf "Ecassioides B: %d (%.1f%% of total)\n\n", ecas_B, ecas_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Erondoui A / Ecassioides B: %d (%.1f%% of total)\n", eron_A_ecas_B, eron_A_ecas_B/total*100
  printf "Erondoui B / Ecassioides A: %d (%.1f%% of total)\n\n", eron_B_ecas_A, eron_B_ecas_A/total*100
  printf "Of Ecassioides A: %.1f%% also A in Erondoui, %.1f%% B in Erondoui\n", both_A/ecas_A*100, eron_B_ecas_A/ecas_A*100
  printf "Of Ecassioides B: %.1f%% also B in Erondoui, %.1f%% A in Erondoui\n", both_B/ecas_B*100, eron_A_ecas_B/ecas_B*100
  printf "Of Erondoui A: %.1f%% also A in Ecassioides, %.1f%% B in Ecassioides\n", both_A/eron_A*100, eron_A_ecas_B/eron_A*100
  printf "Of Erondoui B: %.1f%% also B in Ecassioides, %.1f%% A in Ecassioides\n", both_B/eron_B*100, eron_B_ecas_A/eron_B*100
}' $OUTDIR/Ecas_Eron_single_copy_AB.tsv > $OUTDIR/Ecas_Eron_AB_summary.txt


##### Species pair Enivalis and Ecalcaria
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "X3336\tC0079\tC0079_comp\tX3336_comp" > $OUTDIR/Ecal_Eniv_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $5!="" && $38!="" && $5!~/,/ && $38!~/,/ {print $38"\t"$5}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N24.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/C0079.bed $COMPDIR/C0080/C0080_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/X3336.bed $COMPDIR/X3258/X3258_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Ecal_Eniv_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  Eniv=$3; Ecal=$4
  if (Eniv=="A") Eniv_A++
  if (Eniv=="B") Eniv_B++
  if (Ecal=="A") Ecal_A++
  if (Ecal=="B") Ecal_B++
  if (Eniv=="A" && Ecal=="A") both_A++
  if (Eniv=="B" && Ecal=="B") both_B++
  if (Eniv=="A" && Ecal=="B") Eniv_A_Ecal_B++
  if (Eniv=="B" && Ecal=="A") Eniv_B_Ecal_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Enivalis A: %d (%.1f%% of total)\n", Eniv_A, Eniv_A/total*100
  printf "Enivalis B: %d (%.1f%% of total)\n", Eniv_B, Eniv_B/total*100
  printf "Ecalcaria A: %d (%.1f%% of total)\n", Ecal_A, Ecal_A/total*100
  printf "Ecalcaria B: %d (%.1f%% of total)\n\n", Ecal_B, Ecal_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Enivalis A / Ecalcaria B: %d (%.1f%% of total)\n", Eniv_A_Ecal_B, Eniv_A_Ecal_B/total*100
  printf "Enivalis B / Ecalcaria A: %d (%.1f%% of total)\n\n", Eniv_B_Ecal_A, Eniv_B_Ecal_A/total*100
  printf "Of Ecalcaria A: %.1f%% also A in Enivalis, %.1f%% B in Enivalis\n", both_A/Ecal_A*100, Eniv_B_Ecal_A/Ecal_A*100
  printf "Of Ecalcaria B: %.1f%% also B in Enivalis, %.1f%% A in Enivalis\n", both_B/Ecal_B*100, Eniv_A_Ecal_B/Ecal_B*100
  printf "Of Enivalis A: %.1f%% also A in Ecalcaria, %.1f%% B in Ecalcaria\n", both_A/Eniv_A*100, Eniv_A_Ecal_B/Eniv_A*100
  printf "Of Enivalis B: %.1f%% also B in Ecalcaria, %.1f%% A in Ecalcaria\n", both_B/Eniv_B*100, Eniv_B_Ecal_A/Eniv_B*100
}' $OUTDIR/Ecal_Eniv_single_copy_AB.tsv > $OUTDIR/Ecal_Eniv_AB_summary.txt


##### Species pair Egorge and Epluto 
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "Egorge\tC0096\tC0096_comp\tEgorge_comp" > $OUTDIR/Eplu_Egor_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $8!="" && $19!="" && $8!~/,/ && $19!~/,/ {print $19"\t"$8}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N26.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/C0096.bed $COMPDIR/C0001/C0001_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/Egorge.bed $COMPDIR/Egorge/Egorge_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Eplu_Egor_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  Egor=$4; Eplu=$3
  if (Egor=="A") Egor_A++
  if (Egor=="B") Egor_B++
  if (Eplu=="A") Eplu_A++
  if (Eplu=="B") Eplu_B++
  if (Egor=="A" && Eplu=="A") both_A++
  if (Egor=="B" && Eplu=="B") both_B++
  if (Egor=="A" && Eplu=="B") Egor_A_Eplu_B++
  if (Egor=="B" && Eplu=="A") Egor_B_Eplu_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Egorge A: %d (%.1f%% of total)\n", Egor_A, Egor_A/total*100
  printf "Egorge B: %d (%.1f%% of total)\n", Egor_B, Egor_B/total*100
  printf "Epluto A: %d (%.1f%% of total)\n", Eplu_A, Eplu_A/total*100
  printf "Epluto B: %d (%.1f%% of total)\n\n", Eplu_B, Eplu_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Egorge A / Epluto B: %d (%.1f%% of total)\n", Egor_A_Eplu_B, Egor_A_Eplu_B/total*100
  printf "Egorge B / Epluto A: %d (%.1f%% of total)\n\n", Egor_B_Eplu_A, Egor_B_Eplu_A/total*100
  printf "Of Epluto A: %.1f%% also A in Egorge, %.1f%% B in Egorge\n", both_A/Eplu_A*100, Egor_B_Eplu_A/Eplu_A*100
  printf "Of Epluto B: %.1f%% also B in Egorge, %.1f%% A in Egorge\n", both_B/Eplu_B*100, Egor_A_Eplu_B/Eplu_B*100
  printf "Of Egorge A: %.1f%% also A in Epluto, %.1f%% B in Epluto\n", both_A/Egor_A*100, Egor_A_Eplu_B/Egor_A*100
  printf "Of Egorge B: %.1f%% also B in Epluto, %.1f%% A in Epluto\n", both_B/Egor_B*100, Egor_B_Eplu_A/Egor_B*100
}' $OUTDIR/Eplu_Egor_single_copy_AB.tsv > $OUTDIR/Eplu_Egor_AB_summary.txt

##### Species pair Eepiphron and Epharte
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "Eepiphron\tEpharte\tEpharte_comp\tEepiphron_comp" > $OUTDIR/Epha_Eepi_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $29!="" && $16!="" && $29!~/,/ && $16!~/,/ {print $16"\t"$29}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N35.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/Epharte.bed $COMPDIR/Epharte/Epharte_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/Eepiphron.bed $COMPDIR/Eepiphron/Eepiphron_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Epha_Eepi_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  Eepi=$4; Epha=$3
  if (Eepi=="A") Eepi_A++
  if (Eepi=="B") Eepi_B++
  if (Epha=="A") Epha_A++
  if (Epha=="B") Epha_B++
  if (Eepi=="A" && Epha=="A") both_A++
  if (Eepi=="B" && Epha=="B") both_B++
  if (Eepi=="A" && Epha=="B") Eepi_A_Epha_B++
  if (Eepi=="B" && Epha=="A") Eepi_B_Epha_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Eepiphron A: %d (%.1f%% of total)\n", Eepi_A, Eepi_A/total*100
  printf "Eepiphron B: %d (%.1f%% of total)\n", Eepi_B, Eepi_B/total*100
  printf "Epharte A: %d (%.1f%% of total)\n", Epha_A, Epha_A/total*100
  printf "Epharte B: %d (%.1f%% of total)\n\n", Epha_B, Epha_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Eepiphron A / Epharte B: %d (%.1f%% of total)\n", Eepi_A_Epha_B, Eepi_A_Epha_B/total*100
  printf "Eepiphron B / Epharte A: %d (%.1f%% of total)\n\n", Eepi_B_Epha_A, Eepi_B_Epha_A/total*100
  printf "Of Epharte A: %.1f%% also A in Eepiphron, %.1f%% B in Eepiphron\n", both_A/Epha_A*100, Eepi_B_Epha_A/Epha_A*100
  printf "Of Epharte B: %.1f%% also B in Eepiphron, %.1f%% A in Eepiphron\n", both_B/Epha_B*100, Eepi_A_Epha_B/Epha_B*100
  printf "Of Eepiphron A: %.1f%% also A in Epharte, %.1f%% B in Epharte\n", both_A/Eepi_A*100, Eepi_A_Epha_B/Eepi_A*100
  printf "Of Eepiphron B: %.1f%% also B in Epharte, %.1f%% A in Epharte\n", both_B/Eepi_B*100, Eepi_B_Epha_A/Eepi_B*100
}' $OUTDIR/Epha_Eepi_single_copy_AB.tsv > $OUTDIR/Epha_Eepi_AB_summary.txt


##### Species pair Emedusa and Etriaria
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "Emedusa\tEtriaria\tEtriaria_comp\tEmedusa_comp" > $OUTDIR/Etri_Emed_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $34!="" && $22!="" && $34!~/,/ && $22!~/,/ {print $22"\t"$34}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N31.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/Etriaria.bed $COMPDIR/Etriaria/Etriaria_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/Emedusa.bed $COMPDIR/Emedusa/Emedusa_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Etri_Emed_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  Emed=$4; Etri=$3
  if (Emed=="A") Emed_A++
  if (Emed=="B") Emed_B++
  if (Etri=="A") Etri_A++
  if (Etri=="B") Etri_B++
  if (Emed=="A" && Etri=="A") both_A++
  if (Emed=="B" && Etri=="B") both_B++
  if (Emed=="A" && Etri=="B") Emed_A_Etri_B++
  if (Emed=="B" && Etri=="A") Emed_B_Etri_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Emedusa A: %d (%.1f%% of total)\n", Emed_A, Emed_A/total*100
  printf "Emedusa B: %d (%.1f%% of total)\n", Emed_B, Emed_B/total*100
  printf "Etriaria A: %d (%.1f%% of total)\n", Etri_A, Etri_A/total*100
  printf "Etriaria B: %d (%.1f%% of total)\n\n", Etri_B, Etri_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Emedusa A / Etriaria B: %d (%.1f%% of total)\n", Emed_A_Etri_B, Emed_A_Etri_B/total*100
  printf "Emedusa B / Etriaria A: %d (%.1f%% of total)\n\n", Emed_B_Etri_A, Emed_B_Etri_A/total*100
  printf "Of Etriaria A: %.1f%% also A in Emedusa, %.1f%% B in Emedusa\n", both_A/Etri_A*100, Emed_B_Etri_A/Etri_A*100
  printf "Of Etriaria B: %.1f%% also B in Emedusa, %.1f%% A in Emedusa\n", both_B/Etri_B*100, Emed_A_Etri_B/Etri_B*100
  printf "Of Emedusa A: %.1f%% also A in Etriaria, %.1f%% B in Etriaria\n", both_A/Emed_A*100, Emed_A_Etri_B/Emed_A*100
  printf "Of Emedusa B: %.1f%% also B in Etriaria, %.1f%% A in Etriaria\n", both_B/Emed_B*100, Emed_B_Etri_A/Emed_B*100
}' $OUTDIR/Etri_Emed_single_copy_AB.tsv > $OUTDIR/Etri_Emed_AB_summary.txt


##### Species pair Estirius and Emontana
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "Emontana\tEstirius\tEstirius_comp\tEmontana_comp" > $OUTDIR/Esti_Emon_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $32!="" && $27!="" && $32!~/,/ && $27!~/,/ {print $27"\t"$32}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N37.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/Estirius.bed $COMPDIR/Estirius/Estirius_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/Emontana.bed $COMPDIR/Emontana/Emontana_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Esti_Emon_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  Emon=$4; Esti=$3
  if (Emon=="A") Emon_A++
  if (Emon=="B") Emon_B++
  if (Esti=="A") Esti_A++
  if (Esti=="B") Esti_B++
  if (Emon=="A" && Esti=="A") both_A++
  if (Emon=="B" && Esti=="B") both_B++
  if (Emon=="A" && Esti=="B") Emon_A_Esti_B++
  if (Emon=="B" && Esti=="A") Emon_B_Esti_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Emontana A: %d (%.1f%% of total)\n", Emon_A, Emon_A/total*100
  printf "Emontana B: %d (%.1f%% of total)\n", Emon_B, Emon_B/total*100
  printf "Estirius A: %d (%.1f%% of total)\n", Esti_A, Esti_A/total*100
  printf "Estirius B: %d (%.1f%% of total)\n\n", Esti_B, Esti_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Emontana A / Estirius B: %d (%.1f%% of total)\n", Emon_A_Esti_B, Emon_A_Esti_B/total*100
  printf "Emontana B / Estirius A: %d (%.1f%% of total)\n\n", Emon_B_Esti_A, Emon_B_Esti_A/total*100
  printf "Of Estirius A: %.1f%% also A in Emontana, %.1f%% B in Emontana\n", both_A/Esti_A*100, Emon_B_Esti_A/Esti_A*100
  printf "Of Estirius B: %.1f%% also B in Emontana, %.1f%% A in Emontana\n", both_B/Esti_B*100, Emon_A_Esti_B/Esti_B*100
  printf "Of Emontana A: %.1f%% also A in Estirius, %.1f%% B in Estirius\n", both_A/Emon_A*100, Emon_A_Esti_B/Emon_A*100
  printf "Of Emontana B: %.1f%% also B in Estirius, %.1f%% A in Estirius\n", both_B/Emon_B*100, Emon_B_Esti_A/Emon_B*100
}' $OUTDIR/Esti_Emon_single_copy_AB.tsv > $OUTDIR/Esti_Emon_AB_summary.txt


##### Species pair Etyndarus and Enivalis
# Step 1+2+3: single copy orthologs + compartment assignment + join
echo -e "X3336\tX3738\tX3738_comp\tX3336_comp" > $OUTDIR/Etyn_Eniv_single_copy_AB.tsv
join -t $'\t' -1 1 -2 1 \
  <(join -t $'\t' -1 2 -2 1 \
      <(awk -F'\t' 'NR>1 && $42!="" && $38!="" && $42!~/,/ && $38!~/,/ {print $38"\t"$42}' \
          $ORTHODIR/Phylogenetic_Hierarchical_Orthogroups/N29.tsv | sort -t$'\t' -k2,2) \
      <(get_comp $GENEDIR/X3738.bed $COMPDIR/X3737/X3737_40kb_AB.bedgraph | sort -k1,1) \
    | awk -F'\t' 'OFS="\t" {print $2, $1, $3}' \
    | sort -t$'\t' -k1,1) \
  <(get_comp $GENEDIR/X3336.bed $COMPDIR/X3258/X3258_40kb_AB.bedgraph | sort -k1,1) \
  >> $OUTDIR/Etyn_Eniv_single_copy_AB.tsv
# Step 4: Count % of genes in each compartment in each species
awk -F'\t' '
NR==1 { next }
{
  total++
  Eniv=$3; Etyn=$4
  if (Eniv=="A") Eniv_A++
  if (Eniv=="B") Eniv_B++
  if (Etyn=="A") Etyn_A++
  if (Etyn=="B") Etyn_B++
  if (Eniv=="A" && Etyn=="A") both_A++
  if (Eniv=="B" && Etyn=="B") both_B++
  if (Eniv=="A" && Etyn=="B") Eniv_A_Etyn_B++
  if (Eniv=="B" && Etyn=="A") Eniv_B_Etyn_A++
}
END {
  printf "Total genes: %d\n\n", total
  printf "Enivalis A: %d (%.1f%% of total)\n", Eniv_A, Eniv_A/total*100
  printf "Enivalis B: %d (%.1f%% of total)\n", Eniv_B, Eniv_B/total*100
  printf "Etyndarus A: %d (%.1f%% of total)\n", Etyn_A, Etyn_A/total*100
  printf "Etyndarus B: %d (%.1f%% of total)\n\n", Etyn_B, Etyn_B/total*100
  printf "Both A: %d (%.1f%% of total)\n", both_A, both_A/total*100
  printf "Both B: %d (%.1f%% of total)\n", both_B, both_B/total*100
  printf "Enivalis A / Etyndarus B: %d (%.1f%% of total)\n", Eniv_A_Etyn_B, Eniv_A_Etyn_B/total*100
  printf "Enivalis B / Etyndarus A: %d (%.1f%% of total)\n\n", Eniv_B_Etyn_A, Eniv_B_Etyn_A/total*100
  printf "Of Etyndarus A: %.1f%% also A in Enivalis, %.1f%% B in Enivalis\n", both_A/Etyn_A*100, Eniv_B_Etyn_A/Etyn_A*100
  printf "Of Etyndarus B: %.1f%% also B in Enivalis, %.1f%% A in Enivalis\n", both_B/Etyn_B*100, Eniv_A_Etyn_B/Etyn_B*100
  printf "Of Enivalis A: %.1f%% also A in Etyndarus, %.1f%% B in Etyndarus\n", both_A/Eniv_A*100, Eniv_A_Etyn_B/Eniv_A*100
  printf "Of Enivalis B: %.1f%% also B in Etyndarus, %.1f%% A in Etyndarus\n", both_B/Eniv_B*100, Eniv_B_Etyn_A/Eniv_B*100
}' $OUTDIR/Etyn_Eniv_single_copy_AB.tsv > $OUTDIR/Etyn_Eniv_AB_summary.txt
