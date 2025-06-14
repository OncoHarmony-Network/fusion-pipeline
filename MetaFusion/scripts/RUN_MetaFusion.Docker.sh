#!/bin/bash

#Change date to current date. Can also add tag to this string for multiple runs
date=Jun-08-2024

#DATASETS
sim45_sim52=1
brca_4=1
beers_neg=1
sim_50=1
sim101=1
melanoma=1
ntrk_control=1

fusiontools=/MetaFusion/scripts
#REFERENCE FILES FILES
#runs_dir=/MetaFusion/RUNS
runs_dir=/MetaFusion/Test_Data_RUNS

gene_bed=/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
gene_info=/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
recurrent_bedpe=/MetaFusion/reference_files/blocklist_breakpoints.bedpe


# SIM45.SIM52.combined
if [ $sim45_sim52 -eq 1 ]; then
echo SIM45.SIM52
outdir=$runs_dir/SIM45.SIM52.benchmark.$date
#outdir=$runs_dir/SIM45.SIM52.NAIVE_MERGE_GENES.$date
#outdir=$runs_dir/SIM45.SIM52.NAIVE_MERGE_BREAKPOINTS.$date
echo generating output in $outdir
mkdir $outdir
cff=/MetaFusion/test_data/cff/dream.sim45.sim52.cff
truth_fusions=/MetaFusion/test_data/truth_sets/dream.sim45.sim52.truth_set.dat

#bash MetaFusion.naive_merge.sh --outdir $outdir \
bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools 
fi

#BT474.KPL4.MCF7.SKBR3
if [ $brca_4 -eq 1 ]; then
echo BT474.KPL4.MCF7.SKBR3
outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.$date
#outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.NAIVE_MERGE_GENES.$date
#outdir=$runs_dir/BT474.KPL4.MCF7.SKBR3.NAIVE_MERGE_BREAKPOINTS.$date
echo generating output in $outdir
truth_fusions=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat
cff=/MetaFusion/test_data/cff/BRCA.cff

#bash MetaFusion.naive_merge.sh --outdir $outdir \
bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --gene_info $gene_info \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

#NEGATIVE CONTROL BEERS
if [ $beers_neg -eq 1 ]; then
echo NEGATIVE CONTROL BEERS
outdir=$runs_dir/BEERS.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/beers_neg.cff 
truth_fusions=/MetaFusion/test_data/truth_sets/BRCA.truth_set.dat

#                 --genome_fasta $genome_fasta \
bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --gene_info $gene_info \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

# SIM50 2500 fusions files:
if [ $sim_50 -eq 1 ]; then
echo SIM50
outdir=$runs_dir/SIM50.$date
#outdir=$runs_dir/SIM50.NAIVE_MERGE_GENES.$date
#outdir=$runs_dir/SIM50.NAIVE_MERGE_BREAKPOINTS.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/sim50.cff
truth_fusions=/MetaFusion/test_data/truth_sets/sim50.truth_set.dat

#                 --genome_fasta $genome_fasta \
#bash MetaFusion.naive_merge.sh --outdir $outdir \
bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools

fi

#SIM101 2500 fusions files, same truth set as SIM50
if [ $sim101 -eq 1 ]; then
echo SIM101
outdir=$runs_dir/SIM101.$date
#outdir=$runs_dir/SIM101.NAIVE_MERGE_GENES.$date
#outdir=$runs_dir/SIM101.NAIVE_MERGE_BREAKPOINTS.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/sim101.cff
truth_fusions=/MetaFusion/test_data/truth_sets/sim101.truth_set.dat

#                 --genome_fasta $genome_fasta \
#bash MetaFusion.naive_merge.sh --outdir $outdir \
bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --gene_info $gene_info \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools

fi

# STX16--RAE1 from BT474.KPL4.MCF7.SKBR3 dataset
if [ $stx16_rae1 -eq 1 ]; then
echo STX16--RAE1
outdir=$runs_dir/STX16--RAE1.$date
echo generating output in $outdir
cff=/MetaFusion/test_data/cff/STX16--RAE1.figure_subset.cff


bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools


fi

# Melanoma and CML
if [ $melanoma -eq 1 ]; then
echo MELANOMA and CML 
outdir=$runs_dir/melanoma.CML.$date
#outdir=$runs_dir/melanoma.CML.NAIVE_MERGE_GENES.$date
#outdir=$runs_dir/melanoma.CML.NAIVE_MERGE_BREAKPOINTS.$date
#outdir=$runs_dir/melanoma.CML.no_SRR018269.$date
echo generating output in $outdir
truth_fusions=/MetaFusion/test_data/truth_sets/melanoma.truth_set.dat
#truth_fusions=/MetaFusion/test_data/truth_sets/melanoma.truth_set.no_SRR018269.dat
cff=/MetaFusion/test_data/cff/melanoma.cff
#cff=/MetaFusion/test_data/cff/melanoma.no_SRR018269.cff

#bash MetaFusion.naive_merge.sh --outdir $outdir \
bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --truth_set $truth_fusions \
                 --num_tools=2  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools
fi

# NTRK_control
if [ $ntrk_control -eq 1 ]; then
echo NTRK_control
outdir=$runs_dir/NTRK_control.$date
#outdir=$runs_dir/NTRK.NAIVE_MERGE_BREAKPOINTS.$date
#outdir=$runs_dir/NTRK_control.arriba.star_fusion.star_seqr.$date
echo generating output in $outdir
truth_fusions=/MetaFusion/test_data/truth_sets/NTRK_control.truth_set.dat
cff=/MetaFusion/test_data/cff/NTRK_control.cff
#cff=/MetaFusion/test_data/cff/NTRK_control.arriba.star_fusion.star_seqr.cff

bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --gene_info $gene_info \
                 --genome_fasta $genome_fasta \
                 --truth_set $truth_fusions \
                 --num_tools=2 \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools

fi

