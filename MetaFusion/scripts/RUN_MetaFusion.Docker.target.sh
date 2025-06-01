#!/bin/bash

#Change date to current date. Can also add tag to this string for multiple runs

#date=May-29-2024
# 修改以切换不同的 run
#date=Jul-2024-WSX-t3
#date=Jul-2024-WSX-all
#date=Jul-2024-WSX-metafusion2
date=Jul-2024-WSX-test_final

#DATASETS
target=1

fusiontools=/MetaFusion/scripts
#REFERENCE FILES FILES
runs_dir=/MetaFusion/NAIVE_MERGE_RUNS

gene_bed=/MetaFusion/reference_files/hg38_gene.bed
gene_info=/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/MetaFusion/reference_files/ref_genome.fa
recurrent_bedpe=/MetaFusion/reference_files/blacklist_hg38.bedpe


outdir=$runs_dir/target.$date
echo generating output in $outdir
mkdir -p $outdir
#cff=/MetaFusion/RUNS/RUNS/merged_1e5.cff
cff=/MetaFusion/RUNS/RUNS/merged_test500.cff
#cff=/MetaFusion/RUNS/RUNS/merged.cff
#cff=/MetaFusion/RUNS/RUNS/star_fusion.cff

#bash MetaFusion.naive_merge.sh --outdir $outdir \
bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --num_tools=1  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools

echo "Done for target"
