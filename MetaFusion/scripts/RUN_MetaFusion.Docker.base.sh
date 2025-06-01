#!/bin/bash

echo Start at `date` -----------

# IO settings ---------------
# $1, file path to CFF from RUNS/
# $2, a run name given to this task
cff=/MetaFusion/$1
outdir=/MetaFusion/RUNS_Results/$2

# Reference settings ---------------
fusiontools=/MetaFusion/scripts
runs_dir=/MetaFusion/NAIVE_MERGE_RUNS
gene_bed=/MetaFusion/reference_files/hg38_gene.bed
gene_info=/MetaFusion/reference_files/Homo_sapiens.gene_info
genome_fasta=/MetaFusion/reference_files/ref_genome.fa
recurrent_bedpe=/MetaFusion/reference_files/blacklist_hg38.bedpe

# Checking ---------------
if [ ! -f $cff ]; then
    echo "Warning: the input file does not exist, please check!"
    exit 0
fi

if [ -e $outdir/final*.cluster ]; then
    echo "Warning: the final file existed, please delete it first if you wanna rerun!"
    exit 0
fi

# Running ---------------
echo "Generating output in $outdir"
mkdir -p $outdir

bash MetaFusion.sh --outdir $outdir \
                 --cff $cff  \
                 --gene_bed $gene_bed \
                 --fusion_annotator \
                 --genome_fasta $genome_fasta \
                 --gene_info $gene_info \
                 --num_tools=1  \
                 --recurrent_bedpe $recurrent_bedpe \
                 --scripts $fusiontools

chmod -R 777 $outdir

echo "Done"
echo End at `date` -----------