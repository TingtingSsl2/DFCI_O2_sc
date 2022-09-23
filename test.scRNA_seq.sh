#!/bin/sh

# title: script for single-cell data analysis
# author: Tingting Zhao
# email: tingting_zhao@dfci.harvard.edu
# date: "09/06/2022"
# usage: source scRNA_seq.sh
# usage: sbatch -p short -t 1:00:00 --mem=64G -e %j.err -o %j.out test.scRNA_seq.sh
# platform: HMS O2
# folder structure: project folder (data, output, result, scripts)


# bash setting
#index_file="/data/bioinformatics/projects/hassan2022/data/FC_07189/220118_A01061_0232_BH72NCDMXY/BPF_library_Feany_lab_mod.csv"
reference_folder="/n/data2/dfci/medonc/lindsley/reference/10X/mouse/refdata-gex-mm10-2020-A/genes/"
project_folder="/n/data2/dfci/medonc/lindsley/Tingting/"
data_folder="/n/data2/dfci/medonc/lindsley/Tingting/data/"
#baseSpace_ID =
#bcl_files="/data/bioinformatics/projects/hassan2022/data/FC_07189/220118_A01061_0232_BH72NCDMXY"
#expect_cells=10000
#fastq_folder="01_cellranger_mkfastq"
countTable_folder="BCOR_U2af1_Bcor_10x"


# R setting
pwd="/n/data2/dfci/medonc/lindsley/Tingting/scripts/"
indir="/n/data2/dfci/medonc/lindsley/Tingting/data/BCOR_U2af1_Bcor_10x/"
outdir="/n/data2/dfci/medonc/lindsley/Tingting/output/test/"
scrubletdir="/n/data2/dfci/medonc/lindsley/Tingting/output/scrublet/"
samples="Pool107_23,Pool107_24"
projectName="rahul2022"
marker_link="https://docs.google.com/spreadsheets/d/1gCzAeVe9Ekpyt8XdNOyBrOcw7Letr-WTpBh4-37ySpA/edit#gid=446579886"
marker_sheet="MarkerGenesFiltered"
flag=1 #Options for cell clustering algorithm, 1=louvain, 2=GLMPCA, 3= leiden, 0=louvain and GLMPCA and leiden
mtPattern="^mt"
rbPattern="^Rp[sl]"
qc_cutoff=3
mito_cutoff=10
sex="male,male"
genotypes="U2Bcor,U2Bcor"
refdir="/n/data2/dfci/medonc/lindsley/reference/10X/mouse/refdata-gex-mm10-2020-A/genes/"
scriptdir="/n/data2/dfci/medonc/lindsley/Tingting/scripts/"
geneN=10

# load modules
#conda init /PHShome/tz949/anaconda3/envs/scrnaseq
source ~/.bashrc
source activate sc
module load cellranger/6.0.0 


# step1: download data from BaseSpace
#bs download run -i $baseSpace_ID -o $data_folder


# step2: bcl to fastq
#cd $project_folder/output
#bsub -q big cellranger mkfastq --id =$fastq_folder \
#                   --run =$bcl_files \
#                   --csv =$index_file


# step3: making count table
#cd $project_folder/output/$countTable_folder
#for i in $(cat $index_file | awk -F "," '(NR>1){print $2}'); do
#bsub -q big cellranger count --id = $i \
#--transcriptome = $reference_folder \
#--fastqs =$project_folder/output/$fastq_folder/outs/fastq_path/ \
#--sample =$i \
#--expect-cells =$expect_cells \
#--include-introns
#done


## step4: prepare geneID and geneSymbol table from reference genome
#gtf_file="genes.gtf"
#geneName_file="Mouse_mm10_geneAnnotationTable.txt"
#cd $refdir
## extract Mouse gene ID and gene symbols from gtf file
#cat $gtf_file | awk 'BEGIN{FS="\t"}{split($9,a,";"); if ($3~"gene") print a[1]"\t"a[4]}' | sed 's/gene_id "//' | sed 's/ gene_name "//' | sed 's/"//g' > $geneName_file
## run R script to manipulate gene ID and gene symbols and generate an output file
#cd $scriptdir
#Rscript makeGeneAnnoTable.R $refdir $geneName_file


## step5: run scrublet
#matrixdir="/outs/filtered_feature_bc_matrix"
#~/anaconda3/bin/python scrublet_multi.py $indir $matrixdir $scrubletdir "${samples[@]}" # if run scripts on server


## step6: run Seurat individual, to cell clustering step
##bsub -q big -e seurat_individual.log Rscript seurat_individual.R $pwd $indir $outdir $scrubletdir $samples $projectName $marker_link $marker_sheet $flag $mtPattern $rbPattern $mitoCutoff
Rscript test.seurat_individual.R $pwd $indir $outdir $scrubletdir $samples $projectName $marker_link $marker_sheet $flag $mtPattern $rbPattern $qc_cutoff $mito_cutoff $sex $genotypes $refdir $scriptdir $geneN