#!/bin/bash

#SBATCH -J sc
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q normal
#SBATCH -c 40
#SBATCH --mem=400G

module load cellranger/6.1.1
module load R/4.2.1
module load gcc/11.2.0
module load fftw/3.3.10

### 1.preprocess
### for only single-cell RNA-seq:
rawdata_path="./Data/CleanData/"
sample="G1 G2 G3 G4"
for s in $sample
do
cellranger count --id=$s --fastqs=$rawdata_path/$s --sample=$s --transcriptome=/storage/gaoxiaofeiLab/yaoxingyun/cpu/yaoxingyun/singlecell/refdata-gex-mm10-2020-A/ --localcores=40
done


### 2.quality control for single sample
R_script="Rscript/"
sample="G1 G2 G3 G4"
for s in $sample
do
subpath="${s}/outs/filtered_feature_bc_matrix/"
mkdir $s
Rscript ${R_script}/1.Seurat3_QualityControl.r -i 1.rawdata/${subpath} -w 2.qc/$s/ -o ${s} -s Mus > 2.qc/${s}/1.QualityControl.log
done


### 3.doublets filter and combined samples
sample_names="G1,G2,G4,G3"
rds_files="G1_raw.rds,G2_raw.rds,G4_raw.rds,G3_raw.rds"
filter_fea="7000,7000,7500,7500"
filter_mt="10,15,15,15"
/usr/local/bin/Rscript ${R_script}/2.IntegratedSamples.Seurat3_PCAselection.dedoublet.r -w ./ -l $sample_names \
    -f $rds_files \
    -u $filter_fea -m $filter_mt -o "3.combination.DeDoublets"  > qc_Doublet.output.txt 2>&1

/usr/local/bin/Rscript ${R_script}/3.CCA_tSNEorUMAP.r -w ./3.combination.DeDoublets/ -p 30 \
    -f Origin_Integrated.rds -o combined > ./3.combination.DeDoublets/3.CCA_p30.log