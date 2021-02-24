## **NDID**
## Introduction
**Ch**romatin **I**nteraction **A**nalysis with **P**aired-**E**nd **T**ag (**ChIA-PET**) sequencing is a technology to study genome-wide long-range chromatin interactions bound by protein factors. **NDID** is a statistical technique for the joint normalization and differential chromatin interactions detection from ChIA-PET experiments. 
### The NDID requires the following dependencies:
> 1) &nbsp; R (≥3.4.0) <br />
> 2) &nbsp; fANCOVA (≥ 0.5.1)<br />
> 3) &nbsp; qvalue ( ≥ 2.15.0) <br />
### Data Pre-processing
The data pre-processing has two steps. In the first step, we need to process the two ChIA-PET raw datasets using ChIA-PET Tool V3 (Li et al., 2019). In the second step, new anchors will be defined from the processed results; that is, the two processed data's anchors will be merged and considered the unique anchors. Using the newly defined anchors, we need to re-processed the raw data using an option *--INPUT_ANCHOR_FILE* in ChIA-PET Tool. 

**Example**: for processing the GM12878 versus MCF7 datasets, we will use the following ChIA-PET Tool V3 command lines:
**1)** &nbsp; **For the first-step analysis:**  we will use the following command lines for GM12878 and MCF7 datasets, respectively. <br />

		java -jar ChIA-PET.jar --mode 1 --fastq1 GM12878_1.fastq --fastq2 GM12878_2.fastq --linker ChIA-PET_Tool_V3/linker/linker_long.txt --minimum_linker_alignment_score 14 --GENOME_INDEX hg19.fa --GENOME_LENGTH 3E9 --CHROM_SIZE_INFO ChIA-PET_Tool_V3/chromInfo/hg19.chromSize.txt --CYTOBAND_DATA ChIA-PET_Tool_V3/chromInfo/hg19_cytoBandIdeo.txt --SPECIES 1 --output Output_GM12878 --prefix GM12878 

		java -jar ChIA-PET.jar --mode 1 --fastq1 MCF7_1.fastq --fastq2 MCF7_2.fastq --linker ChIA-PET_Tool_V3/linker/linker_long.txt --minimum_linker_alignment_score 14 --GENOME_INDEX hg19.fa --GENOME_LENGTH 3E9 --CHROM_SIZE_INFO ChIA-PET_Tool_V3/chromInfo/hg19.chromSize.txt --CYTOBAND_DATA ChIA-PET_Tool_V3/chromInfo/hg19_cytoBandIdeo.txt --SPECIES 1 --output Output_MCF7 --prefix MCF7 

**2)** &nbsp; **For the second-step analysis:** we will use the defined anchor from the first step, namely Anchor.bed.<br />

		java -jar ChIA-PET.jar --mode 1 --fastq1 GM12878_1.fastq --fastq2 GM12878_2.fastq --linker ChIA-PET_Tool_V3/linker/linker_long.txt --minimum_linker_alignment_score 14 --GENOME_INDEX hg19.fa --GENOME_LENGTH 3E9 --CHROM_SIZE_INFO ChIA-PET_Tool_V3/chromInfo/hg19.chromSize.txt --CYTOBAND_DATA ChIA-PET_Tool_V3/chromInfo/hg19_cytoBandIdeo.txt --SPECIES 1  --INPUT_ANCHOR_FILE Anchor.bed --output Output_GM12878 --prefix GM12878_MCF7 

		java -jar ChIA-PET.jar --mode 1 --fastq1 MCF7_1.fastq --fastq2 MCF7_2.fastq --linker ChIA-PET_Tool_V3/linker/linker_long.txt --minimum_linker_alignment_score 14 --GENOME_INDEX hg19.fa --GENOME_LENGTH 3E9 --CHROM_SIZE_INFO ChIA-PET_Tool_V3/chromInfo/hg19.chromSize.txt --CYTOBAND_DATA ChIA-PET_Tool_V3/chromInfo/hg19_cytoBandIdeo.txt --SPECIES 1 INPUT_ANCHOR_FILE Anchor.bed --output Output_MCF7 --prefix MCF7_GM12878 


**Remark**: for detailed information on ChIA-PET Tool V3 data analysis, please visit the [ChIA-PET Tool V3](https://github.com/GuoliangLi-HZAU/ChIA-PET_Tool_V3).

From the ChIA-PET Tool output files, the *out.cluster.FDRfiltered.txt* will be used for downstream analysis in NDID. Finally, we will find the overlap results between the two processed results using bedtools pairToPair. We considered the anchors' location, interaction frequency, and self-ligation PETs in each anchor (it used to measure the anchor enrichment) from the overlapped results. When we overlapped the two processed dataset results, we have a chance to find unique interactions only in one dataset. Therefore, we substituted a small interaction frequency (IF=1) for the corresponding dataset. The anchor enrichment is computed from the out.spet file in the ChIA-PET Tool V3 output. It computed for anchor1 (chrom1, start1, end1) and anchor2 (chrom2, start2, end2)  separately and took the average values using the following commands.  
> 1) &nbsp; awk '{if($2<$5){print $1"\t"$2"\t"$5}else{print $1"\t"$5"\t"$2}}' out.spet  > out.spet.bed3 <br />
> 2) &nbsp; bedtools coverage -a Anchor1.bed -b out.spet.bed3|cut -f4  > self1.bed <br />
> 3) &nbsp; bedtools coverage -a Anchor2.bed -b out.spet.bed3|cut -f4  > self2.bed <br />
> 3) &nbsp; then compute the average of it and call the variable name “selfAvg” <br />
### Input files
The input data should have such kind of variables arrangement. The  suffix 1 and 2 indicate the values belogs to sampel-1 and sample-2 respectively. 
|chrom1 |start1 |end1  |chrom2 |start2 |end2  |ipet_1 |selfAvg_1|	ipet_2|	selfAvg_2|
|-------|-------|------|------|--------|------|-------|---------|-------|----------|
|chr1|	27882430|	27885761|	chr1|	27925635|	27937736|	46|	704|	20|	135.5|
### Test data sets
> 1) &nbsp;	RAD21 ChIA-PET data from human MCF7:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127022
> 2) &nbsp;RAD21 ChIA-PET data from human K562:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127027
> 3) &nbsp;RAD21 ChIA-PET data from human H1:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127037
> 4) &nbsp;RAD21 ChIA-PET data from human H9:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127034
> 5) &nbsp;RAD21 ChIA-PET data from human GM12878:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127053
### Usage
	Rscript NDID.R -h
  	Usage: NDID.R [-[-input|i] <character>] [-[-prefix|p] [<character>]] [-[-help|h]]
         -i| --input      input file
         -p| --prefix     output prefix (default "out")
         -h| --help       print help
### Example 
Let us run the NDID on a given dataset, GM12878_MCF7.txt 
-	*Rscript NDID.R  -i GM12878_MCF7.txt
### Result file
We will get the result file named *out_significant_interaction.txt*. <br />

|chrom1 |start1 |end1  |chrom2|	start2|	end2|Normalized_ipet_1|Normalized_ipet_2|	P-value|p.adjust|intensity type|
|-------|-------|------|------|-------|-----|----------------|----------------|--------|--------|--------------|
|chr1	|39648742|39652714|chr1	|39654893|39662163 |3.83360974|	15.20362117|0.000798074	|0.022997704|1|

### Meaning of the columns:
- **chrom1:**  The name of the chromosome on which the cluster anchor 1 exists <br />
- **start1:**  The start coordinates of cluster anchor 1 <br />
- **end1:**  The end coordinate of cluster anchor 1
- **Chrom2:**  The name of the chromosome on which the cluster anchor 2 exists <br />
- **Start2:**  The start coordinates of cluster anchor 2 <br />
- **End2:**  The end coordinate of cluster anchor 2 <br />
- **Normalized_ipet_1:**  normalized number of PETs in sample-1 <br />
- **Normalized_ipet_2:**  normalized number of PETs in sample-2 <br />
- **P-value:** This value represents the statistical significance of the interaction <br />
- **p.adjust:**  p-value adjustment with Benjamini-Hockberg method (1995) <br />
- **intensity type:**  0 and 1 represent decreasing and increasing interactions intensity, respectively <br />


