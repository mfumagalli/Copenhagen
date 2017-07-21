
In this session you will learn how to do:
* basic post-mapping filtering
* genotype calling
* allele frequency estimation
* variant calling

For most of the examples, we will use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen at the University of Copenhagen.
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

According to its website *ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. This is especially useful for low and medium depth data. The software is written in C++ and has been used on large sample sizes. This program is not for manipulating BAM/CRAM files, but solely a tool to perform various kinds of analysis. We recommend the excellent program SAMtools for outputting and modifying bamfiles.*

Please make sure to follow the preparatory instructions on the main page before running these examples.

--------------------------------------------------------------------------

#### Data filtering

First, we will learn how to build a command line in ANGSD, with the specific example of doing a basic post mapping filtering.

To see a full list of options in ANGSD type:
```
$ANGSD/angsd
```
and you should see something like
```
...
Overview of methods:
	-GL		Estimate genotype likelihoods
	-doCounts	Calculate various counts statistics
	-doAsso		Perform association study
	-doMaf		Estimate allele frequencies
	-doError	Estimate the type specific error rates
	-doAncError	Estimate the errorrate based on perfect fastas
	-HWE_pval		Est inbreedning per site or use as filter
	-doGeno		Call genotypes
	-doFasta	Generate a fasta for a BAM file
	-doAbbababa	Perform an ABBA-BABA test
	-sites		Analyse specific sites (can force major/minor)
	-doSaf		Estimate the SFS and/or neutrality tests genotype calling
	-doHetPlas	Estimate hetplasmy by calculating a pooled haploid frequency

	Below are options that can be usefull
	-bam		Options relating to bam reading
	-doMajorMinor	Infer the major/minor using different approaches
	-ref/-anc	Read reference or ancestral genome
	-doSNPstat	Calculate various SNPstat
	-cigstat	Printout CIGAR stat across readlength
	many others

For information of specific options type: 
	./angsd METHODNAME eg 
		./angsd -GL
		./angsd -doMaf
		./angsd -doAsso etc
		./angsd sites for information about indexing -sites files
Examples:
	Estimate MAF for bam files in 'list'
		'./angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1'
```

We will see later of to perform SNP and genotype calling (and many other things) with ANGSD.

ANGSD can accept several input files, as described [here](http://popgen.dk/angsd/index.php/Input):

* BAM/CRAM
* Pileup
* Genotype likelihood/probability files
* VCF

Here we show how ANGSD can also perform some basic filtering of the data.
These filters are based on:

* quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
* SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
* sites, see [here](http://popgen.dk/angsd/index.php/Sites)

Have a look at our list of BAM files:
```
cat $DATA/ALL.bamlist
wc -l $DATA/ALL.bamlist
```

If the input file is in BAM format, the possible options are:
```
$ANGSD/angsd -bam
...
parseArgs_bambi.cpp: bam reader:
	-bam/-b		(null)	(list of BAM/CRAM files)
	-i		(null)	(Single BAM/CRAM file)
	-r		(null)	Supply a single region in commandline (see examples below)
	-rf		(null)	Supply multiple regions in a file (see examples below)
	-remove_bads	1	Discard 'bad' reads, (flag >=256) 
	-uniqueOnly	0	Discards reads that doesnt map uniquely
	-show		0	Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
	-minMapQ	0	Discard reads with mapping quality below
	-minQ		13	Discard bases with base quality below
	-trim		0	Number of based to discard at both ends of the reads
	-trim		0	Number of based to discard at 5 ends of the reads
	-trim		0	Number of based to discard at 3 ends of the reads
	-only_proper_pairs 1	Only use reads where the mate could be mapped
	-C		0	adjust mapQ for excessive mismatches (as SAMtools), supply -ref
	-baq		0	adjust qscores around indels (as SAMtools), supply -ref
	-checkBamHeaders 1	Exit if difference in BAM headers
	-doCheck	1	Keep going even if datafile is not suffixed with .bam/.cram
	-downSample	0.000000	Downsample to the fraction of original data
	-nReads		50	Number of reads to pop from each BAM/CRAMs
	-minChunkSize	250	Minimum size of chunk sent to analyses

Examples for region specification:
		chr:		Use entire chromosome: chr
		chr:start-	Use region from start to end of chr
		chr:-stop	Use region from beginning of chromosome: chr to stop
		chr:start-stop	Use region from start to stop from chromosome: chr
		chr:site	Use single site on chromosome: chr
```

Some basic filtering consists in removing, for instance, read with low quality and/or with multiple hits, and this can be achieved using the parameters ```-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1```.
Let us build such command line.
First we need to define input and output files:
```
# $ANGSD/angsd -b ALL.bamlist -ref $REF -out Results/ALL \
...
```
with `-b` we give the file including paths to all BAM files we need to analyse.
`-ref` specifies the reference sequence.
`-out` states the prefix for all output files that will be generated.

Next we need to define some basic filtering options.
First we define filters based on reads quality.
```
# $ANGSD/angsd -b ALL.bamlist -ref $REF -out Results/ALL \
#        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
...
```
These filters will retain only uniquely mapping reads, not tagged as bad, considering only proper pairs, without trimming, and adjusting for indel/mapping (as in samtools).
`-C 50` reduces the effect of reads with excessive mismatches, while `-baq 1` computes base alignment quality as explained here ([BAQ](http://samtools.sourceforge.net/mpileup.shtml)) used to rule out false SNPs close to INDELS.

Also, you may want to remove reads with low mapping quality and sites with low quality or covered by few reads (low depth).
Under these circumnstances, the assignment of individual genotypes and SNPs is problematic, and can lead to errors. 

However, it is necessary to know the overall distribution of per-site depth, in order to avoid filtering too many (or few) sites.

To make things faster, we are using only 20 samples from European samples (TSI, Italians, of course...).

We first derive the distribution of quality scores and depth on our data set using ```-doQsDist 1 -doDepth 1```.
```
$ANGSD/angsd -b $DATA/TSI.bamlist -ref $REF -out Results/ALL.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 200 -minQ 0 &> /dev/null
```
As an illustration here, ```-maxDepth 200``` corresponds to a per-sample average depth of 10.
This option means that all sites with depth equal or greater than 800 will be binned together.

Have a look at the files generated:
```
ls Results/*
...
-> Output filenames:
		->"Results/ALL.qc.arg"
		->"Results/ALL.qc.qs"
		->"Results/ALL.qc.depthSample"
		->"Results/ALL.qc.depthGlobal"
...
```
```
# counts of quality scores
less -S Results/ALL.qc.qs
# counts of per-sample depth
less -S Results/ALL.qc.depthSample 
wc -l Results/ALL.qc.depthSample # 80
# counts of global depth
less -S Results/ALL.qc.depthGlobal 
```

It is convenient to compute the percentiles of these distributions (and visualize them) in order to make an informative decision on the threshold values we will use for our filtering.
```
Rscript $DIR/Scripts/plotQC.R Results/ALL.qc 2> /dev/null
```
Have a look at the output files:
```
less -S Results/ALL.qc.info
evince Results/ALL.qc.pdf
``` 

Which values would you choose as sensible thresholds on quality score and global depth (min and max)?
We may also want to remove sites where half of the individual have no data. This is achieved by the -minInd option.
A possible command line would contain the following filtering:
```
...
#        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
#        -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
...
```
which corresponds to the following scenario:

Parameter | Meaning |
--- | --- |
-minInd 10 | use only sites with data from at least N individuals |
-setMinDepth 10 | minimum total depth |
-setMaxDepth 100 | maximum total depth |

Specifically here we analyse only reads with a minimum mapping quality of 20, and bases with a minimum quality of 20 (the values are in phred scores).
Also we specify that we are analysing only sites where we have data for half of the individuals and minimum and maximum TOTAL depth of 10 and 100, respectively.
`-doCounts 1` simply forces the calculation of depth.

ANGSD can also compute more sophisticated metrics to filter out SNPs, as described [here](http://popgen.dk/angsd/index.php/SnpFilters), mostly based on:
* strand bias
* deviation from HWE
* quality score bias
The strand bias models are described [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/). Some examples of strand biases, taken from the previously cited study, can be found [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/strand_bias.png).
Different models are implemented, as seen [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/strand_bias_eq.png).
For instance, if we want to remove sites which deviate from HWE genotype frequencies (e.g. due to mapping errors), we can set `-HWE_pval 0.05`.

We have now seen how to build a command line in ANGSD with the example of doing a basic post-mapping filtering.

----------------------------

#### Estimation of allele frequencies and SNP calling

We now want to estimate allele frequencies at each site, for the whole sample.
In other words, at each site we want to to estimate (or count) how many copies of different alleles (two in case of biallelic variants) we observe in our sample (across all sequenced individuals).
However with low depth data direct counting of individually assigned genotypes can lead to biased allele frequencies.

ANGSD has an option to estimate **allele frequencies** taking into account data uncertainty from genotype likelihoods:
```
$ANGSD/angsd -doMaf
...
-doMaf  0 (Calculate persite frequencies '.mafs.gz')
        1: Frequency (fixed major and minor)
        2: Frequency (fixed major unknown minor)
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
        NB. Filedumping is supressed if value is negative
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
Filters:
        -minMaf         -1.000000       (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
Extras:
        -ref    (null)  (Filename for fasta reference)
        -anc    (null)  (Filename for fasta ancestral)
        -eps    0.001000 [Only used for -doMaf &8]
        -beagleProb     0 (Dump beagle style postprobs)
        -indFname       (null) (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor
```

Therefore, the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic).
```
$ANGSD/angsd -doMajorMinor
...
        -doMajorMinor   0
        1: Infer major and minor from GL
        2: Infer major and minor from allele counts
        3: use major and minor from a file (requires -sites file.txt)
        4: Use reference allele as major (requires -ref)
        5: Use ancestral allele as major (requires -anc)
        -rmTrans: remove transitions 0
        -skipTriallelic 0
```

Finally, you need to specify which genotype likelihood model to use.
```
$ANGSD/angsd -GL
...
        -GL=0:
        1: SAMtools
        2: GATK
        3: SOAPsnp
        4: SYK
        5: phys
        -trim           0               (zero means no trimming)
        -tmpdir         angsd_tmpdir/   (used by SOAPsnp)
        -errors         (null)          (used by SYK)
        -minInd         0               (0 indicates no filtering)

Filedumping:
        -doGlf  0
        1: binary glf (10 log likes)    .glf.gz
        2: beagle likelihood file       .beagle.gz
        3: binary 3 times likelihood    .glf.gz
        4: text version (10 log likes)  .glf.gz
```
A description of these different implementation can be found [here](http://www.popgen.dk/angsd/index.php/Genotype_likelihoods).
The GATK model refers to the first GATK paper, SAMtools is somehow more sophisticated (non-independence of errors), SOAPsnp requires a reference sequence for recalibration of quality scores, SYK is error-type specific.
For most applications and data, GATK and SAMtools models should give similar results.

Therefore a possible command line to estimate allele frequencies might be (this may take 1 min to run):
```
$ANGSD/angsd -b $DATA/TSI.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
        -GL 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1 &> /dev/null
```
where we specify:
* -GL 1: genotype likelihood model as in SAMtools
* -doMajorMinor 4: force the major allele to be the reference (the minor is inferred)
* -doMaf 1: major and minor are fixed

What are the output files?
```
->"Results/ALL.arg"
->"Results/ALL.mafs.gz"
```
`.args` file is a summary of all options used, while `.mafs.gz` file shows the allele frequencies computed at each site.

Have a look at this file which contains estimates of allele frequency values.
```
zcat Results/ALL.mafs.gz | head
```
and you may see something like
```
chromo	position	major	minor	ref	knownEM	nInd
2	108999932	G	A	G	0.000003	10
2	108999933	T	A	T	0.000004	10
2	108999934	A	C	A	0.000003	10
2	108999935	T	A	T	0.000004	11
2	108999936	T	A	T	0.000002	11
2	108999937	C	A	C	0.000002	11
2	108999938	T	A	T	0.000005	11
2	108999939	C	A	C	0.000003	11
2	108999940	T	A	T	0.000003	11
```
where `knownEM` specifies the algorithm used to estimate the allele frequency which is given under that column.
Please note that this refers to the allele frequency of the allele labelled as `minor`.
The columns are: chromosome, position, major allele, minor allele, reference allele, allele frequency, p-value for SNP calling (if -SNP-pval was called), number of individuals with data.
The last column gives the number of samples with data (you can see that this never lower than 40 given our filtering).

You can notice that many sites have low allele frequency, probably reflecting the fact that that site is monomorphic.
We may be interested in looking at allele frequencies only for sites that are actually variable in our sample.
Therefore we want to perform a **SNP calling**.
There are two main ways to call SNPs using ANGSD with these options:
```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

As an illustration, let us call SNPs by computing:
 - genotype likelihoods using GATK method;
 - major and minor alleles inferred from genotype likelihoods;
 - frequency from known major allele but unknown minor;
 - SNPs as those having MAF=>0.01.

Try to write down this command by yourself...

Here is the solution (this may take a couple of minutes to run, be patient...):
```
$ANGSD/angsd -b $DATA/TSI.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 \
        -GL 2 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1  \
        -minMaf 0.01 &> /dev/null
```

You can have a look at the results:
```
zcat Results/ALL.mafs.gz | head
...
chromo	position	major	minor	ref	knownEM	nInd
2	108999951	T	G	T	0.048648	13
2	109000064	C	T	C	0.045999	16
2	109000113	C	T	C	0.035418	17
2	109000117	G	T	G	0.036939	17
2	109000122	T	C	T	0.039408	17
2	109000135	T	A	T	0.035312	18
2	109000160	C	A	C	0.037070	18
2	109000181	G	A	G	0.042162	17
2	109000285	A	C	A	0.049312	13
```

How many entries (potential SNPs) you have?
```
zcat Results/ALL.mafs.gz | tail -n+2 | wc -l
```

As a general guidance, `-GL 1`, `-doMaf 1/2` and `-doMajorMinor 1` should be the preferred choice when data uncertainty is high.
If interested in analysing very low frequency SNPs, then `-doMaf 2` should be selected.
When accurate information on reference sequence or outgroup are available, one can use `-doMajorMinor` to 4 or 5.
Also, detecting variable sites based on their probability of being SNPs is generally a better choice than defining a threshold on the allele frequency.
However, various cutoffs and a dedicated filtering should be perform to assess robustness of your called SNPs.

-----------------------------------------------

**EXERCISE**

Try varying the cutoff for SNP calling and record how many sites are predicted to be variable for each scenario.
Identify which sites are not predicted to be variable anymore with a more stringent cutoff (e.g. between a pair of scenario), and plot their allele frequencies.

```
# iterate over some cutoffs (you can change these)
for PV in 0.05 1e-2 1e-4 1e-6
do
        if [ $PV == 0.05 ]; then echo SNP_pval NR_SNPs; fi
        $ANGSD/angsd -b $DATA/TSI.bamlist -ref $REF -out Results/ALL.$PV \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
		-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
		-SNP_pval $PV &> /dev/null
	echo $PV `zcat Results/ALL.$PV.mafs.gz | tail -n+2 | wc -l`
done
```

A possible output is (your numbers may be different):
```
SNP_pval NR_SNPs
0.05 5848
1e-2 3475
1e-4 1825
1e-6 1502
```

Which sites differ from 0.05 and 0.01? What is their frequency?
This script will also print out the first 20 discordant sites (pK.EM is the p-value for the SNP calling test).
```
Rscript -e 'mafs1=read.table(gzfile("Results/ALL.1e-2.mafs.gz"), he=T, strings=F); mafs5=read.table(gzfile("Results/ALL.0.05.mafs.gz"), header=T, stringsAsFact=F); mafs5[!(mafs5[,2] %in% mafs1[,2]),][1:20,]; pdf(file="Results/diff_snpcall.pdf"); par(mfrow=c(1,2)); hist(as.numeric(mafs5[!(mafs5[,2] %in% mafs1[,2]),][,6]), main="Discordant SNPs", xlab="MAF (DAF)"); hist(as.numeric(mafs5[(mafs5[,2] %in% mafs1[,2]),][,6]), main="Concordant SNPs", xlab="MAF"); dev.off();'

evince Results/diff_snpcall.pdf
```

Can you draw some conclusions from these results?
Which frequencies are more difficult to estimate and therefore affect SNP calling?

-----------------------

**EXAMPLE**

We are finally ready to gather some biological insights from the data.
Recalling our research aim, our first goal is to calculate allele frequencies of our target variant in EDAR gene for different human populations.

Our SNP of interest is located in EDAR gene.
According to our reference [paper](http://www.nature.com/ncomms/2016/160519/ncomms11616/full/ncomms11616.html), *The derived G allele at the index SNP in this region (rs3827760) encodes a functional substitution in the intracellular death domain of EDAR (370A) and is associated with reduced chin protrusion*.
The genomic location of this SNP is `chr2:109513601-109513601`.

In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
The file with these positions need to be formatted as (chromosome positions).
```
echo 2 109513601 > snp.txt
```
We need to index this file in order for ANGSD to process it.
```
$ANGSD/angsd sites index snp.txt &> /dev/null
```

We are interested in calculating the derived allele frequencies, so are using the ancestral sequence to polarise the alleles.
We also want to compute the allele frequencies for each population separately.
We need to use a different file for each population, with a different list of BAM files, as provided: 
```
ls $DATA/*.bamlist
# ALL.bamlist  ALL2.bamlist  CHB.bamlist  LWK.bamlist  PEL.bamlist  TSI.bamlist NAM.bamlist
```

Now we can run ANGSD.
Note that we are interested in calculating the derived allele frequency, so we need to specify a putative ancestral sequence.
Finally, since we know analyse single populations, we can filter sites if they deviate from HWE.
Please note that we also add PEL samples, admixed individuals from Peru.
```
for POP in LWK TSI CHB PEL NAM
do
        echo $POP
        $ANGSD/angsd -b $DATA/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 -HWE_pval 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -sites snp.txt &> /dev/null
done
```
An assessment of the deviation from HWE will be print out in files with extension `.hwe.gz`.

You can inspect the results.
```
zcat Results/LWK.mafs.gz Results/TSI.mafs.gz Results/CHB.mafs.gz Results/NAM.mafs.gz Results/PEL.mafs.gz
```

Can you comment on these results?
Do you see any allele frequency differentiation in the derived state?

Indeed, in the original paper, Authors state that *The G allele at rs3827760 is not present in Europeans and Africans but is seen at high frequency in East Asians and is essentially fixed in Native Americans.*

Why is it not fixed in our Peruvians sample? 

------------------------------------------

#### Genotype calling

Here we will explore several ways to call genotypes from sequencing data, once SNPs have been assigned.
We will also calculate genotypes probabilities to each site for each individual.
We will finally compare estimates of allele frequencies from called genotypes to what previously found using a probabilistic framework.

In ANGSD, the option to call genotypes is `-doGeno`:
```
$ANGSD/angsd -doGeno
...
-doGeno	0
	1: write major and minor
	2: write the called genotype encoded as -1,0,1,2, -1=not called
	4: write the called genotype directly: eg AA,AC etc 
	8: write the posterior probability of all possible genotypes
	16: write the posterior probability of called genotype
	32: write the posterior probabilities of the 3 gentypes as binary
	-> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
	-postCutoff=0.333333 (Only genotype to missing if below this threshold)
	-geno_minDepth=-1	(-1 indicates no cutof)
	-geno_maxDepth=-1	(-1 indicates no cutof)
	-geno_minMM=-1.000000	(minimum fraction af major-minor bases)
	-minInd=0	(only keep sites if you call genotypes from this number of individuals)

	NB When writing the posterior the -postCutoff is not used
	NB geno_minDepth requires -doCounts
	NB geno_maxDepth requires -doCounts
```

Therefore, if we set `-doGeno 2`, genotypes are coded as 0,1,2, as the number of alternate alleles.
If we want to print the major and minor alleles as well then we set `-doGeno 3`.

To calculate the posterior probability of genotypes we need to define a model.
```
$ANGSD/angsd -doPost
...
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
...
```
`-doPost 1` uses the estimate per-site allele frequency as a prior for genotype proportions, assuming Hardy Weinberg Equilibrium.
When the assumption of HWE is not valid, you can use an estimate of the inbreeding coefficient, for instance calculated using [ngsF](https://github.com/fgvieira/ngsF).

A typical command for genotype calling assuming HWE is (assuming we analyse our PEL samples):
```
$ANGSD/angsd -b $DATA/PEL.bamlist -ref $REF -out Results/PEL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 15 -setMaxDepth 100 -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-3 \
	-doGeno 3 -doPost 1 -postCutoff 0 &> /dev/null
```

Have a look at the output file:
```
less -S Results/PEL.geno.gz
```
The columns are: chromosome, position, major allele, minor allele, genotypes is 0,1,2 format.

How many sites have at least one missing genotype?
```
zcat Results/PEL.geno.gz | grep -1 - | wc -l
```
Why is that?

You can control how to set missing genotype when their confidence is low with `-postCutoff`.
For instance, we can set as missing genotypes when their (highest) genotype posterior probability is below 0.95:
```
$ANGSD/angsd -b $DATA/PEL.bamlist -ref $REF -out Results/PEL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-3 \
	-doGeno 3 -doPost 1 -postCutoff 0.95 &> /dev/null
```

How many sites do we have in total?
How many sites have at least one missing genotype now?
```
zcat Results/PEL.geno.gz | wc -l
zcat Results/PEL.geno.gz | grep -1 - | wc -l
```

Why are there some many sites with missing genotypes?

The mean depth per sample is around 2-3X, therefore genotypes cannot be assigned with very high confidence.

Setting this threshold depends on the mean sequencing depth of your data, as well as your application.
For some analyses you need to work only with high quality genotypes (e.g. measure of proportion of shared SNPs for gene flow estimate), while for others you can be more relaxed (e.g. estimate of overall nucleotide diversity).
We will show later how to accurately estimate summary statistics with low-depth data.

--------------------------------

**EXAMPLE**

Back to our example of functional variants in EDAR, we want to assign individual genotypes by first computing genotype posterior probabilities for all samples.
Finally we will calculate allele frequencies based on assigned genotypes.

As previously done, let us perform a genotype calling in ANGSD:
```
for POP in LWK TSI CHB NAM PEL
do
        echo $POP
        $ANGSD/angsd -b $DATA/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
		-doGeno 3 -doPost 2 -postCutoff 0.50 \
                -sites snp.txt &> /dev/null
done
```
Note that, as an illustration, here we used a uniform prior (not HWE) with `-doPost 2` with a threshold for setting missing genotypes to 0.50.

Open the output files:
```
for POP in LWK TSI CHB NAM PEL
do
	echo $POP
	zcat Results/$POP.geno.gz
done
```

What is the derived allele frequency for each population?

We have a lot of missing data.
Try to calculate allele frequencies in PEL by using a HWE-prior and comment on the results (e.g. which are the genotypic states more difficult to assign?).

```
# with HWE-prior
for POP in LWK TSI CHB NAM PEL
do
	echo $POP
	$ANGSD/angsd -b $DATA/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -doGeno 3 -doPost 1 -postCutoff 0.50 \
                -sites snp.txt &> /dev/null
done
```

Open the files:
```
for POP in LWK TSI CHB NAM PEL
do
        echo $POP
        zcat Results/$POP.geno.gz
done
```

For instance, you may have 0/40 in LWK and TSI, 40/40 in CHB, and 14/20=0.70 in PEL.
Recall that we previously estimated a minor allele frequency of 0.55 in PEL without assigning individuals.

Why do we observe such a difference in the estimates?

----------------------------

**OPTIONAL**

#### Population structure

**IMPORTANT NOTE**: These commands are given as a mere example as, in practise, such analyses should be performed on larger genomic regions.

We now want to investigate population structure of our sample. 
We perform a principal component analyses (PCA) without relying on called genotypes, but rather by taking their uncertainty into account.
More specifically, the next program we are going to use (ngsTools) takes as input genotype probabilities in binary format, so we need to specify `-doGeno 32` in ANGSD.
Also, we are using a HWE-based prior with `-doPost 1`.

```
$ANGSD/angsd -b $DATA/ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 70 -setMaxDepth 200 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-2 \
        -doGeno 32 -doPost 1 &> /dev/null
```
Unzip the results but you cannot open it since it is in binary format:
```
gunzip Results/ALL.geno.gz
```
We are going to use `ngsCovar`, which estimates the covariance matrix between individuals based on genotype probabilities.
Then this matrix will be decomposed into principal componenets which will be investigated for population structure analyses.
If you type:
```
$NGSTOOLS/ngsPopGen/ngsCovar
```
you will see a list of possible options.

For instance, we need to define how many sites we have.
To retrieve such value, we can inspect the file with allele frequencies:
```
less -S Results/ALL.mafs.gz
```
How many sites do we have?
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```
Now we can perform a PCA by estimating the covariance matrix using:
```
$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/ALL.geno -outfile Results/ALL.covar -nind 80 -nsites $N_SITES -call 0 -norm 0 &> /dev/null
```
with the options `-call 0` meaning that we do not perform genotype calling and `-norm 0` that we are not normalising by allele frequency.
The latter may give more weight to low frequency variants which are harder to estimate.

Look at the output file:
```
less -S Results/ALL.covar
```
which represents a matrix of NxN with N individuals giving the covariance.
Note that this matrix is symmetric.

Finally, we perform an eigenvector decomposition and plot the resulting map:
```
# create a cluster-like file defining the labelling (population) for each sample
Rscript -e 'write.table(cbind(rep(seq(1,20),4),rep(seq(1,20),4),c(rep("CHB",20),rep("LWK",20),rep("NAM",20),rep("TSI",20))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="Results/ALL.clst", quote=F)'
# run and plot
Rscript $DIR/Scripts/plotPCA_ngstools.R Results/ALL.covar Results/ALL.clst 1-2 Results/ALL.pca.pdf
```
where the parameter `1-2` specifies that we are plotting only the first and second component.
On the screen, you will see a series of numbers.
These are the percentage of explained variance for each component.

Finally, you can open the produced image:
```
evince Results/ALL.pca.pdf
```

You can inspect other components:
```
Rscript $DIR/Scripts/plotPCA_ngstools.R Results/ALL.covar Results/ALL.clst 2-3 Results/ALL.pca2.pdf
evince Results/ALL.pca2.pdf
```

Therefore, NAM samples appear very close to EUR but separated.
Indeed, among all Latin American populations present in the 1000G project, Peruvians are the least admixed population.
We can either use all these samples or compute admixture proportions in order to select a subset of putative Native American (unadmixed) samples.
However with such limited data set we cannot refine the structure between NAM and CHB.

-----------------------------------

**OPTIONAL**

#### Admixture

**IMPORTANT NOTE**: These commands are given as a mere example as, in practise, such analyses should be performed on larger genomic regions.

We want to select only a subset of PEL samples with high Native American ancestry, or check the level of admixture in our NAM samples.
We use ngsAdmix, which again works on genotype probabilities and not on individual calls.
This is suitable for low-depth data.

ngsAdmix requires genotype likelihoods in BEAGLE format as input.
We can compute these quantities with ANGSD with `-doGlf 2`.
```
$ANGSD/angsd -b $DATA/ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 50 -setMaxDepth 200 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-4 \
        -doGlf 2 &> /dev/null
```

We assume 4 ancestral populations making up the genetic diversity of our samples.
Therefore we compute admixture proportions with 4 ancestral components.
```
K=4
$NGSADMIX -likes Results/ALL.beagle.gz -K $K -outfiles Results/ALL.admix.K$K -P 4 -minMaf 0.02 &> /dev/null
```

We now combine samples IDs with admixture proportions and inspect the results.
```
paste $DATA/ALL.bamlist Results/ALL.admix.K$K.qopt > Results/ALL.admix.K$K.txt
less -S Results/ALL.admix.K$K.txt
```

From these quantities we can extract how many samples (and which ones) have a high proportion of Native American ancestry (e.g. >0.90). 
------------------

**OPTIONAL**

#### Genetic distances

**IMPORTANT NOTE**: These commands are given as a mere example as, in practise, such analyses should be performed on larger genomic regions.

Genotype probabilities can be used also to infer the structure of your population.
For instance, in our example, they can used to assess whether PEL samples are indeed admixed.

We can compute genetic distances as a basis for population clustering driectly from genotype probabilities, and not from assigned genotypes as we have seen how problematic these latters can be at low-depth.

First, we compute genotype posterior probabilities jointly for all samples:
```
# Assuming HWE, without filtering based on probabilities, with SNP calling
$ANGSD/angsd -b $DATA/ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 50 -setMaxDepth 200 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null
```

Next we record how many sites we retrieve.
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

Then we create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("LWK","TSI","CHB","NAM"),each=20), rep(1:20, 4), sep="_"), sep="\n", file="Data/pops.label")'
cat Data/pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can compute pairwise genetic distances without relying on individual genotype calls.
```
$NGSDIST/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 80 -n_sites $N_SITES -labels Data/pops.label -o Results/ALL.dist -n_threads 4 &> /dev/null
less -S Results/ALL.dist
```

We can visualise the pairwise genetic distances in form of a tree.
```
$FASTME -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
```

Finally, we plot the tree.
```
Rscript -e 'library(ape); library(phangorn); pdf(file="Results/ALL.tree.pdf"); plot(read.tree("Results/ALL.tree"), cex=0.5); dev.off();' &> /dev/null
evince Results/ALL.tree.pdf
```

One can also perform a PCA/MDS from such genetic distances to further explore the population structure.

[HOME](https://github.com/mfumagalli/Copenhagen)

