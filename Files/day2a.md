
Here you will learn how to perform a scan for selection by calculating PBS (population branch statistic) in windows from low-depth data.

As reference, these are the labelling for each population:

- LWK: Africans
- TSI: Europeans
- CHB: East Asians
- NAM: Native Americans

-----------------------------

### Allele frequency differentiation

The joint-SFS can be considered as a summary statistics for the level of genetic differentiation between populations.
We have seen how it can be used as prior information when computing FST (and related metrics) without relying on genotype calling.

Here we see how to compute the 2D-SFS, FST, PBS, and other summary statistics from low-depth data using ANGSD.
Our final goal is to detect signatures of selection in our data, with the specific example of EDAR gene in Native Americans.

To compute FST/PBS we first need to estimate the marginal and joint sites frequency spectra for our populations.

-------------------------------

**IMPORTANT NOTE**: Please skim this part if already covered during other practicals BUT be sure to run the commands as you will need these output files to perform the remaining analyses.

One of the most important aspect of data analysis for population genetics is the estimate of the Site Frequency Spectrum (SFS).
SFS records the proportions of sites at different allele frequencies. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state.
SFS is informative on the demography of the population or on selective events (when estimated at a local scale).

We use ANGSD to estimate SFS using on example dataset, using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
Details on the implementation can be found [here](http://popgen.dk/angsd/index.php/SFS_Estimation).
Briefly, from sequencing data one computes genotype likelihoods (as previously described).
From these quantities ANGSD computes posterior probabilities of Sample Allele Frequency (SAF), for each site.
Finally, an estimate of the SFS is computed.

These steps can be accomplished in ANGSD using `-doSaf 1/2` options and the program `realSFS`.

```
$ANGSD/angsd -doSaf
...
-doSaf		0
	1: perform multisample GL estimation
	2: use an inbreeding version
	3: calculate genotype probabilities (use -doPost 3 instead)
	4: Assume genotype posteriors as input (still beta) 
	-doThetas		0 (calculate thetas)
	-underFlowProtect	0
	-fold			0 (deprecated)
	-anc			(null) (ancestral fasta)
	-noTrans		0 (remove transitions)
	-pest			(null) (prior SFS)
	-isHap			0 (is haploid beta!)
	-doPost			0 (doPost 3,used for accesing saf based variables)
NB:
	  If -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allelefrequency for each site
```

The SFS is typically computed for each population separately.
Also, we want to estimate the unfolded SFS and we use a putative ancestral sequence to polarise our alleles (to ancestral and derived states).

We cycle across all populations and compute SAF files:
```
for POP in LWK TSI CHB NAM
do
        echo $POP
        $ANGSD/angsd -b $DATA/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 15 -setMaxDepth 150 -doCounts 1 \
                -GL 1 -doSaf 1 &> /dev/null
done
```

Have a look at the output file.
```
$ANGSD/misc/realSFS print Results/NAM.saf.idx | less -S
```
These values represent the sample allele frequency likelihoods at each site, as seen during the lecture.
So the first value (after the chromosome and position columns) is the likelihood of having 0 copies of the derived allele, the second indicates the probability of having 1 copy and so on.
Note that these values are in log format and scaled so that the maximum is 0.

Can you spot any site which is likely to be variable?
What does this mean? It means that you should look for sites where the highest likelihood does not correspond to allele frequencies of 0 or 1.

The next step would be to use these likelihoods and estimate the overall SFS.
This is achieved by the program `realSFS`.
```
$ANGSD/misc/realSFS
-> ---./realSFS------
	-> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:

	-> Estimate the SFS for entire genome??
	-> ./realSFS afile.saf.idx 

	-> 1) Estimate the SFS for entire chromosome 22 ??
	-> ./realSFS afile.saf.idx -r chr22 

	-> 2) Estimate the 2d-SFS for entire chromosome 22 ??
	-> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22 

	-> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??
	-> ./realSFS afile.saf.idx -nSites 500000000 

	-> 4) Estimate the SFS around a gene ??
	-> ./realSFS afile.saf.idx -r chr2:135000000-140000000 

	-> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications]

	-> See realSFS print for possible print options
	-> Use realSFS print_header for printing the header
	-> Use realSFS cat for concatenating saf files

	->------------------
	-> NB: Output is now counts of sites instead of log probs!!
	-> NB: You can print data with ./realSFS print afile.saf.idx !!
	-> NB: Higher order SFSs can be estimated by simply supplying multiple .saf.idx files!!
	-> NB: Program uses accelerated EM, to use standard EM supply -m 0 
```

Therefore, this command will estimate the SFS for each population separately:
```
for POP in LWK TSI CHB NAM
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx 2> /dev/null > Results/$POP.sfs
done
```
The output will be saved in `Results/POP.sfs` files.

You can now have a look at the output file, for instance for the African (LWK) samples:
```
cat Results/LWK.sfs
```
The first value represent the expected number of sites with derived allele frequency equal to 0, the second column the expected number of sites with frequency equal to 1 and so on.

How many values do you expect?
```
awk -F' ' '{print NF; exit}' Results/LWK.sfs
```
Indeed this represents the unfolded spectrum, so it has 2N+1 values with N diploid individuals.

Why is it so bumpy?

The maximum likelihood estimation of the SFS should be performed at the whole-genome level to have enough information.
However, for practical reasons, here we could not use large genomic regions.
Also, as we will see later, this region is not really a proxy for neutral evolution so the SFS is not expected to behave neutrally for some populations.
Nevertheless, these SFS should be a reasonable prior to be used for estimation of summary statistics.

Let us plot the SFS for each pop using this simple R script.
```
Rscript $DIR/Scripts/plotSFS.R Results/LWK.sfs Results/TSI.sfs Results/CHB.sfs Results/NAM.sfs
evince Results/LWK_TSI_CHB_NAM.pdf
```

Do they behave like expected?
Which population has more SNPs?
Which population has a higher proportion of common (not rare) variants?

---------------------------------------

**VERY OPTIONAL** (which means you should ignore this)

It is sometimes convenient to generate bootstrapped replicates of the SFS, by sampling with replacements genomic segments.
This could be used for instance to get confidence intervals when using the SFS for demographic inferences.
This can be achieved in ANGSD using:
```
$ANGSD/misc/realSFS Results/NAM.saf.idx -bootstrap 10  2> /dev/null > Results/NAM.boots.sfs
cat Results/NAM.boots.sfs
```
This command may take some time.
The output file has one line for each boostrapped replicate.

More examples on how to estimate the SFS with ANGSD can be found [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/sfs.md).

---------------------------------------

Secondly, we need to estimate a **multi-dimensional SFS**, for instance the joint SFS between 2 populations (2D).
This can be used for making inferences on their divergence process (time, migration rate and so on).
However, here we are interested in estimating the 2D-SFS as prior information for our FST/PBS.

An important issue when doing this is to be sure that we are comparing the exactly same corresponding sites between populations.
ANGSD does that automatically and considers only a set of overlapping sites.

We are performing PBS assuming NAM (Native Americans) being the targeted population.
The 2D-SFS between all populations and NAM are computed with:
```
POP2=NAM
for POP in LWK TSI CHB
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx Results/$POP2.saf.idx 2> /dev/null > Results/$POP.$POP2.sfs
done
# we also need the comparison between LWK and TSI as we will see later 
$ANGSD/misc/realSFS Results/LWK.saf.idx Results/TSI.saf.idx 2> /dev/null > Results/LWK.TSI.sfs
```

The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.
```
less -S Results/LWK.NAM.sfs
```
You can plot it, but you need to define how many samples you have per population.
```
Rscript $DIR/Scripts/plot2DSFS.R Results/LWK.NAM.sfs 20 20
evince Results/LWK.NAM.sfs.pdf
```

You can even estimate SFS with higher order of magnitude.
This command may take some time and you should skip it if not interested.
```
# $ANGSD/misc/realSFS Results/LWK.saf.idx Results/TSI.saf.idx Results/NAM.saf.idx 2> /dev/null > Results/LWK.TSI.NAM.sfs
```

------------------------------------

Here we are going to calculate **allele frequency differentiation** using the PBS (population branch statistic) metric.
Again, we can achieve this by avoid genotype calling using ANGSD.
From the sample allele frequencies likelihoods (.saf files) we can estimate PBS using the following pipeline.

Note that here we use the previously calculated SFS as prior information.
Also, PEL is our target population, while CHB and TSI are reference populations.
If not already done, you should calculate .saf.idx files for each population, as explained in the section above.

Therefore, we need to use the 2D-SFS between TSI and CHB and NAM (already done).

If we also assume CHB being the target population, as a possible separate analysis is (you don't have to run this if not interested):
```
POP2=CHB
for POP in LWK TSI
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx Results/$POP2.saf.idx 2> /dev/null > Results/$POP.$POP2.sfs
done
```

The 2D-SFS will be used as prior information for the joint allele frequency probabilities at each site.
From these probabilities we will calculate the population branch statistic (PBS) using the NAM (and/or CHB) as target population and LWK and TSI as reference populations.
Our goal is to detect selection in NAM (and/or CHB) in terms of allele frequency differentiation.

Specifically, we are computing a slinding windows scan, with windows of 50kbp and a step of 10kbp.
This can be achieved using the following commands.

1) This command will compute per-site FST indexes (please note the order of files):
```
# NAM
$ANGSD/misc/realSFS fst index Results/LWK.saf.idx Results/TSI.saf.idx Results/NAM.saf.idx -sfs Results/LWK.TSI.sfs -sfs Results/LWK.NAM.sfs -sfs Results/TSI.NAM.sfs -fstout Results/NAM.pbs &> /dev/null
# CHB
$ANGSD/misc/realSFS fst index Results/LWK.saf.idx Results/TSI.saf.idx Results/CHB.saf.idx -sfs Results/LWK.TSI.sfs -sfs Results/LWK.CHB.sfs -sfs Results/TSI.CHB.sfs -fstout Results/CHB.pbs &> /dev/null
```
and you can have a look at their values:
```
# NAM
$ANGSD/misc/realSFS fst print Results/NAM.pbs.fst.idx | less -S
# CHB
$ANGSD/misc/realSFS fst print Results/CHB.pbs.fst.idx | less -S
```
where columns are: chromosome, position, (a), (a+b) values for the three FST comparisons, where FST is defined as a/(a+b).
Note that FST on multiple SNPs is calculated as sum(a)/sum(a+b).

2) The next command will perform a sliding-window analysis:
```
# NAM
$ANGSD/misc/realSFS fst stats2 Results/NAM.pbs.fst.idx -win 50000 -step 10000 > Results/NAM.pbs.txt 2> /dev/null
# CHB
$ANGSD/misc/realSFS fst stats2 Results/CHB.pbs.fst.idx -win 50000 -step 10000 > Results/CHB.pbs.txt 2> /dev/null
```

Have a look at the output file:
```
# NAM
less -S Results/NAM.pbs.txt
# CHB
less -S Results/CHB.pbs.txt
```
The header is:
```
region  chr     midPos  Nsites  Fst01   Fst02   Fst12   PBS0    PBS1    PBS2
```
Where are interested in the column `PB2` which gives the PBS values assuming our population (coded here as 2) being the target population.
Note that negative PBS and FST values are equivalent to 0.

We are also provided with the individual FST values.
You can see that high values of PBS2 are indeed associated with high values of both Fst02 and Fst12 but not Fst01.
We can plot the results along with the gene annotation.
```
# NAM
Rscript $DIR/Scripts/plotPBS.R Results/NAM.pbs.txt Results/NAM.pbs.pdf
# CHB
Rscript $DIR/Scripts/plotPBS.R Results/CHB.pbs.txt Results/CHB.pbs.pdf
```

It will also print out the maximum PBS value observed as this value will be used in the next part.
This script will also plot the PBS variation in LWK as a control comparison.
```
# NAM
evince Results/NAM.pbs.pdf
# CHB
evince Results/CHB.pbs.pdf
```
Comment on the results.

-------------------------

**OPTIONAL**

We are also interested in assessing whether an increase in allele frequency differentiation is also associated with a change of **nucleotide diversity** in CHB.
Again, we can achieve this using ANGSD by estimating levels of diversity without relying on called genotypes.

The procedure is similar to what done for PBS, and the SFS is again used as a prior to compute allele frequencies probabilities.
From these quantities, expectations of various diversity indexes are compute.
This can be achieved using the following pipeline.

First we compute the allele frequency posterior probabilities and associated statistics (-doThetas) using the SFS as prior information (-pest)
```
POP=CHB
echo $POP
$ANGSD/angsd -b $DATA/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
	-GL 1 -doSaf 1 \
	-doThetas 1 -pest Results/$POP.sfs &> /dev/null
```
Then we need to index thess file and perform a sliding windows analysis using a window length of 50kbp and a step size of 10kbp.
```
POP=CHB
echo $POP
# index files
$ANGSD/misc/thetaStat make_bed Results/$POP.thetas.gz &> /dev/null
# perform a sliding-window analysis
$ANGSD/misc/thetaStat do_stat Results/$POP.thetas.gz -nChr 1 -win 50000 -step 10000 -outnames Results/$POP.thetas &> /dev/null
```

Look at the results:
```
less -S Results/CHB.thetas.pestPG
```
and plot the sliding windows scan for nucleotide diversity:
```
Rscript $DIR/Scripts/plotSS.R
evince Results/CHB.ss.pdf
```

------------------------

[HOME](https://github.com/mfumagalli/Copenhagen)



