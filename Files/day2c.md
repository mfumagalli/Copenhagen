
We have some suggestive results pointing towards selection acting in EDAR based on patterns of genetic differentiation.
However, we employed low-depth data and our analyses were not based on called genotypes.
Therefore, we did not have the power, for instance, to identify the causal variant (if any) and patterns of haplotype distribution.

It is usual to then perform a targeted deep resequencing of our region of interest, to further refine our selection analysis and highlight putative causal variants.
Here we can even include more samples (as the experimental cost will be lower anyway).

Therefore, we are now using high-depth phased data (in VCF format) for 80 samples from LWK (Africans), CHB (Chinese, East Asians), NAM and CLM (Colombians) to assess whether selection signatures are shared across other populations.

Please make sure to follow the preparatory instructions on the main page before running these examples.
```
SS=/ricco/data/matteo/Software/selscan/bin/linux

DIR=/home/matteo/Copenhagen
DATA=/ricco/data/matteo/Data
```

--------------------------

As we are using phased data, we are able to perform selection tests based on haplotype diversity and/or homozygosity, as seen during the lecture.
We are going to use the software [selscan](https://github.com/szpiech/selscan), which implements several selection tests (iHS, nSL, XP-EHH and so on).
According to its manual *[these] statistics are designed to use phased genotypes to identify putative regions of recent or ongoing positive selection in genomes. They are all based on the model of a hard selective sweep, where a de novo adaptive mutation arises on a haplotype that quickly sweeps toward fixation, reducing diversity around the locus. If selection is strong enough, this occurs faster than recombination or mutation can act to break up the haplotype, and thus a signal of high haplotype homozygosity can be observed extending from an adaptive locus.*

A powerful statistic to delect selection from hard sweeps is XP-EHH, which measures the differential decay of haplotype homozygosity. 
We assume that our target population is NAM (or CHB, East Asians) and our reference population is TSI (Europeans).

First, look at the options in selscan:
```
$SS/selscan --help
```
The basic usage of selscan to compute nSL from VCF file is the following. Please note that this command is commented so it won't run but it's given just to
show the usage of selscan.
```
# $SS/selscan --xpehh --vcf Data/NAM.chr2.vcf --ref Data/TSI.chr2.vcf --map Data/genetic_map_chrom2.map --out Results/NAM
```
This statistic can be computed using selscan, provided that we have a genetic map file.
A genetic map for chromosome 2 (based on 1000 Genomes data and PEL samples) is provided here:
```
less -S $DATA/genetic_map_chrom2.map
```
However we need to extract only the sites that correspond in our VCF file.
We also need to interpolate over the sites that are not recorded in our genetic map.
A simple R script to do that is here:
```
Rscript $DIR/Scripts/getGenMap.R $DATA/genetic_map_chrom2.map $DATA/NAM.chr2.vcf > Results/genetic.map
```

Now we can run XP-EHH giving the resulting file as input.
This may take some time...choose one population!...

If you really cannot wait, then I provided you with the output file for CHB and NAM in Files, so you can run something like:
```
cp $DIR/Files/*.xpehh.out Results/.
```
or run it yourself (choose one population, this takes a while! Launch it and take a break):
```
# NAM
$SS/selscan --xpehh --vcf $DATA/NAM.chr2.vcf --vcf-ref $DATA/TSI.chr2.vcf --map Results/genetic.map --out Results/NAM --threads 1
# CHB
$SS/selscan --xpehh --vcf $DATA/CHB.chr2.vcf --vcf-ref $DATA/TSI.chr2.vcf --map Results/genetic.map --out Results/CHB --threads 1
# CLM
$SS/selscan --xpehh --vcf $DATA/CLM.chr2.vcf --vcf-ref $DATA/TSI.chr2.vcf --map Results/genetic.map --out Results/CLM --threads 1
# LWK
$SS/selscan --xpehh --vcf $DATA/LWK.chr2.vcf --vcf-ref $DATA/TSI.chr2.vcf --map Results/genetic.map --out Results/LWK --threads 1
```
The output file has the header:
`< locusID > < physicalPos > < geneticPos > < popA ’1 ’ freq > < ihhA > < popB ’1 ’ freq > < ihhB > < unstandardized XPEHH >`
```
# NAM
less -S Results/NAM.xpehh.out
# CHB
less -S Results/CHB.xpehh.out
# ...
```

Finally, these unstandardized scores are normalized across the entire genome.
Here such normalisation makes little sense since we are just analysing a small portion of the chromsome, but let us see how to do it anyway.
This can be achieved by the extra command `norm`:
```
$SS/norm --help
```
We normalise the results (knowing that this should be done genome-wide):
```
# NAM
$SS/norm --xpehh --files Results/NAM.xpehh.out
# CHB
$SS/norm --xpehh --files Results/CHB.xpehh.out
# CLM
# $SS/norm --xpehh --files Results/CLM.xpehh.out
# LWK
# $SS/norm --xpehh --files Results/LWK.xpehh.out
```
Have a look at the results:
```
# NAM
less -S Results/NAM.xpehh.out.norm
# CHB
less -S Results/CHB.xpehh.out.norm
# ...
```

In order to perform a sliding window scan we can assign the maximum XP-EHH score to each window and plot it.
Please note that here we may be better off using unstandardised values here as we are focusing on one small region.
```
# NAM
Rscript $DIR/Scripts/plotXPEHH.R Results/NAM.xpehh.out.norm Results/NAM.xpehh.pdf
evince Results/NAM.xpehh.pdf
# CHB
Rscript $DIR/Scripts/plotXPEHH.R Results/CHB.xpehh.out.norm Results/CHB.xpehh.pdf
evince Results/CHB.xpehh.pdf
# CLM
# Rscript $DIR/Scripts/plotXPEHH.R Results/CLM.xpehh.out.norm Results/CLM.xpehh.pdf
# evince Results/CLM.xpehh.pdf
# LWK
# Rscript $DIR/Scripts/plotXPEHH.R Results/LWK.xpehh.out.norm Results/LWK.xpehh.pdf
# evince Results/LWK.xpehh.pdf
```

**QUESTIONS**

What conclusions can you make from these plots?

----------------------------------------------------------

**VERY OPTIONAL**

As a further illustration, we are here calculating nSL statistic.
Note that we have already performed a filtering on these VCF files.
First, look at the options in selscan:
```
$SS/selscan --help
```
The basic usage of selscan to compute nSL from VCF file is:
```
# $SS/selscan --nsl --vcf Data/NAM.chr2.vcf --out Results/NAM
```
As you can see, unlike for iHS and EHH, we are not required to provide a genetic map since values are computed in windows based on physical distances.

For nSL, the integration is not cut after a certain threshold, but we need to set a value of maximum length for the window (in number of SNPs units) for building the haplotype.
This is set with the option `--max-extend-nsl`.
It is also common to filter out variant with very low frequency.
Therefore our command line might be:
```
$SS/selscan --nsl --vcf $DATA/NAM.chr2.vcf --out Results/NAM --max-extend-nsl 200 --maf 0.02
```
Have a look at the output file, knowning that the header is:
`< locusID > < physicalPos > < ’1 ’ freq > <sl1 > <sl0 > < unstandardized nSL >`
```
less -S Results/NAM.nsl.out
```

Finally, these unstandardized scores are normalized in allele frequency bins across the entire genome.
This can be achieved by the extra command `norm`:
```
$SS/norm --help
...
--bins <int>: The number of frequency bins in [0,1] for score normalization.
        Default: 100
...
```
Thus, our command would be (note tha ihs and nsl normalisation are equivalent):
```
$SS/norm --ihs --files Results/NAM.nsl.out --bins 20
```
The output file is called `Results/NAM.nsl.out.20bins.norm`.
```
less -S Results/NAM.nsl.out.20bins.norm
```
Please note that, again, we should perform the normalization on genome-wide data (or at least chromosome-wide).

One strategy to plot such statistic is to record the number of SNPs above a certain threshold, e.g. |nSL|>2.
We can do this using the following script:
```
Rscript $DIR/Scripts/plotnSL.R Results/NAM.nsl.out.20bins.norm Results/NAM.nsl.pdf
evince Results/NAM.nsl.pdf
```

What is happening here? What is wrong?
Note that in the VCF alleles have been polarised according to the reference sequence, and not the ancestral one!
Therefore the standardisation into bins of "non-reference" allele frequency is wrong.

To overcome this issue, either you should use absolute values of nSL and perform a normalisation on the minor allele frequency rather than the non-reference allele frequency.
This exercise stresses the importance of knowing that the reference sequence is just an arbitrary choice of alleles and should NOT be used for evolutionary inferences.

-------------------------------

Finally, we are interested in investigating the **haplotype distribution** in this region.
Specifically, we want to draw a haplotype network, where all (unique) haplotypes are clustered based on their mutual genetic distance.

First, we convert our VCF files into FASTA files.
```
> Results/EDAR.fa
Rscript $DIR/Scripts/vcf2fasta.R $DATA/NAM.edar.vcf NAM Results/NAM.edar.snp >> Results/EDAR.fa
Rscript $DIR/Scripts/vcf2fasta.R $DATA/TSI.edar.vcf TSI NULL >> Results/EDAR.fa
```
Have a look at the resulting file:
```
less -S Results/EDAR.fa
```
We are using [pegas](https://bioinformatics.oxfordjournals.org/content/26/3/419.full) package in R to draw haplotype network.
```
Rscript $DIR/Scripts/plotNet.R Results/EDAR.fa Results/NAM.edar.snp Results/EDAR.pdf > Results/EDAR.diff
```
Open the plot:
```
evince Results/EDAR.pdf
```
Each unique haplotype is represented as a circle whose size is proportional to its frequency.

Our variant of interest is rs3827760:
```
grep rs3827760 Results/NAM.edar.snp 
# 2	109513601	rs3827760	A	G
```

------------------------

[HOME](https://github.com/mfumagalli/Copenhagen)




