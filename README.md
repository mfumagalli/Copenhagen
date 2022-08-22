
# Summer course in analysis of high throughput data for population genetics


Tuesday - Analysis of NGS data and population structure
 - 9am - 12pm: Estimation of allele frequencies, SNP calling and genotype calling from NGS data (theory 9-10.15 + practical 10.30-12.00)

Wednesday - Research lecture
 - 4.15pm - 5pm: Research lecture

Thursday - Demography and selection
 - 1pm - 4pm: Selection scans (theory 1-2.15 + practical 2.15-4)

## Material

Slides can be found [here](https://github.com/mfumagalli/Copenhagen/tree/master/Slides).

The data for practicals has been already downloaded and it is provided in `/ricco/data/matteo/Data`.
These instructions, including all relevant files and scripts, can be found at `/home/matteo/Copenhagen`.
In short, you don't have to worry about anything for the practicals.

### Data

As an illustration, we will use 40 BAM files of human samples (of African, European, East Asian, and Native American descent), a reference genome, and putative ancestral sequence.
We will also use 10 BAM files of Latinos in one example.
*To make things more interesting, we have downsampled our data to an average mean depth of 2X!*.

We will also use VCF files for 120 individuals from the same populations.
The human data represents a small genomic region (1MB on chromosome 2) extracted from the [1000 Genomes Project](http://www.internationalgenome.org/home).

## Preparation

For each day there will be indications on which software and scripts you will be using.
However, before doing anything else, please create a folder where you will put all the results and some temporary data.
```
mkdir Results
mkdir Data
```
That's it.

## Case study

*MOTIVATION*

Detecting signatures of natural selection in the genome has the twofold meaning of (i) understanding which adaptive processes shaped genetic variation and (ii) identifying putative functional variants.
In case of humans, biological pathways enriched with selection signatures include pigmentation, immune-system regulation and metabolic processes.
The latter may be related to human adaptation to different diet regimes, depending on local food availability (e.g. the case of lactase persistence in dairy-practicing populations).

The human Ectodysplasin A receptor gene, or EDAR, is part of the EDA signalling pathway which specifies prenatally the location, size and shape of ectodermal appendages (such as hair follicles, teeth and glands).
EDAR is a textbook example of positive selection in East Asians (Sabeti et al. 2007 Nature) with tested phenotypic effects (using transgenic mice).

Recently, a genome-wide association study found the same functional variant in EDAR associated to several human facial traits (ear shape, chin protusion, ...) in Native American populations (Adhikari et al. Nat Commun 2016).

![](Slides/NGS_analysis/Pics/practical.png)


*HYPOTHESIS*

- Is the functional allele in East Asian at high frequency in other human populations (e.g. Native Americans)?
- Can we identify signatures of natural selection on EDAR in Native Americans?
- Is selection targeting the same functional variant?

*CHALLENGES*
- Admixed population
- Low-depth sequencing data
- Effect of genetif drift
- ...

*PLAN OF ACTION*

Goal day 1:

- Estimate allele frequencies for tested variant for African, European, East Asian and Native American samples from low-depth sequencing data

Optional:
- Investigate population structure of American samples related to Europeans and Africans
- Select individuals with high Native American ancestry

Goal day 2:

- Perfom a sliding windows scan based on allele frequency differentiation
- Assess statistical significance of selection signatures through simulations
- Test for extended haplotype homozygosity on high-depth sequencing data

## Agenda

### Tuesday morning -  introduction to NGS data

#### Lecture

* Maximum likelihood and Bayesian estimation
* Genotype likelihoods
* Allele frequencies, SNPs and genotypes calling

#### [Practical](Files/day1.md)

* Basic data filtering
* Estimation of allele frequencies and SNP calling
* Genotype calling
* Example: estimation of allele frequencies from low-depth sequencing data: the case of _EDAR_ genetic variation in Native Americans

#### Research lecture

* TBA

### Thursday afternoon - detecting selection

#### Lecture

* The effect of selection on the genome
* Methods to detect selection signals
* The problem of assessing significance
* Bias introduced by NGS data
* Summary statistics from low-depth data

#### Practical

* Selection scan based on genetic [differentiation and diversity](Files/day2a.md) from low-depth data
* Assessing significance through [simulations](Files/day2b.md)
* Selection test based on [haplotype](Files/day2c.md) diversity
* Example: detection of natural selection from low-depth sequencing data and haplotype data: the case of EDAR genetic variation in Native Americans

## Credits

Thanks to Thorfinn Korneliussen, Anders Albrechtsen, Tyler Linderoth, Filipe G. Vieira, Dean Ousby, Javier Mendoza, Ryan Waples, and possibly many others I forgot to mention.


