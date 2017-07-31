
## Pipeline to download and process the data to be used for this workshop

## Ignore it unless you know what you are doing!

# set path
SAMTOOLS=/data/data/Software/samtools-1.4.1/samtools
#SAMTOOLS=samtools
echo Is this your path to samtools? $SAMTOOLS
BGZIP=bgzip
echo Is this your path to bgzip? $BGZIP

VCFLIB=/data/data/Software/vcflib/bin
echo Is this your path to VCFlib? $VCFLIB

# chrom start end
CHROM=2
START=109000000
END=110000000

mkdir Data
mkdir Results

# get unrelated samples iDs
Rscript Scripts/getUnrelated.R

# download BAM files
bash Scripts/getBams.sh $CHROM $START $END $SAMTOOLS
# this creates files and folders in Data/PEL.BAMs/* and TSI and LWK and CHB

# create file with list of BAMs
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/CHB.BAMs/*.bam Data/NAM.BAMs/*.bam > Data/ALL.bamlist
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/CHB.BAMs/*.bam Data/PEL.BAMs/*.bam > Data/ALL2.bamlist
ls Data/LWK.BAMs/*.bam > Data/LWK.bamlist
ls Data/TSI.BAMs/*.bam > Data/TSI.bamlist
ls Data/CHB.BAMs/*.bam > Data/CHB.bamlist
ls Data/PEL.BAMs/*.bam > Data/PEL.bamlist
ls Data/NAM.BAMs/*.bam > Data/NAM.bamlist
mv *.txt Data/. # bam and names lists

# perhaps use only 10 samples per population (these files are called .bams)
# also change their labels to human-readable
# also remove PEL from ALL
ls Data/LWK.BAMs/*.bam | head -n 10 > Data/AFR.bams
ls Data/TSI.BAMs/*.bam | head -n 10 > Data/EUR.bams
ls Data/CHB.BAMs/*.bam | head -n 10 > Data/EAS.bams
ls Data/PEL.BAMs/*.bam | head -n 10 > Data/LAT.bams
ls Data/NAM.BAMs/*.bam | head -n 10 > Data/NAM.bams
cat Data/AFR.bams Data/EUR.bams Data/EAS.bams Data/NAM.bams > Data/ALL.bams

# download ancestral sequence
echo Downloading and processing ancestral sequence...
wget http://popgen.dk/software/download/angsd/hg19ancNoChr.fa.gz
zcat hg19ancNoChr.fa.gz > Data/anc.fa
$BGZIP Data/anc.fa
$SAMTOOLS faidx Data/anc.fa.gz
rm hg19ancNoChr.fa.gz

# download reference sequence
echo Downloading and processing reference sequence...
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz &> /dev/null
zcat human_g1k_v37.fasta.gz > Data/ref.fa 2> /dev/null
$BGZIP Data/ref.fa
$SAMTOOLS faidx Data/ref.fa.gz
rm human_g1k_v37.fasta.gz

# download VCF files
echo Downloading VCF file...
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr$CHROM.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz &> /dev/null
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr$CHROM.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi &> /dev/null

echo Processing VCF file...
VCF=ALL.chr$CHROM.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# chrom start end
CHROM=2
START=109000000
END=110000000

# whole region for selscan
$VCFLIB/vcffilter -f "VT = SNP" -f "AC > 1" -f "QUAL > 20" -f "DP > 40" -f "DP < 40000" -r $CHROM:$START-$END $VCF > ALL.chr$CHROM.tmp.vcf
grep -v MULTI_ALLELIC ALL.chr$CHROM.tmp.vcf > Data/ALL.chr$CHROM.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr$CHROM.vcf HG01926 HG01961 HG01938 HG01920 HG02259 HG02271 HG02291 HG01968 HG01951 HG01572 HG01974 HG02275 HG01923 HG01927 HG02105 HG02266 HG02265 HG02147 HG01954 HG02278 HG01992 HG02104 HG01953 HG02146 HG01997 HG02292 HG02299 HG01917 HG01950 HG01942 > Data/NAM.chr$CHROM.vcf
sed -i -e 's/\t\./\t0|0/g' Data/NAM.chr$CHROM.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr$CHROM.vcf NA20502 NA20503 NA20504 NA20505 NA20506 NA20507 NA20508 NA20509 NA20510 NA20511 NA20512 NA20513 NA20514 NA20515 NA20516 NA20517 NA20518 NA20519 NA20520 NA20521 NA20522 NA20524 NA20525 NA20526 NA20527 NA20528 NA20529 NA20530 NA20531 NA20532 > Data/TSI.chr$CHROM.vcf
sed -i -e 's/\t\./\t0|0/g' Data/TSI.chr$CHROM.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr$CHROM.vcf NA18524 NA18525 NA18526 NA18527 NA18528 NA18529 NA18530 NA18531 NA18532 NA18533 NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18540 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 > Data/CHB.chr$CHROM.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CHB.chr$CHROM.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr$CHROM.vcf HG01112 HG01113 HG01119 HG01121 HG01122 HG01124 HG01125 HG01130 HG01131 HG01133 HG01134 HG01136 HG01137 HG01139 HG01140 HG01142 HG01148 HG01149 HG01250 HG01251 HG01253 HG01254 HG01256 HG01257 HG01259 HG01260 HG01269 HG01271 HG01272 HG01274 > Data/CLM.chr$CHROM.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CLM.chr$CHROM.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr$CHROM.vcf NA19313 NA19331 NA19381 NA19382 NA19432 NA19434 NA19444 NA19445 NA19453 NA19469 NA19470 NA19017 NA19019 NA19020 NA19023 NA19024 NA19025 NA19026 NA19027 NA19028 NA19030 NA19031 NA19035 NA19036 NA19037 NA19038 NA19041 NA19042 NA19043 NA19044 > Data/LWK.chr$CHROM.vcf
sed -i -e 's/\t\./\t0|0/g' Data/LWK.chr$CHROM.vcf


# EDAR is our target gene
# chr2:109,510,927-109,605,828
START=109510927
END=109605828
GENE=edar

$VCFLIB/vcffilter -f "VT = SNP" -f "AC > 1" -f "QUAL > 30" -f "DP > 50" -f "DP < 40000" -r $CHROM:$START-$END $VCF > ALL.$GENE.tmp.vcf
grep -v MULTI_ALLELIC ALL.$GENE.tmp.vcf > Data/ALL.$GENE.vcf

$VCFLIB/vcfkeepsamples Data/ALL.$GENE.vcf HG01926 HG01961 HG01938 HG01920 HG02259 HG02271 HG02291 HG01968 HG01951 HG01572 HG01974 HG02275 HG01923 HG01927 HG02105 HG02266 HG02265 HG02147 HG01954 HG02278 HG01992 HG02104 HG01953 HG02146 HG01997 HG02292 HG02299 HG01917 HG01950 HG01942 > Data/NAM.$GENE.vcf
sed -i -e 's/\t\./\t0|0/g' Data/NAM.$GENE.vcf

$VCFLIB/vcfkeepsamples Data/ALL.$GENE.vcf NA20502 NA20503 NA20504 NA20505 NA20506 NA20507 NA20508 NA20509 NA20510 NA20511 NA20512 NA20513 NA20514 NA20515 NA20516 NA20517 NA20518 NA20519 NA20520 NA20521 NA20522 NA20524 NA20525 NA20526 NA20527 NA20528 NA20529 NA20530 NA20531 NA20532 > Data/TSI.$GENE.vcf
sed -i -e 's/\t\./\t0|0/g' Data/TSI.$GENE.vcf

$VCFLIB/vcfkeepsamples Data/ALL.$GENE.vcf NA18524 NA18525 NA18526 NA18527 NA18528 NA18529 NA18530 NA18531 NA18532 NA18533 NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18540 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 > Data/CHB.$GENE.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CHB.$GENE.vcf

$VCFLIB/vcfkeepsamples Data/ALL.$GENE.vcf HG01112 HG01113 HG01119 HG01121 HG01122 HG01124 HG01125 HG01130 HG01131 HG01133 HG01134 HG01136 HG01137 HG01139 HG01140 HG01142 HG01148 HG01149 HG01250 HG01251 HG01253 HG01254 HG01256 HG01257 HG01259 HG01260 HG01269 HG01271 HG01272 HG01274 > Data/CLM.$GENE.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CLM.$GENE.vcf

$VCFLIB/vcfkeepsamples Data/ALL.$GENE.vcf NA19313 NA19331 NA19381 NA19382 NA19432 NA19434 NA19444 NA19445 NA19453 NA19469 NA19470 NA19017 NA19019 NA19020 NA19023 NA19024 NA19025 NA19026 NA19027 NA19028 NA19030 NA19031 NA19035 NA19036 NA19037 NA19038 NA19041 NA19042 NA19043 NA19044 > Data/LWK.$GENE.vcf
sed -i -e 's/\t\./\t0|0/g' Data/LWK.$GENE.vcf

# clean up
rm $VCF $VCF.tbi ALL.chr$CHROM.tmp.vcf ALL.$GENE.tmp.vcf

echo Downloading recombination map
mkdir Map
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/PEL_omni_recombination_20130507.tar &> /dev/null
tar -C Map -xvf PEL_omni_recombination_20130507.tar &> /dev/null
zcat Map/PEL/PEL-$CHROM-final.txt.gz > Data/genetic_map_chrom$CHROM.map
rm -rf Map
rm PEL_omni_recombination_20130507.tar

echo Done!
ls -lh Data/* > Data/download.log
echo Open Data/download.log to see what files have been generated.


