
# ex1
```
for POP in AFR EUR EAS LAT NAM
do
        echo $POP
        $ANGSD/angsd -b $DATA/$POP.bams -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -doGeno 3 -doPost 1 -postCutoff 0.50 \
                -sites Data/snp.txt
done
```

```
for POP in AFR EUR EAS LAT NAM
do
        echo $POP
        zcat Results/$POP.geno.gz
done
```

For instance, you may have 0/20 in AFR and EUR, 20/20 in EAS, while there are only 4 called genotypes in NAM.
Recall that we previously estimated a minor allele frequency of 0.84% in NAM without assigning individuals.

# quick

```
$ANGSD/angsd -b $DATA/EUR.bams -ref $REF -out Results/EUR \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 2 -doMajorMinor 1 -doMaf 1 \
        -minMaf 0.05
```

```
zcat Results/EUR.mafs.gz | head
```
```
chromo  position        major   minor   ref     knownEM nInd
2       109000113       C       T       C       0.060517        10
2       109000117       G       T       G       0.062740        10
2       109000122       T       C       T       0.064025        10
2       109000160       C       A       C       0.074944        9
2       109000285       A       C       A       0.153066        5
2       109000483       C       T       C       0.120045        6
2       109000505       A       T       A       0.109368        7
2       109000539       G       T       G       0.110355        7
2       109000688       A       T       A       0.401350        7
```

How many entries (potential SNPs) you have?
```
zcat Results/EUR.mafs.gz | tail -n +2 | wc -l
```

# ex2
```
for POP in AFR EUR EAS LAT NAM
do
        echo $POP
        $ANGSD/angsd -b $DATA/$POP.bams -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 1 -doHWE 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 \
                -sites Data/snp.txt
done
```
An assessment of the deviation from HWE will be print out in files with extension `.hwe.gz`.

You can inspect the results.
```
zcat Results/AFR.mafs.gz Results/EUR.mafs.gz Results/EAS.mafs.gz Results/LAT.mafs.gz Results/NAM.mafs.gz
```








