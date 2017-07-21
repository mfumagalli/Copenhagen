
CHROM=$1
START=$2
END=$3
SAMTOOLS=$4

echo Retrieving file names...
NS=20
# get IDs 
for POP in LWK TSI CHB PEL;
do

	INPUT=$POP.txt
	INDS=`cat $INPUT`
	OUTPUT=$POP.BAMs.txt
	> $OUTPUT
	> tmp

	for i in $INDS;
	do
		grep $i Files/phase3_bamlist.txt >> tmp
	done
	# wc -l tmp
	head -n $NS tmp > $OUTPUT 
done
cp Files/NAM.BAMs.txt .

echo Downloading BAM files...
# download and index bams
for POP in LWK TSI CHB PEL NAM;
do
	mkdir Data/$POP.BAMs
	echo $POP
	INDLIST=`cat $POP.BAMs.txt`
	for i in $INDLIST;
	do
		NAME=`echo -n $i | tail -c 58`
		echo $NAME
		$SAMTOOLS view -s 1 -h -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/$i $CHROM:$START-$END > Data/$POP.BAMs$NAME 2> /dev/null
	done
done

echo Removing index files...
rm *.bai

exit



