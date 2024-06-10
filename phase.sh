# The vcf files should be split by chromosome

MY_VCF = "PATH_TO_VCF"
RECOMBINATION_MAP = "PATH_TO_RECOMBINATIONMAP"
OUTPUT_BEAGLE = 'PATH_TO_BEAGLE_OUTPUT'
OUTPUT_SHAPEIT = "PATH_TO_SHAPEIT_OUTPUT"

# phase with beagle 4
java -Xmx64000m -jar Software/beagle.r1399.jar \
gt=$MY_VCF \
out=$OUTPUT_BEAGLE \
nthreads=20 \
ped=helperfiles/pedfile_long.txt


# phase with shapeit4
shapeit4 \
--input $MY_VCF \
--map $RECOMBINATION_MAP \
--region chr22 \ # if you are using chromosome 22 for instance
--thread 20 \
--output $OUTPUT_SHAPEIT
