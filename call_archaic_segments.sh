# a step-by-step guide on how to run hmmix can be found here:
# https://github.com/LauritsSkov/Introgression-detection

# here we show the steps to infer archaic admixture


# STEP 1 - download relavant files from: https://zenodo.org/records/11212339 (hg38 reference genome build). 
# remember to download archaic variant information from here hg38 (https://zenodo.org/records/10806726)
# Install hmmix
pip install hmmix

# STEP 2 - remove variants found in outgroup
hmmix create_ingroup  -ind=NAME_OF_IND \
	-vcf=PATH_TO_LASIDAD*.bcf \
	-weights=hg38_strick_callability_mask.bed \
	-out=obs \
	-outgroup=hg38_Outgroup_1000g_HGDP.txt \
	-ancestral=hg38_ancestral/homo_sapiens_ancestor_*.fa

# STEP 3 - train on each individuals
hmmix train  \
	-obs=obs.NAME_OF_IND.txt \
	-weights=hg38_strick_callability_mask.bed \
	-mutrates=hg38_mutationrate.bed \
	-out=trained.NAME_OF_IND.json \
	-haploid

# STEP 4 - decode on each individual
hmmix decode \
	-obs=obs.NAME_OF_IND.txt \
	-weights=hg38_strick_callability_mask.bed \
	-mutrates=hg38_mutationrate.bed \
	-param=trained.NAME_OF_IND.json \
	-admixpop=archaicvar/*bcf \
	-haploid
