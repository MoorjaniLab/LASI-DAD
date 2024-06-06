#!/bin/sh
##Script to create the parfile and to run ALDER, to download : https://groups.csail.mit.edu/cb/alder/
#Need data in eigentrat format, more about the format : https://reich.hms.harvard.edu/software/InputFileFormats

studied_pop=$1 ##the mix pop to study
ref_pop=$2 ##the reference population

        par="par${studied_pop}"

        # create parfile
        echo "DIR: /your_directory" >${par}
        echo "SSS: eigenstratfile_name" >> ${par}
        echo "indivname:    DIR/SSS.ind" >> ${par}
        echo "snpname:      DIR/SSS.snp" >> ${par}
        echo "genotypename:      DIR/SSS.geno" >> ${par}
        echo "admixpop: ${studied_pop}" >> ${par}
        echo "refpops:  ${ref_pop}" >> ${par}
        echo "binsize:      .001" >> ${par}
        echo "raw_outname:    ${studied_pop}_${anc1}.out" >> ${par}
        echo "maxdis:       1.0" >> ${par}
        echo "seed:         77" >> ${par}
        echo "runmode:      1" >> ${par}
        echo "chithresh: 0.0" >> ${par}
        echo "zdipcorrmode: YES" >> ${par}
        echo "jackknife: YES" >> ${par}
        echo "mincount: 4" >> ${par}

##run ALDER
/bin/alder -p ${par}

echo "done"
