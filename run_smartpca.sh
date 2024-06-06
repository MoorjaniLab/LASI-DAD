##Script to create a parfile and run SMARTPCA, part of the EIGENSOFT package, can be download here: https://github.com/DReichLab/EIG
#Need data in eigentrat format, more about the format : https://reich.hms.harvard.edu/software/InputFileFormats


eigenstrat_file=$1

list_of_population=$2 ##List of population you want to based the PCA on

DIR="/your_directory"
pop="${list_of_population}"
        # create parfile
        echo "DIR: ${DIR} " >par_${pop}
        echo "SSS: ${eigenstrat_file}   " >>par_${pop}
        echo "S1:  out_${eigenstrat_file}_${pop} " >>par_${pop}
        echo "genotypename: DIR/SSS.geno" >>par_${pop} 
        echo "snpname: DIR/SSS.snp" >>par_${pop}
        echo "indivname: DIR/SSS.ind" >>par_${pop}
        echo "evecoutname: S1.evec" >>par_${pop}
        echo "evaloutname: S1.eval" >>par_${pop}
        echo "nsnpldregress: 2       " >>par_${pop}
        echo "snpweightoutname: S1_SNPWeightout.txt " >>par_${pop}
        echo "poplistname: DIR/${list_of_population}.txt " >>par_${pop}
        echo "numoutliter: 0" >>par_${pop}
        echo "numoutleigs: 0" >>par_${pop}
        echo "numoutevec: 10" >>par_${pop}
        echo "pubmean: YES" >>par_${pop}
        echo "fastmode: YES" >>par_${pop}
        echo "fstonly: NO" >>par_${pop}

        # Run smartpca
        EIG/bin/smartpca -p par_${pop}

echo "done"
