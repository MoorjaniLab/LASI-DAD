##Script to create a parfile and run qpAdm with our without the "allsnps:YES" option
## Software part of the ADMIXTOOLS package, can be download here: https://github.com/DReichLab/AdmixTools
#Need data in eigentrat format, more about the format : https://reich.hms.harvard.edu/software/InputFileFormats

##population or individual of interest
pop=$1

#file to use (name without .ind, .geno or .snp)
eigenstrat_file=$2

        par="par_${pop}"


        ##List of the population/individual to test and the source populations, in this example Sarazm_EN, Central_Steppe_MLBA and Onge, can be done separately
        list="list_${pop}.txt"        
        echo -e "${ind}\nSarazm_EN\nCentral_Steppe_MLBA\nOnge" > ${list}
        
##With "allsnps: NO"
        echo "DIR: /your_directory" >${par}
        echo "S1: ${eigenstrat_file}" >>${par}
        echo "indivname: DIR/S1.ind" >>${par}
        echo "snpname: DIR/S1.snp" >>${par}
        echo "genotypename: DIR/S1.geno" >>${par}
        echo "popright: rigthpop_list.txt ">>${par} ## list of outgroup populations (rigth populations)
        echo "popleft:  ${list}" >>${par}
        echo "hashcheck: NO" >> ${par}
        echo "details:       YES " >> ${par}    
        echo "allsnps: NO" >> ${par} 
        echo "inbreed: YES" >> ${par}   ##Need if use of ancient DNA 
        
        AdmixTools/bin/qpAdm -p ${par}
        echo "done" ${par}

##With "allsnps: YES"
        echo "DIR: /your_directory" >${par}
        echo "S1: ${eigenstrat_file}" >>${par}
        echo "indivname: DIR/S1.ind" >>${par}
        echo "snpname: DIR/S1.snp" >>${par}
        echo "genotypename: DIR/S1.geno" >>${par}
        echo "popright: rigthpop_list.txt ">>${par} ## list of outgroup populations (rigth populations)
        echo "popleft:  ${list}" >>${par}
        echo "hashcheck: NO" >> ${par}
        echo "details:       YES " >> ${par}    
        echo "allsnps: YES" >> ${par} 
        echo "inbreed: NO" >> ${par}
        echo "inbreedlistname: rigthpop_list_inbreed.txt" >> ${par}       ##list of ancient DNA population with more than one individual
        echo "basepop: popname"  >> ${par} ## Population to compute qpfstats on, then use these fstats to use the allspns: YES option.

        qpAdm -p ${par}
        echo "done" ${par}

