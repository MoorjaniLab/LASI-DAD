
for chr in {1..22}
do
        #genetic map in PLINK format, can be download here: https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
        map_plinkformat="/plink.GRCh38.map/chr_plink.chr${chr}.GRCh38.map" 

        ##Need phased vcf
        vcf="/PhasedVcf/chr${chr}.myvcf.gz"

        #hap-ibd can be download here : https://github.com/browning-lab/hap-ibd      
        java -jar hap-ibd.jar  gt=${vcf} map=${map_plinkformat} out=hapIBD_out_chr${chr}        min-seed=0.5 min-output=1
done
