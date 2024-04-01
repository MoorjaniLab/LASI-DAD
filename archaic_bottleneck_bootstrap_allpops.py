#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 16:40:14 2023

@author: elisekerdoncuff
"""

#Simple demography + Bottenecks to compare the amount of common variants In LASI-DAD dataset: 95% of genome with Neanderthal frequency < 0.06679

##Archaic simulations 

import tskit
import msprime
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import demes
import demesdraw
import time
import random


##Introgression functions
##from : https://tskit-dev.github.io/tutorials/introgression.html

def get_migrating_tracts(ts):
    neanderthal_id = [p.id for p in ts.populations() if p.metadata['name']=='N'][0] ##find the id of Neanderthal
    migrating_tracts = []
    # Get all tracts that migrated into the neanderthal population (backwards)
    for migration in ts.migrations():
        if migration.dest == neanderthal_id:
            migrating_tracts.append((migration.left, migration.right))
    return np.array(migrating_tracts)

##In the migrating tracts, get the ones that coalesce in Neanderthal pop
def get_coalescing_tracts(ts):
    neanderthal_id = [p.id for p in ts.populations() if p.metadata['name']=='N'][0]
    coalescing_tracts = []
    tract_left = None
    for tree in ts.trees():
        # 0 is the samplepop sample and 2 is the Neanderthal
        mrca_pop = ts.node(tree.mrca(0, 2)).population
        left = tree.interval[0]
        if mrca_pop == neanderthal_id and tract_left is None:
            # Start a new tract
            tract_left = left      
        elif mrca_pop != neanderthal_id and tract_left is not None:
            # End the last tract
            coalescing_tracts.append((tract_left, left))
            tract_left = None
    if tract_left is not None:
        coalescing_tracts.append((tract_left, ts.sequence_length))
    return np.array(coalescing_tracts)

##get tracts for which the sample individual is closer to Neanderthal than to the AMH (by comparing to the root)
##May need to change that part if considering mutliple individuals
def get_samplepop_nea_tracts(ts):
    tracts = []
    tract_left = None
    for tree in ts.trees():
        # 0 is the samplepop sample and 2 is the Neanderthal
        mrca = tree.mrca(0, 2)
        left = tree.interval[0]
        if mrca != tree.root and tract_left is None:
            # Start a new tract
            tract_left = left      
        elif mrca == tree.root and tract_left is not None:
            # End the last tract
            tracts.append((tract_left, left))
            tract_left = None
    if tract_left is not None:
        tracts.append((tract_left, ts.sequence_length))
    return np.array(tracts)

                                    
##My function, similar to get_migrating_tracts, but using tree.sequence (useful for get_migrating_tracts_count)
def get_migrating_tracts_per_indiv(ts,indiv_i):
    neanderthal_id = [p.id for p in ts.populations() if p.metadata['name']=='N'][0] ##find the id of Neanderthal
    migrating_tracts = []
    # Get all tracts that migrated into the neanderthal population (backwards)
    for migration in ts.migrations():
        if migration.dest == neanderthal_id:
            mig_node=migration.node
            for tree in ts.trees():
                if tree.interval.left>migration.right:
                    break
                elif tree.interval.left>=migration.left and tree.interval.right<=migration.right and tree.is_descendant(indiv_i,mig_node):
                    migrating_tracts.append((tree.interval.left, tree.interval.right))
    return np.array(migrating_tracts)


##My function, input: msprime simulation with only sample from introgressed population + 1 Archaic (Neanderthal in this case)
##ouput: for each migrating segment, the length and the number of Individuals who carried it (number of leaves below the migration node)
def get_migrating_tracts_count(ts):
    neanderthal_id = [p.id for p in ts.populations() if p.metadata['name']=='N'][0] ##find the id of Neanderthal
    migrating_tracts = []
    n_leaves=[]
    # Get all tracts that migrated into the neanderthal population (backwards)
    for migration in ts.migrations():
        if migration.dest == neanderthal_id:
            mig_node=migration.node
            #find the number of leaves below the node
            for tree in ts.trees():
                if tree.interval.left>migration.right:
                    break
                elif tree.interval.left>=migration.left and tree.interval.right<=migration.right:
                    migrating_tracts.append((tree.interval.left, tree.interval.right))
                    n_leaves.append(len(list(tree.leaves(mig_node))))
    return np.array(migrating_tracts),np.array(n_leaves)

#Get ancestry tract for one individual when multiples are samples
def get_migrating_tracts_oneind(ts,ind1):
    neanderthal_id = [p.id for p in ts.populations() if p.metadata['name']=='N'][0] ##find the id of Neanderthal
    migrating_tracts = []
    # Get all tracts that migrated into the neanderthal population (backwards)
    for migration in ts.migrations():
        if migration.dest == neanderthal_id:
            mig_node=migration.node
            ##find if ind is leave of node
            for tree in ts.trees():
                if tree.interval.left>=migration.left and tree.interval.right<=migration.right and ind1 in list(tree.leaves(mig_node)):
                    migrating_tracts.append((tree.interval.left, tree.interval.right))
    return np.array(migrating_tracts)

###Count % genome above thresold, from output of get_migrating_tracts_count():
def inarchaic_percentage_above_limit(freqs,length_tracts,lim):
    den=0
    num=0
    for i in range(len(freqs)):
        if (np.isnan(length_tracts[i])==False):
            den+=length_tracts[i]
            if (freqs[i]>lim):
                num+=length_tracts[i]
            
    return num/den

##compute value of percentile from count_per_bin (in frequence) and bin_width:
def percentile_value(count_per_bin, bin_width, percentile):
    cumm=0
    i=0
    while(cumm+count_per_bin[i]<percentile):
        cumm+=count_per_bin[i]
        i+=1     
    
    return max((i-1)*bin_width+((percentile-cumm)/count_per_bin[i])*bin_width,0)


####
## LASI-DAD info
####

genome_wide_95=0.06679
archaic_95= 0.08426

####
## INIT
####
n_simu=50
i_simu=50
n_boot=50
bin_width=0.01

sum_results="results_metrics_per_scenario_bt"
outfile_gw95="gw_95_bt"
outfile_gw20="gw_above20_bt"
outfile_seg="seg_simu_bt"




####
## Sscenario
####

# Simple demography (modified from Laurits file)
# -----------------------------------------------------------------------------------------------------
# Parameters for demography
# -----------------------------------------------------------------------------------------------------

# Times are provided in years, so we convert into generations.
gen_time = 29.0 #generation time
rec_rate = 1.2e-8 #recombination rate
mutation_rate = 1.2e-8 #mutation rate

# Population sizes
N_Neanderthals = 1_500
N_Africa = 20_000
N_Europe = 10_000
N_Asia = 10_000

N_ancestral_eurasians = 2_000 
N_ancestral_humans = 7000 
#N_bottleneck_nonAfricans = 250

# Split Times - converting years BP into generations BP
T_Human_Neandertal = 575000 / gen_time 
T_out_of_Africa = 60_000 / gen_time 
T_Eurasian_bottleneck = 50_000 / gen_time ##end of Bottleneck ##before 58
T_europa_asia = 50_000 / gen_time


# Admixture times and proportions from archaic humans
admixtureproportion = 0.013 #(0.013 = 1.3% in this example)
T_geneflow_Neanderthal = 50_000 / gen_time 


sequence_length=10_000_000
length_chr14=106880170 ##as in chr14
sequence_length=length_chr14

n_sample_A=500
samples_A = msprime.SampleSet(n_sample_A,population="Asia")
samples_N = msprime.SampleSet(1,population="N",time=T_geneflow_Neanderthal)
sample_list_simple = [samples_A,samples_N]


##constant

demography_simple = msprime.Demography() #initializing an msprime "demography" object
# This is the "trunk" population that we merge other populations into
demography_simple.add_population(
    name="Africa",
    description="Africa",
    initial_size=N_Africa,
    initially_active=True,
);
demography_simple.add_population(
    name="Europe",
    description="European",
    initial_size=N_Europe,
    initially_active=True,
)
demography_simple.add_population(
    name="Asia",
    description="Asian",
    initial_size=N_Asia,
    initially_active=True,
)    
demography_simple.add_population(
    name="N",
    description="Introgressing-Neanderthal population",
    initial_size=N_Neanderthals,
    initially_active=True,
)


# Split Europeans and Asians
demography_simple.add_population_split(
    time=T_europa_asia, derived=["Asia"], ancestral="Europe"
)

# Add neanderthal geneflow for one generation (pulse model)
demography_simple.add_migration_rate_change(time = T_geneflow_Neanderthal, rate = admixtureproportion, dest = "N", source = "Europe")
demography_simple.add_migration_rate_change(time = T_geneflow_Neanderthal+1, rate = 0, dest = "N", source = "Europe")

# Add split out of Africa
demography_simple.add_population_split(
    time=T_out_of_Africa, derived=["Europe"], ancestral='Africa'
)
demography_simple.add_population_parameters_change(
    time=T_Eurasian_bottleneck, population="Europe", initial_size=N_ancestral_eurasians
)
demography_simple.add_population_parameters_change(
    time=T_out_of_Africa, population="Africa", initial_size=N_ancestral_humans
)
demography_simple.add_population_split(
    time=T_Human_Neandertal, derived=["N"], ancestral="Africa" 
)

demography_simple.sort_events()

print(demography_simple)

#plot the demography
graph_simple = demography_simple.to_demes()
fig, ax = plt.subplots(figsize=(10,6))  # use plt.rcParams["figure.figsize"]
demesdraw.tubes(graph_simple, ax=ax, seed=1)
plt.savefig("demo_simple.png")
plt.show()

##Recombination map 
#to model chromosome 14
full_map_chr14=msprime.RateMap.read_hapmap("formsprime_genetic_map_chr14.txt" )
print(full_map_chr14)
#full_map_chr1=msprime.RateMap.read_hapmap("formsprime_genetic_map_chr1.txt" )
#print(full_map_chr1)


###
## Bootstrap constant
###

##simulate and keep info
for i_simu in range(n_simu):
    sum_mig=0
    if (i_simu%5==0):
        print(i_simu)
    ts = msprime.sim_ancestry(
        ploidy=1,
        samples=sample_list_simple, 
        demography=demography_simple,
        record_migrations=True, 
        sequence_length=sequence_length,
        #recombination_rate=rec_rate,
        recombination_rate=full_map_chr14,
        record_provenance=False,
        random_seed=None)

    migrating_tracts, counts = get_migrating_tracts_count(ts)
    if (len(migrating_tracts)>0):
        length_tracts=migrating_tracts[:,1]-migrating_tracts[:,0]
        sum_mig=sum(length_tracts)
        freqs=counts/n_sample_A
        
    ##for global on the genome (remove last variable to only consider archaic segments)
    no_archaic=sequence_length-sum_mig
    length_tracts=np.append(length_tracts,no_archaic)
    freqs=np.append(freqs,0)
    if (i_simu==0):
        df_for_boot=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
        df_for_boot.insert(1,"freqs"+str(i_simu),freqs)
    else:
        df_for_boot_add=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
        df_for_boot_add.insert(1,"freqs"+str(i_simu),freqs)
        df_for_boot = pd.concat([df_for_boot,df_for_boot_add],axis=1)
            
###Bootstrap
#sample i_simu simulation, n_boot times
##collect gw_95 and above20 param
gw_95_vec=np.zeros(n_boot)
gw_above20_vec=np.zeros(n_boot)
#archaic_95_vec=np.zeros(n_boot)

for ite_n_boot in range(n_boot):
    count_per_bin=np.zeros(int(1/bin_width)+1)
    count_per_bin_nozero=np.zeros(int(1/bin_width)+1)

    gw_above_limit=0
    gw_above20_limit=0
    archaic_above_limit=0

    for ite_i_simu in range(i_simu):
        bt_i=random.randint(0,n_simu-1)
        bin_for_distri=[math.floor(float(x)) for x in df_for_boot["freqs"+str(bt_i)]/bin_width if str(x) != 'nan']
        for i_bin in range(len(bin_for_distri)):
            count_per_bin[bin_for_distri[i_bin]]+=df_for_boot["length_tracts"+str(bt_i)][i_bin]
            if (df_for_boot["freqs"+str(bt_i)][i_bin]!=0):
                count_per_bin_nozero[bin_for_distri[i_bin]]+=df_for_boot["length_tracts"+str(bt_i)][i_bin]

        #gw_above_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],genome_wide_95)/i_simu
        gw_above20_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],0.2)/i_simu
        #archaic_above_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)][:len(df_for_boot["freqs"+str(bt_i)])-1],df_for_boot["length_tracts"+str(bt_i)][:len(df_for_boot["freqs"+str(bt_i)])-1],archaic_95)/i_simu

    count_per_bin=count_per_bin/sum(count_per_bin)
    count_per_bin_nozero=count_per_bin_nozero/sum(count_per_bin_nozero)

    gw_95_vec[ite_n_boot]=percentile_value(count_per_bin, bin_width, 0.95)
    #archaic_95_vec[ite_n_boot]=percentile_value(count_per_bin_nozero, bin_width, 0.95)
    gw_above20_vec[ite_n_boot]=gw_above20_limit
    
    
all_results = pd.DataFrame([["Constant",N_Asia,0,np.mean(gw_95_vec),np.var(gw_95_vec),np.mean(gw_above20_vec),np.var(gw_above20_vec)]],columns=['Scenario','Ne', 't', 'mean_gw95','var_gw_95','mean_gw_above20','var_gw_above20'])
gw_95= pd.DataFrame(gw_95_vec, columns=["constant"])
gw_above20= pd.DataFrame(gw_above20_vec, columns=["constant"])


####
## Bottleneck simple
####
size_bt=5000
size_bt_vec=[5000,1000,100,90,70,50,30,10]

##Bottleneck
start_bt=100
start_bt_vec=[100,30]


for size_bt in size_bt_vec:
    print("size_bt=",size_bt)
    for start_bt in start_bt_vec:
        print("start_bt=",start_bt)
        end_bt=start_bt-10
        
        ##scenario
        demography_bt = msprime.Demography() #initializing an msprime "demography" object
        # This is the "trunk" population that we merge other populations into
        demography_bt.add_population(
            name="Africa",
            description="Africa",
            initial_size=N_Africa,
            initially_active=True,
        );
        demography_bt.add_population(
            name="Europe",
            description="European",
            initial_size=N_Europe,
            initially_active=True,
        )
        demography_bt.add_population(
            name="Asia",
            description="Asian",
            initial_size=N_Asia,
            initially_active=True,
        )    
        demography_bt.add_population(
            name="N",
            description="Introgressing-Neanderthal population",
            initial_size=N_Neanderthals,
            initially_active=True,
        )
        
        
        # Split Europeans and Asians
        demography_bt.add_population_split(
            time=T_europa_asia, derived=["Asia"], ancestral="Europe"
        )
        
        # Add neanderthal geneflow for one generation (pulse model)
        demography_bt.add_migration_rate_change(time = T_geneflow_Neanderthal, rate = admixtureproportion, dest = "N", source = "Europe")
        demography_bt.add_migration_rate_change(time = T_geneflow_Neanderthal+1, rate = 0, dest = "N", source = "Europe")
        
        # Add split out of Africa
        demography_bt.add_population_split(
            time=T_out_of_Africa, derived=["Europe"], ancestral='Africa'
        )
        demography_bt.add_population_parameters_change(
            time=T_Eurasian_bottleneck, population="Europe", initial_size=N_ancestral_eurasians
        )
        demography_bt.add_population_parameters_change(
            time=T_out_of_Africa, population="Africa", initial_size=N_ancestral_humans
        )
        demography_bt.add_population_split(
            time=T_Human_Neandertal, derived=["N"], ancestral="Africa" 
        )
        
        # Add Bottleneck in Asian population
        demography_bt.add_population_parameters_change(
            time=start_bt, population="Asia", initial_size=N_Asia
        )
        
        demography_bt.add_population_parameters_change(
            time=end_bt, population="Asia", initial_size=size_bt
        )
        
        demography_bt.sort_events()
        
        ##simulations
        ##simulate and keep info
        for i_simu in range(n_simu):
            sum_mig=0
            #if (i_simu%5==0):
            print(i_simu)
            ts = msprime.sim_ancestry(
                ploidy=1,
                samples=sample_list_simple, 
                demography=demography_bt,
                record_migrations=True, 
                sequence_length=sequence_length,
                #recombination_rate=rec_rate,
                recombination_rate=full_map_chr14,
                record_provenance=False,
                random_seed=None)
        
            migrating_tracts, counts = get_migrating_tracts_count(ts)
            if (len(migrating_tracts)>0):
                length_tracts=migrating_tracts[:,1]-migrating_tracts[:,0]
                sum_mig=sum(length_tracts)
                freqs=counts/n_sample_A
                
            ##for global on the genome (remove last variable to only consider archaic segments)
            no_archaic=sequence_length-sum_mig
            length_tracts=np.append(length_tracts,no_archaic)
            freqs=np.append(freqs,0)
            if (i_simu==0):
                df_for_boot=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
                df_for_boot.insert(1,"freqs"+str(i_simu),freqs)
            else:
                df_for_boot_add=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
                df_for_boot_add.insert(1,"freqs"+str(i_simu),freqs)
                df_for_boot = pd.concat([df_for_boot,df_for_boot_add],axis=1)
                    
        ##Bootstrap
        ###Bootstrap
        #sample i_simu simulation, n_boot times
        ##collect gw_95 and above20 param
        gw_95_vec=np.zeros(n_boot)
        gw_above20_vec=np.zeros(n_boot)
        archaic_95_vec=np.zeros(n_boot)
        
        for ite_n_boot in range(n_boot):
            print("n_boot=",ite_n_boot)
            count_per_bin=np.zeros(int(1/bin_width)+1)
            count_per_bin_nozero=np.zeros(int(1/bin_width)+1)
        
            gw_above_limit=0
            gw_above20_limit=0
            archaic_above_limit=0
        
            for ite_i_simu in range(i_simu):
                bt_i=random.randint(0,n_simu-1)
                bin_for_distri=[math.floor(float(x)) for x in df_for_boot["freqs"+str(bt_i)]/bin_width if str(x) != 'nan']
                print("ite=",ite_i_simu,"bt=",bt_i)
                for i_bin in range(len(bin_for_distri)):
                    count_per_bin[bin_for_distri[i_bin]]+=df_for_boot["length_tracts"+str(bt_i)][i_bin]
                    if (df_for_boot["freqs"+str(bt_i)][i_bin]!=0):
                        count_per_bin_nozero[bin_for_distri[i_bin]]+=df_for_boot["length_tracts"+str(bt_i)][i_bin]
        
                gw_above_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],genome_wide_95)/i_simu
                gw_above20_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],0.2)/i_simu
                archaic_above_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)][:len(df_for_boot["freqs"+str(bt_i)])-1],df_for_boot["length_tracts"+str(bt_i)][:len(df_for_boot["freqs"+str(bt_i)])-1],archaic_95)/i_simu
        
            count_per_bin=count_per_bin/sum(count_per_bin)
            count_per_bin_nozero=count_per_bin_nozero/sum(count_per_bin_nozero)
        
            gw_95_vec[ite_n_boot]=percentile_value(count_per_bin, bin_width, 0.95)
            archaic_95_vec[ite_n_boot]=percentile_value(count_per_bin_nozero, bin_width, 0.95)
            gw_above20_vec[ite_n_boot]=gw_above20_limit
            
          
        all_results = pd.concat([all_results, pd.DataFrame([str(size_bt)+"bt"+str(start_bt),size_bt,start_bt,np.mean(gw_95_vec),np.var(gw_95_vec),np.mean(gw_above20_vec),np.var(gw_above20_vec)])], ignore_index=True)
        gw_95.insert(1,str(size_bt)+"bt"+str(start_bt),gw_95_vec)
        gw_above20.insert(1,str(size_bt)+"bt"+str(start_bt),gw_above20_vec)


####
## Double bottleneck
####
n_sample_A1=250
n_sample_A2=250
samples_A1 = msprime.SampleSet(n_sample_A1,population="Asia1")
samples_A2 = msprime.SampleSet(n_sample_A2,population="Asia2")
samples_N = msprime.SampleSet(1,population="N",time=T_geneflow_Neanderthal)
sample_list_b2t = [samples_A1,samples_A2,samples_N]

##TWo bottlenecks
size_bt1=100
size_bt2=100

##Bottleneck
start_b2t_vec=[100,30]
for start_b2t in start_b2t_vec:
    end_bt=start_b2t-10
    print("start:",start_b2t)
    
    demography_b2t = msprime.Demography() #initializing an msprime "demography" object
    # This is the "trunk" population that we merge other populations into
    demography_b2t.add_population(
        name="Africa",
        description="Africa",
        initial_size=N_Africa,
        initially_active=True,
    );
    demography_b2t.add_population(
        name="Europe",
        description="European",
        initial_size=N_Europe,
        initially_active=True,
    )
    demography_b2t.add_population(
        name="Asia1",
        description="Asian_bt1",
        initial_size=N_Asia/2,
        initially_active=True,
    )    
    demography_b2t.add_population(
        name="Asia2",
        description="Asian_bt2",
        initial_size=N_Asia/2,
        initially_active=True,
    )  
    demography_b2t.add_population(
        name="N",
        description="Introgressing-Neanderthal population",
        initial_size=N_Neanderthals,
        initially_active=True,
    )
    
    
    # Split Europeans and Asians
    demography_b2t.add_population_split(
        time=T_europa_asia, derived=["Asia1"], ancestral="Europe"
    )
    
    
    # Add neanderthal geneflow for one generation (pulse model)
    demography_b2t.add_migration_rate_change(time = T_geneflow_Neanderthal, rate = admixtureproportion, dest = "N", source = "Europe")
    demography_b2t.add_migration_rate_change(time = T_geneflow_Neanderthal+1, rate = 0, dest = "N", source = "Europe")
    
    # Add split out of Africa
    demography_b2t.add_population_split(
        time=T_out_of_Africa, derived=["Europe"], ancestral='Africa'
    )
    demography_b2t.add_population_parameters_change(
        time=T_Eurasian_bottleneck, population="Europe", initial_size=N_ancestral_eurasians
    )
    demography_b2t.add_population_parameters_change(
        time=T_out_of_Africa, population="Africa", initial_size=N_ancestral_humans
    )
    demography_b2t.add_population_split(
        time=T_Human_Neandertal, derived=["N"], ancestral="Africa" 
    )
    
    # Add split in Asian population
    demography_b2t.add_population_split(
        time=start_b2t, derived=["Asia2"], ancestral="Asia1"
    )
    ## Bottleneck
    demography_b2t.add_population_parameters_change(
        time=start_b2t, population="Asia1", initial_size=N_Asia
    )
    
    demography_b2t.add_population_parameters_change(
        time=end_bt, population="Asia1", initial_size=size_bt1
    )
    
    demography_b2t.add_population_parameters_change(
        time=end_bt, population="Asia2", initial_size=size_bt2
    )
    
    demography_b2t.sort_events()
    
    ###simulation
    ##simulate and keep info
    for i_simu in range(n_simu):
        sum_mig=0
        #if (i_simu%5==0):
        print(i_simu)
        ts = msprime.sim_ancestry(
            ploidy=1,
            samples=sample_list_b2t, 
            demography=demography_b2t,
            record_migrations=True, 
            sequence_length=sequence_length,
            #recombination_rate=rec_rate,
            recombination_rate=full_map_chr14,
            record_provenance=False,
            random_seed=None)
    
        migrating_tracts, counts = get_migrating_tracts_count(ts)
        if (len(migrating_tracts)>0):
            length_tracts=migrating_tracts[:,1]-migrating_tracts[:,0]
            sum_mig=sum(length_tracts)
            freqs=counts/n_sample_A
            
        ##for global on the genome (remove last variable to only consider archaic segments)
        no_archaic=sequence_length-sum_mig
        length_tracts=np.append(length_tracts,no_archaic)
        freqs=np.append(freqs,0)
        if (i_simu==0):
            df_for_boot=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
            df_for_boot.insert(1,"freqs"+str(i_simu),freqs)
        else:
            df_for_boot_add=pd.DataFrame(length_tracts,columns =["length_tracts"+str(i_simu)])
            df_for_boot_add.insert(1,"freqs"+str(i_simu),freqs)
            df_for_boot = pd.concat([df_for_boot,df_for_boot_add],axis=1)
                
     ##bootstrap
     #sample i_simu simulation, n_boot times
    ##collect gw_95 and above20 param
    gw_95_vec=np.zeros(n_boot)
    gw_above20_vec=np.zeros(n_boot)
    archaic_95_vec=np.zeros(n_boot)
    
    for ite_n_boot in range(n_boot):
        print("n_boot=",ite_n_boot)
        count_per_bin=np.zeros(int(1/bin_width)+1)
        count_per_bin_nozero=np.zeros(int(1/bin_width)+1)
    
        gw_above_limit=0
        gw_above20_limit=0
        archaic_above_limit=0
    
        for ite_i_simu in range(i_simu):
            bt_i=random.randint(0,n_simu-1)
            bin_for_distri=[math.floor(float(x)) for x in df_for_boot["freqs"+str(bt_i)]/bin_width if str(x) != 'nan']
            print("ite=",ite_i_simu,"bt=",bt_i)
            for i_bin in range(len(bin_for_distri)):
                count_per_bin[bin_for_distri[i_bin]]+=df_for_boot["length_tracts"+str(bt_i)][i_bin]
                if (df_for_boot["freqs"+str(bt_i)][i_bin]!=0):
                    count_per_bin_nozero[bin_for_distri[i_bin]]+=df_for_boot["length_tracts"+str(bt_i)][i_bin]
    
            gw_above_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],genome_wide_95)/i_simu
            gw_above20_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)],df_for_boot["length_tracts"+str(bt_i)],0.2)/i_simu
            archaic_above_limit+=inarchaic_percentage_above_limit(df_for_boot["freqs"+str(bt_i)][:len(df_for_boot["freqs"+str(bt_i)])-1],df_for_boot["length_tracts"+str(bt_i)][:len(df_for_boot["freqs"+str(bt_i)])-1],archaic_95)/i_simu
    
        count_per_bin=count_per_bin/sum(count_per_bin)
        count_per_bin_nozero=count_per_bin_nozero/sum(count_per_bin_nozero)
    
        gw_95_vec[ite_n_boot]=percentile_value(count_per_bin, bin_width, 0.95)
        archaic_95_vec[ite_n_boot]=percentile_value(count_per_bin_nozero, bin_width, 0.95)
        gw_above20_vec[ite_n_boot]=gw_above20_limit  
    
    all_results = pd.concat([all_results, pd.DataFrame([str(size_bt1)+"b2t"+str(start_b2t),"2*"+str(size_bt1),start_bt,np.mean(gw_95_vec),np.var(gw_95_vec),np.mean(gw_above20_vec),np.var(gw_above20_vec)])], ignore_index=True)
    gw_95.insert(1,str(size_bt1)+"b2t"+str(start_b2t),gw_95_vec)
    gw_above20.insert(1,str(size_bt1)+"b2t"+str(start_b2t),gw_above20_vec)
 

####
##write FINAL
####
all_results.to_csv(sum_results+str(n_simu)+"_"+str(n_boot)+".txt", index=None, sep=' ', mode='w')
gw_95.to_csv(outfile_gw95+str(n_simu)+"_"+str(n_boot)+".txt", index=None, sep=' ', mode='w')
gw_above20.to_csv(outfile_gw20+str(n_simu)+"_"+str(n_boot)+".txt", index=None, sep=' ', mode='w')
