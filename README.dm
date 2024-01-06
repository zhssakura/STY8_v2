######## Basic information of the symbiotic microorganisms associated with the sponge Stylissa sp. ######## 

## Taxonomic classification of the eight microbial species, their abbreviations (i.e., OTU + a figure), and the functional population they represent:
# OTU1: Spongiicolonus versatilis     -- SOB
# OTU2: Parasynechococcus marinus     -- photosymbiont
# OTU3: Deltaporibacter stylissae     -- Denitrifier
# OTU4: Spongiihabitans thiooxidans   -- SOB
# OTU5: Gammaporibacter thiotrophicus -- SOB
# OTU6: Nitrospongiibacter stylissae  -- NOB
# OTU7: Oxydemutator thiolithooxidans -- SOB
# OTU8: Cenoporarchaeum stylissae     -- AOA


####### Command lines of using BacArena to conduct network simulations #######

# 1. To simulate an in silico growth test for a single microbial species, here is an example using the model of NOB:
~/Gh_Sflabelliformis_8_MAGs/10_A_genome-scale_metabolic_network_for_S.flabelliformis_publication/Cmds_for_running_R_scripts/Cmds_Single-species_simulations/Cmds_for_Rscript_BacArena_SingleSpecies__simuDiet_optparse_pFBA_ka.sh


# 2. To simulate a network of the two-species microbial consortium, which includes both the NOB and AOA models, execute the following command line:
~/Gh_Sflabelliformis_8_MAGs/10_A_genome-scale_metabolic_network_for_S.flabelliformis_publication/Cmds_for_running_R_scripts/Cmds_2-species_network/Cmds_for_Rscript_BacArena_otu6n8__simuDiet_optparse_FBA_ka.sh


# 3. To simulate a network of the three-species microbial consortium, which includes the NOB, AOA and photosymbiont models, run the following command line:
~/Gh_Sflabelliformis_8_MAGs/10_A_genome-scale_metabolic_network_for_S.flabelliformis_publication/Cmds_for_running_R_scripts/Cmds_3-species_network/Cmds_for_Rscript_BacArena_otu2n6n8__simuDiet_optparse_FBA_ka.sh


# 4. To simulate a network of the eight-species microbiome, here is an example using hypotaurine as the sole organic carbon at a concentration of 0.2 mM:
~/Gh_Sflabelliformis_8_MAGs/10_A_genome-scale_metabolic_network_for_S.flabelliformis_publication/Cmds_for_running_R_scripts/Cmds_8-species_network/Cmds_for_Rscript_BacArena_8species__simuDiet_optparse_pFBA_ka.sh


######## Citation ####### 

Shan Zhang, Weizhi Song, Geogios Marinos, Silvio Waschina, Johannes Zimmermann, Christoph Kaleta, Torsten Thomas. Insights into the metabolic interactions in sponge microbiome: a genome-scale metabolic network for the microbiome of the marine sponge Stylissa sp. (2024)