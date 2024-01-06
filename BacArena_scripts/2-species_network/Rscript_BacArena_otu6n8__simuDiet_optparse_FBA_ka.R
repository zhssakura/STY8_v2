################################################################# katana ####################################################################################################
############ 
library(lattice)
library(ReacTran)
library(rootSolve)
library(deSolve)
library(shape)
library(sybil)
library(Matrix)
library(BacArena)
library(parallel)
library(optparse)

######################## argument ########################

option_list = list(
  optparse::make_option(c("-n", "--nh3_con"),   type="double",       default=0.5,       help="NH3 concentration, defalt: 0.5"),
  optparse::make_option(c("-r", "--rm_rate"),   type="double",       default=0,         help="removeM value, defalt: 0"),
  optparse::make_option(c("-a", "--arena_mn"),  type="double",       default=20,        help="integer indicating the length of an arena, defalt: 20"),
  optparse::make_option(c("-d", "--death_r"),   type="double",       default=0,         help="A percentage of biomass reduce due to the nutrient limitation, defalt: 0"),
  optparse::make_option(c("-s", "--tstep"),     type="double",       default=1,         help="Time step, defalt: 1h per iteration"),
  optparse::make_option(c("-t", "--iter"),      type="double",       default=200,       help="Number of iteration, defalt:200"));

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

NH3_con =  opt$nh3_con
rm_rate =  opt$rm_rate
grid_no =  opt$arena_mn
death_r =  opt$death_r
time_step =  opt$tstep
iter_no   =  opt$iter
#####################################################################################################################################################################

############ Katana ############ 
diet <- readRDS('~/Gh_Sflabelliformis_8_MAGs/10_A_genome-scale_metabolic_network_for_S.flabelliformis_publication/Growth_media/2-species_network/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified_9999_NH3only.RDS')

getwd <- getwd()
getwd

#####################################################################################################################################################################
dir_my_model = '~/Gh_Sflabelliformis_8_MAGs/10_A_genome-scale_metabolic_network_for_S.flabelliformis_publication/GEMs/'

# the NOB:
STYtaxon_6='/STY_Merged_OTU06.RDS'

# the AOA - NO forming deactivate!!:
STYtaxon_8='/STY_Merged_OTU08_MgfM.RDS'

model6 <- readRDS(paste(dir_my_model,STYtaxon_6,sep = "")) #d__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__UBA8639;g__bin75;s__
model8 <- readRDS(paste(dir_my_model,STYtaxon_8,sep = "")) #d__Archaea;p__Crenarchaeota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Cenarchaeum;s__
#################################################################################

combination_prefix = paste(NH3_con, rm_rate, sep = '_')

print(combination_prefix)


##### add main code here
sink(file = paste(getwd,'/NH3_',NH3_con,'_mM__removeM_',rm_rate,'.txt',sep = '')) # get the sink function started

print("============================ OTU6: the NOB ========================================")
print("============================ OTU8: the modified AOA (NO-forming deactivated!!!!) ========================================")
print("============================ lyse = F,deathrate=0.0,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=FALSE ========================================")
print("============================ OTU6:OTU8(3.690:28.570); space = 400; [NH3] = see in the file name; nitrite=0; removeM=see in the file name; iter = 10 ========================================")
print("============================ Continous culture (with replenishment of NH3): yes ========================================")
print("============================ Batch culture: no ========================================")
print("============================ Diet: mineral_sw_210514_OTU08_3_BacArena_realistic_simplified_9999_NH3_n_light.RDS ========================================")

bacterium6 <- BacArena::Bac(model6,lyse = F,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=FALSE)
bacterium8 <- BacArena::Bac(model8,lyse = F,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=FALSE)

arena <- BacArena::Arena(n=grid_no,m=grid_no)
inocc <- arena@n * arena@m * 0.03099814 # 0.1% inocculum
arena <- BacArena::addOrg(object = arena, specI = bacterium6, amount = inocc*3.690) # amount can be added by fraction.
arena <- BacArena::addOrg(object = arena, specI = bacterium8, amount = inocc*28.57) # amount can be added by fraction.

arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange, smax = diet$Input_mM, unit = "mM", add = T)
arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[35], smax = NH3_con, unit = "mM", add = F) #[1] EX_cpd00013_e0	N as ammonia
arena@tstep <- time_step #hour
BacArena::chemotaxis(bacterium6, arena, 1, chemo='EX_cpd00075_e0', arena@occupyM)#nitrite
BacArena::chemotaxis(bacterium8, arena, 1, chemo='EX_cpd00013_e0', arena@occupyM)#ammonia


######## iteration 1:
simulation <- BacArena::simEnv(object = arena, time = 1, continue = T) #FBA; iter 1
# Get time points:
tp0=simulation@simlist[[1]]
tp1=simulation@simlist[[2]]

# 1) to calculate number of individuals: from simlist.
tp0_spe1_cell_no = length(which(tp0[["type"]] == '1')) # number of individuals of 1st model (OTU6) at time point 0.
tp0_spe2_cell_no = length(which(tp0[["type"]] == '2')) # number of individuals of 2nd model (OTU8) at time point 0.
tp0_all_cell_no = length(tp0[["type"]]) # number of individuals in total at time point 0.

tp1_spe1_cell_no = length(which(tp1[["type"]] == '1')) # number of individuals of 1st model (OTU6) at time point 1.
tp1_spe2_cell_no = length(which(tp1[["type"]] == '2')) # number of individuals of 2nd model (OTU8) at time point 1.
tp1_all_cell_no = length(tp1[["type"]]) # number of individuals in total at time point 1.

# 2) to calculate the biomass: from simlist.
# OTU6:
tp0_spe1_ord<-which(tp0[["type"]] == '1') # Find the 1st model (OTU6) in the orgdat at time point 0.
tp1_spe1_ord<-which(tp1[["type"]] == '1') # Find the 1st model (OTU6) in the orgdat at time point 1.
# OTU8:
tp0_spe2_ord<-which(tp0[["type"]] == '2') # Find the 1st model (OTU8) in the orgdat at time point 0.
tp1_spe2_ord<-which(tp1[["type"]] == '2') # Find the 1st model (OTU8) in the orgdat at time point 1.

current_biomass_OTU6_tp0 =0 # assign a variable with integer of 0.
for (i in tp0_spe1_ord) {
  bb<-tp0[["biomass"]][i] # find the biomass of OTU6 at time point 0
  current_biomass_OTU6_tp0= current_biomass_OTU6_tp0+bb # calculate the total amount of bb
}

current_biomass_OTU8_tp0 =0 # assign a variable with integer of 0.
for (i in tp0_spe2_ord) {
  bb<-tp0[["biomass"]][i] # find the biomass of OTU6 at time point 0
  current_biomass_OTU8_tp0= current_biomass_OTU8_tp0+bb # calculate the total amount of bb
}

print(paste('↑↑↑Total biomass of simlist befor the 1st iteration (2)    = ', round(current_biomass_OTU6_tp0,5)+round(current_biomass_OTU8_tp0,5),', OTU6=',round(current_biomass_OTU6_tp0,5),', OTU8=',round(current_biomass_OTU8_tp0,5), sep = ''))# equal to arena 
print(paste('↑↑↑Total cells   of simlist befor the 1st iteration        = ', tp0_all_cell_no,', OTU6=',tp0_spe1_cell_no,', OTU8=',tp0_spe2_cell_no, sep = ''))# equal to arena 

current_biomass_OTU6_tp1 =0 # assign a variable with integer of 0.
for (i in tp1_spe1_ord) {
  bb<-tp1[["biomass"]][i] # find the biomass of OTU6 at time point 0
  current_biomass_OTU6_tp1= current_biomass_OTU6_tp1+bb # calculate the total amount of bb
}

current_biomass_OTU8_tp1 =0 # assign a variable with integer of 0.
for (i in tp1_spe2_ord) {
  bb<-tp1[["biomass"]][i] # find the biomass of OTU6 at time point 0
  current_biomass_OTU8_tp1= current_biomass_OTU8_tp1+bb # calculate the total amount of bb
}

print(paste('↑↑↑Total biomass of simlist after the 1st iteration (2)    = ', round(current_biomass_OTU6_tp1,5)+round(current_biomass_OTU8_tp1,5),', OTU6=',round(current_biomass_OTU6_tp1,5),', OTU8=',round(current_biomass_OTU8_tp1,5), sep = ''))
print(paste('↑↑↑Total cells   of simlist after the 1st iteration        = ', tp1_all_cell_no,', OTU6=',tp1_spe1_cell_no,', OTU8=',tp1_spe2_cell_no, sep = ''))# equal to arena 
print('                          ')
print('                          ')

######## iteration of 2 and on:
for (i in 2:iter_no) {
  arena2 <- BacArena::getArena(simulation,1) # Add the arena in the simulist of 1st iteration.
  arena2@orgdat <- arena2@orgdat[-sample(nrow(arena2@orgdat), round(nrow(arena2@orgdat) * rm_rate/100)), ] #remove ~10% of individuals randomly
  
  arena2 <- BacArena::addSubs(object = arena2, mediac = diet$Exchange, smax = diet$Input_mM, unit = "mM", add = T)
  arena2 <- BacArena::addSubs(object = arena2, mediac = diet$Exchange[35], smax = NH3_con, unit = "mM", add = F) #replenishement of nutrient:[1] EX_cpd00013_e0	N as ammonia
  simulation <- BacArena::simEnv(object = arena2, time = 1, continue = T) #FBA for 1 iter
  
  # An old way to print biomass:       
  bimass_plotGrothCurve = plotGrowthCurve(simulation,use_biomass = T,ret_data = T)
  print(paste('↑↑↑Total biomass of plotGrowthCurve befor the ',i,'th iteration = ', round(bimass_plotGrothCurve[1,3]+bimass_plotGrothCurve[2,3]+bimass_plotGrothCurve[3,3],5),'. OTU06:',round(bimass_plotGrothCurve[1,3],5),', OTU08:',round(bimass_plotGrothCurve[2,3],5),', OTU02:',round(bimass_plotGrothCurve[3,3],5), sep = ''))
  print(paste('↑↑↑Total biomass of plotGrowthCurve after the ',i,'th iteration = ', round(bimass_plotGrothCurve[4,3]+bimass_plotGrothCurve[5,3]+bimass_plotGrothCurve[6,3],5),'. OTU06:',round(bimass_plotGrothCurve[4,3],5),', OTU08:',round(bimass_plotGrothCurve[5,3],5),', OTU02:',round(bimass_plotGrothCurve[6,3],5), sep = ''))
  print('                          ')
  print('                          ')
  
  name <- paste("biomass_bf_tp", i, sep = "")
  assign(name, round(bimass_plotGrothCurve[1,3]+bimass_plotGrothCurve[2,3],5))

  name <- paste("biomass_aft_tp", i, sep = "")
  assign(name, round(bimass_plotGrothCurve[3,3]+bimass_plotGrothCurve[4,3],5))

  name <- paste("biomass_OTU6_bf_tp", i, sep = "")
  assign(name, round(bimass_plotGrothCurve[1,3],5))

  name <- paste("biomass_OTU6_aft_tp", i, sep = "")
  assign(name, round(bimass_plotGrothCurve[3,3],5))

  name <- paste("biomass_OTU8_bf_tp", i, sep = "")
  assign(name, round(bimass_plotGrothCurve[2,3],5))

  name <- paste("biomass_OTU8_aft_tp", i, sep = "")
  assign(name, round(bimass_plotGrothCurve[4,3],5))


  # Get time points:
  tp0=simulation@simlist[[1]]
  tp1=simulation@simlist[[2]]
  
  # # 1) to calculate number of individuals: from simlist.
  name <- paste("Total_cell_no_bf_tp", i, sep = "")
  assign(name, length(tp0[["type"]]))

  name <- paste("Total_cell_no_af_tp", i, sep = "")
  assign(name, length(tp1[["type"]]))

  name <- paste("spe1_cell_no_bf_tp", i, sep = "")
  assign(name, length(which(tp0[["type"]] == '1')))

  name <- paste("spe1_cell_no_af_tp", i, sep = "")
  assign(name, length(which(tp1[["type"]] == '1')))

  name <- paste("spe2_cell_no_bf_tp", i, sep = "")
  assign(name, length(which(tp0[["type"]] == '2')))

  name <- paste("spe2_cell_no_af_tp", i, sep = "")
  assign(name, length(which(tp1[["type"]] == '2')))

  print(paste('↑↑↑Total cells of simlist befor the ',i,'th iteration = ', length(tp0[["type"]]),'. OTU06:',length(which(tp0[["type"]] == '1')),', OTU08:',length(which(tp0[["type"]] == '2')), sep = ''))
  print(paste('↑↑↑Total cells of simlist after the ',i,'th iteration = ', length(tp1[["type"]]),'. OTU06:',length(which(tp1[["type"]] == '1')),', OTU08:',length(which(tp1[["type"]] == '2')), sep = ''))
  print('                          ')
  print('                          ')
}

sink()

########################### get three summary files including info of biomass and cell number ########################### 
# A) export a file including variance of Total biomass + cell of all species:
sink(file = paste(getwd,'/NH3_',NH3_con,'_mM__removeM_',rm_rate,'_','_summary.txt',sep = '')) # get the sink function started
print(paste('NH3_con','rm_rate','iter_no',
            'Total biomass of simlist befor the 1st iter','Total biomass of simlist after the 1st iter',
            'biomass_bf_tp2','biomass_aft_tp2','biomass_bf_tp3','biomass_aft_tp3','biomass_bf_tp4','biomass_aft_tp4','biomass_bf_tp5','biomass_aft_tp5','biomass_bf_tp6','biomass_aft_tp6','biomass_bf_tp7','biomass_aft_tp7','biomass_bf_tp8','biomass_aft_tp8','biomass_bf_tp9','biomass_aft_tp9','biomass_bf_tp10','biomass_aft_tp10',
            'biomass_bf_tp41,biomass_aft_tp41,biomass_bf_tp42,biomass_aft_tp42,biomass_bf_tp43,biomass_aft_tp43,biomass_bf_tp44,biomass_aft_tp44,biomass_bf_tp45,biomass_aft_tp45,biomass_bf_tp46,biomass_aft_tp46,biomass_bf_tp47,biomass_aft_tp47,biomass_bf_tp48,biomass_aft_tp48,biomass_bf_tp49,biomass_aft_tp49,biomass_bf_tp50,biomass_aft_tp50',
            
            'tp0_all_cell_no','tp1_all_cell_no','Total_cell_no_bf_tp2','Total_cell_no_af_tp2','Total_cell_no_bf_tp3','Total_cell_no_af_tp3','Total_cell_no_bf_tp4','Total_cell_no_af_tp4','Total_cell_no_bf_tp5','Total_cell_no_af_tp5','Total_cell_no_bf_tp6','Total_cell_no_af_tp6','Total_cell_no_bf_tp7','Total_cell_no_af_tp7','Total_cell_no_bf_tp8','Total_cell_no_af_tp8','Total_cell_no_bf_tp9','Total_cell_no_af_tp9','Total_cell_no_bf_tp10','Total_cell_no_af_tp10',
            'ratio of sp2/sp1 bf tp1','ratio of sp2/sp1 af tp1','ratio of sp2/sp1 bf tp2','ratio of sp2/sp1 af tp2','ratio of sp2/sp1 bf tp3','ratio of sp2/sp1 af tp3','ratio of sp2/sp1 bf tp4','ratio of sp2/sp1 af tp4','ratio of sp2/sp1 bf tp5','ratio of sp2/sp1 af tp5','ratio of sp2/sp1 bf tp6','ratio of sp2/sp1 af tp6','ratio of sp2/sp1 bf tp7','ratio of sp2/sp1 af tp7','ratio of sp2/sp1 bf tp8','ratio of sp2/sp1 af tp8','ratio of sp2/sp1 bf tp9','ratio of sp2/sp1 af tp9','ratio of sp2/sp1 bf tp10','ratio of sp2/sp1 af tp10',
            
            'tp0_spe1_cell_no,tp0_spe2_cell_no','tp1_spe1_cell_no,tp1_spe2_cell_no',
            'spe1_cell_no_bf_tp41,spe2_cell_no_bf_tp41,spe1_cell_no_af_tp41,spe2_cell_no_af_tp41',
            'spe1_cell_no_bf_tp42,spe2_cell_no_bf_tp42,spe1_cell_no_af_tp42,spe2_cell_no_af_tp42',
            'spe1_cell_no_bf_tp43,spe2_cell_no_bf_tp43,spe1_cell_no_af_tp43,spe2_cell_no_af_tp43',
            'spe1_cell_no_bf_tp44,spe2_cell_no_bf_tp44,spe1_cell_no_af_tp44,spe2_cell_no_af_tp44',
            'spe1_cell_no_bf_tp45,spe2_cell_no_bf_tp45,spe1_cell_no_af_tp45,spe2_cell_no_af_tp45',
            'spe1_cell_no_bf_tp46,spe2_cell_no_bf_tp46,spe1_cell_no_af_tp46,spe2_cell_no_af_tp46',
            'spe1_cell_no_bf_tp47,spe2_cell_no_bf_tp47,spe1_cell_no_af_tp47,spe2_cell_no_af_tp47',
            'spe1_cell_no_bf_tp48,spe2_cell_no_bf_tp48,spe1_cell_no_af_tp48,spe2_cell_no_af_tp48',
            'spe1_cell_no_bf_tp49,spe2_cell_no_bf_tp49,spe1_cell_no_af_tp49,spe2_cell_no_af_tp49',
            'spe1_cell_no_bf_tp50,spe2_cell_no_bf_tp50,spe1_cell_no_af_tp50,spe2_cell_no_af_tp50',
            'ratio of sp2/sp1 bf tp1','ratio of sp2/sp1 af tp1','ratio of sp2/sp1 af tp41','ratio of sp2/sp1 af tp42','ratio of sp2/sp1 af tp43','ratio of sp2/sp1 af tp44','ratio of sp2/sp1 af tp45','ratio of sp2/sp1 af tp46','ratio of sp2/sp1 af tp47','ratio of sp2/sp1 af tp48','ratio of sp2/sp1 af tp49','ratio of sp2/sp1 af tp50',
            sep=','))

print(paste(NH3_con,rm_rate,iter_no,
            round(current_biomass_OTU6_tp0,5)+round(current_biomass_OTU8_tp0,5),round(current_biomass_OTU6_tp1,5)+round(current_biomass_OTU8_tp1,5),#biomass of (bf/af) the 1st iter.
            biomass_bf_tp2,biomass_aft_tp2,biomass_bf_tp3,biomass_aft_tp3,biomass_bf_tp4,biomass_aft_tp4,biomass_bf_tp5,biomass_aft_tp5,biomass_bf_tp6,biomass_aft_tp6,biomass_bf_tp7,biomass_aft_tp7,biomass_bf_tp8,biomass_aft_tp8,biomass_bf_tp9,biomass_aft_tp9,biomass_bf_tp10,biomass_aft_tp10,
            biomass_bf_tp41,biomass_aft_tp41,biomass_bf_tp42,biomass_aft_tp42,biomass_bf_tp43,biomass_aft_tp43,biomass_bf_tp44,biomass_aft_tp44,biomass_bf_tp45,biomass_aft_tp45,biomass_bf_tp46,biomass_aft_tp46,biomass_bf_tp47,biomass_aft_tp47,biomass_bf_tp48,biomass_aft_tp48,biomass_bf_tp49,biomass_aft_tp49,biomass_bf_tp50,biomass_aft_tp50,
            
            tp0_all_cell_no,tp1_all_cell_no,Total_cell_no_bf_tp2,Total_cell_no_af_tp2,Total_cell_no_bf_tp3,Total_cell_no_af_tp3,Total_cell_no_bf_tp4,Total_cell_no_af_tp4,Total_cell_no_bf_tp5,Total_cell_no_af_tp5,Total_cell_no_bf_tp6,Total_cell_no_af_tp6,Total_cell_no_bf_tp7,Total_cell_no_af_tp7,Total_cell_no_bf_tp8,Total_cell_no_af_tp8,Total_cell_no_bf_tp9,Total_cell_no_af_tp9,Total_cell_no_bf_tp10,Total_cell_no_af_tp10,
            tp0_spe2_cell_no/tp0_spe1_cell_no,tp1_spe2_cell_no/tp1_spe1_cell_no,spe2_cell_no_bf_tp2/spe1_cell_no_bf_tp2,spe2_cell_no_af_tp2/spe1_cell_no_af_tp2,spe2_cell_no_bf_tp3/spe1_cell_no_bf_tp3,spe2_cell_no_af_tp3/spe1_cell_no_af_tp3,spe2_cell_no_bf_tp4/spe1_cell_no_bf_tp4,spe2_cell_no_af_tp4/spe1_cell_no_af_tp4,spe2_cell_no_bf_tp5/spe1_cell_no_bf_tp5,spe2_cell_no_af_tp5/spe1_cell_no_af_tp5,spe2_cell_no_bf_tp6/spe1_cell_no_bf_tp6,spe2_cell_no_af_tp6/spe1_cell_no_af_tp6,spe2_cell_no_bf_tp7/spe1_cell_no_bf_tp7,spe2_cell_no_af_tp7/spe1_cell_no_af_tp7,spe2_cell_no_bf_tp8/spe1_cell_no_bf_tp8,spe2_cell_no_af_tp8/spe1_cell_no_af_tp8,spe2_cell_no_bf_tp9/spe1_cell_no_bf_tp9,spe2_cell_no_af_tp9/spe1_cell_no_af_tp9,spe2_cell_no_bf_tp10/spe1_cell_no_bf_tp10,spe2_cell_no_af_tp10/spe1_cell_no_af_tp10,
            
            tp0_spe1_cell_no,tp0_spe2_cell_no,tp1_spe1_cell_no,tp1_spe2_cell_no,
            spe1_cell_no_bf_tp41,spe2_cell_no_bf_tp41,spe1_cell_no_af_tp41,spe2_cell_no_af_tp41,
            spe1_cell_no_bf_tp42,spe2_cell_no_bf_tp42,spe1_cell_no_af_tp42,spe2_cell_no_af_tp42,
            spe1_cell_no_bf_tp43,spe2_cell_no_bf_tp43,spe1_cell_no_af_tp43,spe2_cell_no_af_tp43,
            spe1_cell_no_bf_tp44,spe2_cell_no_bf_tp44,spe1_cell_no_af_tp44,spe2_cell_no_af_tp44,
            spe1_cell_no_bf_tp45,spe2_cell_no_bf_tp45,spe1_cell_no_af_tp45,spe2_cell_no_af_tp45,
            spe1_cell_no_bf_tp46,spe2_cell_no_bf_tp46,spe1_cell_no_af_tp46,spe2_cell_no_af_tp46,
            spe1_cell_no_bf_tp47,spe2_cell_no_bf_tp47,spe1_cell_no_af_tp47,spe2_cell_no_af_tp47,
            spe1_cell_no_bf_tp48,spe2_cell_no_bf_tp48,spe1_cell_no_af_tp48,spe2_cell_no_af_tp48,
            spe1_cell_no_bf_tp49,spe2_cell_no_bf_tp49,spe1_cell_no_af_tp49,spe2_cell_no_af_tp49,
            spe1_cell_no_bf_tp50,spe2_cell_no_bf_tp50,spe1_cell_no_af_tp50,spe2_cell_no_af_tp50,
            tp0_spe2_cell_no/tp0_spe1_cell_no,tp1_spe2_cell_no/tp1_spe1_cell_no,spe2_cell_no_af_tp41/spe1_cell_no_af_tp41,spe2_cell_no_af_tp42/spe1_cell_no_af_tp42,spe2_cell_no_af_tp43/spe1_cell_no_af_tp43,spe2_cell_no_af_tp44/spe1_cell_no_af_tp44,spe2_cell_no_af_tp45/spe1_cell_no_af_tp45,spe2_cell_no_af_tp46/spe1_cell_no_af_tp46,spe2_cell_no_af_tp47/spe1_cell_no_af_tp47,spe2_cell_no_af_tp48/spe1_cell_no_af_tp48,spe2_cell_no_af_tp49/spe1_cell_no_af_tp49,spe2_cell_no_af_tp50/spe1_cell_no_af_tp50,
            sep=','))
sink()
