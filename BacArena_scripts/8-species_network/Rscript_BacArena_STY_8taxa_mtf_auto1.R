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
library(cplexAPI)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
######################## argument ########################
option_list = list(
  optparse::make_option(c("-o", "--otu_name"),                 type="character",      default='8taxa',           help="OTU_name, defalt: 8taxa"),
  optparse::make_option(c("-D", "--infile_BacarenaRDS_diet"),  type="character",      default="",                help="Diet RDS file for bacarena with full path"),
  optparse::make_option(c("-k", "--keywd"),                    type="character",      default="sw_v20210716",    help="Keyword (Nutirent/date) used in the folder name of output files, defalt: sw_v20210716"),
  optparse::make_option(c("-r", "--rm_rate"),                  type="double",         default=0,                 help="RemoveM value, defalt: 0"),
  optparse::make_option(c("-p", "--tstep"),                    type="double",         default=1,                 help="Time step, defalt: 1h per iteration"),
  optparse::make_option(c("-a", "--arena_mn"),                 type="double",         default=100,               help="An integer indicating the length of an arena, defalt: 20"),
  optparse::make_option(c("-d", "--death_r"),                  type="double",         default=0,                 help="A percentage of biomass reduce due to the nutrient limitation, defalt: 0"),
  optparse::make_option(c("-i", "--inocc_no"),                 type="double",         default=1,                 help="inocculum, defalt: 1"),
  optparse::make_option(c("-c", "--cl_no"),                    type="double",         default=1,                 help="Number of replicates, defalt: 1"),
  optparse::make_option(c("-s", "--setAllExInf_value"),        type="character",      default='FALSE',           help="setAllExInf, defalt: FALSE"),
  optparse::make_option(c("-n", "--auto_num"),                 type="double",         default=1,                 help="Surfix in the simulation RDS file name."),
  optparse::make_option(c("-N", "--NH3_con"),                  type="double",         default=0.5,               help="NH3 concentration, defalt: 0.5"),
  optparse::make_option(c("-H", "--HT_con"),                   type="double",         default=0.5,               help="HT concentration, defalt: 0.5"),
  optparse::make_option(c("-P", "--photon_con"),               type="double",         default=1000,              help="Number indicating the photon concentration, defalt: 1000 mM"),
  optparse::make_option(c("-G", "--ddca_con"),                 type="double",         default=5,                 help="Number indicating the ddca concentration, defalt: 5 mM"),
  optparse::make_option(c("-L", "--lxlycm"),                   type="double",         default=0.01,              help="Number giving the horizontal grid size in cm, defalt: 0.01 cm. For 10000 cells in a 0.01^3 cm3, cell concentrition is 1*10^10 cell/mL, each grid unit is 1 um^2"),
  optparse::make_option(c("-C", "--cconc"),                    type="double",         default=1E+10,             help="Number of cell density in cell/mL, defalt: 1E+10"),
  optparse::make_option(c("-m", "--diameter_um"),              type="double",         default=0.82,              help="Number of cell diameter in um, defalt: 0.82"),
  optparse::make_option(c("-t", "--iter"),                     type="double",         default=1,                 help="Number of iteration, defalt:1"));

opt_parser = optparse::OptionParser(option_list=option_list, add_help_option=FALSE);
opt = optparse::parse_args(opt_parser);

otu_name        = opt$otu_name
infile_diet     = opt$infile_BacarenaRDS_diet
keywd           = opt$keywd
arena_mn        = opt$arena_mn
death_r         = opt$death_r
Inocc_no        = opt$inocc_no
Cl_no           = opt$cl_no
setAllExInf_value = opt$setAllExInf_value
iter_no         = opt$iter
rm_rate         = opt$rm_rate
time_step       = opt$tstep
auto_num        = opt$auto_num
NH3_con         = opt$NH3_con
HT_con          = opt$HT_con
photon_con      = opt$photon_con
ddca_con        = opt$ddca_con
lxlycm          = opt$lxlycm
Cconc           = opt$ccon
diameter_um     = opt$diameter_um

#####################################################################################################################################################################

############ Katana ############ 
diet <- readRDS(infile_diet)
#####################################################################################################################################################################
getwd <- getwd()
getwd

# create a new folder:
new_folder = paste(otu_name,'_',keywd,'_death_',death_r,'_rmRate_',rm_rate,'_gridno_',arena_mn,'_lxlycm_',lxlycm,'_Ccon_',Cconc,'_Cd_',diameter_um,'_ammo',NH3_con,'_HT',HT_con,'_ddca',ddca_con,'/',sep = '')
dir.create(paste(getwd,'/',new_folder,sep = ''))

dir_gs_models='/srv/scratch/z5265700/Shan_z5095298/z5095298/sponge_modeling/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/'

# OTU1 v20210902:
STYtaxon_1='/OTU01_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU01_MOM_20210825/20210831__OTU01_MOM_20210825__TrainDiet_O2_bothN_no0_sulfite_NoHT/STY_Merged_OTU01_MOM.RDS'
# OTU2 v20210909 (Photosymbiont use nitrate -> NH3. This was used in 3-taxa co-culture)
STYtaxon_2='/OTU02_renamed/20210505_diet_sw_no_oxygen_STY_Merged_OTU02-draft__12_NO3/STY_Merged_OTU02.RDS'
# OTU3 v20210902:
STYtaxon_3='/OTU03_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU03_MOM_20210802hypotau_sulfitetransport_nitritetransport_adh/20210802_diet_sw_STY_Merged_OTU03-MOM_20210802hypotau_sulfitetransport_nitritetransport_adh__12_hypotaurine_nitrite_NoO2_noVB_noNO3/STY_Merged_OTU03_MOM.RDS'
# OTU4 v20210902:
STYtaxon_4='/OTU04_renamed/2021-05-11_ReactionPool_to_ObjectModel/07_adapt_addReact_R/04_Rscp/OTU04_MOM_20210825/20210825__OTU04_MOM_20210825__TrainDiet_O2_bothN_no0_sulfite_NoHT/STY_Merged_OTU04_MOM.RDS'
# OTU5 v20210902:
STYtaxon_5='/OTU05_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU05_MOM_20210813_sulftTrans/20210815_diet_sw_STY_Merged_OTU05-MOM_20210813_sulftTrans__12_SOB_NoHT_O2_bothN_Sulfite_no0/STY_Merged_OTU05_MOM.RDS'
# OTU6 v20210902:
STYtaxon_6='/OTU06_renamed/20210515_STY_Merged_OTU06_mineral_sw_3_NoH2S/STY_Merged_OTU06.RDS'
# OTU7 v20210902:
STYtaxon_7='/OTU07_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU07_MOM_20210813/20210815_diet_sw_STY_Merged_OTU07-MOM_20210813__12_SOB_noHT_O2_bothN_sulfite_no0/STY_Merged_OTU07_MOM.RDS'
# OTU8 v20210909 (The AOA - NO forming deactivate!! This was used in 3-taxa co-culture):
STYtaxon_8='/OTU08_renamed_arc/Step_protocals_using_R/07_gf_model/20210515_STY_Merged_OTU08_mineral_sw_3_unlimited_CO2/STY_Merged_OTU08_MgfM/STY_Merged_OTU08_MgfM.RDS'


model1 <- readRDS(paste(dir_gs_models,STYtaxon_1,sep = "")) #d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__UBA10353;f__LS-SOB;g__;s__
model2 <- readRDS(paste(dir_gs_models,STYtaxon_2,sep = "")) #d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__Synechococcales_A;f__Cyanobiaceae;g__Synechococcus_C;s__GCA_001628295.1
model3 <- readRDS(paste(dir_gs_models,STYtaxon_3,sep = "")) #d__Bacteria;p__Myxococcota;c__UBA9160;o__UBA9160;f__UBA6930;g__;s__
model4 <- readRDS(paste(dir_gs_models,STYtaxon_4,sep = "")) #d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__UBA10353;f__LS-SOB;g__;s__
model5 <- readRDS(paste(dir_gs_models,STYtaxon_5,sep = "")) #d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__UBA6729;f__;g__;s__
model6 <- readRDS(paste(dir_gs_models,STYtaxon_6,sep = "")) #d__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__UBA8639;g__bin75;s__
model7 <- readRDS(paste(dir_gs_models,STYtaxon_7,sep = "")) #d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__;f__;g__;s__
model8 <- readRDS(paste(dir_gs_models,STYtaxon_8,sep = "")) #d__Archaea;p__Crenarchaeota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Cenarchaeum;s__


replicates <- Cl_no
cores <- Cl_no
cl <- makeCluster(cores, type="PSOCK")
clusterExport(cl, c("diet","model1","model2","model3","model4","model5","model6","model7","model8",
                    "getwd","new_folder","arena_mn","death_r","Inocc_no","iter_no","replicates",
                    "cores","setAllExInf_value","rm_rate","time_step","auto_num","NH3_con","HT_con","photon_con","ddca_con",
                    "Cconc", "diameter_um","lxlycm"))
clusterEvalQ(cl, sink(paste0(getwd,'/',new_folder, Sys.time(), ".txt")))

simlist <- parLapply(cl, 1:replicates, function(i){
  sybil::SYBIL_SETTINGS("SOLVER", "cplexAPI") #https://github.com/euba/BacArena/issues/152
  
  print("====================================================================================================")
  print("============================= Best_models_STY_8taxa_20230920: auto step1 =======================================")
  print(paste("=============================", Sys.time(), '=======================================', sep = ' '))
  print(paste("==Parysynechococcus use nitrate -> NH3. This was used in 3-taxa co-culture ===================================", sep = ' '))
  print(paste("==AOA - NO forming deactivate!! This was used in 3-taxa co-culture =======================================", sep = ' '))
  print(paste("==OTU1 & 4: SOB and NOB =======================================", sep = ' '))
  print(paste("============================= auto_num=", auto_num, '=======================================', sep = ' '))
  print("====================================================================================================")
  # bacterium_SOB <- BacArena::Bac(my_model,lyse = F,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  # lxlycm = (arena_mn * arena_mn * Inocc_no / ( Cconc * diameter_um / 10000) )^0.5
  arena <- BacArena::Arena(n=arena_mn,m=arena_mn,Lx=lxlycm,Ly=lxlycm)
  inocc <- arena@n * arena@m * Inocc_no # Inocc_no = 1
  bacterium1        <- BacArena::Bac(model1,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, cellarea = 3.14*(diameter_um/2)*2, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium2        <- BacArena::Bac(model2,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, cellarea = 3.14*(diameter_um/2)*2, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium3        <- BacArena::Bac(model3,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, cellarea = 3.14*(diameter_um/2)*2, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium4        <- BacArena::Bac(model4,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, cellarea = 3.14*(diameter_um/2)*2, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium5        <- BacArena::Bac(model5,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, cellarea = 3.14*(diameter_um/2)*2, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium6        <- BacArena::Bac(model6,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, cellarea = 3.14*(diameter_um/2)*2, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium7        <- BacArena::Bac(model7,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, cellarea = 3.14*(diameter_um/2)*2, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium8        <- BacArena::Bac(model8,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, cellarea = 3.14*(diameter_um/2)*2, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  
  arena <- BacArena::addOrg(object = arena, specI = bacterium1, amount = inocc*0.18220)
  arena <- BacArena::addOrg(object = arena, specI = bacterium2, amount = inocc*0.02740)
  arena <- BacArena::addOrg(object = arena, specI = bacterium3, amount = inocc*0.11580)
  arena <- BacArena::addOrg(object = arena, specI = bacterium4, amount = inocc*0.20940)
  arena <- BacArena::addOrg(object = arena, specI = bacterium5, amount = inocc*0.05150)
  arena <- BacArena::addOrg(object = arena, specI = bacterium6, amount = inocc*0.03690)
  arena <- BacArena::addOrg(object = arena, specI = bacterium7, amount = inocc*0.09090)
  arena <- BacArena::addOrg(object = arena, specI = bacterium8, amount = inocc*0.28570)
  
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[1:20], difunc = "pde",
                             pde = "Diff2d", smax = diet$Input_mM[1:20], unit = "mM", add = F) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[21:29], difunc = "pde",
                             pde = "Diff2d", smax = diet$Input_mM[21:29], unit = "mM", add = F) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[22], difunc = "pde",
                             pde = "Diff2d", smax = NH3_con, unit = "mM", add = F) #[1] replace EX_cpd00013_e0	N as ammonia
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[17], difunc = "pde",
                             pde = "Diff2d", smax = HT_con, unit = "mM", add = F) #[1] replace EX_cpd00406_e0	S as HT
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[15], difunc = "pde",
                             pde = "Diff2d", smax = photon_con, unit = "mM", add = F) #[1] replace EX_cpd11632_e0	photon (hn)
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[61], difunc = "pde",
                             pde = "Diff2d", smax = ddca_con, unit = "mM", add = F) #[1] replace EX_cpd00222_e0	Gluconate
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[62], difunc = "pde",
                             pde = "Diff2d", smax = diet$Input_mM[62], unit = "mM", add = F) #[1] replace EX_cpd00222_e0	Gluconate

  arena@tstep <- time_step
  
  BacArena::chemotaxis(arena@specs[["STY_Merged_OTU06"]], arena, 1, chemo='EX_cpd00075_e0', arena@occupyM)#N as nitrite

  simulation <- BacArena::simEnv(object = arena, time = 1, sec_obj = "mtf", continue = T) #pFBA; iter 1

})
stopCluster(cl)
saveRDS(simlist, file = paste(getwd,'/',new_folder,'/BacArena_',otu_name,'_',as.character(arena_mn*arena_mn),'_grids_mineral_sw_',auto_num,'.RDS',sep = ''))

# ############################ Download RDS file and run on MacOS ############################
getwd <- getwd()
getwd
simulation_loop <- readRDS(paste(getwd,'/',new_folder,'/BacArena_',otu_name,'_',as.character(arena_mn*arena_mn),'_grids_mineral_sw_',auto_num,'.RDS',sep = ''))
