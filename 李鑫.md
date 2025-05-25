# create leaf trait dataset for analyses

#########################################
# packages necessary
#########################################
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)

#################### load data containing leaf traits (Jen Firn) ####################################################################
#####################################################################################################################################

leaf_traits = read.csv("./input/NutNet-foliar-traits-7JAN2017.csv")
names(leaf_traits)
nrow(leaf_traits) ## 2764
length (unique(leaf_traits$site_code)) # 27 sites
length (unique(leaf_traits$Taxon)) # 244 species 
length (unique(leaf_traits$trt)) # 10 trt
table(leaf_traits$site_code, leaf_traits$block)

block_counts <- leaf_traits %>% 
  group_by(site_code) %>% 
  summarise(num_blocks = n_distinct(block))
View(block_counts)
## each site contains at least three blocks (except bldr.us which has 2 blocks), each block contains different plots corresponding to different trt.
## in each trt we have leaf data for different species, the number and type of species vary between plots 

##################### load data containing SLA version 2 ##################### 
leaf_v2 =  read.csv("./input/foliar_cover_updated_2.csv")
names(leaf_v2)
nrow(leaf_v2) ## 2664
length (unique(leaf_v2$site_code)) # 27 sites
length (unique(leaf_v2$Taxon)) # 243 species 
length (unique(leaf_v2$trt)) # 10 trt

leaf_new_sla = leaf_v2[, c(1:5, 13,29:35)]  # selection of site_code, plot, year, Taxon, block, SLA_v2 
leaf = left_join(leaf_traits, leaf_new_sla) ### join by Taxon, site_code, year, block, plot 
names(leaf)
names(leaf_new_sla)


leaf[leaf == 'NULL'] <- NA ## replace cells with NULL by NA

leaf[, 8:38] <- sapply(leaf[, 8:38], as.numeric) ## set columns with measurements to numeric 
sapply(leaf, class) ## verify the data types for each column
summary(leaf)

leaf$Ntrt = 0 ## creation of a column with Ntrt= 0 
leaf$Ntrt[leaf$trt == 'N' | leaf$trt == 'NP' | leaf$trt == 'NK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1 ## set Ntrt=1 if any point receives N
### repeat the same for the other treatments 
leaf$Ptrt = 0
leaf$Ptrt[leaf$trt == 'P' | leaf$trt == 'NP' | leaf$trt == 'PK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1
leaf$Ktrt = 0
leaf$Ktrt[leaf$trt == 'K' | leaf$trt == 'NK' | leaf$trt == 'PK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1

leaf$Nfix = 'no' ## creation of a column with Nfix as no
leaf$Nfix[leaf$Family == 'Fabaceae'] = 'yes' ## when the families are Fabaceae, set them as N fixers, the other families as non fixers

write.csv(leaf, "./output/leaf.csv")
#######################################################################################################################################

#################### load data containing biomass and cover percentage ####################################################################
#####################################################################################################################################
#### This file contains information on species richness, percentage cover, biomass for each plot (not per species but per plot)

core =  read.csv("./input/comb-by-plot-09-April-2018.csv")
names(core) 
nrow(core) ## 14939
length (unique(core$site_code)) # 100 sites

core[core == 'NULL'] <- NA ## replace NULL by NA

core[, 26:45] <- sapply(core[, 26:45], as.numeric) # set columns with measurements as numeric
sapply(core, class)

core$par_rat = core$Ground_PAR / core$Ambient_PAR ### calculate the ratio of ground PAR over Ambiant PAR 
core$LAI = -log(core$Ground_PAR / core$Ambient_PAR)/0.5 ### calculate the ratio of ground PAR over Ambiant PAR 

#### selection of the columns: site_code, year, block, plot, trt, from sum_INT_cover to the end 
core_data = core[, c(2, 25, 15, 16, 22, 26:47)]  

#### selection of the columns: site_code, block, plot, trt, site_name, lat,long and other columns containing information
### about each site
core_info = core[, c(2, 15, 16, 22, 1, 3:14, 17:21, 23, 24)]


############## merge leaf with core data for each site/year/block/plot/trt  ##########################
leaf_core_data = join(leaf, core_data, by = c("site_code", "year", "block", "plot", "trt"), type = 'left', match = 'first')

leaf_core = join(leaf_core_data, core_info, by = c("site_code", "block", "plot", "trt"), type = 'left', match = 'first')

write.csv (leaf_core, "./output/leaf_core.csv")
##### here for each species in the same plot: same data of core for the whole plot, no data of biomass or PAR per species 
##########################################################################################################################


########################################### load and add climate data to leaf_core #########################################################

spei = read.csv("./input/CRU-annual_2018-07-06.csv")
names(spei)   ## this file contains for each site the annual data from 1903 to 2016 of precip, PET and SPEI-6,SPEI-12, SPEI-24
table(spei$year)
spei$p_pet = spei$precip / spei$PET  ### Calculate aridity Index 
sapply(spei, class)

## precip_mean and pet_mean in the file SPEI are  averages for the period 1975-2016
leaf_core_spei = join(leaf_core, spei, by = c("site_code", "year"), type = 'left', match = 'first')
names(leaf_core_spei)

write.csv(leaf_core_spei, "./output/leaf_core_spei.csv")

# check
nrow(leaf) # 2788
nrow(leaf_core) ## 2788
nrow(leaf_core_spei) ## 2788

######################################## load growing season climate data: temperature, par, vpd, z  average for the period 1901-2015 ###################
## these data are for each longitude and latitude, .... ########################### 

tmp_globe = read.csv("./input/cru_tmp_climExtract_growingseason_globe.csv")
par_globe = read.csv("./input/cru_par_climExtract_growingseason_globe.csv")
vpd_globe = read.csv("./input/cru_vpd_climExtract_growingseason_globe.csv")
z_globe =  read.csv("./input/z_globe.csv")

leaf_core_spei$tmp = NA
leaf_core_spei$par = NA
leaf_core_spei$vpd = NA
leaf_core_spei$z = NA

climate_df = c()
for (i in 1:nrow(leaf_core_spei)){
  
  currentLat = leaf_core_spei$latitude[i]
  currentLon = leaf_core_spei$longitude[i]
  
  clim_comb = tmp_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  tmp = subset(clim_comb, lat==best_lat & lon==best_lon)$tmp
  
  clim_comb = par_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  par = subset(clim_comb, lat==best_lat & lon==best_lon)$par
  
  clim_comb = vpd_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  vpd = subset(clim_comb, lat==best_lat & lon==best_lon)$vpd
  
  clim_comb = z_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  z = subset(clim_comb, lat==best_lat & lon==best_lon)$z
  
  climate_df = rbind(climate_df, c(tmp, par, vpd, z, currentLat, currentLon, best_lat, best_lon, i))
  
  
}

plot(climate_df[,6] ~ climate_df[,8]) # longitude match verification
plot(climate_df[,5] ~ climate_df[,7]) # latitude match verification

leaf_core_spei$tmp = climate_df[, 1] # tmp extraction from climate_df + addition to leaf_core_spei df
leaf_core_spei$par = climate_df[, 2] # par extraction from climate_df + addition to leaf_core_spei df
leaf_core_spei$vpd = climate_df[, 3] # vpd extraction from climate_df + addition to leaf_core_spei df
leaf_core_spei$z = climate_df[, 4] # elevation extraction from climate_df + addition to leaf_core_spei df
names(leaf_core_spei)

write.csv(leaf_core_spei, "./output/leaf_core_spei.csv")

########################################################
# add species information
########################################################
info = read.csv("./input/full-species-info-25-January-2019.csv")
nrow(info) ## 5897
length(unique(info$site_code)) # 108 sites

### this file contains information about species whether they are grass or Forb, C3 or C4, perennial or annual 

code_taxon_data = paste(toupper(leaf_core_spei$site_code), toupper(leaf_core_spei$Taxon), sep = ' ')

info_caps = mutate_all(info, .funs = toupper) ### transform all lower cases to upper cases 
names(info_caps)
info_caps$code_taxon = paste(info_caps$site_code, info_caps$standard_taxon, sep = ' ') #bind site_code to standard_taxon with space between them


n_info = NULL

for (i in 1:length(code_taxon_data)){
  
  ancillary_data = subset(info_caps, code_taxon == code_taxon_data[i])[, c(6:10)]
  
  ancillary_data$n = i
  
  n_info = rbind(n_info, ancillary_data)
  
}

leaf_core_spei_info = cbind(leaf_core_spei, n_info)
leaf_core_spei_info


#############################
## add in new growing season climate data
#############################

gs_climate <- read.csv("./input/cru_gs_climate.csv")
gs_climate

leaf_core_spei_info_gsclimate <- left_join(leaf_core_spei_info, gs_climate, by = c("site_code"))

write.csv(leaf_core_spei_info_gsclimate, "./output/leaf_core_spei_info_gsclimate.csv")
nrow(leaf_core_spei_info_gsclimate) #2788

###################################### add in species composition data ###############################################
######################################################################################################################
cover = read.csv("./input/full-cover-31-August-2020.csv")
names(cover)
nrow(cover) ## 239683

cover_select = select(cover, site_code, plot, year, Taxon, max_cover)

core_leaf_spei_cover = left_join(leaf_core_spei_info_gsclimate, cover_select)
nrow(core_leaf_spei_cover) # 2788
names(core_leaf_spei_cover)
names(leaf_core_spei_info_gsclimate)
####################################### load biomass data ##############################
full_biomass<-read.csv("./input/full-biomass-nutrients-06-December-2017.csv")
names(full_biomass)
length(unique(full_biomass$site_code)) ## 27
full_biomass$trt
########### for each site, block and plot, we have in this file data per category
#### forb, graminoid etc... not per species 

full_biomass[full_biomass == 'NULL'] <- NA # replace empty values by NA
full_biomass$pct_N <- as.numeric(full_biomass$pct_N)
sapply(full_biomass, class)

full_biomass$Ntrt = 0
full_biomass$Ntrt[full_biomass$trt == 'N' | full_biomass$trt == 'NP' | full_biomass$trt == 'NK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1 # when N is added: Ntrt=1

full_biomass$Ptrt = 0
full_biomass$Ptrt[full_biomass$trt == 'P' | full_biomass$trt == 'NP' | full_biomass$trt == 'PK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1 # when P is added: Ptrt=1

full_biomass$Ktrt = 0
full_biomass$Ktrt[full_biomass$trt == 'K' | full_biomass$trt == 'NK' | full_biomass$trt == 'PK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1 # when K is added: Ktrt=1

full_biomass$Ntrt_fac = as.factor(full_biomass$Ntrt)
full_biomass$Ptrt_fac = as.factor(full_biomass$Ptrt)
full_biomass$Ktrt_fac = as.factor(full_biomass$Ktrt)


full_biomass$pct_C = as.numeric(as.character(full_biomass$pct_C))
full_biomass$pct_N = as.numeric(as.character(full_biomass$pct_N))
full_biomass$pct_P = as.numeric(as.character(full_biomass$pct_P))
full_biomass$pct_K = as.numeric(as.character(full_biomass$pct_K))
full_biomass$pct_Mg = as.numeric(as.character(full_biomass$pct_Mg))
full_biomass$pct_Ca = as.numeric(as.character(full_biomass$pct_Ca))


## AGN = above ground N!, mass is leaf mass
full_biomass$AGN <- full_biomass$mass*(full_biomass$pct_N*0.01) # calculation of aboveground N quatity

full_biomass$Nfix = 'no'
full_biomass$Nfix[full_biomass$category == 'LEGUME'] = 'yes' # set legumes as N fixers 

nrow(full_biomass) # 2102

biomass_core_leaf_spei = join(core_leaf_spei_cover, full_biomass, by = c("site_code", "year"), type = 'left', match = 'first')
names(biomass_core_leaf_spei)
nrow(biomass_core_leaf_spei) # 2788

core_leaf_spei_cover$Ntrt_fac
### add growing season length data#######################
gs_information <- read.csv("./input/Weather_site_inventory_20200921.csv")
head(gs_information)
colnames(gs_information)

gs_length <- gs_information[,c(2, 20)] ## selection of the colum site_code and gs_len

biomass_core_leaf_spei <- join(biomass_core_leaf_spei, gs_length, by = c("site_code"), type = 'left', match = 'first')

biomass_core_leaf_spei$site_code[is.na(biomass_core_leaf_spei$gs_len)]
biomass_core_leaf_spei$gs_len[biomass_core_leaf_spei$site_code == 'gilb.za'] <- 6
biomass_core_leaf_spei$gs_len[biomass_core_leaf_spei$site_code == 'summ.za'] <- 8

biomass_core_leaf_spei$gs_frac <- biomass_core_leaf_spei$gs_len / 12

###### add soil moisture data#######################
soil <- read.csv("./input/soil_data.csv")
biomass_core_leaf_spei <- join(biomass_core_leaf_spei, soil, by = c("site_code"), type = 'left', match = 'first')

write.csv(biomass_core_leaf_spei, "./output/biomass_core_leaf_spei_final.csv")
biomass_core_leaf_spei$max_cover

#################################### get PET and aridity index : mean annual variables between 1970 and 2000  ###############
#######################################################################################################################

########## Load necessary librairies ########
library(sf)
library(raster)
library(tmap)
library(mapview)

getwd()
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Ecology_Letters/input")
l = list.files(pattern = "v3_yr") #reading world maps containing the string mask in the crus_rasters (maps from jeff work)
l
################# loop to extract all the pixels (raster) from the world maps of PET data ############################# 
######## https://cgiarcsi.community/data/global-aridity-and-pet-database/################

l = list.files(pattern = "v3_yr")
l = c(l, list.files(pattern = "proj_elev", full.names = T)) # bring in new, projected elev data

r = list()
for(i in l){
  r1 = raster(i)
  r[[i]] = r1
  print(r1)
}

r = raster::stack(l)
crs(r) = 4326


records =  biomass_core_leaf_spei %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = crs(r))

records = st_transform(records, crs = st_crs(r))
st_crs(records)==st_crs(r)


traits = cbind(records,
            raster::extract(r, st_coordinates(records), method = 'simple'))
traits

traits$Long = st_coordinates(traits)[,1]
traits$Lat = st_coordinates(traits)[,2]
st_geometry(traits) = NULL
traits$crs = 4326
names(traits)

colnames(traits)[colnames(traits) == "awi_pm_sr_yr"] <- "aridity"
traits$aridity= traits$aridity * 0.0001
hist(traits$aridity)
plot(traits$aridity, traits$SPEI_12)

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Ecology_Letters")
write.csv(traits, "./output/traits.csv")
hist(traits$aridity)
