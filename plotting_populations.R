### Get worldclim data from point locality information
### Locality information came from the Boechera_stricta_latslong.xlsx file.

#skip to line 245 if just comparing climate variables for the populations


library(raster)
library(sp)
library(rasterVis)
library(dplyr)
library(geodata) #to download the tif files straight from world clim. 
#library(rgdal) #load for the AI 
library(ncdf4) #for extracting the snowclim data (they are in .nc format)



#set the working directory to the location of climate data and files. this is on the 24 TB. make sure you are on UGA network
setwd("/Volumes/JTA_24TB/Climate_variables/") 

# read in datafile with the longs/lats
latlong <-read.table("boechera_population_data/lat_long.txt", header=T, sep="\t",row.names=1)
# I learned that the longitude needs to be before the latitude to get extract to work

latlongs <-select(latlong,-population) #remove the pop column for the extract command. it expects only the lat and long columns

#### extracting data from World clim ####
#use  worldclim_tile to extract tiles, the resolution (minutes of a degree) is specified in "res". Have to specify a path to save the tif files, so i saved them originally on my computer. they are now on the 24 tb so no need to run this again unless you are updating

# read more here: https://rdrr.io/cran/geodata/man/worldclim.html

#avg_temp<-worldclim_tile(var="tavg", res=0.5,lon=-107, lat=38, path="~/Desktop/Anderson_data/other/worldclim/", version="2.1")
#min_temp<-worldclim_tile(var="tmin", res=0.5,lon=-107, lat=38, path="~/Desktop/Anderson_data/other/worldclim/", version="2.1")
#max_temp<-worldclim_tile(var="tmax", res=0.5,lon=-107, lat=38, path="~/Desktop/Anderson_data/other/worldclim/", version="2.1")
#precip<-worldclim_tile(var="prec", res=0.5,lon=-107, lat=38, path="~/Desktop/Anderson_data/other/worldclim/", version="2.1")
#bioclim_elev<-worldclim_tile(var="elev", res=0.5,lon=-107, lat=38, path="~/Desktop/Anderson_data/other/worldclim/", version="2.1")
#bioclim_var<-worldclim_tile(var="bio", res=0.5,lon=-107, lat=38, path="~/Desktop/Anderson_data/other/worldclim/", version="2.1")
#For a definition of the bioclimatic variables: http://www.worldclim.org/bioclim

#reading in the tif files from the 24tb, 
#avg_temp <- raster("wc2.1_tiles/tile_15_wc2.1_30s_tavg.tif") 
#min_temp <- raster("wc2.1_tiles/tile_15_wc2.1_30s_tmin.tif")
#max_temp <- raster("wc2.1_tiles/tile_15_wc2.1_30s_tmax.tif")
#precip <- raster("wc2.1_tiles/tile_15_wc2.1_30s_prec.tif")
elev <-raster("wc2.1_tiles/tile_15_wc2.1_30s_elev.tif")
#bioclim_elev <- raster("wc2.1_tiles/tile_15_wc2.1_30s_bio.tif") #i think some of these values are incorrect. for example, bio clim 12 is supposed to be the precip of the driest month, but is in the 800s

##Plotting
plot(avg_temp)
plot(min_temp)
plot(max_temp)
plot(precip)
plot(bioclim_var)
#For a definition of the bioclimatic variables: http://www.worldclim.org/bioclim
#For help with extracting data: http://creativemorphometrics.co.vu/blog/2014/03/27/extracting-climate-data-in-r/


#Here is another resource:https://ecologicaconciencia.wordpress.com/tag/worldclim/

# read in datafile with the longs/lats


setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper/") 
latlong <-read.csv("grasshopper_population_lats_longs.csv")



setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/field/") 
latlong <-read.csv("field_population_lats_longs.csv")


newext<-c(-107.10,-106.79,38.69,39.05) #coordinates of gothic valley

elev.new<-crop(elev,newext)




plot(elev.new$tile_15_wc2.1_30s_elev)
#plotting the points onto the map
points(latlong $Longitude, latlong $Latitude, col='black', pch=20)

library(dplyr)

#schofield <-filter(latlong, population == 250)

#points(schofield $Longitude, schofield $Latitude, col='red', pch=8, cex=2)


gothic <-filter(latlong, population == 283)

points(gothic $Longitude, gothic $Latitude, col='red', pch=18, cex=2)

#pcoordinates(latlong)<-c("Longitude", "Latitude")


tmin.new<-crop(min_temp,newext)
plot(tmin.new$tile_15_wc2.1_30s_tmin_1)
plot(tmin.new$tile_15_wc2.1_30s_tmin_2)
plot(tmin.new$tile_15_wc2.1_30s_tmin_3)
plot(tmin.new$tile_15_wc2.1_30s_tmin_4)
plot(tmin.new$tile_15_wc2.1_30s_tmin_5)
plot(tmin.new$tile_15_wc2.1_30s_tmin_6)
plot(tmin.new$tile_15_wc2.1_30s_tmin_7)
plot(tmin.new$tile_15_wc2.1_30s_tmin_8)
plot(tmin.new$tile_15_wc2.1_30s_tmin_9)
plot(tmin.new$tile_15_wc2.1_30s_tmin_10)
plot(tmin.new$tile_15_wc2.1_30s_tmin_11)
plot(tmin.new$tile_15_wc2.1_30s_tmin_12)

tmax.new<-crop(max_temp,newext)
plot(tmax.new$tile_15_wc2.1_30s_tmax_1)
plot(tmax.new$tile_15_wc2.1_30s_tmax_2)
plot(tmax.new$tile_15_wc2.1_30s_tmax_3)
plot(tmax.new$tile_15_wc2.1_30s_tmax_4)
plot(tmax.new$tile_15_wc2.1_30s_tmax_5)
plot(tmax.new$tile_15_wc2.1_30s_tmax_6)
plot(tmax.new$tile_15_wc2.1_30s_tmax_7)
plot(tmax.new$tile_15_wc2.1_30s_tmax_8)
plot(tmax.new$tile_15_wc2.1_30s_tmax_9)
plot(tmax.new$tile_15_wc2.1_30s_tmax_10)
plot(tmax.new$tile_15_wc2.1_30s_tmax_11)
plot(tmax.new$tile_15_wc2.1_30s_tmax_12)


#plotting the points onto the map
points(latlong $Longitude, latlong $Latitude, col='red', pch=19)
#pcoordinates(latlongs)<-c("Longitude", "Latitude")


#Let's move forward using bilinear
mintemp.bil<- extract(min_temp, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
avgtemp.bil<-extract(avg_temp, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
maxtemp.bil<-extract(max_temp, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
precip.bil<-extract(precip, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
elev.bil<-extract(bioclim_elev, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
bioclim.bil<-extract(bioclim_var, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)


##combining back with lat long data.
climate<-cbind(latlong,mintemp.bil, avgtemp.bil, maxtemp.bil, precip.bil, elev.bil,bioclim.bil)

write.table(climate, file="bioclimdata.txt")
#Make a new column with latitude in radians

latlong$latrad<-latlong$Latitude*(pi/180)


#Calculating extraterrestrial solar radiation. Units here are MJ/m-squared per day. We need mm/day. For that conversion, see http://www.fao.org/docrep/X0490E/x0490e07.htm we need to multiply by 0.408
library(sirad)
extrat(dayOfYear("2015-01-15"),0.6803669)
extrat(dayOfYear("2015-01-15"), latlong$latrad[1])
extrat(dayOfYear("2015-01-15"), latlong$latrad[2])
#Creating necessary vectors
population<-sort(unique(latlong $population))
RAJan<-rep(NA, 207) #the 207 represents the total number of populations
RAFeb<-rep(NA, 207)
RAMarch<-rep(NA, 207)
RAApril<-rep(NA, 207)
RAMay<-rep(NA, 207)
RAJune<-rep(NA, 207)
RAJuly<-rep(NA, 207)
RAAug<-rep(NA, 207)
RASept<-rep(NA, 207)
RAOct<-rep(NA, 207)
RANov<-rep(NA, 207)
RADec<-rep(NA, 207)

RA_Jan<-extrat(dayOfYear("2015-01-15"), latlong[,4])
#Converting to the proper units
RAJandaily<-RA_Jan$ExtraTerrestrialSolarRadiationDaily*0.408
RAJan<-RAJandaily*31
RAJan

RA_feb<-extrat(dayOfYear("2015-02-15"), latlong[,4])
#Converting to the proper units
RAfebdaily<-RA_feb$ExtraTerrestrialSolarRadiationDaily*0.408
RAfebdaily
RAFeb<-RAfebdaily*28
RAFeb

RA_m<-extrat(dayOfYear("2015-03-15"), latlong[,4])
#Converting to the proper units
RAmdaily<-RA_m $ExtraTerrestrialSolarRadiationDaily*0.408
RAmdaily
RAMarch <-RAmdaily*31
RAMarch

RA_a<-extrat(dayOfYear("2015-04-15"), latlong[,4])
#Converting to the proper units
RAadaily<-RA_a $ExtraTerrestrialSolarRadiationDaily*0.408
RAadaily
RAApril <-RAadaily*30
RAApril

RA_may<-extrat(dayOfYear("2015-05-15"), latlong[,4])
#Converting to the proper units
RAmaydaily<-RA_may $ExtraTerrestrialSolarRadiationDaily*0.408
RAmaydaily
RAMay <-RAmaydaily*31
RAMay

RA_june<-extrat(dayOfYear("2015-06-15"), latlong[,4])
#Converting to the proper units
RAjunedaily<-RA_june $ExtraTerrestrialSolarRadiationDaily*0.408
RAjunedaily
RAJune <-RAjunedaily*30
RAJune

RA_july<-extrat(dayOfYear("2015-07-15"), latlong[,4])
#Converting to the proper units
RAjulydaily<-RA_july $ExtraTerrestrialSolarRadiationDaily*0.408
RAjulydaily
RAJuly <-RAjulydaily*31
RAJuly

RA_aug<-extrat(dayOfYear("2015-08-15"), latlong[,4])
#Converting to the proper units
RAaugdaily<-RA_aug $ExtraTerrestrialSolarRadiationDaily*0.408
RAaugdaily
RAAug <-RAaugdaily*31
RAAug

RA_sept<-extrat(dayOfYear("2015-09-15"), latlong[,4])
#Converting to the proper units
RAseptdaily<-RA_sept $ExtraTerrestrialSolarRadiationDaily*0.408
RAseptdaily
RASept <-RAseptdaily*30
RASept

RA_o<-extrat(dayOfYear("2015-10-15"), latlong[,4])
#Converting to the proper units
RAodaily<-RA_o$ExtraTerrestrialSolarRadiationDaily*0.408
RAodaily
RAOct <-RAodaily*31
RAOct

RA_n<-extrat(dayOfYear("2015-11-15"), latlong[,4])
#Converting to the proper units
RAndaily<-RA_n $ExtraTerrestrialSolarRadiationDaily*0.408
RAndaily
RANov <-RAndaily*30
RANov

RA_d<-extrat(dayOfYear("2015-12-15"), latlong[,4])
#Converting to the proper units
RAddaily<-RA_d $ExtraTerrestrialSolarRadiationDaily*0.408
RAddaily
RADec <-RAddaily*31
RADec


 
RA<-data.frame(population, RAJan, RAFeb, RAMarch, RAApril, RAMay, RAJune, RAJuly, RAAug, RASept, RAOct, RANov, RADec)

write.table(RA, file="RA.txt")




#### extracting data from projected future climates, CMIP6 climate model data #### 
#use  worldclim_tile to extract tiles, the resolution (minutes of a degree) is specified in "res". Have to specify a path to save the tif files, so i saved them originally on my computer. they are now on the 24 tb so no need to run this again unless you are updating



# read more here: https://rdrr.io/cran/geodata/man/cmip6.html

future_temp<- cmip6_tile(lon=-107, lat=38, "EC-Earth3-Veg-LR", "585", "2081-2100", var="tmax", path="~/Desktop/Anderson_data/other/worldclim/",res=0.5)

bio10 <- cmip6_world("CNRM-CM6-1", "585", "2061-2080", var="bioc", res=10, path=tempdir())




#### extracting aridity from global aridity index ####

setwd("Global_Aridity_Index_and_Potential_Evapotranspiration")#location of the Global AI and PET data. It is on the 24TB, folder named "Global_Aridity_Index_and_Potential_Evapotranspiration"

AI_annual_average <- raster("Global-AI_ET0_v3_annual/ai_v3_yr.tif") #specifying the annual average AI

AI_pops<- extract(AI_annual_average, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)# using extract to get AI for all the population lats and longs. 

AI_pops <- AI_pops$awi_pm_sr_yr * 0.0001 #the AI values need to be multiplied by 0.0001

AI_dataframe<-cbind(latlong,AI_pops)
write.table(AI_dataframe, file="~/Desktop/Anderson_data/other/worldclim/AI.txt")



#to plot AI
newext<-c(-107.2,-106.7,38.69,39.05) #coordinates for gothic valley
AI_new<-crop(AI_annual_average,newext)
AI_new <- AI_new * 0.0001

plot(AI_new)
#plotting the points onto the map
points(latlong $Longitude, latlong $Latitude, col='red', pch=18, lwd=0.2)


#### extracting snowmelt parameters from snowclim ####

#there are more variables that i did not extract here. check the read me in the SnowClim folder in the 24TB

setwd("SnowClim/CTRL_current")#location of the snowclim data for contemporary snow data. It is on the 24TB, folder named "Snowclim"


snowendct<-stack("date_of_snowcover_end_CTRL.nc")
plot(snowendct)

snowend_pops<- extract(snowendct, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)# using extract to get Julian day of end of snow cover, averaged across years for all the population lats and longs. The day of the end of snow cover is defined as the last day of the last period of 5 consecutive days with snow depth greater than 10 mm. 

minswe<-stack("min_swe_CTRL.nc")
plot(minswe)

minswe_pops<- extract(minswe, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)# using extract to get Monthly minimum snow water equivalent (m), averaged across years for all the population lats and longs. 


maxswe<-stack("max_swe_CTRL.nc")
maxswe_pops<- extract(maxswe, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)# using extract to get Monthly maximum snow water equivalent (m), averaged across years for all the population lats and longs.

meanswe<-stack("mean_swe_CTRL.nc")
meanswe_pops<- extract(meanswe, latlongs, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)# using extract to get Monthly mean snow water equivalent (m), averaged across years for all the population lats and longs.

#compile the snow metrics to one file with the lats/longs
snow_dataframe<-cbind(latlong,snowend_pops,minswe_pops,maxswe_pops,meanswe_pops)

#write to file
write.table(snow_dataframe, file="/Volumes/JTA_24TB/Climate_variables/snow.txt")


#to plot snowmelt timing 
newext<-c(-107.08,-106.79,38.69,39.05) #coordinates for gothic valley
snowend_new<-crop(snowendct,newext)

plot(snowend_new)
#plotting the points onto the map
points(latlong $Longitude, latlong $Latitude, col='red', pch=18, lwd=0.2)




#### comparing the climate variables, i combined the files manually ####

climate <-read.csv("boechera_population_data/Boechera_stricta_latslongs052323.csv")


sapply(climate,class)
climate$elev2<-climate$elevation^2
climate$W_Elevation2<-climate$elevation_worldclim^2
climateNONA<-na.omit(climate)

#elevation
Elev<-lm(elevation_worldclim ~ elevation, data=climate)
summary(Elev)
plot(elevation_worldclim ~ elevation, data=climate, xlab="original elevation", ylab="world clim elevation", cex=1.3)
abline(Elev)
legend("topleft",legend=paste("R2=", format(summary(Elev)$r.squared,digits=3)))


#Aridity Index
AI<-lm(AI_GAI~ AI_Population, data=climate)
summary(AI)
plot(AI_GAI~ AI_Population, data=climate, xlab="original Aridity Index", ylab="GAI Aridity Index (MAP/MAE)", cex=1.3)
abline(AI)
legend("topleft",legend=paste("R2=", format(summary(AI)$r.squared,digits=3)))



elevation_AI2<-lm(AI_GAI~elevation+elev2, data= climate)
summary(elevation_AI2)
plot(AI_GAI~elevation, data= climate, xlab="Elevation (m above sea level)", ylab="Aridity Index (Mean annual precipitation/Mean annual potential evapotranspiration)", cex=1.3)
curve((8.653e-01-4.414e-04*x+9.912e-08*x^2), add=T,  from=2499, to=sort(climate $elevation, TRUE)[2])


elevation_AI<-lm(AI_GAI~elevation, data=climate)
summary(elevation_AI)
plot(AI_GAI~elevation, data=climate, xlab="Elevation (m above sea level)", ylab="Aridity Index (MAP/MAE)", cex=1.3)
abline(elevation_AI)


Welevation_AI2<-lm(AI_GAI~ elevation_worldclim + W_Elevation2, data=climate)
summary(Welevation_AI2)
plot(AI_GAI~elevation_worldclim, data=climate, cex=1.3)

curve((1.054-5.690e-04*x+1.199e-07*x^2), add=T,  from=2614, to=sort(climate$elevation_worldclim, TRUE)[2])

Welevation_AI<-lm(AI_GAI~ elevation_worldclim, data=climate)
summary(Welevation_AI)
plot(AI_GAI~ elevation_worldclim, data=climate, xlab="Elevation (m above sea level)", ylab="Aridity Index (MAP/MAE)", cex=1.3)
abline(Welevation_AI)





#percipitation

elevation_MAP<-lm(MAP_Population~elevation+elev2, data= climate)
summary(elevation_MAP)

plot(MAP_Population~elevation, data= climate, xlab="Elevation (m above sea level)", ylab="Mean Annual Precipitation (mm)", cex=1.3)
curve(( 1.199e+03 -6.476e-01*x+1.428e-04 *x^2), add=T,  from=min(climateNONA $elevation), to=sort(climateNONA $elevation, TRUE)[2])


elevation_MAP2<-lm(MAP~elevation, data=climate)
summary(elevation_MAP2)
plot(MAP~elevation, data=climate, xlab="Elevation (m above sea level)", ylab="Mean Annual Precipitation (mm)", cex=1.3)
abline(elevation_MAP2)

elevation_MAT<-lm(MAT~elevation, data=climate)
summary(elevation_MAT)


#mean annual precipitation
plot(MAP_Population~MAP_worldclim, data=climate, xlab="Worldclim Mean Annual Precipitation (degrees C)", ylab="Original Mean Annual Precipitation (degrees C)", cex=1.3)



#mean annual temp
plot(MAT_Population~MAT_worldclim, data=climate, xlab="Worldclim Mean Annual Temperature (degrees C)", ylab="Original Mean Annual Temperature (degrees C)", cex=1.3)



#### plots of the variables on elevation

p1 <- plot(MAT_Population~elevation, data=climate, xlab="Original Elevation (m above sea level)", ylab="Original Mean Annual Temperature (degrees C)", cex=1.3)

p2 <- plot(MAT_worldclim~elevation_worldclim, data=climate, xlab="Worldclim Elevation (m above sea level)", ylab="Worldclim Mean Annual Temperature (degrees C)", cex=1.3)


p3 <- plot(MAP_Population~elevation, data=climate, xlab="Original Elevation (m above sea level)", ylab="Original Mean Annual Precipitation (degrees C)", cex=1.3)


p4 <- plot(MAP_worldclim~elevation_worldclim, data=climate, xlab="Worldclim Elevation (m above sea level)", ylab="Worldclim Mean Annual Precipitation (degrees C)", cex=1.3)


p5 <- plot(AI_Population~elevation, data=climate, xlab="Original Elevation (m above sea level)", ylab="Original Aridity Index", cex=1.3)

p6 <- plot(AI_GAI~elevation_worldclim, data=climate, xlab="Original Elevation (m above sea level)", ylab="Global Aridity Index", cex=1.3)

p7 <- plot(DoY_snowmelt_bioclim~elevation, data=climate, xlab="Original Elevation (m above sea level)", ylab="Original Aridity Index", cex=1.3)


#snowmelt plots


elevation_snow<-lm(DoY_snowmelt_bioclim~elevation+elev2, data= climate)
summary(elevation_snow)

plot(MAP~elevation, data=climate, xlab="Elevation (m above sea level)", ylab="Mean Annual Precipitation (mm)", cex=1.3)
abline(elevation_MAP2)


p7 <- plot(DoY_snowmelt_bioclim~elevation, data=climate, xlab="Source Elevation (m above sea level)", ylab="Day of snowmelt, averaged across years", cex=1.3)






