
################################################################################
# Name        : ecological_modelling/s1_collect_data.R
# description : Yamin's MS project
# script      : Collect data and save to disk
################################################################################

rm(list = ls()); gc() # or restart R from the menu to free up more space
# Species     : Chromolaena odorata (Siam weed)
# Read        : https://en.wikipedia.org/wiki/Chromolaena_odorata
# Read        : http://www.floraofbangladesh.com/2020/09/boro-shialmuti-or-siam-weed-chromolaena.html

# Article 1   : https://doi.org/10.1016/j.gecco.2020.e01196
# Article 2   : https://www.researchgate.net/publication/348002595_Invasive_Alien_Species_of_Bangladesh
# Article 3   : https://doi.org/10.3759/tropics.SINT04
# Article 4   : https://doi.org/10.3390/agronomy12071592
################################################################################
# Step 1: Install and load required libraries

# install libraries, run only first time 
# install.packages("terra",dependencies = T)
# install.packages("geodata",dependencies = T)
# install.packages("spocc",dependencies = T)

# load the librariers
library(terra)
library(sf)
library(geodata)
library(spocc)
library(dplyr)
library(ggplot2)
################################################################################
# Step 2: Get occurrence data of Siam weed from "GBIF" database

# Get Bangladesh country boundary
temp <-tempdir()
bd_0 <- gadm(country = "Bangladesh",level = 0,path = "./data") %>% st_as_sf()
ind_0 <- gadm(country = "India",level = 0,path = temp)

y_max <- 26.63
y_min <- 20.60
x_max <- 92.7
x_min <- 87.96

# Calculate the height and width of the rectangle
height_km <- y_max - y_min
width_km <- x_max - x_min
# Calculate the area in square kilometers
area_km2 <- height_km * width_km
# Print the result
cat("The area is", area_km2, "KM^2\n")

boundary <- ext(x_min,x_max,y_min,y_max)
plot(ind_0)
plot(boundary,add=T,lty=2)


# Download Siam weed occurrence data 
df <- spocc::occ(
  query = 'Chromolaena odorata',           # the species data to collect
  from = 'gbif',                           # grom which database
  gbifopts = list(hasCoordinate = TRUE),   # we collect only these points that have Longitude and Latitude information
  geometry =  st_bbox(bd_0),               # collect point only within Bangladesh
  limit = 100) %>%  occ2df()               # limit the size to 100 and finally convert to data.frame object


# Download predictors from the WorldClim website
bio<- geodata::worldclim_country(country = "Bangladesh",res=0.5,var="bio",path = "./data")
names(bio) <- paste0("Bio",1:19) # rename the layers

################################################################################
# Step 3:  Save the data to disk to use later
# write.csv(df,row.names = F,file = "./data/occ.csv") # save the occurrence data
# writeRaster(bio,filename = "./data/bio.grd",overwrite=T)        # save the predictors (image) data



# Test to see if the csv file was saved correctly
plot(read.csv("./data/occ.csv")[,c("longitude","latitude")])
# Test if the saving process was okay or not.
plot(rast("./data/bio.grd"))

################################################################################
# Step 4: Make a study area map
ggplot()+
  geom_sf(data = bd_0)+
  geom_point(data = df,aes(longitude,latitude))+
  labs(x="Longitude",y="Latitude")

ggsave("./result/Fig1_study_area.jpg",dpi=300)  

