########################################
####### MAXNET MODELING ##########
########################################

# clean R session
rm(list = ls()); gc() 
library(terra)
library(geodata)
library(sf)
library(ggplot2)

############################ Data Prepare Start ##############################
# the bio data sets are 19 predictors, each one image for one predictors
bio<- rast("./data/bio.grd")
bd_0 <- gadm(country = "Bangladesh",level = 0,path = "./data") %>% st_as_sf()

# read occurrence data 
df <- read.csv("./data/occ.csv")[,c("longitude","latitude")]
names(df) <- c("x","y") # Rename long and lat as X,Y

# We have 19 variables, but we need to choose variables that are not correlated
library(usdm)
(cor <- vifcor(bio,th=0.70,size = 10000))
bio <- bio[[c(2,3,9,15,18,19)]]  #Selected Bio

# Collect random absence points within Bangladesh
set.seed(11)
bg_df <- spatSample( bio, 1000, "random", cells=F, xy=TRUE, values=FALSE, as.df=TRUE, na.rm=T)

# Prepare a SWD object for model building
data <- prepareSWD(species="Chromolaena odorata",
                   p = df, #data.frame. presence locations.
                   a = bg_df, #data.frame. absence/background locations.
                   env= bio, #rast environmental variables ext locations.
)
############################ Data Prepare END ##############################

# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(data, test = 0.2, only_presence = TRUE, seed = 25)
train <- datasets[[1]]
test <- datasets[[2]]


############ @@@ Model Traing Start @@@ #####################
# Train a MAXNET "model"
model <- train(method = "Maxnet", data = train)
# cat("Training auc: ", auc(model))
# cat("Testing auc: ", auc(model, test = test))
pred <- predict(model, data = train,type = "cloglog")
map <- predict(model, data = bio,type = "cloglog")
plotPred(map)
############ @@@ Model Traing End @@@ #####################


#### @@@ K-FOld cv_model Start @@@ #####################################################

# Create the folds from the training dataset
folds <- randomFolds(test,k = 4, only_presence = TRUE,seed = 25)

# Train "cv_model" model
cv_model <- train(method = "Maxnet", data = train, folds = folds)
cat("Training auc: ", auc(cv_model))
cat("Testing auc: ", auc(cv_model, test = test))
m <- combineCV(cv_model)
plotROC(m, test=test)

#### @@@ K-FOld cv_model End @@@ #####################################################


############ @@@ hyperparameters Tuning Start @@@ #####################

# Define the hyperparameters to test for Maxnet
h <- list(reg = seq(0.1, 3, 0.1), fc = c("lq", "lh", "lqp", "lqph", "lqpht"))

# Test all the possible combinations with gridSearch
hytune <- gridSearch(cv_model, hypers = h, metric = "auc", test = test)
head(hytune@results[order(-hytune@results$test_AUC), ])  # Best combinations

# Combine cross validation models
hytune@models[[5]]
hytuneModel <- combineCV(hytune@models[[5]])
plotROC(hytuneModel, test = test)

# Train Maxnet
predRF <- predict(hytune@models[[5]], data = bio, type = "cloglog")
plotPred(predRF, lt = "Habitat\nsuitability",colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
############ @@@ hyperparameters Tuning End @@@ #####################


####### @@@ Genetic algorithm instead with optimizeModel start @@@ ####### 
om <- optimizeModel(hytune@models[[5]], hypers = h, metric = "auc", test = test, seed = 4)
head(om@results[order(-om@results$test_AUC), ])  # Best combinations

# Combine cross validation models for output with highest test AUC model
YAMiN <- combineCV(om@models[[1]])
plotROC(YAMiN, test = test)
ggsave("./result/maxnet/plotROC_optimize.jpg",dpi=300,width = 8,height = 5)

map <- predict(YAMiN, data = bio, type = "cloglog")
plotPred(map)
ggsave("./result/maxnet/map__optimizeModel.png",dpi=300,width = 8,height = 5)

####### @@@ Genetic algorithm instead with optimizeModel End @@@ ####### 




# ############ EXTRA TEST CODES #####################
# (ths <- thresholds(YAMiN, type = "cloglog"))
# plotPA(map, th = ths[3, 2])
# ggsave("./result/maxnet/map_Present_optimizeModel.png",dpi=300,width = 8,height = 5)
# 
# plotVarImp(varImp(YAMiN, permut = 5),color = "green")
# ggsave("./result/maxnet/plotVarImp_optimizeModel.png",dpi=300,width = 8,height = 5)
# 
# 
# 
# #-----------------------------------------------------------------#
# utm_crs <- crs("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs")
# # Transform the SpatRaster to the target UTM CRS
# map_transformed <- project(map, utm_crs)
# 
# #88, 93, 20.5, 27 
# # Access the extent information
# # Extract specific values from the extent object
# xmin <- map_transformed@ptr$extent$vector[1]
# xmax <- map_transformed@ptr$extent$vector[2]
# ymin <- map_transformed@ptr$extent$vector[3]
# ymax <- map_transformed@ptr$extent$vector[4]
# # Calculate the height and width of the rectangle
# height_km <- ymax - ymin
# width_km <- xmax - xmin
# # Calculate the area in square kilometers
# area_km2 <- height_km * width_km
# # Print the result
# cat("The area is", area_km2, "KM^2\n")
# #-----------------------------------------------------------------#
# 
# 
