# clean R session
rm(list = ls()); gc() 

library(SDMtune)
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

print(data)
attributes(data)
str(data@data)
hist(bio)
############################ Data Prepare END ##############################





# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(data, test = 0.2, only_presence = TRUE, seed = 25)
train <- datasets[[1]]
test <- datasets[[2]]

# Train a RF model
model <- train(method = "RF", data = train)
cat("Training auc: ", auc(model))
cat("Testing auc: ", auc(model, test = test))
plotVarImp(varImp(model, permut = 5),color = "green")
#ggsave("./result/RF/RF_plotVarImp.jpg",dpi=300,width = 8,height = 5)
plotROC(model, test = test)
#ggsave("./result/RF/RF_train_plotROC.jpg",dpi=300,width = 8,height = 5)

pred <- predict(model, data = data,type = "cloglog")
map <- predict(model, data = bio,type = "cloglog")
plotPred(map)
#ggsave("./result/RF/RF_Predict_Map.jpg",dpi=300,width = 8,height = 5)

(ths <- thresholds(model, type = "logistic", test = test))
plotPA(map, th = ths[3, 2])
#ggsave("./result/RF/RF_Thresholds.jpg",dpi=300,width = 8,height = 5)


#### K-FOld ####################################################################################
# Create the folds from the training dataset
folds <- randomFolds(test,k = 4, only_presence = TRUE,seed = 25)
# Train the model
cv_model <- train(method = "RF", data = test, folds = folds)
cat("Training auc: ", auc(cv_model))
cat("Testing auc: ", auc(cv_model, test = test))
h_rf <- list(mtry= seq(2,8,1),ntree= seq(100,1000,10),nodesize= 1)
exp_8 <- randomSearch(cv_model,hypers = h_rf,metric = "auc",pop = 50,seed = 65466)
head(exp_8@results[order(-exp_8@results$test_AUC), ])  # Best combinations

############################################################################
##Select First Model
exp_8@models[[1]]
# Train for Random Forest
predRF <- predict(exp_8@models[[1]], data = bio, type = "cloglog")
plotPred(predRF, lt = "Habitat\nsuitability",colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
############################################################################
#### K-FOld ####################################################################################

# Define the hyperparameters to test for Maxnet
#h <- list(reg = seq(0.1, 3, 0.1), fc = c("lq", "lh", "lqp", "lqph", "lqpht"))

##FOR R
h <- list(mtry = seq(2,8,1), ntree= seq(100,1000,10),nodesize  = 1)

# Test all the possible combinations with gridSearch
gs <- gridSearch(model, hypers = h, metric = "auc", test = test)
head(gs@results[order(-gs@results$test_AUC), ])  # Best combinations

############################################################################
##Select First Model
gs@models[[1]]
plotROC(gs@models[[1]], test = test)
# Train for Random Forest
predRF <- predict(gs@models[[1]], data = bio, type = "cloglog")
plotPred(predRF, lt = "Habitat\nsuitability",colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
############################################################################








# Use the genetic algorithm instead with optimizeModel
om <- optimizeModel(model, hypers = h, metric = "auc", test = test, seed = 4)
head(om@results)  # Best combinations




pred <- predict(model, data = data, type = "cloglog")
head(pred)
p <- data@data[data@pa == 1, ]
pred <- predict(model,data = p,type = "cloglog")
tail(pred)
map <- predict(model, data = bio, type = "cloglog")
plotPred(map)
#ggsave("./result/maxnet/Maxnet_without_train.png",dpi=300,width = 8,height = 5)
(ths <- thresholds(model, type = "cloglog"))
plotPA(map, th = ths[3, 2])
#ggsave("./result/maxnet/Maxnet_presence_absence.png",dpi=300,width = 8,height = 5)