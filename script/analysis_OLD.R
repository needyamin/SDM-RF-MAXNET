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
cor <- vifcor(bio,th=0.70,size = 10000)
cor
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

slotNames(data)

############################ Data Prepare END ##############################




#Maxnet model Train
maxnet_model <- train(method = "Maxnet", data = data)
auc(maxnet_model)
plotROC(maxnet_model)
#ggsave("./result/maxnet/Maxnet_ROC_curve.jpg",dpi=300,width = 8,height = 5)

pred <- predict(maxnet_model,
                data = data,
                type = "cloglog")

map <- predict(maxnet_model,
               data = bio,
               type = "cloglog")
plotPred(map)
#ggsave("./result/maxnet/Maxnet_without_train.png",dpi=300,width = 8,height = 5)

ths <- thresholds(maxnet_model, type = "cloglog")
ths
plotPA(map, th = ths[3, 2])
#ggsave("./result/maxnet/Maxnet_presence_absence.png",dpi=300,width = 8,height = 5)
## maxnet end



#RF model Train
RF_model <- train(method = "RF", data = data)
auc(RF_model)
plotROC(RF_model)
#ggsave("./result/RF/RF_ROC_curve.jpg",dpi=300,width = 8,height = 5)


pred <- predict(RF_model,
                data = data,
                type = "cloglog")

map <- predict(RF_model,
               data = bio,
               type = "cloglog")
plotPred(map)
#ggsave("./result/RF/RF__train.png",dpi=300,width = 8,height = 5)

ths <- thresholds(RF_model, type = "cloglog")
ths
plotPA(map, th = ths[3, 2])
#ggsave("./result/RF/RF_presence_absence.png",dpi=300,width = 8,height = 5)
## RF end





### Split DataSet for Maxnet Training and testing datasets ###
library(zeallot)  # For unpacking assignment
c(train, test) %<-% trainValTest(data, #data from SWD object 
                                 test = 0.2, 
                                 only_presence = TRUE, 
                                 seed = 25)

maxnet_model <- train("Maxnet", data = train)
cat("Training auc: ", auc(maxnet_model))
cat("Testing auc: ", auc(maxnet_model, test = test))

plotROC(maxnet_model, test = test)
#ggsave("./result/maxnet/ROC_Training_testing.jpg",dpi=300,width = 8,height = 5)

#Permutation importance
vi_maxnet <- varImp(maxnet_model, permut = 5)
vi_maxnet
plotVarImp(vi_maxnet)
#ggsave("./result/maxnet/permutation_importance.jpg",dpi=300,width = 8,height = 5)

#Jackknife test for variable importance 
jk <- doJk(maxnet_model, metric = "auc", test = test)
jk
plotJk(jk, type = "train", ref = auc(maxnet_model))
#ggsave("./result/maxnet/Jackknife_train_variable_importance.jpg",dpi=300,width = 8,height = 5)


plotJk(jk, type = "test", ref = auc(maxnet_model, test = test))
#ggsave("./result/maxnet/Jackknife_test_variable_importance.jpg",dpi=300,width = 8,height = 5)

#Response curves for Bio19
plotResponse(maxnet_model, 
             var = "Bio19", 
             type = "cloglog", 
             only_presence = TRUE, 
             marginal = FALSE, 
             rug = TRUE)
#ggsave("./result/maxnet/response_curves_Bio19.jpg",dpi=300,width = 8,height = 5)


#logistic marginal response curve of biome that is a categorical variable, 
plotResponse(maxnet_model, 
             var = "Bio19", 
             type = "logistic", 
             only_presence = TRUE, 
             marginal = TRUE, 
             fun = mean, 
             color = "blue")

#ggsave("./result/maxnet/logistic_marginal_response_curve_Bio19.jpg",dpi=300,width = 8,height = 5)



plotResponse(maxnet_model, 
             var = "Bio19", 
             type = "cloglog", 
             only_presence = TRUE, 
             marginal = TRUE, 
             fun = mean, 
             rug = TRUE)





str(data)

output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
colnames(output) <- c("seed", "trainAUC", "testAUC")

set.seed(25)
seeds <- sample.int(1000, 10) # Create 10 different random seeds

for (i in seq_along(seeds)) { # Loop through the seeds
  c(train, test) %<-% trainValTest(data, 
                                   test = 0.2, 
                                   seed = seeds[i], 
                                   only_presence = TRUE) # Make the train/test split
  
  m <- train("Maxnet", 
             data = train) # train the model
  
  # Populate the output data.frame
  output[i, 1] <- seeds[i]
  output[i, 2] <- auc(m)
  output[i, 3] <- auc(m, test = test)
}
output
order(output[,"trainAUC"])


## Cross validation start ###
folds <- randomFolds(data, 
                     k = 4, 
                     only_presence = TRUE, 
                     seed = 25)

cv_model <- train("Maxnet", 
                  data = data, 
                  folds = folds)
cv_model

cat("Training AUC: ", auc(cv_model))
cat("Testing AUC: ", auc(cv_model, test = TRUE))

## Cross validation end ###




### Split DataSet for RF Training and testing datasets ###
library(zeallot)  # For unpacking assignment
c(train, test) %<-% trainValTest(data, #data from SWD object 
                                 test = 0.2, 
                                 only_presence = TRUE, 
                                 seed = 25)

RF_model <- train("RF", data = train)
cat("Training auc: ", auc(RF_model))
cat("Testing auc: ", auc(RF_model, test = test))

plotROC(RF_model, test = test)
ggsave("./result/RF/ROC_Training_testing.jpg",dpi=300,width = 8,height = 5)

str(data)

output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
colnames(output) <- c("seed", "trainAUC", "testAUC")

set.seed(25)
seeds <- sample.int(1000, 10) # Create 10 different random seeds

for (i in seq_along(seeds)) { # Loop through the seeds
  c(train, test) %<-% trainValTest(data, 
                                   test = 0.2, 
                                   seed = seeds[i], 
                                   only_presence = TRUE) # Make the train/test split
  
  m <- train("RF", 
             data = train) # train the model
  
  # Populate the output data.frame
  output[i, 1] <- seeds[i]
  output[i, 2] <- auc(m)
  output[i, 3] <- auc(m, test = test)
}
output
order(output[,"trainAUC"])


## Cross validation start ###
folds <- randomFolds(data, 
                     k = 4, 
                     only_presence = TRUE, 
                     seed = 25)

cv_model <- train("RF", 
                  data = data, 
                  folds = folds)
cv_model

cat("Training AUC: ", auc(cv_model))
cat("Testing AUC: ", auc(cv_model, test = TRUE))

## Cross validation end ###


























#########################################################
#########################################################
#################### Random Forest ######################
#########################################################
#########################################################
# Create the folds from the training dataset
folds <- randomFolds(data,
                     k = 4,
                     only_presence = TRUE,
                     seed = 25)

# Train the model
cv_model <- train("RF", data = data, folds = folds,
                  # mtry     = 4,
                  # ntree    = 500,
                  # nodesize = 1
)

# Check initial accuracy
auc(cv_model, test = TRUE)



#Permutation importance
vi_maxnet <- varImp(cv_model, permut = 5)
vi_maxnet
plotVarImp(vi_maxnet)
#ggsave("./result/maxnet/permutation_importance.jpg",dpi=300,width = 8,height = 5)

#Jackknife test for variable importance 
jk <- doJk(cv_model, metric = "auc", test = test)
jk
plotJk(jk, type = "train", ref = auc(cv_model))
#ggsave("./result/maxnet/Jackknife_train_variable_importance.jpg",dpi=300,width = 8,height = 5)


plotJk(jk, type = "test", ref = auc(cv_model, test = test))
#ggsave("./result/maxnet/Jackknife_test_variable_importance.jpg",dpi=300,width = 8,height = 5)

#Response curves for Bio19
plotResponse(cv_model, 
             var = "Bio19", 
             type = "cloglog", 
             only_presence = TRUE, 
             marginal = FALSE, 
             rug = TRUE)
#ggsave("./result/maxnet/response_curves_Bio19.jpg",dpi=300,width = 8,height = 5)


#logistic marginal response curve of biome that is a categorical variable, 
plotResponse(cv_model, 
             var = "Bio19", 
             type = "logistic", 
             only_presence = TRUE, 
             marginal = TRUE, 
             fun = mean, 
             color = "blue")

#ggsave("./result/maxnet/logistic_marginal_response_curve_Bio19.jpg",dpi=300,width = 8,height = 5)



plotResponse(cv_model, 
             var = "Bio19", 
             type = "cloglog", 
             only_presence = TRUE, 
             marginal = TRUE, 
             fun = mean, 
             rug = TRUE)













## hyperparameters Tune
# For random forest, we can tune three parameters, namely "mtry", "ntree" and "nodesize"

h_rf <- list(
  mtry      = seq(2,8,1),
  ntree     = seq(100,1000,10),
  nodesize  = 1
)

exp_8 <- randomSearch(
  cv_model,
  hypers = h_rf,
  metric = "auc",
  pop = 50,
  seed = 65466
)


# random search result

(res <- exp_8@results)
library(ggplot2)
ggplot()+
  geom_point(data = res,aes(x=1:50,y=test_AUC))+
  labs(x="model",y="AUC")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("./result/random_search_res.tiff",dpi=300,width = 8,height = 5)


exp_8@models[[1]]
# Train for Random Forest
predRF <- predict(exp_8@models[[1]], data = bio, type = "cloglog")
plotPred(predRF, lt = "Habitat\nsuitability",
         colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
ggsave("./result/predict_randomForest.jpg",height = 6, width = 6, dpi = 300)
