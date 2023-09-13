################################################################################
# Script      : s2_build_model.R
# Description : Run a maxent model  
################################################################################

# clean R session
rm(list = ls()); gc() 

# We will use SDMtune package for building our models; so intall the package
#install.packages("SDMtune",dependencies = T)

# Image data
library(SDMtune)
library(terra)
library(geodata)

################################################################################
# Step 1: Read data
################################################################################

# the bio data sets are 19 predictors, each one image for one predictors
bio<- rast("./data/bio.grd")
bd_0 <- gadm(country = "Bangladesh",level = 0,path = "./data") %>% st_as_sf()

# we can make them as a data.frame file
# bio_df<- as.data.frame(bio,xy=T)

# read occurrence data but this time keep only longitude and latitude column
df <- read.csv("./data/occ.csv")[,c("longitude","latitude")]
names(df) <- c("x","y")


# Collect random absence points within Bangladesh
set.seed(11)
bg_df <- spatSample(
  bio, 
  1000, 
  "random", 
  cells=F, 
  xy=TRUE, 
  values=FALSE,
  as.df=TRUE,
  na.rm=T
  )

################################################################################
# Step 2: Prepare data for model
################################################################################

# We have 19 variables, but we need to choose variables that are not correlated
library(usdm)
cor <- vifcor(bio,th=0.70,size = 10000)
cor
bio <- bio[[c(2,3,9,15,18,19)]]


# SDMtune has very detailed vignette (instructions) on how to use
vignette("basic-use")

# Prepare a SWD object for model building
data <- prepareSWD(species="Chromolaena odorata",
                   p = df, #data.frame. the presence locations.
                   a = bg_df, #data.frame. absence/background locations.
                   env= bio, 
                   #rast environmental variables extract locations.
                   )

## save the presence and absense data                   
#swd2csv(data, file_name = "data.csv")

#df <- data.frame(data@data, weed_y = data@pa,data@data)
#write.csv(x = df,file = "./data/SWDObject.csv",row.names = TRUE)


#df <- data.frame(data@coords, weed_y = data@pa,data@data)
#write.csv(x = df,
#          file = "./data/my_model_data.csv",
#          row.names = FALSE)


################################################################################
# Step 3: Tuning model parameters
################################################################################

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

## We need to furthur configure the hyperparameters of our model
# For random forest, we can tune three parameters, namely "mtry", "ntree" and
# "nodesize"


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


################################################################################
# Resuts of random search
################################################################################
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

# ##Trail for hyperparameter tuning for RandomForest Model
# library(zeallot)  # For unpacking assignment
# c(train, val, test) %<-% trainValTest(data, 
#                                       val = 0.2, 
#                                       test = 0.2, 
#                                       only_presence = TRUE, 
#                                       seed = 61516)
# 
# RandomForest_tune <-train(method = "RF",data = train)
# 
# # Define the values for the regularization multiplier
# h <- list(mtry = c(1,2,3,4,5), 
#           nodesize = c(1,2),
#           ntree= c(100,200,300,400,500,600,700,800,900, 1000))
# 
# # Call the gridSearch function
# exp_1 <- gridSearch(RandomForest_tune, 
#                     hypers = h, 
#                     metric = "auc", 
#                     test = test)
# 
# results <- exp_1@results
# 

library(MVN)
data(iris)
setosa <- iris[1:50, 1:4]
mvn(setosa, subset = NULL, mvnTest = c("energy"))
qqplot(setosa$v1)

