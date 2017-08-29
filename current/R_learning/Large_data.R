#Returns double: Takes matrix and computes accuracy
my_accuracy<-function(mat){
  size = dim(mat)[1]
  total = sum(sum(mat))
  diagonal = 0
  for(i in 1:size){
    diagonal = diagonal + mat[i,i]
  }
  return(diagonal/total)
}

#Function to load all csv from path
load_data <- function(path) { 
  files <- dir(path, pattern = '*.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

load_data_chunck <- function(files) { 
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

#load trainging data
sim <- load_data("/Volumes/LaCiE/physics/e1d/v2/sim_csv/")

#Change some rows to factors for classification
sim[,1]<-as.factor(sim[,1])
sim[,2]<-as.factor(sim[,2])
sim[,3]<-as.factor(sim[,3])
sim[,4]<-as.factor(sim[,4])
sim[,5]<-as.factor(sim[,5])
sim[,6]<-as.factor(sim[,6])
sim[,7]<-as.factor(sim[,7])
sim[,8]<-as.factor(sim[,8])
sim[,13]<-as.factor(sim[,13])
sim[,14]<-as.factor(sim[,14])
sim[,15]<-as.factor(sim[,15])
sim[,16]<-as.factor(sim[,16])
sim[,21]<-as.factor(sim[,21])

#Create training and test sets of the data
#train_set<-sim[,c(1,9,10,11,12,19,20)]
train_set<-sim[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,24)]

#### Bayes Model
library("e1071")
bayes_model<-naiveBayes(train_set[,-1], train_set[,1])
bayes_train<-predict(bayes_model,train_set[,-1])
confusion_mat_bayes_train<-table(bayes_train,train_set[,1])
my_accuracy(confusion_mat_bayes_train)


### Decision Tree
library("rpart")
library("rpart.plot")
tree1<-rpart(id_0 ~ ., method = "class", data=train_set, control=rpart.control(minsplit=2, cp=0.02))
tree_train<-predict(tree1,train_set[,-1],type='class')
confusion_mat_tree_train<-table(tree_train,train_set[,1])
my_accuracy(confusion_mat_tree_train)
library(rattle)
fancyRpartPlot(tree1, main="Classification Tree Simulated Data")

### Get Predictions from Models

#load a chunck of the real data
ana_files <- dir("/Volumes/LaCiE/physics/e1d/v2/csv", pattern = '*.csv', full.names = TRUE)
###TODO loop here for each in ana_files
ana <- load_data_chunck(ana_files[1:10])
ana <- ana[ana$ec_0 < 7,]
ana <- ana[ana$gpart < 14,]
ana <- ana[ana$stat_0 >= 0,]
ana <- ana[ana$sc_0 < 10,]
ana <- ana[ana$dc_0 < 10,]
ana <- ana[ana$dc_stat < 6,]
ana <- ana[ana$cc_0 < 7,]

ana[,1]<-as.factor(ana[,1])
ana[,2]<-as.factor(ana[,2])
ana[,3]<-as.factor(ana[,3])
ana[,4]<-as.factor(ana[,4])
ana[,5]<-as.factor(ana[,5])
ana[,6]<-as.factor(ana[,6])
ana[,7]<-as.factor(ana[,7])
ana[,8]<-as.factor(ana[,8])
ana[,13]<-as.factor(ana[,13])
ana[,14]<-as.factor(ana[,14])
ana[,15]<-as.factor(ana[,15])
ana[,16]<-as.factor(ana[,16])
ana[,21]<-as.factor(ana[,21])

#test_set<-ana[,c(1,9,10,11,12,19,20,22,23,24)]
test_set<-ana[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)]

system.time(tree_test<-suppressWarnings(predict(tree1,test_set[,-1],type='class')))
system.time(bayes_test<-suppressWarnings(predict(bayes_model,test_set[,-1])))

temp <- test_set[,-1]
temp$id_0 <- bayes_test
bayes_elec <- subset(temp, id_0 == 11)

temp <- test_set[,-1]
temp$id_0 <- tree_test
tree_elec <- subset(temp, id_0 == 11)

library("ggplot2")
wq2 <- ggplot(bayes_elec, aes(W, Q2)) + 
  xlim(0, 3.5) +
  ylim(0,10) +
  stat_bin2d(bins = 500)

print(wq2)

wq2 <- ggplot(tree_elec, aes(W, Q2)) + 
  xlim(0, 3.5) +
  ylim(0,10) +
  stat_bin2d(bins = 500)

print(wq2)
