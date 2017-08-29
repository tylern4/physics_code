## My Functions

#Returns double: Takes confusion matrix and computes the sensativity
my_sensitivity<-function(mat){
  tp = mat[1,1]
  fn = mat[1,2]
  return(tp/(tp+fn))
}

#Returns double: Takes confusion matrix and computes the specificity
my_specificity<-function(mat){
  tn = mat[2,2]
  fp = mat[2,1]
  return(tn/(tn+fp))
}

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

load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

sim <- load_data("/Volumes/LaCiE/physics/e1d/v2/sim_csv/")
ana <- load_data("/Volumes/LaCiE/physics/e1d/v2/csv/")
## Part 1
#library(data.table)  
#files <- list.files(path = "/mnt/Bolden_data2/e1d/sim/singularity/",pattern = "sim_sing_07_18_2017.000*.csv")
#sim <- do.call(rbind,lapply(file,read.csv))

#Copy training data
#sim <- read.csv(file="/Volumes/LaCiE/physics/e1d/v2/sim_csv/sim_sing_07_18_2017.0000.csv", header=TRUE, sep=",",stringsAsFactors = TRUE)

#Convert column 2 to factor
sim[,1]<-as.factor(sim[,1])
sim[,2]<-as.factor(sim[,2])
sim[,3]<-as.factor(sim[,3])
sim[,4]<-as.factor(sim[,4])
sim[,5]<-as.factor(sim[,5])
sim[,6]<-as.factor(sim[,6])
sim[,7]<-as.factor(sim[,7])
sim[,8]<-as.factor(sim[,8])

ana[,1]<-as.factor(ana[,1])
ana[,2]<-as.factor(ana[,2])
ana[,3]<-as.factor(ana[,3])
ana[,4]<-as.factor(ana[,4])
ana[,5]<-as.factor(ana[,5])
ana[,6]<-as.factor(ana[,6])
ana[,7]<-as.factor(ana[,7])
ana[,8]<-as.factor(ana[,8])

#Create training and test sets of the data
test_set<-ana[,c(1,9,10,11,12,19,20,22,23,24)]
train_set<-sim[,c(1,9,10,11,12,19,20,22,23,24)]


## Part 2


#load library
library("e1071")

#Part 2.1
#train with training data
bayes_model<-naiveBayes(train_set[,-1], train_set[,1])
bayes_model


### Part 2.2

#predict and create table confusion matrix
bayes_train<-predict(bayes_model,train_set[,-1])
#confusion_mat_bayes_train<-table(bayes_train,train_set[,1])
#confusion_mat_bayes_train

#my_sensitivity(confusion_mat_bayes_train)
#my_specificity(confusion_mat_bayes_train)

### Part 2.3

#my_accuracy(confusion_mat_bayes_train)


### Part 2.4

#predict and create table confusion matrix
bayes_test<-predict(bayes_model,test_set[,-1])
#confusion_mat_bayes_test<-table(bayes_test,test_set[,1])
#confusion_mat_bayes_test

#my_sensitivity(confusion_mat_bayes_test)
#my_specificity(confusion_mat_bayes_test)
#print(paste("bayes:",my_accuracy(confusion_mat_bayes_test)))


## Part 3

#Load rpart and run classification function
# Part 3.1
library("rpart")
library("rpart.plot")
tree1<-rpart(id_0 ~ ., method = "class", data=train_set, control=rpart.control(minsplit=2, cp=0.02))

### Part 3.2
tree_train<-predict(tree1,train_set[,-1],type='class')
#confusion_mat_tree_train<-table(tree_train,train_set[,1])
#confusion_mat_tree_train

#my_sensitivity(confusion_mat_tree_train)
#my_specificity(confusion_mat_tree_train)

### Part 3.3

#my_accuracy(confusion_mat_tree_train)


### Part 3.4

tree_test<-suppressWarnings(predict(tree1,test_set[,-1],type='class'))
#confusion_mat_tree_test<-table(tree_test,test_set[,1])
#confusion_mat_tree_test

#my_sensitivity(confusion_mat_tree_test)
#my_specificity(confusion_mat_tree_test)
#print(paste("tree:",my_accuracy(confusion_mat_tree_test)))


## Part 4

plot(tree1, uniform=TRUE, main="Classification Tree Simulated Data")
text(tree1, use.n=TRUE, all=TRUE, cex=0.7)

library(rattle)
t <- fancyRpartPlot(tree1, main="Classification Tree Simulated Data")


library("ggplot2")

p_plot<-ggplot(data=sim, aes(sim$p)) + 
  geom_histogram(aes(y =..density..), alpha = .5, binwidth=5) + 
  geom_density() +
  xlim(0,6) +
  labs(title="Histogram for Momentum") +
  labs(x="Momentum", y="Count")

print(p_plot)

temp <- test_set[,-1]
temp$id_0 <- bayes_test
bayes_elec <- subset(temp, id_0 == 11)

temp <- test_set[,-1]
temp$id_0 <- tree_test
tree_elec <- subset(temp, id_0 == 11)

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
