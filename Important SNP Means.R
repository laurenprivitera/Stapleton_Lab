dat <- read.csv(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/data/tidied/Random2.csv"))
set.seed(12345)
keepc <- sort(c(1,2,3,sample(4:dim(dat)[2],100)))
keepr <- sort(c(1,2,sample(3:dim(dat)[1],100)))
dat <- dat[keepr,keepc]

outp <- read.csv(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/RandomvQTL_LOD%2CPvals%2CEffectSizes-2-22-18.csv"))
write.csv(outp,file = "C://Users/Lauren/Documents/Github/Stapleton_Lab/randoutput.csv")
#this is reading in the data file
importantsnps <- read.csv(file = "C://Users/Lauren/Documents/Github/Stapleton_Lab/FamilyQTLoverlap.csv")
dim(importantsnps)
colnames(importantsnps)
importantsnps <- importantsnps[,-1]
vq <- importantsnps[which(importantsnps$vQTL_ID != ""),]
keepv <- sapply(unique(vq$vQTL_ID), function(x){
  grow <- which(vq$vQTL_ID == x)
  estgrow <- mean(vq$incre_new[grow])
  diff <- abs(vq$incre_new - estgrow)
  keep <- which(diff == min(diff))
  return(vq$markername[keep])
})

mq <- importantsnps[which(importantsnps$meanQTL_ID != ""),]

keepm <- sapply(unique(mq$meanQTL_ID), function(x){
  grow <- which(mq$meanQTL_ID == x)
  estgrow <- mean(mq$incre_new[grow])
  diff <- abs(mq$incre_new - estgrow)
  keep <- which(diff == min(diff))
  return(mq$markername[keep])
})
ggplot2::ggplot(data = )