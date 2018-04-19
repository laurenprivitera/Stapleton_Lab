install.packages("vqtl")
install.packages("qtl")
library("qtl")
library("vqtl")
#we also got rid of "(" in SNP ID rz44bd and rz574bc
dat <- read.csv(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/data/tidied/Random2.csv"))
set.seed(12345)
keepc <- sort(c(1,2,3,sample(4:dim(dat)[2],100)))
keepr <- sort(c(1,2,sample(3:dim(dat)[1],100)))
dat <- dat[keepr,keepc]
dat[1:2,1:3] <- ""
write.csv(dat,file = "C://Users/Lauren/Documents/Github/Stapleton_Lab/randcross.csv", row.names = F)
random <-read.cross(file = "C://Users/Lauren/Documents/Github/Stapleton_Lab/randcross.csv")
random <- drop.nullmarkers(random)
#scan with variance
random <- calc.genoprob(random)
routv <- scanonevar(cross = random,
                    mean.formula = height.in. ~ mean.QTL.add + mean.QTL.dom,
                    var.formula = 
                      ~ var.QTL.add + var.QTL.dom)
mean_var_plot_model_free(cross = random, phenotype.name = "height.in.", grouping.factor.names = "gpm705a")
library("dplyr")
effect.sizes = function (cross, phenotype.name, focal.groups = NULL, nuisance.groups = NULL, 
                         genotype.names = c("AA", "AB", "BB"), xlim = NULL, ylim = NULL, 
                         title = paste(phenotype.name, "by", paste(focal.groups, 
                                                                   collapse = ", ")), draw_ribbons = TRUE, se_line_size = 1, 
                         point_size = 1) 
{
  indiv.mean.estim <- indiv.mean.lb <- indiv.mean.ub <- "fake_global_for_CRAN"
  indiv.sd.estim <- indiv.sd.lb <- indiv.sd.ub <- "fake_global_for_CRAN"
  group.mean.estim <- group.mean.ub <- group.mean.lb <- "fake_global_for_CRAN"
  group.sd.estim <- group.sd.ub <- group.sd.lb <- "fake_global_for_CRAN"
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  modeling.df[[phenotype.name]] <- cross[["pheno"]][[phenotype.name]]
  marker.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.geno(cross = cross))], 
                    nuisance.groups[nuisance.groups %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.pheno(cross = cross))], 
                  nuisance.groups[nuisance.groups %in% colnames(qtl::pull.pheno(cross = cross))])
  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(x = qtl::pull.geno(cross = cross)[, 
                                                                           marker.name], labels = genotype.names)
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[[phen.name]])
  }
  modeling.df[["placeholder"]] <- NULL
  covar.form <- paste(focal.groups, collapse = "+")
  if (!is.null(nuisance.groups)) {
    covar.form <- paste(covar.form, "+", paste(nuisance.groups, 
                                               collapse = "+"))
  }
  mean.form <- paste(phenotype.name, "~", covar.form)
  var.form <- paste("~", covar.form)
  dglm.fit <- dglm::dglm(formula = stats::formula(mean.form), 
                         dformula = stats::formula(var.form), data = modeling.df)
  mean.pred <- stats::predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit
  sd.pred <- stats::predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  sd.estim <- sd.pred$fit/sd.pred$residual.scale
  sd.se <- sd.pred$se.fit
  indiv.prediction.tbl <- dplyr::bind_cols(stats::na.omit(modeling.df), 
                                           dplyr::data_frame(indiv.mean.estim = mean.estim, indiv.mean.lb = mean.estim - 
                                                               mean.se, indiv.mean.ub = mean.estim + mean.se, indiv.sd.estim = exp(sd.estim), 
                                                             indiv.sd.lb = exp(sd.estim - sd.se), indiv.sd.ub = exp(sd.estim + 
                                                                                                                      sd.se)))
  group.prediction.tbl <- indiv.prediction.tbl %>% dplyr::group_by_(.dots = c(focal.groups)) %>% 
    dplyr::summarise(group.mean.estim = mean(indiv.mean.estim), 
                     group.mean.lb = mean(indiv.mean.lb), group.mean.ub = mean(indiv.mean.ub), 
                     group.sd.estim = mean(indiv.sd.estim), group.sd.lb = mean(indiv.sd.lb), 
                     group.sd.ub = mean(indiv.sd.ub))
  return(group.prediction.tbl)
}

y = 1:length(routv$result$loc.name)
#effect sizes can not be computed for these 3 SNPs
rsizedf = sapply(y, function(x){
  tryCatch({
    print(x)
    tempm =  effect.sizes(cross = random,
                          phenotype.name = "height.in.",
                          genotype.names = c("AA","BB"),
                          focal.groups = routv$result$loc.name[x])
  }, error = function(e) message(e),
  finally = function(tempm){
    tempv = c(tempm[1,2:7],tempm[2,2:7])
    return(unlist(tempv))
  }
  )
})
rsizedf1 <- as.data.frame(matrix(rep(0,length(y)*12), ncol= 12))
sapply(1:length(rsizedf), function(x){
  print(x)
  if(!is.null(rsizedf[[x]])){
    rsizedf1[x,] <<- c(rsizedf[[x]][1,2:7],rsizedf[[x]][2:7])
  }
})
nall0 <-sapply(1:dim(rsizedf1)[1], function(x){
  !all(rsizedf1[x,] == 0)
})  
ditch <- which(nall0 == F)
rsizedf1 <- rsizedf1[-ditch,]
keep <- 1:3235; keep<- keep[-ditch]
routvdf<- data.frame(routv$result$loc.name[keep],
                     routv$result$pos[keep],
                     routv$result$mean.lod[keep],
                     routv$result$mean.asymp.p[keep],
                     routv$result$var.lod[keep],
                     routv$result$var.asymp.p[keep],
                     routv$result$joint.lod[keep],
                     routv$result$joint.asymp.p[keep])
routvdf = cbind(routvdf,rsizedf1)
colnames(routvdf) = c("SNP Name",
                      "Position (cM)",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value",
                      "A Mean Est",
                      "A Mean Lower Bound",
                      "A Mean Upper Bound",
                      "A Standard Deviation Est",
                      "A Standard Deviation Lower Bound",
                      "A Standard Deviation Upper Bound",
                      "B Mean Est",
                      "B Mean Lower Bound",
                      "B Mean Upper Bound",
                      "B Standard Deviation Est",
                      "B Standard Deviation Lower Bound",
                      "B Standard Deviation Upper Bound")

write.csv(routvdf, file = "RandomvQTL_LOD,Pvals,EffectSizes-2-22-18.csv")

#####now with the Family dataset####

family <-read.cross(file = "data/tidied/Family.csv")
family <- drop.nullmarkers(family)
#scan with variance
family <- calc.genoprob(family)
foutv <- scanonevar(cross = family,
                    mean.formula = PlantHeight ~ mean.QTL.add + mean.QTL.dom,
                    var.formula = ~ var.QTL.add + var.QTL.dom)
y = 1:length(foutv$result$loc.name)
#effect sizes can not be computed for these 3 SNPs
fsizedf = sapply(y, function(x){
  tryCatch({
    print(x)
    tempm =  effect.sizes(cross = family,
                          phenotype.name = "PlantHeight",
                          genotype.names = c("AA","BB"),
                          focal.groups = foutv$result$loc.name[x])
  }, error = function(e) message(e),
  finally = function(tempm){
    tempv = c(tempm[1,2:7],tempm[2,2:7])
    return(unlist(tempv))
  }
  )
})
fsizedf1 <- as.data.frame(matrix(rep(0,length(y)*12), ncol= 12))
sapply(1:length(fsizedf), function(x){
  print(x)
  if(!is.null(fsizedf[[x]])){
    fsizedf1[x,] <<- c(fsizedf[[x]][1,2:7],fsizedf[[x]][2:7])
  }
})
nall0 <-sapply(1:dim(fsizedf1)[1], function(x){
  !all(fsizedf1[x,] == 0)
})  
ditch <- which(nall0 == F)
fsizedf1 <- fsizedf1[-ditch,]
keep <- 1:3235; keep<- keep[-ditch]
foutvdf<- data.frame(foutv$result$loc.name[keep],
                     foutv$result$pos[keep],
                     foutv$result$mean.lod[keep],
                     foutv$result$mean.asymp.p[keep],
                     foutv$result$var.lod[keep],
                     foutv$result$var.asymp.p[keep],
                     foutv$result$joint.lod[keep],
                     foutv$result$joint.asymp.p[keep])
foutvdf = cbind(foutvdf,fsizedf1)
colnames(foutvdf) = c("SNP Name",
                      "Position (cM)",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value",
                      "A Mean Est",
                      "A Mean Lower Bound",
                      "A Mean Upper Bound",
                      "A Standard Deviation Est",
                      "A Standard Deviation Lower Bound",
                      "A Standard Deviation Upper Bound",
                      "B Mean Est",
                      "B Mean Lower Bound",
                      "B Mean Upper Bound",
                      "B Standard Deviation Est",
                      "B Standard Deviation Lower Bound",
                      "B Standard Deviation Upper Bound")

write.csv(foutvdf, file = "FamilyvQTL_LOD,Pvals,EffectSizes-2-22-18.csv")
# #trying to script fitqtls for all the markers
# varlist<- paste("c",1:10, sep = "")
# values <- rep(1,10)
# df<-data.frame(varlist,values)
# pval <- rep(2,3235)
# names <- rep("empty",3235)
# names[1:531] <- colnames(random$geno$`1`$data)
# names[532:879] <- colnames(random$geno$`2`$data)
# names[880:1263] <- colnames(random$geno$`3`$data)
# names[1264:1583] <- colnames(random$geno$`4`$data)
# names[1584:1908] <- colnames(random$geno$`5`$data)
# names[1909:2205] <- colnames(random$geno$`6`$data)
# names[2206:2475] <- colnames(random$geno$`7`$data)
# names[2476:2734] <- colnames(random$geno$`8`$data)
# names[2735:3003] <- colnames(random$geno$`9`$data)
# names[3004:3235] <- colnames(random$geno$`10`$data)
# for (x in 1:531) {
#   for (y in 1:10) {
#     if (length(random$geno[[y]]$data[1,]) >= x) {
#       df[y,2] <- x
#     }
#     else {
#       df[y,2] <- 1
#     }
#   }
#   poslist<- c(routv[colnames(random$geno$`1`$data)[df[1,2]],2],
#               routv[colnames(random$geno$`2`$data)[df[2,2]],2],
#               routv[colnames(random$geno$`3`$data)[df[3,2]],2],
#               routv[colnames(random$geno$`4`$data)[df[4,2]],2],
#               routv[colnames(random$geno$`5`$data)[df[5,2]],2],
#               routv[colnames(random$geno$`6`$data)[df[6,2]],2],
#               routv[colnames(random$geno$`7`$data)[df[7,2]],2],
#               routv[colnames(random$geno$`8`$data)[df[8,2]],2],
#               routv[colnames(random$geno$`9`$data)[df[9,2]],2],
#               routv[colnames(random$geno$`10`$data)[df[10,2]],2])
#   qtl <- makeqtl(random, chr = 1:10, pos = poslist, what = "prob")
#   fit <- fitqtl(random, qtl = qtl, get.ests = TRUE)
#   sum <- summary(fit)
#   pval[x] <- sum$result.drop[1,7]
#   pval[df[2,2]+531] <- sum$result.drop[2,7]
#   pval[df[3,2]+879] <- sum$result.drop[3,7]
#   pval[df[4,2]+1263] <- sum$result.drop[4,7]
#   pval[df[5,2]+1583] <- sum$result.drop[5,7]
#   pval[df[6,2]+1908] <- sum$result.drop[6,7]
#   pval[df[7,2]+2205] <- sum$result.drop[7,7]
#   pval[df[8,2]+2475] <- sum$result.drop[8,7]
#   pval[df[9,2]+2734] <- sum$result.drop[9,7]
#   pval[df[10,2]+3003] <- sum$result.drop[10,7]
# }
# allpos<-routv[,2]
# pval <- data.frame(names,allpos,pval)
# pval
# write.csv(pval, file = "RandomPvalues.csv")
# max(pval[,3])



#Creating a Manhattan plot
#   Specifically using the qqman package
library(qqman)
CHR = rep(1,length(pval[,1]))
plotframe <- data.frame(SNP = pval[,1], BP = pval[,2], P = pval[,3], CHR = CHR)
manhattan(plotframe)

#trying to actually run this correctly
random <-read.cross(dir = "~/PCR\ data/", file = "Random2.csv")
random <- drop.nullmarkers(random)
scan <- scanonevar(cross = random,
                   mean.formula = height.in. ~ mean.QTL.add + mean.QTL.dom,
                   var.formula = ~ var.QTL.add + var.QTL.dom,
                   chrs = 1)
perm = scanonevar.perm(scan,100, n.cores = 2)

#now doing it will all 10

scanall <- scanonevar(cross = random,
                      mean.formula = height.in. ~ mean.QTL.add + mean.QTL.dom,
                      var.formula = ~ var.QTL.add + var.QTL.dom)
perm = scanonevar.perm(scanall,100, n.cores = 2)

#finding important SNPs for mean, variance or both
#current functionality:
#Inputs, perm, threshold, type
#type legend

# For one scanonevar
# mean = mean.asymp.p
# var = var.asymp.p
# both = joint.asymp.p


# For perm
# mean = mean.empir.p
# var = var.empir.p
# both = joint.empir.p
SNP <- function(perm, type, thresh) {
  values <- perm$result[[type]]
  thresh = thresh
  newv = NULL
  for (x in values) {
    if (x <= thresh) {
      newv = c(newv,x)
    }
  } 
  return(newv)
}

#finding significant joint SNPS for different thresholds
sigj<-SNP(perm = perm, type = "joint.empir.p", thresh = .001)
sigom<- SNP(perm = routv, type = "mean.asymp.p", thresh = .005)

#matching this for SNP   index
SNP.match <- function(perm, sig, type) {
  markers = NULL
  for (x in sig) {
    markers = c(markers,match(x, perm$result[[type]]))
  }
  return(markers)
} 
#continuing example
matchj <- SNP.match(perm = perm, sig = sigj, type = "joint.empir.p")
matchom <- SNP.match(perm = routv, sig = sigom, type = "mean.asymp.p")
#matching indeces to names
SNP.name <- function(perm, match, type) {
  snp.names = NULL
  for (x in match) {
    snp.names = c(snp.names, perm$result$loc.name[x])
  }
  return(snp.names)
}
#example
namesj <- SNP.name(perm = perm, match = matchj, type = "joint.empir.p")
namesom <- SNP.name(perm = routv, match = matchom, type = "mean.empir.p")
#there seems to be issues with it adding markers we don't have
#let's try getting rid of these extra ones
SNP.scrub <- function(names) {
  scrubbed.names = names
  bad = grep("_",names)
  if (length(bad) == 0){
    return(scrubbed.names)
  } else{
    return(scrubbed.names[-bad])
  }
}
scrub <- SNP.scrub(names = namesj)
scrubom <- SNP.scrub(names = namesom)

#outputting it all as a nice list

#trying to plot these SNPs
mean_var_plot_model_based(cross = random,
                          phenotype.name = "height.in.",
                          genotype.names = c("AA","BB"),
                          focal.groups = scrub[2])
#exporting plots
x = 1 
for (x in 1:length(scrub)) {
  pdf(file = paste("VQTLplots/",scrub[x],".pdf",sep = ""))
  print(mean_var_plot_model_based(cross = random,
                                  phenotype.name = "height.in.",
                                  genotype.names = c("AA","BB"),
                                  focal.groups = scrub[x])
  )
  dev.off()
}
dev.cur()


#specifically looking for mean and variance effect size estimates
#looking to get it to work once by debugging mean_var_plot_based
mean_var_plot_model_based(cross = random,
                          phenotype.name = "height.in.",
                          genotype.names = c("AA","BB"),
                          focal.groups = scrub[1])
#using source code for mean_var_plot_model_based, we try to cut out the graph, and extract the estimates
library("dplyr")
effect.sizes = function (cross, phenotype.name, focal.groups = NULL, nuisance.groups = NULL, 
                         genotype.names = c("AA", "AB", "BB"), xlim = NULL, ylim = NULL, 
                         title = paste(phenotype.name, "by", paste(focal.groups, 
                                                                   collapse = ", ")), draw_ribbons = TRUE, se_line_size = 1, 
                         point_size = 1) 
{
  indiv.mean.estim <- indiv.mean.lb <- indiv.mean.ub <- "fake_global_for_CRAN"
  indiv.sd.estim <- indiv.sd.lb <- indiv.sd.ub <- "fake_global_for_CRAN"
  group.mean.estim <- group.mean.ub <- group.mean.lb <- "fake_global_for_CRAN"
  group.sd.estim <- group.sd.ub <- group.sd.lb <- "fake_global_for_CRAN"
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  modeling.df[[phenotype.name]] <- cross[["pheno"]][[phenotype.name]]
  marker.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.geno(cross = cross))], 
                    nuisance.groups[nuisance.groups %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.pheno(cross = cross))], 
                  nuisance.groups[nuisance.groups %in% colnames(qtl::pull.pheno(cross = cross))])
  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(x = qtl::pull.geno(cross = cross)[, 
                                                                           marker.name], labels = genotype.names)
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[[phen.name]])
  }
  modeling.df[["placeholder"]] <- NULL
  covar.form <- paste(focal.groups, collapse = "+")
  if (!is.null(nuisance.groups)) {
    covar.form <- paste(covar.form, "+", paste(nuisance.groups, 
                                               collapse = "+"))
  }
  mean.form <- paste(phenotype.name, "~", covar.form)
  var.form <- paste("~", covar.form)
  dglm.fit <- dglm::dglm(formula = stats::formula(mean.form), 
                         dformula = stats::formula(var.form), data = modeling.df)
  mean.pred <- stats::predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit
  sd.pred <- stats::predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  sd.estim <- sd.pred$fit/sd.pred$residual.scale
  sd.se <- sd.pred$se.fit
  indiv.prediction.tbl <- dplyr::bind_cols(stats::na.omit(modeling.df), 
                                           dplyr::data_frame(indiv.mean.estim = mean.estim, indiv.mean.lb = mean.estim - 
                                                               mean.se, indiv.mean.ub = mean.estim + mean.se, indiv.sd.estim = exp(sd.estim), 
                                                             indiv.sd.lb = exp(sd.estim - sd.se), indiv.sd.ub = exp(sd.estim + 
                                                                                                                      sd.se)))
  group.prediction.tbl <- indiv.prediction.tbl %>% dplyr::group_by_(.dots = c(focal.groups)) %>% 
    dplyr::summarise(group.mean.estim = mean(indiv.mean.estim), 
                     group.mean.lb = mean(indiv.mean.lb), group.mean.ub = mean(indiv.mean.ub), 
                     group.sd.estim = mean(indiv.sd.estim), group.sd.lb = mean(indiv.sd.lb), 
                     group.sd.ub = mean(indiv.sd.ub))
  return(group.prediction.tbl)
}
#this works
#this is an example of the first one
sizes = effect.sizes(cross = random,
                     phenotype.name = "height.in.",
                     genotype.names = c("AA","BB"),
                     focal.groups = scrub[1])
sizes

rsizedf <- data.frame(NULL)

y = 1:length(routv$result$loc.name)
y = y[-c(458,2482,2483)]
#effect sizes can not be computed for these 3 SNPs
rsizedf = sapply(y, function(x){
  tempm =  effect.sizes(cross = random,
                        phenotype.name = "height.in.",
                        genotype.names = c("AA","BB"),
                        focal.groups = routv$result$loc.name[x])
  tempv = c(tempm[1,2:7],tempm[2,2:7])
  return(unlist(tempv))
})
undebug(effect.sizes)
effect.sizes(cross = random,
             phenotype.name = "height.in.",
             genotype.names = c("AA","BB"),
             focal.groups = routv$result$loc.name[470])


routvdf<- data.frame(routv$result$loc.name,
                     routv$result$pos,
                     routv$result$mean.lod,
                     routv$result$mean.asymp.p,
                     routv$result$var.lod,
                     routv$result$var.asymp.p,
                     routv$result$joint.lod,
                     routv$result$joint.asymp.p)
#dropping the SNPs whose effect sizes could not be computed
routvdf = routvdf[-c(458,2482,2483),]
routvdf = cbind(routvdf,t(rsizedf))
colnames(routvdf) = c("SNP Name",
                      "Position (cM)",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value",
                      "A Mean Est",
                      "A Mean Lower Bound",
                      "A Mean Upper Bound",
                      "A Standard Deviation Est",
                      "A Standard Deviation Lower Bound",
                      "A Standard Deviation Upper Bound",
                      "B Mean Est",
                      "B Mean Lower Bound",
                      "B Mean Upper Bound",
                      "B Standard Deviation Est",
                      "B Standard Deviation Lower Bound",
                      "B Standard Deviation Upper Bound")

write.csv(routvdf, file = "RandomvQTL_LOD,Pvals,EffectSizes1.csv")