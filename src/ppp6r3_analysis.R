library(emmeans)
library(ggplot2)
library(ggsignif)
library(ggpubr)#mainly for stat_pvalue_manual
library(car)
library(patchwork)

#
options(stringsAsFactors = FALSE)
set.seed(8675309)
options("na.action")
#$na.action
#[1] "na.omit"

#read in the combined phenotype file

ppp6r3_data = read.csv("./data/ppp6r3_cleaned.csv")

ppp6r3_data$GENOTYPE = factor(ppp6r3_data$GENOTYPE, levels=c("WT","HET","MUT"))


# define plotting function
plot_lsmeans = function(lm, sex="A", y.lab = "Phenotype", main.pheno = "Phenotype", mult = c(0.002, 0.007, 0.005)){
  lff<-lsmeans(lm,"GENOTYPE")
  pval = as.data.frame(lsmeans(lff,pairwise~GENOTYPE, adj="tukey")$contrasts)
  pval = as.data.frame(pval$p.value)
  lsmeans = as.data.frame(lsmeans(lff,pairwise~GENOTYPE, adj="tukey")$lsmeans)
  
  colnames(pval)[1] = "pval"
  pval$group1 = c("WT","WT","HET")
  pval$group2 = c("HET","MUT","MUT")
  ciel = max(lsmeans$lsmean+lsmeans$SE)
  y.pos = c(ciel+mult[1]*ciel, ciel+mult[2]*ciel, ciel+mult[3]*ciel)
  pval$y.position = y.pos
  pval$sig = apply(pval, 1, function(x) ifelse(as.numeric(x[1])<=0.05, return("*"), return("")))
  pval$pval = formatC(pval$pval,format="g",digits = 3,flag="#")
  pval$pval_comp = paste0(pval$pval, pval$sig)
  #start plotting
  p = ggplot(lsmeans,aes(x=GENOTYPE)) + geom_errorbar(data = lsmeans, aes(ymin=lsmean-SE, ymax=lsmean+SE),color="blue",size=3,width=0.11)
  p = p + geom_point(aes(y=lsmean),size=5)
  p = p + theme_bw()
  if(sex == "A"){
    samp = table(lm$model$SEX, lm$model$GENOTYPE)
    lsmeans$samp_size = colSums(samp)
    p = p + scale_x_discrete(labels=c(paste0("WT \n(n=",lsmeans$samp_size[1],")"),paste0("HET \n(n=",lsmeans$samp_size[2],")"), paste0("MUT \n(n=",lsmeans$samp_size[3],")")), expand = c(0,0.1)) 
    p = p + ylab(y.lab)+ggtitle(paste0(main.pheno, ", sex-combined"))
  } else if(sex == "F"){
    samp = table(lm$model$GENOTYPE)
    lsmeans$samp_size = samp
    p = p + scale_x_discrete(labels=c(paste0("WT \n(n=",lsmeans$samp_size[1],")"),paste0("HET \n(n=",lsmeans$samp_size[2],")"), paste0("MUT \n(n=",lsmeans$samp_size[3],")")), expand = c(0,0.1)) 
    p = p + ylab(y.lab)+ggtitle(paste0(main.pheno, ", females"))
  } else if(sex == "M"){
    samp = table(lm$model$GENOTYPE)
    lsmeans$samp_size = samp
    p = p + scale_x_discrete(labels=c(paste0("WT \n(n=",lsmeans$samp_size[1],")"),paste0("HET \n(n=",lsmeans$samp_size[2],")"), paste0("MUT \n(n=",lsmeans$samp_size[3],")")), expand = c(0,0.1)) 
    p = p + ylab(y.lab)+ggtitle(paste0(main.pheno, ", males"))
  }
  
  p = p + stat_pvalue_manual(data = pval, label="pval_comp", y.position = "y.position", face="bold", bracket.size=1.1, tip.length = 0.02, label.size = 4.5)
  p = p + theme(plot.title = element_text(hjust = 0.5, size=12))+ theme(axis.text = element_text(size=12, color = "black"), axis.title = element_text(size=12))
  p = p + xlab("Genotype")
  
  
  
  return(p)
}

###
###
###
#####LENGTH AND WEIGHT ANALYSIS####
###
###

# Length
lm.length <-lm(length~GENOTYPE+SEX+sac.weight+age,data=ppp6r3_data)
Anova(lm.length)
lff<-lsmeans(lm.length,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
pval = as.data.frame(lsmeans(lff,pairwise~GENOTYPE, adj="tukey")$contrasts)
lff = as.data.frame(lff)
table(lm.length$model$SEX, lm.length$model$GENOTYPE)
p1 = plot_lsmeans(lm.length, sex = "A", y.lab = "Nose-Anus Length (cm)", main.pheno = "Nose-Anus Length")
ggsave(plot = p,filename = "weight_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)
###


# Weight

lm.weight <-lm(sac.weight~GENOTYPE+SEX+age,data=ppp6r3_data)
Anova(lm.weight)
lff<-lsmeans(lm.weight,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.weight$model$SEX, lm.weight$model$GENOTYPE)
p = plot_lsmeans(lm.weight, sex = "A", y.lab = "Weight (g)", main.pheno = "Weight", mult = c(0.004, 0.01, 0.02))
ggsave(plot = p,filename = "weight_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)

####

####Caliper-based analyses####
# ML, FL, AP
###
###
###

#removed FL with no condyles, and samples with no FL/AP/ML
caliper_data = read.csv(file = "./data/ppp6r3_caliper_cleaned.csv")

caliper_data$GENOTYPE = factor(caliper_data$GENOTYPE, levels=c("WT","HET","MUT"))


###
lm.AP <-lm(AP~GENOTYPE+SEX+sac.weight+age,data=caliper_data)
Anova(lm.AP)
lff<-lsmeans(lm.AP,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.AP$model$SEX, lm.AP$model$GENOTYPE)
#car::vif(lm.AP)
p7 = plot_lsmeans(lm.AP, sex = "A", y.lab = "Anterior-Posterior Width (mm)", main.pheno = "Anterior-Posterior Femoral Width", mult = c(0.003, 0.015, 0.008))

ggsave(plot = p7,filename = "AP_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)


#################################################################
lm.ML <-lm(ML~GENOTYPE+SEX+sac.weight+age,data=caliper_data)
Anova(lm.ML)
lff<-lsmeans(lm.ML,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.ML$model$SEX, lm.ML$model$GENOTYPE)
p10 = plot_lsmeans(lm.ML, sex = "A", y.lab = "Medial-Lateral Width (mm)", main.pheno = "Medial-Lateral Femoral Width", mult = c(0.003, 0.015, 0.008))
ggsave(plot = p10,filename = "ML_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)

#################################################################
lm.FL <-lm(FL~GENOTYPE+SEX+sac.weight+age,data=caliper_data_cond)
Anova(lm.FL)
lff<-lsmeans(lm.FL,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.FL$model$SEX, lm.FL$model$GENOTYPE)
p13 = plot_lsmeans(lm.FL, sex = "A", y.lab = "Femoral Length (mm)", main.pheno = "Femoral Length", mult = c(0.003, 0.015, 0.008))
ggsave(plot = p13,filename = "FL_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)
####DXA analyses####



#NOTE: THREE OF THE DATAPOINTS ARE MISSING GENOTYPES, SO REMOVED
ls_bmd = read.csv("./data/ppp6r3_DXA_spine.csv")
fn_bmd = read.csv("./data/ppp6r3_DXA_femur.csv")

ls_bmd$GENOTYPE = factor(ls_bmd$GENOTYPE, levels=c("WT","HET","MUT"))

fn_bmd$GENOTYPE = factor(fn_bmd$GENOTYPE, levels=c("WT","HET","MUT"))


#################################################################
lm.BMD <-lm(BMD~GENOTYPE+SEX+sac.weight+age+CenterRectX+CenterRectY,data=ls_bmd)
Anova(lm.BMD)
lff<-lsmeans(lm.BMD,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.BMD$model$SEX, lm.BMD$model$GENOTYPE)
(p1 = plot_lsmeans(lm.BMD, sex = "A", y.lab = "areal Bone Mineral Density - Lumbar Spine (g HA/cm2)", main.pheno = "areal Bone Mineral Density - Lumbar Spine (g HA/cm2)", mult = c(0.003, 0.015, 0.009)))
ggsave(plot = p1,filename = "ls_aBMD_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)


########################################### REPEAT FOR FEMURS ###################################################
lm.BMD <-lm(BMD~GENOTYPE+SEX+sac.weight+age+CenterRectX+CenterRectY,data=fn_bmd)
Anova(lm.BMD)
lff<-lsmeans(lm.BMD,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.BMD$model$SEX, lm.BMD$model$GENOTYPE)
(p1 = plot_lsmeans(lm.BMD, sex = "A", y.lab = "areal Bone Mineral Density - Femur (g HA/cm2)", main.pheno = "areal Bone Mineral Density - Femur (g HA/cm2)", mult = c(0.003, 0.015, 0.009)))
ggsave(plot = p1,filename = "fn_aBMD_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)

#################################################################

####microCT analyses####
#analyze: "Tb.BV.TV", "Tb.vBMD", "Tb.vTMD", "Tb.Th", "Tb.Sp", "Tb.N"


uCT = read.csv("./data/ppp6r3_uCT_spine.csv")

#9 uCT samples thrown out because no weight info recorded.

uCT$GENOTYPE = factor(uCT$GENOTYPE, levels=c("WT","HET","MUT"))

lm.Tb.BV.TV <-lm(Tb.BV.TV~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.BV.TV)
lff<-lsmeans(lm.Tb.BV.TV,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.BV.TV$model$SEX, lm.Tb.BV.TV$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.BV.TV, sex = "A", y.lab = "Trabecular BV/TV", main.pheno = "Trabecular Bone Volume Fraction - Lumbar Spine", mult = c(0.007, 0.025, 0.015)))
ggsave(plot = p1,filename = "tbbvtv_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)

#################################################################
lm.Tb.vBMD <-lm(Tb.vBMD~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.vBMD)
lff<-lsmeans(lm.Tb.vBMD,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.vBMD$model$SEX, lm.Tb.vBMD$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.vBMD, sex = "A", y.lab = "Tb.vBMD (mgHA/cm3)", main.pheno = "Trabecular vBMD - Lumbar Spine (mgHA/cm3)", mult = c(0.007, 0.025, 0.015)))
ggsave(plot = p1,filename = "vBMD_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)


####
lm.Tb.vTMD <-lm(Tb.vTMD~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.vTMD)
lff<-lsmeans(lm.Tb.vTMD,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.vTMD$model$SEX, lm.Tb.vTMD$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.vTMD, sex = "A", y.lab = "Tb.TMD (mgHA/cm3)", main.pheno = "Trabecular TMD - Lumbar Spine (mgHA/cm3)", mult = c(0.002, 0.007, 0.004)))
ggsave(plot = p1,filename = "TMD_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)

####
lm.Tb.Sp <-lm(Tb.Sp~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.Sp)
lff<-lsmeans(lm.Tb.Sp,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.Sp$model$SEX, lm.Tb.Sp$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.Sp, sex = "A", y.lab = "Trabecular Separation", main.pheno = "Trabecular Separation - Lumbar Spine", mult = c(0.01, 0.03, 0.02)))
ggsave(plot = p1,filename = "TbSp_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)


####
lm.Tb.Th <-lm(Tb.Th~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.Th)
lff<-lsmeans(lm.Tb.Th,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.Th$model$SEX, lm.Tb.Th$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.Th, sex = "A", y.lab = "Trabecular Thickness", main.pheno = "Trabecular Thickness - Lumbar Spine", mult = c(0.01, 0.035, 0.02)))
ggsave(plot = p1,filename = "TbTh_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)

####
lm.Tb.N <-lm(Tb.N~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.N)
lff<-lsmeans(lm.Tb.N,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.N$model$SEX, lm.Tb.N$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.N, sex = "A", y.lab = "Trabecular Number", main.pheno = "Trabecular Number - Lumbar Spine", mult = c(0.01, 0.035, 0.02)))
ggsave(plot = p1,filename = "TbN_c.pdf", device = "pdf",path = "./",dpi = 300,width = 5,height = 5)






##Add in the Raman Data###
ppp6r3_s1 = read.csv("./data//PPP6R3_mice_s1.csv")
ppp6r3_s1 = ppp6r3_s1[!is.na(ppp6r3_s1$MOUSE.ID),]
#experimental mouse sheet
ppp6r3_s2 = read.csv("./data/PPP6R3_mice_s2.csv")
ppp6r3_s2 = ppp6r3_s2[!is.na(ppp6r3_s2$Mouse),]

ppp6r3_merged = merge(ppp6r3_s1,ppp6r3_s2, by.x ="MOUSE.ID", by.y="Mouse")



#TAKE RELEVANT COLS
ppp6r3_merged = ppp6r3_merged[,c("MOUSE.ID","SEX","GENOTYPE","length","sac.weight")]


#
#fix from HET? to HET


ppp6r3_merged[which(ppp6r3_merged$MOUSE.ID=="163"),"GENOTYPE"] = "HET"

#4 samples have unknown genotypes
table(ppp6r3_merged$GENOTYPE)

#remove mice with no genotype
ppp6r3_merged = ppp6r3_merged[-which((ppp6r3_merged$GENOTYPE == "")),]
ppp6r3_merged = ppp6r3_merged[-which((ppp6r3_merged$GENOTYPE == "?")),]




raman_femur = read.csv("data/Raman/raman_femur.csv", stringsAsFactors = F)

min.mat = aggregate(Min.Mat ~ Mouse.., data = raman_femur, FUN=mean) 
carb.phos = aggregate(Carb.Phos ~ Mouse.., data = raman_femur, FUN=mean) 
crystallinity = aggregate(Crystallinity ~ Mouse.., data = raman_femur, FUN=mean) 

raman_femur_means = cbind(min.mat, carb.phos, crystallinity)
raman_femur_means = raman_femur_means[,c(1,2,4,6)]

raman_femur_means = merge(raman_femur_means, ppp6r3_merged, by.x = "Mouse..", by.y = "MOUSE.ID")
raman_femur_means_F = raman_femur_means[which(raman_femur_means$SEX == "F"),]
raman_femur_means_M = raman_femur_means[which(raman_femur_means$SEX == "M"),]

raman_femur_means$GENOTYPE = factor(raman_femur_means$GENOTYPE, levels=c("WT","HET","MUT"))
raman_femur_means_F$GENOTYPE = factor(raman_femur_means_F$GENOTYPE, levels=c("WT","HET","MUT"))
raman_femur_means_M$GENOTYPE = factor(raman_femur_means_M$GENOTYPE, levels=c("WT","HET","MUT"))

raman_spine = read.csv("data/Raman/raman_spine.csv", stringsAsFactors = F)

min.mat = aggregate(Min.Mat ~ Mouse.., data = raman_spine, FUN=mean) 
carb.phos = aggregate(Carb.Phos ~ Mouse.., data = raman_spine, FUN=mean) 
crystallinity = aggregate(Crystallinity ~ Mouse.., data = raman_spine, FUN=mean) 

raman_spine_means = cbind(min.mat, carb.phos, crystallinity)
raman_spine_means = raman_spine_means[,c(1,2,4,6)]

raman_spine_means = merge(raman_spine_means, ppp6r3_merged, by.x = "Mouse..", by.y = "MOUSE.ID")
raman_spine_means_F = raman_spine_means[which(raman_spine_means$SEX == "F"),]
raman_spine_means_M = raman_spine_means[which(raman_spine_means$SEX == "M"),]

raman_spine_means$GENOTYPE = factor(raman_spine_means$GENOTYPE, levels=c("WT","HET","MUT"))
raman_spine_means_F$GENOTYPE = factor(raman_spine_means_F$GENOTYPE, levels=c("WT","HET","MUT"))
raman_spine_means_M$GENOTYPE = factor(raman_spine_means_M$GENOTYPE, levels=c("WT","HET","MUT"))

####################################################################################################
min.mat <-lm(Min.Mat ~ GENOTYPE + SEX, data=raman_femur_means)
Anova(min.mat) 
lff<-lsmeans(min.mat,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p1 = plot_lsmeans(min.mat, sex = "A", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Femur", mult = c(0.01, 0.035, 0.02)))

min.mat_F <-lm(Min.Mat ~ GENOTYPE, data=raman_femur_means_F)
Anova(min.mat_F) 
lff<-lsmeans(min.mat_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p2 = plot_lsmeans(min.mat_F, sex = "F", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Femur", mult = c(0.009, 0.03, 0.022)))

min.mat_M <-lm(Min.Mat ~ GENOTYPE, data=raman_femur_means_M)
Anova(min.mat_M) 
lff<-lsmeans(min.mat_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p3 = plot_lsmeans(min.mat_M, sex = "M", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Femur", mult = c(0.01, 0.035, 0.02)))


###
carb.phos <-lm(Carb.Phos ~ GENOTYPE + SEX, data=raman_femur_means)
Anova(carb.phos) 
lff<-lsmeans(carb.phos,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p4 = plot_lsmeans(carb.phos, sex = "A", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Femur", mult = c(0.01, 0.035, 0.02)))

carb.phos_F <-lm(Carb.Phos ~ GENOTYPE, data=raman_femur_means_F)
Anova(carb.phos_F) 
lff<-lsmeans(carb.phos_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p5 = plot_lsmeans(carb.phos_F, sex = "F", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Femur", mult = c(0.01, 0.035, 0.02)))

carb.phos_M <-lm(Carb.Phos ~ GENOTYPE, data=raman_femur_means_M)
Anova(carb.phos_M) 
lff<-lsmeans(carb.phos_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p6 = plot_lsmeans(carb.phos_M, sex = "M", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Femur", mult = c(0.015, 0.025, 0.005)))

##
crys <-lm(Crystallinity ~ GENOTYPE + SEX, data=raman_femur_means)
Anova(crys) 
lff<-lsmeans(crys,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p7 = plot_lsmeans(crys, sex = "A", y.lab = "Crystallinity", main.pheno = "Crystallinity - Femur", mult = c(0.001, 0.003, 0.005)))

crys_F <-lm(Crystallinity ~ GENOTYPE, data=raman_femur_means_F)
Anova(crys_F) 
lff<-lsmeans(crys_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p8 = plot_lsmeans(crys_F, sex = "F", y.lab = "Crystallinity", main.pheno = "Crystallinity - Femur", mult = c(0.0007, 0.003, 0.007)))

crys_M <-lm(Crystallinity ~ GENOTYPE, data=raman_femur_means_M)
Anova(crys_M) 
lff<-lsmeans(crys_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p9 = plot_lsmeans(crys_M, sex = "M", y.lab = "Crystallinity", main.pheno = "Crystallinity - Femur", mult = c(0.004, 0.013, 0.008)))

p = (p1|p2|p3)/(p4|p5|p6)/(p7|p8|p9)

ggsave(plot = p,filename = "mean_raman_femur.pdf", device = "pdf",path = "C://Users/basel/Desktop/",dpi = 300,width = 15,height = 20)


## SPINES ##

min.mat <-lm(Min.Mat ~ GENOTYPE + SEX, data=raman_spine_means)
Anova(min.mat) 
lff<-lsmeans(min.mat,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p10 = plot_lsmeans(min.mat, sex = "A", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Spine", mult = c(0.01, 0.035, 0.02)))

min.mat_F <-lm(Min.Mat ~ GENOTYPE, data=raman_spine_means_F)
Anova(min.mat_F) 
lff<-lsmeans(min.mat_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p11 = plot_lsmeans(min.mat_F, sex = "F", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Spine", mult = c(0.03, 0.07, 0.045)))

min.mat_M <-lm(Min.Mat ~ GENOTYPE, data=raman_spine_means_M)
Anova(min.mat_M) 
lff<-lsmeans(min.mat_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p12 = plot_lsmeans(min.mat_M, sex = "M", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Spine", mult = c(0.01, 0.035, 0.02)))

######
carb.phos <-lm(Carb.Phos ~ GENOTYPE + SEX, data=raman_spine_means)
Anova(carb.phos) 
lff<-lsmeans(carb.phos,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p13 = plot_lsmeans(carb.phos, sex = "A", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Spine", mult = c(0.00005, 0.012, 0.005)))

carb.phos_F <-lm(Carb.Phos ~ GENOTYPE, data=raman_spine_means_F)
Anova(carb.phos_F) 
lff<-lsmeans(carb.phos_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p14 = plot_lsmeans(carb.phos_F, sex = "F", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Spine", mult = c(0.0001, 0.035, 0.02)))

carb.phos_M <-lm(Carb.Phos ~ GENOTYPE, data=raman_spine_means_M)
Anova(carb.phos_M) 
lff<-lsmeans(carb.phos_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p15 = plot_lsmeans(carb.phos_M, sex = "M", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Spine", mult = c(0.012, 0.035, 0.02)))

######
crys <-lm(Crystallinity ~ GENOTYPE + SEX, data=raman_spine_means)
Anova(crys) 
lff<-lsmeans(crys,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p16 = plot_lsmeans(crys, sex = "A", y.lab = "Crystallinity", main.pheno = "Crystallinity - Spine", mult = c(0.003, 0.007, 0.005)))

crys_F <-lm(Crystallinity ~ GENOTYPE, data=raman_spine_means_F)
Anova(crys_F) 
lff<-lsmeans(crys_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p17 = plot_lsmeans(crys_F, sex = "F", y.lab = "Crystallinity", main.pheno = "Crystallinity - Spine", mult = c(0.001, 0.005, 0.003)))

crys_M <-lm(Crystallinity ~ GENOTYPE, data=raman_spine_means_M)
Anova(crys_M) 
lff<-lsmeans(crys_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p18 = plot_lsmeans(crys_M, sex = "M", y.lab = "Crystallinity", main.pheno = "Crystallinity - Spine", mult = c(0.005, 0.015, 0.01)))

p = (p10|p11|p12)/(p13|p14|p15)/(p16|p17|p18)

ggsave(plot = p,filename = "mean_raman_spine.pdf", device = "pdf",path = "C://Users/basel/Desktop/",dpi = 300,width = 15,height = 20)


##################################### Raman SDs ###

raman_femur = read.csv("data/Raman/raman_femur.csv", stringsAsFactors = F)

min.mat = aggregate(Min.Mat ~ Mouse.., data = raman_femur, FUN=sd) 
carb.phos = aggregate(Carb.Phos ~ Mouse.., data = raman_femur, FUN=sd) 
crystallinity = aggregate(Crystallinity ~ Mouse.., data = raman_femur, FUN=sd) 

raman_femur_sd = cbind(min.mat, carb.phos, crystallinity)
raman_femur_sd = raman_femur_sd[,c(1,2,4,6)]

raman_femur_sd = merge(raman_femur_sd, ppp6r3_merged, by.x = "Mouse..", by.y = "MOUSE.ID")
raman_femur_sd_F = raman_femur_sd[which(raman_femur_sd$SEX == "F"),]
raman_femur_sd_M = raman_femur_sd[which(raman_femur_sd$SEX == "M"),]

raman_femur_sd$GENOTYPE = factor(raman_femur_sd$GENOTYPE, levels=c("WT","HET","MUT"))
raman_femur_sd_F$GENOTYPE = factor(raman_femur_sd_F$GENOTYPE, levels=c("WT","HET","MUT"))
raman_femur_sd_M$GENOTYPE = factor(raman_femur_sd_M$GENOTYPE, levels=c("WT","HET","MUT"))


raman_spine = read.csv("data/Raman/raman_spine.csv", stringsAsFactors = F)

min.mat = aggregate(Min.Mat ~ Mouse.., data = raman_spine, FUN=sd) 
carb.phos = aggregate(Carb.Phos ~ Mouse.., data = raman_spine, FUN=sd) 
crystallinity = aggregate(Crystallinity ~ Mouse.., data = raman_spine, FUN=sd) 

raman_spine_sd = cbind(min.mat, carb.phos, crystallinity)
raman_spine_sd = raman_spine_sd[,c(1,2,4,6)]

raman_spine_sd = merge(raman_spine_sd, ppp6r3_merged, by.x = "Mouse..", by.y = "MOUSE.ID")
raman_spine_sd_F = raman_spine_sd[which(raman_spine_sd$SEX == "F"),]
raman_spine_sd_M = raman_spine_sd[which(raman_spine_sd$SEX == "M"),]

raman_spine_sd$GENOTYPE = factor(raman_spine_sd$GENOTYPE, levels=c("WT","HET","MUT"))
raman_spine_sd_F$GENOTYPE = factor(raman_spine_sd_F$GENOTYPE, levels=c("WT","HET","MUT"))
raman_spine_sd_M$GENOTYPE = factor(raman_spine_sd_M$GENOTYPE, levels=c("WT","HET","MUT"))



####################################################################################################
min.mat <-lm(Min.Mat ~ GENOTYPE + SEX, data=raman_femur_sd)
Anova(min.mat) 
lff<-lsmeans(min.mat,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p1 = plot_lsmeans(min.mat, sex = "A", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Femur", mult = c(0.011, 0.038, 0.025)))

min.mat_F <-lm(Min.Mat ~ GENOTYPE, data=raman_femur_sd_F)
Anova(min.mat_F) 
lff<-lsmeans(min.mat_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p2 = plot_lsmeans(min.mat_F, sex = "F", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Femur", mult = c(0.05, 0.11, 0.07)))

min.mat_M <-lm(Min.Mat ~ GENOTYPE, data=raman_femur_sd_M)
Anova(min.mat_M) 
lff<-lsmeans(min.mat_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p3 = plot_lsmeans(min.mat_M, sex = "M", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Femur", mult = c(0.05, 0.11, 0.07)))

###
carb.phos <-lm(Carb.Phos ~ GENOTYPE + SEX, data=raman_femur_sd)
Anova(carb.phos) 
lff<-lsmeans(carb.phos,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p4 = plot_lsmeans(carb.phos, sex = "A", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Femur",mult = c(0.05, 0.11, 0.07)))

carb.phos_F <-lm(Carb.Phos ~ GENOTYPE, data=raman_femur_sd_F)
Anova(carb.phos_F) 
lff<-lsmeans(carb.phos_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p5 = plot_lsmeans(carb.phos_F, sex = "F", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Femur", mult = c(0.055, 0.15, 0.1)))

carb.phos_M <-lm(Carb.Phos ~ GENOTYPE, data=raman_femur_sd_M)
Anova(carb.phos_M) 
lff<-lsmeans(carb.phos_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p6 = plot_lsmeans(carb.phos_M, sex = "M", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Femur", mult = c(0.055, 0.15, 0.1)))

##
crys <-lm(Crystallinity ~ GENOTYPE + SEX, data=raman_femur_sd)
Anova(crys) 
lff<-lsmeans(crys,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p7 = plot_lsmeans(crys, sex = "A", y.lab = "Crystallinity", main.pheno = "Crystallinity - Femur", mult = c(0.025, 0.1, 0.07)))

crys_F <-lm(Crystallinity ~ GENOTYPE, data=raman_femur_sd_F)
Anova(crys_F) 
lff<-lsmeans(crys_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p8 = plot_lsmeans(crys_F, sex = "F", y.lab = "Crystallinity", main.pheno = "Crystallinity - Femur", mult = c(0.025, 0.1, 0.07)))

crys_M <-lm(Crystallinity ~ GENOTYPE, data=raman_femur_sd_M)
Anova(crys_M) 
lff<-lsmeans(crys_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p9 = plot_lsmeans(crys_M, sex = "M", y.lab = "Crystallinity", main.pheno = "Crystallinity - Femur", mult = c(0.025, 0.1, 0.07)))

p = (p1|p2|p3)/(p4|p5|p6)/(p7|p8|p9)

ggsave(plot = p,filename = "var_raman_femur.pdf", device = "pdf",path = "C://Users/basel/Desktop/",dpi = 300,width = 15,height = 20)


## SPINES ##

min.mat <-lm(Min.Mat ~ GENOTYPE + SEX, data=raman_spine_sd)
Anova(min.mat) 
lff<-lsmeans(min.mat,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p10 = plot_lsmeans(min.mat, sex = "A", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Spine", mult = c(0.01, 0.035, 0.02)))

min.mat_F <-lm(Min.Mat ~ GENOTYPE, data=raman_spine_sd_F)
Anova(min.mat_F) 
lff<-lsmeans(min.mat_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p11 = plot_lsmeans(min.mat_F, sex = "F", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Spine", mult = c(0.03, 0.07, 0.045)))

min.mat_M <-lm(Min.Mat ~ GENOTYPE, data=raman_spine_sd_M)
Anova(min.mat_M) 
lff<-lsmeans(min.mat_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p12 = plot_lsmeans(min.mat_M, sex = "M", y.lab = "Mineral:Matrix Ratio", main.pheno = "Mineral:Matrix Ratio - Spine", mult = c(0.01, 0.035, 0.02)))

######
carb.phos <-lm(Carb.Phos ~ GENOTYPE + SEX, data=raman_spine_sd)
Anova(carb.phos) 
lff<-lsmeans(carb.phos,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p13 = plot_lsmeans(carb.phos, sex = "A", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Spine", mult = c(0.025, 0.1, 0.07)))

carb.phos_F <-lm(Carb.Phos ~ GENOTYPE, data=raman_spine_sd_F)
Anova(carb.phos_F) 
lff<-lsmeans(carb.phos_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p14 = plot_lsmeans(carb.phos_F, sex = "F", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Spine", mult = c(0.025, 0.1, 0.05)))

carb.phos_M <-lm(Carb.Phos ~ GENOTYPE, data=raman_spine_sd_M)
Anova(carb.phos_M) 
lff<-lsmeans(carb.phos_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p15 = plot_lsmeans(carb.phos_M, sex = "M", y.lab = "Carbonate:Phosphate Ratio", main.pheno = "Carbonate:Phosphate Ratio - Spine", mult = c(0.025, 0.1, 0.07)))

######
crys <-lm(Crystallinity ~ GENOTYPE + SEX, data=raman_spine_sd)
Anova(crys) 
lff<-lsmeans(crys,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p16 = plot_lsmeans(crys, sex = "A", y.lab = "Crystallinity", main.pheno = "Crystallinity - Spine", mult = c(0.025, 0.1, 0.07)))

crys_F <-lm(Crystallinity ~ GENOTYPE, data=raman_spine_sd_F)
Anova(crys_F) 
lff<-lsmeans(crys_F,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p17 = plot_lsmeans(crys_F, sex = "F", y.lab = "Crystallinity", main.pheno = "Crystallinity - Spine",mult = c(0.025, 0.1, 0.07)))

crys_M <-lm(Crystallinity ~ GENOTYPE, data=raman_spine_sd_M)
Anova(crys_M) 
lff<-lsmeans(crys_M,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
(p18 = plot_lsmeans(crys_M, sex = "M", y.lab = "Crystallinity", main.pheno = "Crystallinity - Spine", mult = c(0.025, 0.1, 0.06)))

p = (p10|p11|p12)/(p13|p14|p15)/(p16|p17|p18)

ggsave(plot = p,filename = "var_raman_spine.pdf", device = "pdf",path = "C://Users/basel/Desktop/",dpi = 300,width = 15,height = 20)

