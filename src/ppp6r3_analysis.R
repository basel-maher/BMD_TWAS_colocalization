library(emmeans)
library(ggplot2)
library(ggsignif)
library(ggpubr)#mainly for stat_pvalue_manual
library(car)

#
options(stringsAsFactors = FALSE)
set.seed(8675309)
options("na.action")
#$na.action
#[1] "na.omit"

#read in the combined phenotype file

ppp6r3_data = read.csv("./data/ppp6r3_cleaned.csv")

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
ppp6r3_caliper_cleaned = read.csv(file = "./data/ppp6r3_caliper_cleaned.csv")



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

ggsave(plot = p7,filename = "AP_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)


#################################################################
lm.ML <-lm(ML~GENOTYPE+SEX+sac.weight+age,data=caliper_data)
Anova(lm.ML)
lff<-lsmeans(lm.ML,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.ML$model$SEX, lm.ML$model$GENOTYPE)
p10 = plot_lsmeans(lm.ML, sex = "A", y.lab = "Medial-Lateral Width (mm)", main.pheno = "Medial-Lateral Femoral Width", mult = c(0.003, 0.015, 0.008))
ggsave(plot = p10,filename = "ML_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)

#################################################################
lm.FL <-lm(FL~GENOTYPE+SEX+sac.weight+age,data=caliper_data_cond)
Anova(lm.FL)
lff<-lsmeans(lm.FL,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.FL$model$SEX, lm.FL$model$GENOTYPE)
p13 = plot_lsmeans(lm.FL, sex = "A", y.lab = "Femoral Length (mm)", main.pheno = "Femoral Length", mult = c(0.003, 0.015, 0.008))
ggsave(plot = p13,filename = "FL_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)
####DXA analyses####
#analyze: BMD, BMC, B.Area, T.Area, RST?, Fat?


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
ggsave(plot = p1,filename = "ls_aBMD_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)


########################################### REPEAT FOR FEMURS ###################################################
lm.BMD <-lm(BMD~GENOTYPE+SEX+sac.weight+age+CenterRectX+CenterRectY,data=fn_bmd)
Anova(lm.BMD)
lff<-lsmeans(lm.BMD,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.BMD$model$SEX, lm.BMD$model$GENOTYPE)
(p1 = plot_lsmeans(lm.BMD, sex = "A", y.lab = "areal Bone Mineral Density - Femur (g HA/cm2)", main.pheno = "areal Bone Mineral Density - Femur (g HA/cm2)", mult = c(0.003, 0.015, 0.009)))
ggsave(plot = p1,filename = "fn_aBMD_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)

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
ggsave(plot = p1,filename = "tbbvtv_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)

#################################################################
lm.Tb.vBMD <-lm(Tb.vBMD~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.vBMD)
lff<-lsmeans(lm.Tb.vBMD,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.vBMD$model$SEX, lm.Tb.vBMD$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.vBMD, sex = "A", y.lab = "Tb.vBMD (mgHA/cm3)", main.pheno = "Trabecular vBMD - Lumbar Spine (mgHA/cm3)", mult = c(0.007, 0.025, 0.015)))
ggsave(plot = p1,filename = "vBMD_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)


####
lm.Tb.vTMD <-lm(Tb.vTMD~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.vTMD)
lff<-lsmeans(lm.Tb.vTMD,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.vTMD$model$SEX, lm.Tb.vTMD$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.vTMD, sex = "A", y.lab = "Tb.TMD (mgHA/cm3)", main.pheno = "Trabecular TMD - Lumbar Spine (mgHA/cm3)", mult = c(0.002, 0.007, 0.004)))
ggsave(plot = p1,filename = "TMD_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)

####
lm.Tb.Sp <-lm(Tb.Sp~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.Sp)
lff<-lsmeans(lm.Tb.Sp,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.Sp$model$SEX, lm.Tb.Sp$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.Sp, sex = "A", y.lab = "Trabecular Separation", main.pheno = "Trabecular Separation - Lumbar Spine", mult = c(0.01, 0.03, 0.02)))
ggsave(plot = p1,filename = "TbSp_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)


####
lm.Tb.Th <-lm(Tb.Th~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.Th)
lff<-lsmeans(lm.Tb.Th,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.Th$model$SEX, lm.Tb.Th$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.Th, sex = "A", y.lab = "Trabecular Thickness", main.pheno = "Trabecular Thickness - Lumbar Spine", mult = c(0.01, 0.035, 0.02)))
ggsave(plot = p1,filename = "TbTh_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)

####
lm.Tb.N <-lm(Tb.N~GENOTYPE+SEX+sac.weight+age, data=uCT)
Anova(lm.Tb.N)
lff<-lsmeans(lm.Tb.N,"GENOTYPE")
plot(lff,horizontal=FALSE)
lsmeans(lff,pairwise~GENOTYPE, adj="tukey")
lff = as.data.frame(lff)
table(lm.Tb.N$model$SEX, lm.Tb.N$model$GENOTYPE)
(p1 = plot_lsmeans(lm.Tb.N, sex = "A", y.lab = "Trabecular Number", main.pheno = "Trabecular Number - Lumbar Spine", mult = c(0.01, 0.035, 0.02)))
ggsave(plot = p1,filename = "TbN_c.pdf", device = "pdf",path = "E:/Users/Basel/Desktop/",dpi = 300,width = 5,height = 5)


