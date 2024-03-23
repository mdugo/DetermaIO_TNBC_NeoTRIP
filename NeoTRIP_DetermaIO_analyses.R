
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(grid)
library(cowplot)
library(vcd)
library(egg)
library(grDevices)
library(rcartocolor)
library(epiR)
library(openxlsx)

# put your path here
inputPath<-""
outputPath<-gsub("\\/data", "\\/results", inputPath)
dir.create(path = outputPath)

mycols<-carto_pal(12,"Bold")[c(1,2,12,12)]

df<-read.table(paste0(inputPath, "/NeoTRIP_data_IOScore.txt"),sep="\t", header=T, as.is=T)
df$sTILs_IHC<-factor(ifelse(df$sTILs_IHC=="High","High","Low/Mid"), levels=c("Low/Mid","High"))
df$TNBCtypes<-factor(df$TNBCtypes, levels=c("LAR","BL1","BL2","M","MSL","ND"))
df$pCR.num<-ifelse(df$pCR=="pCR",1,0)
df2<-df[!is.na(df$IO_cat_qPCR),]


#********************************************************
# Figure 1 - Pre-treatment RT-qPCR IO score
#********************************************************

# overall pie chart

df.temp<-df2
freq<-as.data.frame(round(table(df.temp$IO_cat_qPCR)/nrow(df.temp)*100,1))
freq<-freq[order(freq$Freq,decreasing=T),]
freq$lab.ypos = rev(cumsum(rev(freq$Freq)) - 0.5*rev(freq$Freq))
freq$Var1<-factor(freq$Var1,levels=freq$Var1)
freq$Label<-paste0(round(freq$Freq,1),"%")
colnames(freq)[1]<-"IOScore"
g1<-ggplot(freq, aes(x="", y=Freq,fill=IOScore)) +
  geom_bar(width=1, stat = "identity") +
  coord_polar(theta="y", start=0) +
  geom_text(aes(y = lab.ypos, label = Label), color = "white",size=5, fontface="plain") +
  scale_fill_manual(name="", labels=c("IO-","IO+"), values=mycols[1:2]) +
  theme_void() +
  theme(plot.margin=margin(t=0,r=0,b=0,l=0),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))

# barplot frequency vs arms

conta<-table(df.temp$IO_cat_qPCR,df.temp$Arm)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="CT"] <- cumsum(freq$Freq[freq$Var1=="CT"])-(freq$Freq[freq$Var1=="CT"]/2)
freq$lab.pos[freq$Var1=="CT/A"] <- cumsum(freq$Freq[freq$Var1=="CT/A"])-(freq$Freq[freq$Var1=="CT/A"]/2)
freq$Label<-paste0(round(freq$Freq,0),"%")
freq$Label[freq$Label=="0%"]<-""
freq$Var1<-factor(freq$Var1, levels=c("CT/A","CT"))
test<-round(chisq.test(conta)$p.value, digits=2)

g2<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=mycols[1:2], labels=c("IO-","IO+")) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="right") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "white",size=3) +
  coord_cartesian(ylim = c(0, 102), clip = "off") +
  annotate('text', 1.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))


# barplot IOScore vs PD-L1

conta<-table(df.temp$PDL1_IHC,df.temp$IO_cat_qPCR)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="IO+"] = cumsum(freq$Freq[freq$Var1=="IO+"])-(freq$Freq[freq$Var1=="IO+"]/2)
freq$lab.pos[freq$Var1=="IO-"] = cumsum(freq$Freq[freq$Var1=="IO-"])-(freq$Freq[freq$Var1=="IO-"]/2)
freq$Label<-paste0(round(freq$Freq,1),"%")
freq$Label[freq$Label=="0%"]<-""
test<-format(fisher.test(conta,simulate.p.value=FALSE, workspace = 2000000)$p.value,scientific = T, digits=2)

g3<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=pal_d3("category10", alpha = 0.8)(4)[c(2,4)], labels=c("PD-L1-","PD-L1+")) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="bottom") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "black",size=2.5) +
  annotate('text', 1.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))


# barplot IOScore vs sTILs

conta<-table(df.temp$sTILs_IHC,df.temp$IO_cat_qPCR)
freq<-round(t(conta)/colSums(conta)*100,2)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="IO+"] = cumsum(freq$Freq[freq$Var1=="IO+"])-(freq$Freq[freq$Var1=="IO+"]/2)
freq$lab.pos[freq$Var1=="IO-"] = cumsum(freq$Freq[freq$Var1=="IO-"])-(freq$Freq[freq$Var1=="IO-"]/2)
freq$Label<-paste0(round(freq$Freq,1),"%")
freq$Label[freq$Label=="0%"]<-""
test<-format(fisher.test(conta,simulate.p.value=F)$p.value,scientific = T, digits=2)

g4<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=brewer.pal(9,"Blues")[c(4,6)]) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="bottom") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "black",size=2.5) +
  annotate('text', 1.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0)) 


# pCR by Arm and IOScore (categorical)

df.temp2<-data.frame(GeneSet="IOScore",
                     IOScore=df.temp$IO_cat_qPCR,
                     pCR=df.temp$pCR,
                     pCRnum=df.temp$pCR.num,
                     Arm=df.temp$Arm)
df.temp2$Group<-factor(paste(df.temp2$Arm,df.temp2$IOScore),levels=c("CT/A IO-","CT/A IO+","CT IO-","CT IO+"))

t2x2<-table(df.temp2$Group,df.temp2$pCR)
perc<-round(t2x2/rowSums(t2x2)*100,1)
t2x2.overall<-table(df.temp2$IOScore,df.temp2$pCR)
perc.overall<-round(t2x2.overall/rowSums(t2x2.overall)*100,1)

df.perc<-data.frame(Perc=c(perc[,1],perc.overall[,1]),
                    N=c(t2x2[,1],t2x2.overall[,1]),
                    IOscore=c(gsub(".* ","",rownames(t2x2)),rownames(t2x2.overall)),
                    Arm=c(gsub(" .*","",rownames(t2x2)),rep("CT/A+CT",2)),
                    Group=c(rownames(t2x2),rownames(t2x2.overall)))
df.perc$Arm<-factor(df.perc$Arm, levels=c("CT/A+CT","CT/A","CT"))
df.perc$Group<-factor(df.perc$Group,levels=c("IO-","IO+","CT/A IO-","CT/A IO+","CT IO-","CT IO+"))
df.perc$Position<-c(4,5,7,8,1,2)

fit.overall<-glm(pCRnum~IOScore, family =binomial, data=df.temp2)
fit.cta<-glm(pCRnum~IOScore, family =binomial, data=df.temp2[df.temp2$Arm=="CT/A",])
fit.ct<-glm(pCRnum~IOScore, family =binomial, data=df.temp2[df.temp2$Arm=="CT",])
fit.interact<-glm(pCRnum~IOScore*Arm, family =binomial, data=df.temp2)

g5<-ggplot(data=df.perc, aes(x=Position, y=Perc, fill=IOscore)) +
  geom_bar(stat="identity", width=0.99, position = position_dodge(width=0.9)) +
  theme_pubr(border=F, legend="bottom") +
  xlab("") +
  ylab("pCR rate (%)") +
  scale_x_continuous(breaks=seq(1.5,7.5,3), labels=c("CT/A+CT\n(N = 220)","CT/A\n(N = 108)","CT\n(N = 112)")) +
  scale_y_continuous(breaks=seq(0,80,20),limits = c(0,90)) +
  geom_text(aes(x=Position,y=Perc+2,label=paste0(round(Perc,1),"%")), color="black", size=3) +
  scale_fill_manual(name="",values = mycols[1:2]) +
  theme(plot.margin=margin(t=10),
        axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=4)) +
  annotate('text', 1.5, 80, 
           label=paste("italic(P)==", round(summary(fit.overall)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', 4.5, 80, 
           label=paste("italic(P)==", round(summary(fit.cta)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', 7.5, 80, 
           label=paste("italic(P)==", round(summary(fit.ct)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', x= 6, y=89.5, 
           label=paste("italic(P)[interaction]==", round(summary(fit.interact)$coefficients[4,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate("rect", xmin = 4.5, xmax = 7.5, ymin = 87, ymax =87, alpha=1,colour = "black") +
  annotate("rect", xmin = 4.5, xmax = 4.5, ymin = 85, ymax =87, alpha=1, colour = "black") +
  annotate("rect", xmin = 7.5, xmax = 7.5, ymin = 85, ymax =87, alpha=1, colour = "black")


layout.matrix<-matrix(c(1,2,5,3,4,5),ncol=3,byrow=T)

myplot1 <- arrangeGrob(g1,top = textGrob("A",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot2 <- arrangeGrob(g2,top = textGrob("",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot3 <- arrangeGrob(g3,top = textGrob("B",
                                         x = unit(0, "npc"),
                                         y = unit(0.5, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot4 <- arrangeGrob(g4,top = textGrob("",
                                         x = unit(0, "npc"),
                                         y = unit(0.5, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot5 <- arrangeGrob(g5,top = textGrob("C",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))

pdf(paste0(outputPath, "/Figure_1.pdf"),width=7.48,height=5)
grid.arrange(myplot1,myplot2,myplot3,myplot4,myplot5,
             layout_matrix=layout.matrix,
             widths=unit(c(4.6,4.6,9.5), c("cm")))
dev.off()

#********************************************************
# Figure 2 - Pre-treatment TNBC types
#********************************************************

# overall pie chart

df.temp<-df
freq<-as.data.frame(round(table(df.temp$TNBCtypes)/nrow(df.temp)*100,1))
freq<-freq[order(freq$Freq,decreasing=T),]
freq$lab.ypos = rev(cumsum(rev(freq$Freq)) - 0.5*rev(freq$Freq))
freq$Var1<-factor(freq$Var1,levels=freq$Var1)
freq$Label<-paste0(round(freq$Freq,0),"%")
colnames(freq)[1]<-"TNBC"

subcol<-structure(c(pal_npg("nrc")(5),"grey60"), names=c("BL1", "BL2", "LAR", "M", "MSL", "ND"))
subcol<-subcol[match(freq$TNBC, names(subcol))]
g1<-ggplot(freq, aes(x="", y=Freq,fill=TNBC)) +
  geom_bar(width=1, stat = "identity") +
  coord_polar(theta="y", start=0) +
  geom_text(aes(x=1.2,y = lab.ypos, label = Label), color = "white",size=2.5, fontface="plain") +
  scale_fill_manual(name="",values=subcol, labels=names(subcol))  +
  theme_void() +
  theme(plot.margin=margin(t=0,r=0,b=0,l=0),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-10,b=10))

# barplot frequency vs arms

conta<-table(df.temp$TNBCtypes,df.temp$Arm)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq$Var2<-as.vector(freq$Var2)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="CT"] <- cumsum(freq$Freq[freq$Var1=="CT"])-(freq$Freq[freq$Var1=="CT"]/2)
freq$lab.pos[freq$Var1=="CT/A"] <- cumsum(freq$Freq[freq$Var1=="CT/A"])-(freq$Freq[freq$Var1=="CT/A"]/2)
freq$Label<-paste0(round(freq$Freq,0),"%")
freq$Label[freq$Label=="0%"]<-""
test<-round(chisq.test(conta)$p.value, digits=2)

g2<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=c(pal_npg("nrc")(5),"grey60"), labels=c("BL1", "BL2", "LAR", "M", "MSL", "ND")) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="right") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "white",size=2) +
  coord_cartesian(ylim = c(0, 102), clip = "off") +
  annotate('text', 1.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))

# barplot TNBC types vs PD-L1

conta<-table(df.temp$PDL1_IHC,df.temp$TNBCtypes)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="ND"] = cumsum(freq$Freq[freq$Var1=="ND"])-(freq$Freq[freq$Var1=="ND"]/2)
freq$lab.pos[freq$Var1=="MSL"] = cumsum(freq$Freq[freq$Var1=="MSL"])-(freq$Freq[freq$Var1=="MSL"]/2)
freq$lab.pos[freq$Var1=="M"] = cumsum(freq$Freq[freq$Var1=="M"])-(freq$Freq[freq$Var1=="M"]/2)
freq$lab.pos[freq$Var1=="LAR"] = cumsum(freq$Freq[freq$Var1=="LAR"])-(freq$Freq[freq$Var1=="LAR"]/2)
freq$lab.pos[freq$Var1=="BL2"] = cumsum(freq$Freq[freq$Var1=="BL2"])-(freq$Freq[freq$Var1=="BL2"]/2)
freq$lab.pos[freq$Var1=="BL1"] = cumsum(freq$Freq[freq$Var1=="BL1"])-(freq$Freq[freq$Var1=="BL1"]/2)
freq$Label<-paste0(round(freq$Freq,0),"%")
freq$Label[freq$Label=="0%"]<-""
test<-format(fisher.test(conta,simulate.p.value=FALSE, workspace = 2000000)$p.value,scientific = T, digits=2)

g3<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=pal_d3("category10", alpha = 0.8)(4)[c(2,4)], labels=c("PD-L1-","PD-L1+")) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="bottom") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "black",size=2.2) +
  annotate('text', 3.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))

# barplot TNBC types vs sTILs

conta<-table(df.temp$sTILs_IHC,df.temp$TNBCtypes)
freq<-round(t(conta)/colSums(conta)*100,2)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="ND"] = cumsum(freq$Freq[freq$Var1=="ND"])-(freq$Freq[freq$Var1=="ND"]/2)
freq$lab.pos[freq$Var1=="MSL"] = cumsum(freq$Freq[freq$Var1=="MSL"])-(freq$Freq[freq$Var1=="MSL"]/2)
freq$lab.pos[freq$Var1=="M"] = cumsum(freq$Freq[freq$Var1=="M"])-(freq$Freq[freq$Var1=="M"]/2)
freq$lab.pos[freq$Var1=="LAR"] = cumsum(freq$Freq[freq$Var1=="LAR"])-(freq$Freq[freq$Var1=="LAR"]/2)
freq$lab.pos[freq$Var1=="BL2"] = cumsum(freq$Freq[freq$Var1=="BL2"])-(freq$Freq[freq$Var1=="BL2"]/2)
freq$lab.pos[freq$Var1=="BL1"] = cumsum(freq$Freq[freq$Var1=="BL1"])-(freq$Freq[freq$Var1=="BL1"]/2)
freq$Label<-paste0(round(freq$Freq,0),"%")
freq$Label[freq$Label=="0%"]<-""
test<-format(fisher.test(conta,simulate.p.value=T, B=20000)$p.value,scientific = T, digits=2)

g4<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=brewer.pal(9,"Blues")[c(4,6)]) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="bottom") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "black",size=2.2) +
  annotate('text', 3.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0)) 


# pCR by Arm and TNBC types

df.temp2<-df.temp[df.temp$TNBCtypes!="ND",]
df2<-data.frame(GeneSet="TNBC",
                TNBC=droplevels(df.temp2$TNBCtypes),
                pCR=df.temp2$pCR,
                Arm=df.temp2$Arm)
df2$Group<-factor(paste(df2$Arm,df2$TNBC),levels=c("CT/A BL1","CT/A BL2","CT/A LAR","CT/A M","CT/A MSL",
                                                   "CT BL1","CT BL2","CT LAR","CT M","CT MSL"))
t2x2<-table(df2$Group,df2$pCR)
perc<-round(t2x2/rowSums(t2x2)*100,1)
t2x2.overall<-table(df2$TNBC,df2$pCR)
perc.overall<-round(t2x2.overall/rowSums(t2x2.overall)*100,1)

df.perc<-data.frame(Perc=c(perc[,1],perc.overall[,1]),
                    N=c(t2x2[,1],t2x2.overall[,1]),
                    TNBC=c(gsub(".* ","",rownames(t2x2)),rownames(t2x2.overall)),
                    Arm=c(gsub(" .*","",rownames(t2x2)),rep("CT/A+CT",5)),
                    Group=c(rownames(t2x2),rownames(t2x2.overall)))
df.perc$Arm<-factor(df.perc$Arm, levels=c("CT/A+CT","CT/A","CT"))
df.perc$Group<-factor(df.perc$Group,levels=c("CT/A BL1","CT/A BL2","CT/A LAR","CT/A M","CT/A MSL",
                                             "CT BL1","CT BL2","CT LAR","CT M","CT MSL",
                                             "BL1","BL2","LAR","M","MSL"))
df.perc$TNBC<-factor(df.perc$TNBC, levels=c("LAR","BL1","BL2","M","MSL"))
df.perc<-df.perc[order(df.perc$Arm, df.perc$TNBC),]
df.perc$Position<-c(1:5,7:11,13:17)

g5<-ggplot(data=df.perc, aes(x=Position, y=Perc, fill=TNBC)) +
  geom_bar(stat="identity", width=0.99, position = position_dodge(width=0.9)) +
  theme_pubr(border=F, legend="bottom") +
  xlab("") +
  ylab("pCR rate (%)") +
  scale_x_continuous(breaks=seq(3.5,15.5,6), labels=c("CT/A+CT","CT/A","CT")) +
  scale_y_continuous(breaks=seq(0,100,20),limits = c(0,100)) +
  geom_text(aes(x=Position,y=Perc+1,label=paste0(round(Perc,0),"%")), color="black", size=2.5, angle=90, hjust=0) +
  scale_fill_manual(name="",values = pal_npg("nrc")(5)[c(3,1,2,4,5)]) +
  theme(plot.margin=margin(t=10),
        axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=4))


layout.matrix<-matrix(c(1,2,5,3,4,5),ncol=3,byrow=T)

myplot1 <- arrangeGrob(g1,top = textGrob("A",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot2 <- arrangeGrob(g2,top = textGrob("",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot3 <- arrangeGrob(g3,top = textGrob("B",
                                         x = unit(0, "npc"),
                                         y = unit(0.5, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot4 <- arrangeGrob(g4,top = textGrob("",
                                         x = unit(0, "npc"),
                                         y = unit(0.5, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot5 <- arrangeGrob(g5,top = textGrob("C",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))

pdf(paste0(outputPath, "/Figure_2.pdf"),width=7.48,height=5)
grid.arrange(myplot1,myplot2,myplot3,myplot4,myplot5,
             layout_matrix=layout.matrix,
             widths=unit(c(4.6,4.6,9.5), c("cm")))
dev.off()

#***************************************************
# Figure 3 - Validation of IO score in I-SPY2 trial
#***************************************************

ispy2<-read.table(paste0(inputPath, "/ISPY2_data_IOScore.txt"),sep="\t", header=T, as.is=T)

# pie chart

freq<-as.data.frame(round(table(ispy2$IO_ScoreBin)/nrow(ispy2)*100,1))
freq<-freq[order(freq$Freq,decreasing=T),]
freq$lab.ypos = rev(cumsum(rev(freq$Freq)) - 0.5*rev(freq$Freq))
freq$Var1<-factor(freq$Var1,levels=freq$Var1)
freq$Label<-paste0(round(freq$Freq,1),"%")
colnames(freq)[1]<-"IOScore"
g1<-ggplot(freq, aes(x="", y=Freq,fill=IOScore)) +
  geom_bar(width=1, stat = "identity") +
  coord_polar(theta="y", start=0) +
  geom_text(aes(y = lab.ypos, label = Label), color = "white",size=5, fontface="plain") +
  scale_fill_manual(name="", labels=c("IO-","IO+"), values=mycols[1:2]) +
  theme_void() +
  theme(plot.margin=margin(t=0,r=0,b=0,l=0),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))


# barplot frequency vs arms

conta<-table(ispy2$IO_ScoreBin,ispy2$Therapy)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="Paclitaxel"] <- cumsum(freq$Freq[freq$Var1=="Paclitaxel"])-(freq$Freq[freq$Var1=="Paclitaxel"]/2)
freq$lab.pos[freq$Var1=="Paclitaxel + Pembrolizumab"] <- cumsum(freq$Freq[freq$Var1=="Paclitaxel + Pembrolizumab"])-(freq$Freq[freq$Var1=="Paclitaxel + Pembrolizumab"]/2)
freq$Var1<-factor(freq$Var1, levels=c("Paclitaxel + Pembrolizumab","Paclitaxel"))
freq$Label<-paste0(round(freq$Freq,0),"%")
freq$Label[freq$Label=="0%"]<-""
test<-round(chisq.test(conta)$p.value, digits=2)

g2<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=mycols[1:2], labels=c("IO-","IO+")) +
  scale_x_discrete(labels=c("Paclitaxel +\nPembrolizumab","Paclitaxel")) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="right") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "white",size=4) +
  coord_cartesian(ylim = c(0, 102), clip = "off") +
  annotate('text', 1.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  theme(axis.title = element_text(size=8),
        axis.text.x=element_text(size=9, angle=45, hjust=1),
        axis.text.y=element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))

# pCR by Arm and IOScore

ispy2$Therapy<-factor(ispy2$Therapy, levels=c("Paclitaxel + Pembrolizumab","Paclitaxel"))
fit.overall<-glm(PCR~IO_ScoreBin, family =binomial, data=ispy2)
fit1<-glm(PCR~IO_ScoreBin, family =binomial, data=ispy2[ispy2$Therapy=="Paclitaxel + Pembrolizumab",])
c(round(exp(coef(fit1))[2],3),
  round(exp(confint.default(fit1, level = 0.95))[2,1],3),
  round(exp(confint.default(fit1, level = 0.95))[2,2],3),
  summary(fit1)$coefficients[2,3],
  summary(fit1)$coefficients[2,4])

fit2<-glm(PCR~IO_ScoreBin, family =binomial, data=ispy2[ispy2$Therapy=="Paclitaxel",])
c(round(exp(coef(fit2))[2],3),
  round(exp(confint.default(fit2, level = 0.95))[2,1],3),
  round(exp(confint.default(fit2, level = 0.95))[2,2],3),
  summary(fit2)$coefficients[2,3],
  summary(fit2)$coefficients[2,4])
fit.interact<-glm(PCR~IO_ScoreBin*Therapy, family =binomial, data=ispy2)

pvals<-data.frame(Therapy=c("Overall",levels(ispy2$Therapy)),
                  pval=c(round(summary(fit.overall)$coefficients[2,4],3),
                         round(summary(fit1)$coefficients[2,4],3),
                         round(summary(fit2)$coefficients[2,4],3)),
                  PosLab=seq(1.5,4.5,1.5))

df3<-data.frame(GeneSet="IOScore",
                IOScore=ispy2$IO_ScoreBin,
                pCR=ifelse(ispy2$PCR=="1","pCR","RD"),
                Arm=ispy2$Therapy)
df3$Group<-factor(paste(df3$Arm,df3$IOScore),
                  levels=paste(rep(c("Paclitaxel + Pembrolizumab","Paclitaxel"), each=2), rep(c("IO-","IO+"),2)))

t2x2<-table(df3$Group,df3$pCR)
perc<-round(t2x2/rowSums(t2x2)*100,1)
t2x2.overall<-table(df3$IOScore,df3$pCR)
perc.overall<-round(t2x2.overall/rowSums(t2x2.overall)*100,1)

df.perc<-data.frame(Perc=c(perc[,1],perc.overall[,1]),
                    N=c(t2x2[,1],t2x2.overall[,1]),
                    IOscore=c(gsub(".* ","",rownames(t2x2)),rownames(t2x2.overall)),
                    Arm=c(gsub(" .*","",rownames(t2x2)),rep("CT/A+CT",2)),
                    Group=c(rownames(t2x2),rownames(t2x2.overall)))
df.perc$Arm<-factor(df.perc$Arm, levels=c("Paclitaxel + Pembrolizumab","Paclitaxel"))
df.perc$Group<-factor(df.perc$Group,levels=paste(rep(c("Paclitaxel + Pembrolizumab","Paclitaxel"), each=2), rep(c("IO-","IO+"),2)))
df.perc$Position<-c(3,3,4.5,4.5,1.5,1.5)
df.perc$NTot<-rep(c(table(ispy2$Therapy), nrow(ispy2)), each=2)

g3<-ggplot(data=df.perc, aes(x=Position, y=Perc, fill=IOscore)) +
  geom_bar(stat="identity", width=0.99, position = position_dodge(width=0.99)) +
  theme_pubr(border=F, legend="bottom") +
  xlab("") +
  ylab("pCR rate (%)") +
  scale_x_continuous(breaks=c(1.5, 3, 4.5), labels=paste0(c("Overall","Paclitaxel + Pembrolizumab","Paclitaxel"),"\n(N = ",df.perc$NTot[c(6,2,4)],")")) +
  scale_y_continuous(breaks=seq(0,100,20),limits = c(0,134)) +
  geom_text(aes(x=Position,y=Perc+2,label=paste0(round(Perc,1),"%")), color="black", size=3.3, angle=90, hjust=0,position = position_dodge(width = 0.99)) +
  scale_fill_manual(name="",values = mycols[1:2]) +
  theme(plot.margin=margin(t=10),
        axis.title = element_text(size=9),
        axis.text.x=element_text(size=9, angle=0, hjust=0.5),
        axis.text.y=element_text(size=9),
        legend.text = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=4)) +
  annotate('text', x= c(1.5,3,4.5), y=111, 
           label=paste("italic(P)==", pvals$pval), 
           parse=TRUE, 
           hjust=0.5, size=3.3) +
  annotate('text', x= 4.5-((4.5-3)/2), y=131, 
           label=paste("italic(P)[interaction]==", round(summary(fit.interact)$coefficients[4,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3.3) +
  annotate("rect", xmin = 3, xmax = 4.5, ymin = 124, ymax =124, alpha=1,colour = "black") +
  annotate("rect", xmin = 3, xmax = 3, ymin = 119, ymax =124, alpha=1, colour = "black") +
  annotate("rect", xmin = 4.5, xmax = 4.5, ymin = 119, ymax =124, alpha=1, colour = "black")

# continuous IO score by arm and pCR

df3<-ispy2[ispy2$Therapy%in%c("Paclitaxel","Paclitaxel + Pembrolizumab"),]
df3b<-df3
df3b$Therapy<-"Overall"
df3<-rbind(df3,df3b)
df3$Therapy<-factor(df3$Therapy, levels=c("Overall","Paclitaxel + Pembrolizumab","Paclitaxel"))

fit.overall<-glm(PCR~IO_ScoreCont, family =binomial, data=df3[df3$Therapy=="Overall",])
fit.ct<-glm(PCR~IO_ScoreCont, family =binomial, data=df3[df3$Therapy=="Paclitaxel",])
fit.cta<-glm(PCR~IO_ScoreCont, family =binomial, data=df3[df3$Therapy=="Paclitaxel + Pembrolizumab",])
fit.interact<-glm(PCR~IO_ScoreCont*Therapy, family =binomial, data=df3[df3$Therapy!="Overall",])
df3$pCR<-ifelse(df3$PCR==1,"pCR","RD")
g4<-ggplot(data=df3, aes(x=Therapy, y=IO_ScoreCont, fill=pCR)) +
  geom_point(shape=21, 
             position=position_jitterdodge(
               jitter.width = 0.2,
               jitter.height = 0,
               dodge.width = 0.5,
               seed = NA),
             cex=2,alpha = 0.7,show.legend=T)   +
  geom_boxplot(alpha=0.1,show.legend=F,width=0.5) +
  #stat_compare_means(method = "t.test",paired=F,size=7,label.y =1) +
  theme_pubr() +
  xlab("") +
  ylab("IO score") +
  ylim(-1,1.5) +
  annotate('text', 1, 1.1, 
           label=paste("italic(P)==",round(summary(fit.overall)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', 2, 1.1, 
           label=paste("italic(P)==",round(summary(fit.cta)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', 3, 1.1, 
           label=paste("italic(P)==",round(summary(fit.ct)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', x= 2.5, y=1.49, 
           label=paste("italic(P)[interaction]==",round(summary(fit.interact)$coefficients[4,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate("rect", xmin = 2, xmax = 3, ymin = 1.35, ymax =1.35, alpha=1,colour = "black") +
  annotate("rect", xmin = 2, xmax = 2, ymin = 1.25, ymax =1.35, alpha=1, colour = "black") +
  annotate("rect", xmin = 3, xmax = 3, ymin = 1.25, ymax =1.35, alpha=1, colour = "black") +
  scale_fill_manual(name="", values=carto_pal(3,"Vivid")[1:2]) +
  theme(plot.margin=margin(t=10),
        axis.title = element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        legend.text = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=4)) +
  scale_x_discrete(labels=paste0(c("Overall","Paclitaxel + Pembrolizumab","Paclitaxel"),"\n(N = ",df.perc$NTot[c(6,2,4)],")"))

layout.matrix<-matrix(c(1,3,2,4),ncol=2,byrow=T)

myplot1 <- arrangeGrob(g1,top = textGrob("A",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot2 <- arrangeGrob(g2,top = textGrob("",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot3 <- arrangeGrob(g3,top = textGrob("B",
                                         x = unit(0, "npc"),
                                         y = unit(0.5, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot4 <- arrangeGrob(g4,top = textGrob("C",
                                         x = unit(0, "npc"),
                                         y = unit(0.5, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))

pdf(paste0(outputPath, "/Figure_3.pdf"),width=8,height=7.5)
grid.arrange(myplot1,myplot2,myplot3,myplot4,
             layout_matrix=layout.matrix,
             widths=unit(c(6.5,11.5), c("cm")),
             heights=unit(c(10,8), c("cm")))
dev.off()


#*************************************************************
# Supplementary Figure 1 - QC
#*************************************************************

df<-df[order(df$Arm,df$pCR, df$Total_mapped_reads),]
df$PatientID<-factor(df$PatientID, levels=df$PatientID)
df$Group<-paste(df$Arm, df$pCR)
df$xPos<-1:nrow(df)
g1<-ggplot(data=df, aes(x=xPos, y=Total_mapped_reads)) +
  geom_point(size=1) +
  geom_hline(yintercept = 4e07, color="grey70", linetype="dashed") +
  theme_pubr() +
  xlab("") +
  ylim(0,2e+08) +
  ylab("Million reads aligned to coding genes") + 
  scale_x_continuous(breaks=c(32.5,94,151,214), labels=c("CT, RD","CT, pCR","CT/A, RD","CT/A, pCR")) +
  geom_boxplot(data=df[df$Arm=="CT" & df$pCR=="RD",], aes(x=32.5,y=Total_mapped_reads), inherit.aes = F, coef=NULL, width=10, color="brown", alpha=0) +
  geom_boxplot(data=df[df$Arm=="CT" & df$pCR=="pCR",], aes(x=94,y=Total_mapped_reads), inherit.aes = F, coef=NULL, width=10, color="brown", alpha=0) +
  geom_boxplot(data=df[df$Arm=="CT/A" & df$pCR=="RD",], aes(x=151,y=Total_mapped_reads), inherit.aes = F, coef=NULL, width=10, color="brown", alpha=0) +
  geom_boxplot(data=df[df$Arm=="CT/A" & df$pCR=="pCR",], aes(x=214,y=Total_mapped_reads), inherit.aes = F, coef=NULL, width=10, color="brown", alpha=0)
  
pdf(paste0(outputPath, "/Supplementary_Figure_1.pdf"),width=8,height=5)
g1
dev.off()


#****************************************************************************
# Supplementary Figure 2 - Concordance RT-qPCR IOScore vs RNA-Seq IOScore
#****************************************************************************

df.temp<-df[!is.na(df$IO_cat_qPCR),]
df.temp$IOscore<-ifelse(df.temp$IO_cat_NGS=="IO-","NGS IO-","NGS IO+")
df.temp$IOscore_PCR<-ifelse(df.temp$IO_cat_qPCR=="IO-","RT-qPCR IO-","RT-qPCR IO+")
df.temp$IOscoreConcordance<-paste(df.temp$IOscore,df.temp$IOscore_PCR,sep="/")
df.temp$IOscoreConcordance<-factor(df.temp$IOscoreConcordance, levels=sort(unique(df.temp$IOscoreConcordance))[c(1,4,3,2)])
ccc<-round(epi.ccc(df.temp$IO_cont_NGS, df.temp$IO_cont_qPCR)$rho.c,2)

pearson<-round(cor(df.temp$IO_cont_NGS, df.temp$IO_cont_qPCR),2)
pearson.p<-format(cor.test(df.temp$IO_cont_NGS, df.temp$IO_cont_qPCR)$p.value,digits=3,scientific = T)

t2x2<-table(df.temp$IO_cat_qPCR,df.temp$IO_cat_NGS)
res.k <- Kappa(t2x2)
ci<-round(confint(res.k),2)
capture <- capture.output(result <- print(res.k))
pvalue<-strsplit(capture," ")[[2]][5]
rownames(t2x2)<-c("IO neg (RT-qPCR)","IO pos (RT-qPCR)")
colnames(t2x2)<-c("IO neg (NGS)","IO pos (NGS)")

g1<-ggplot(data=df.temp, aes(x=IO_cont_qPCR, y=IO_cont_NGS, fill=IOscoreConcordance)) +
  geom_hline(yintercept = 0.09, linetype="dashed",color="grey20") +
  geom_vline(xintercept = 0.09, linetype="dashed",color="grey20") +
  geom_point(size=1.3, pch=21) +
  scale_fill_manual(values=carto_pal(12,"Bold")[c(1,2,12,12)]) +
  xlab("IO Score (RT-qPCR)") +
  ylab("IO Score (RNA-Seq)") +
  theme_pubr(legend="none") +
  theme(axis.text = element_text(size=9), axis.title = element_text(size=12), axis.text.x = element_text(angle=90,hjust=1)) +
  scale_y_continuous(breaks=c(seq(-1,0,0.5),0.09,seq(0.5,1,0.5)),limits = c(-1,1)) +
  scale_x_continuous(breaks=c(seq(-1,0,0.5),0.09,seq(0.5,1,0.5)),limits = c(-1,1))

g2<-ggplot(data=df.temp, aes(x=IO_cont_qPCR, y=IO_cont_NGS)) +
  annotate("rect",xmin=-Inf, xmax=0, ymin=-Inf, ymax=0, fill=mycols[1], alpha=1) +
  annotate("rect",xmin=0, xmax=Inf, ymin=0, ymax=Inf, fill=mycols[2], alpha=1) +
  annotate("rect",xmin=-Inf, xmax=0, ymin=0, ymax=Inf, fill=mycols[3], alpha=1) +
  annotate("rect",xmin=0, xmax=Inf, ymin=-Inf, ymax=0, fill=mycols[4], alpha=1) +
  geom_hline(yintercept = 0, color="grey20", linewidth=1) +
  geom_vline(xintercept = 0, color="grey20", linewidth=1) +
  ylim(-1,1) +
  xlim(-1,1) +
  theme_void() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  annotate('text', -0.55, -0.5, 
           label=paste0("RT-qPCR IO-\nRNA-Seq IO-\n(N = ", t2x2[1,1],")"), 
           parse=FALSE, 
           hjust=0.5, size=4,color="white", fontface=2) +
  annotate('text', 0.55, 0.5, 
           label=paste0("RT-qPCR IO+\nRNA-Seq IO+\n(N = ", t2x2[2,2],")"), 
           parse=FALSE, 
           hjust=0.5, size=4,color="white", fontface=2) +
  annotate('text', -0.55, 0.5, 
           label=paste0("RT-qPCR IO-\nRNA-Seq IO+\n(N = ", t2x2[1,2],")"), 
           parse=FALSE, 
           hjust=0.5, size=4,color="white", fontface=2) +
  annotate('text', 0.55, -0.5, 
           label=paste0("RT-qPCR IO+\nRNA-Seq IO-\n(N = ", t2x2[2,1],")"), 
           parse=FALSE, 
           hjust=0.5, size=4,color="white", fontface=2)

pdf(paste0(outputPath, "/Supplementary_Figure_2.pdf"),width=7.08,height=3.54)
plot_grid(g1,g2,ncol=2, labels = c("A","B"))
dev.off()

#*************************************************************
# Supplementary Figure 3 - Continuous RT-qPCR IOScore vs pCR
#*************************************************************

df.temp<-df[!is.na(df$IO_cat_qPCR),]
df3<-data.frame(GeneSet="IOScore",
                IOScore=df.temp$IO_cont_qPCR,
                pCR=df.temp$pCR,
                pCRnum=df.temp$pCR.num,
                Arm=df.temp$Arm)
df3b<-df3
df3b$Arm<-"CT/A+CT"
df3<-rbind(df3,df3b)
df3$Arm<-factor(df3$Arm, levels=c("CT/A+CT","CT/A","CT"))

fit.overall<-glm(pCRnum~IOScore, family =binomial, data=df3[df3$Arm=="CT/A+CT",])
fit.cta<-glm(pCRnum~IOScore, family =binomial, data=df3[df3$Arm=="CT/A",])
fit.ct<-glm(pCRnum~IOScore, family =binomial, data=df3[df3$Arm=="CT",])
fit.interact<-glm(pCRnum~IOScore*Arm, family =binomial, data=df3[df3$Arm!="CT/A+CT",])

g1<-ggplot(data=df3, aes(x=Arm, y=IOScore, fill=pCR)) +
  geom_point(shape=21, 
             position=position_jitterdodge(
               jitter.width = 0.2,
               jitter.height = 0,
               dodge.width = 0.5,
               seed = NA),
             cex=3,alpha = 0.7,show.legend=T)   +
  geom_boxplot(alpha=0.1,show.legend=F,width=0.5) +
  #stat_compare_means(method = "t.test",paired=F,size=7,label.y =1) +
  theme_pubr() +
  xlab("") +
  ylab("IO score") +
  ylim(-1,1.35) +
  annotate('text', 1, 1, 
           label=paste("italic(P)==",round(summary(fit.overall)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=4) +
  annotate('text', 2, 1, 
           label=paste("italic(P)==",round(summary(fit.cta)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=4) +
  annotate('text', 3, 1, 
           label=paste("italic(P)==",round(summary(fit.ct)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=4) +
  annotate('text', x= 2.5, y=1.3, 
           label=paste("italic(P)[interaction]==",round(summary(fit.interact)$coefficients[4,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=4) +
  annotate("rect", xmin = 2, xmax = 3, ymin = 1.2, ymax =1.2, alpha=1,colour = "black") +
  annotate("rect", xmin = 2, xmax = 2, ymin = 1.14, ymax =1.2, alpha=1, colour = "black") +
  annotate("rect", xmin = 3, xmax = 3, ymin = 1.14, ymax =1.2, alpha=1, colour = "black") +
  scale_fill_manual(name="", values=carto_pal(3,"Vivid")[1:2])

pdf(paste0(outputPath, "/Supplementary_Figure_3.pdf"),width=7.48,height=5.5)
g1
dev.off()

#********************************************************
# Supplementary Figure 4 - Pre-treatment RNA-Seq IO score
#********************************************************

# overall pie chart

df.temp<-df
freq<-as.data.frame(round(table(df.temp$IO_cat_NGS)/nrow(df.temp)*100,1))
freq<-freq[order(freq$Freq,decreasing=T),]
freq$lab.ypos = rev(cumsum(rev(freq$Freq)) - 0.5*rev(freq$Freq))
freq$Var1<-factor(freq$Var1,levels=freq$Var1)
freq$Label<-paste0(round(freq$Freq,1),"%")
colnames(freq)[1]<-"IOScore"
g1<-ggplot(freq, aes(x="", y=Freq,fill=IOScore)) +
  geom_bar(width=1, stat = "identity") +
  coord_polar(theta="y", start=0) +
  geom_text(aes(y = lab.ypos, label = Label), color = "white",size=5, fontface="plain") +
  scale_fill_manual(name="", labels=c("IO-","IO+"), values=mycols[1:2]) +
  theme_void() +
  theme(plot.margin=margin(t=0,r=0,b=0,l=0),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))

# barplot frequency vs arms

conta<-table(df.temp$IO_cat_NGS,df.temp$Arm)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="CT"] <- cumsum(freq$Freq[freq$Var1=="CT"])-(freq$Freq[freq$Var1=="CT"]/2)
freq$lab.pos[freq$Var1=="CT/A"] <- cumsum(freq$Freq[freq$Var1=="CT/A"])-(freq$Freq[freq$Var1=="CT/A"]/2)
freq$Label<-paste0(round(freq$Freq,0),"%")
freq$Label[freq$Label=="0%"]<-""
test<-round(chisq.test(conta)$p.value, digits=2)

g2<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=mycols[1:2], labels=c("IO-","IO+")) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="right") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "white",size=3) +
  coord_cartesian(ylim = c(0, 102), clip = "off") +
  annotate('text', 1.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))

# barplot IOScore vs PD-L1

conta<-table(df.temp$PDL1_IHC,df.temp$IO_cat_NGS)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="IO+"] = cumsum(freq$Freq[freq$Var1=="IO+"])-(freq$Freq[freq$Var1=="IO+"]/2)
freq$lab.pos[freq$Var1=="IO-"] = cumsum(freq$Freq[freq$Var1=="IO-"])-(freq$Freq[freq$Var1=="IO-"]/2)
freq$Label<-paste0(round(freq$Freq,1),"%")
freq$Label[freq$Label=="0%"]<-""
test<-format(fisher.test(conta,simulate.p.value=FALSE, workspace = 2000000)$p.value,scientific = T, digits=2)

g3<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=pal_d3("category10", alpha = 0.8)(4)[c(2,4)], labels=c("PD-L1-","PD-L1+")) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="bottom") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "black",size=2.5) +
  annotate('text', 1.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))


# barplot IOScore vs sTILs

conta<-table(df.temp$sTILs_IHC,df.temp$IO_cat_NGS)
freq<-round(t(conta)/colSums(conta)*100,2)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="IO+"] = cumsum(freq$Freq[freq$Var1=="IO+"])-(freq$Freq[freq$Var1=="IO+"]/2)
freq$lab.pos[freq$Var1=="IO-"] = cumsum(freq$Freq[freq$Var1=="IO-"])-(freq$Freq[freq$Var1=="IO-"]/2)
freq$Label<-paste0(round(freq$Freq,1),"%")
freq$Label[freq$Label=="0%"]<-""
test<-format(fisher.test(conta,simulate.p.value=F)$p.value,scientific = T, digits=2)

g4<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=brewer.pal(9,"Blues")[c(4,6)]) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="bottom") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "black",size=2.5) +
  annotate('text', 1.5, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=0.5, size=2.5) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme(axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0)) 


# pCR by Arm and IOScore

df3<-data.frame(GeneSet="IOScore",
                IOScore=df.temp$IO_cat_NGS,
                pCR=df.temp$pCR,
                pCRnum=df.temp$pCR.num,
                Arm=df.temp$Arm)
df3$Group<-factor(paste(df3$Arm,df3$IOScore),levels=c("CT/A IO-","CT/A IO+","CT IO-","CT IO+"))

t2x2<-table(df3$Group,df3$pCR)
perc<-round(t2x2/rowSums(t2x2)*100,1)
t2x2.overall<-table(df3$IOScore,df3$pCR)
perc.overall<-round(t2x2.overall/rowSums(t2x2.overall)*100,1)

df.perc<-data.frame(Perc=c(perc[,1],perc.overall[,1]),
                    N=c(t2x2[,1],t2x2.overall[,1]),
                    IOscore=c(gsub(".* ","",rownames(t2x2)),rownames(t2x2.overall)),
                    Arm=c(gsub(" .*","",rownames(t2x2)),rep("CT/A+CT",2)),
                    Group=c(rownames(t2x2),rownames(t2x2.overall)))
df.perc$Arm<-factor(df.perc$Arm, levels=c("CT/A+CT","CT/A","CT"))
df.perc$Group<-factor(df.perc$Group,levels=c("IO-","IO+","CT/A IO-","CT/A IO+","CT IO-","CT IO+"))
df.perc$Position<-c(4,5,7,8,1,2)

fit.overall<-summary(glm(pCRnum ~ IOScore, data=df3, family = "binomial"))
fit.cta<-summary(glm(pCRnum ~ IOScore, data=df3[df3$Arm=="CT/A",], family = "binomial"))
fit.ct<-summary(glm(pCRnum ~ IOScore, data=df3[df3$Arm=="CT",], family = "binomial"))
fit.int<-summary(glm(pCRnum ~ IOScore*Arm, data=df3, family = "binomial"))

g5<-ggplot(data=df.perc, aes(x=Position, y=Perc, fill=IOscore)) +
  geom_bar(stat="identity", width=0.99, position = position_dodge(width=0.9)) +
  theme_pubr(border=F, legend="bottom") +
  xlab("") +
  ylab("pCR rate (%)") +
  scale_x_continuous(breaks=seq(1.5,7.5,3), labels=c("CT/A+CT","CT/A","CT")) +
  scale_y_continuous(breaks=seq(0,100,20),limits = c(0,100)) +
  geom_text(aes(x=Position,y=Perc+5,label=paste0(round(Perc,1),"%")), color="black", size=3) +
  scale_fill_manual(name="",values = mycols[1:2]) +
  theme(plot.margin=margin(t=10),
        axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=4)) +
  annotate('text', 1.5, 85, 
           label=paste0("~italic(P)== ", round(fit.overall$coefficients[2,4],3)), #y~`=`~
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', 4.5, 85, 
           label=paste("~italic(P)== ", round(fit.cta$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', 7.5, 85, 
           label=paste("~italic(P)== ", round(fit.ct$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', x= 6, y=99, 
           label=paste("~italic(P)[interaction]== ", round(fit.int$coefficients[4,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate("rect", xmin = 4.5, xmax = 7.5, ymin = 93, ymax =93, alpha=1,colour = "black") +
  annotate("rect", xmin = 4.5, xmax = 4.5, ymin = 90, ymax =93, alpha=1, colour = "black") +
  annotate("rect", xmin = 7.5, xmax = 7.5, ymin = 90, ymax =93, alpha=1, colour = "black")

# IOScore (continuous) by arm and pCR

df3<-data.frame(GeneSet="IOScore",
                IOScore=df.temp$IO_cont_NGS,
                pCR=df.temp$pCR,
                pCRnum=df.temp$pCR.num,
                Arm=df.temp$Arm)
df3b<-df3
df3b$Arm<-"CT/A+CT"
df3<-rbind(df3,df3b)
df3$Arm<-factor(df3$Arm, levels=c("CT/A+CT","CT/A","CT"))

fit.overall<-glm(pCRnum~IOScore, family =binomial, data=df3)
fit.cta<-glm(pCRnum~IOScore, family =binomial, data=df3[df3$Arm=="CT/A",])
fit.ct<-glm(pCRnum~IOScore, family =binomial, data=df3[df3$Arm=="CT",])
fit.interact<-glm(pCRnum~IOScore*Arm, family =binomial, data=df3)

g6<-ggplot(data=df3, aes(x=Arm, y=IOScore, fill=pCR)) +
  geom_point(shape=21, 
             position=position_jitterdodge(
               jitter.width = 0.2,
               jitter.height = 0,
               dodge.width = 0.5,
               seed = NA),
             cex=1.3,alpha = 0.7,show.legend=T)   +
  geom_boxplot(alpha=0.1,show.legend=F,width=0.5) +
  #stat_compare_means(method = "t.test",paired=F,size=7,label.y =1) +
  theme_pubr() +
  xlab("") +
  ylab("IO score") +
  ylim(-1,1.5) +
  annotate('text', 1, 1.1, 
           label=paste("italic(P)==",round(summary(fit.overall)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', 2, 1.1, 
           label=paste("italic(P)==",round(summary(fit.cta)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', 3, 1.1, 
           label=paste("italic(P)==",round(summary(fit.ct)$coefficients[2,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate('text', x= 2.5, y=1.49, 
           label=paste("italic(P)[interaction]==",round(summary(fit.interact)$coefficients[4,4],3)), 
           parse=TRUE, 
           hjust=0.5, size=3) +
  annotate("rect", xmin = 2, xmax = 3, ymin = 1.35, ymax =1.35, alpha=1,colour = "black") +
  annotate("rect", xmin = 2, xmax = 2, ymin = 1.25, ymax =1.35, alpha=1, colour = "black") +
  annotate("rect", xmin = 3, xmax = 3, ymin = 1.25, ymax =1.35, alpha=1, colour = "black") +
  scale_fill_manual(name="", values=carto_pal(3,"Vivid")[1:2]) +
  theme(plot.margin=margin(t=10),
        axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=4))


layout.matrix<-matrix(c(1,2,5,7,3,4,6,7),ncol=4,byrow=T)

myplot1 <- arrangeGrob(g1,top = textGrob("A",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot2 <- arrangeGrob(g2,top = textGrob("",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot3 <- arrangeGrob(g3,top = textGrob("B",
                                         x = unit(0, "npc"),
                                         y = unit(0.5, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot4 <- arrangeGrob(g4,top = textGrob("",
                                         x = unit(0, "npc"),
                                         y = unit(0.5, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot5 <- arrangeGrob(g5,top = textGrob("C",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot6 <- arrangeGrob(g6,top = textGrob("D",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))

df.temp2<-df.temp[!is.na(df.temp$IO_cat_qPCR),]

fit.overall.pcr<-glm(pCR.num~IO_cat_qPCR, family =binomial, data=df.temp2)
fit.overall.ngs<-glm(pCR.num~IO_cat_NGS, family =binomial, data=df.temp2)

fit.cta.pcr<-glm(pCR.num~IO_cat_qPCR, family =binomial, data=df.temp2[df.temp2$Arm=="CT/A",])
fit.cta.ngs<-glm(pCR.num~IO_cat_NGS, family =binomial, data=df.temp2[df.temp2$Arm=="CT/A",])

fit.ct.pcr<-glm(pCR.num~IO_cat_qPCR, family =binomial, data=df.temp2[df.temp2$Arm=="CT",])
fit.ct.ngs<-glm(pCR.num~IO_cat_NGS, family =binomial, data=df.temp2[df.temp2$Arm=="CT",])

aic<-matrix(0,nrow=3,ncol=2)
rownames(aic)<-c("CT/A+CT","CT/A","CT")
colnames(aic)<-c("RT-qPCR","RNA-Seq")
aic[1,1]<-fit.overall.pcr$aic
aic[1,2]<-fit.overall.ngs$aic
aic[2,1]<-fit.cta.pcr$aic
aic[2,2]<-fit.cta.ngs$aic
aic[3,1]<-fit.ct.pcr$aic
aic[3,2]<-fit.ct.ngs$aic

dfplot<-data.frame(Arm=factor(rep(rownames(aic),2),levels=c("CT/A+CT","CT/A","CT")),
                   Platform=factor(rep(colnames(aic),each=3),levels=c("RT-qPCR","RNA-Seq")),
                   AIC=c(aic[,1],aic[,2]))
dfplot<-dfplot[order(dfplot$Arm,dfplot$Platform),]
dfplot$Position<-c(1,2,4,5,7,8)

g7<-ggplot(data=dfplot, aes(x=Position, y=AIC, fill=Platform)) +
  geom_bar(stat="identity", width=0.99, position = position_dodge(width=0.9)) +
  theme_pubr(border=F, legend="bottom") +
  xlab("") +
  ylab("AIC") +
  scale_x_continuous(breaks=seq(1.5,7.5,3), labels=c("CT/A+CT","CT/A","CT")) +
  scale_y_continuous(breaks=seq(0,300,50),limits = c(0,310)) +
  scale_fill_manual(name="",values = c("grey60","grey80")) +
  theme(plot.margin=margin(t=10),
        axis.title = element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=4))

myplot7 <- arrangeGrob(g7,top = textGrob("E",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))

pdf(paste0(outputPath, "/Supplementary_Figure_4.pdf"),width=10.8,height=5)
grid.arrange(myplot1,myplot2,myplot3,myplot4,myplot5,myplot6,myplot7,
             layout_matrix=layout.matrix,
             widths=unit(c(4.6,4.6,9.5,8), c("cm")))
dev.off()

#********************************************************
# Figure S5 - TNBC types vs IO score
#********************************************************

# RT-qPCR IO score

df.temp<-df[!is.na(df$IO_cat_qPCR) & df$TNBCtypes!="ND",]
df.temp$TNBCtypes<-droplevels(df.temp$TNBCtypes)
conta<-table(df.temp$IO_cat_qPCR,df.temp$TNBCtypes)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="MSL"] = cumsum(freq$Freq[freq$Var1=="MSL"])-(freq$Freq[freq$Var1=="MSL"]/2)
freq$lab.pos[freq$Var1=="M"] = cumsum(freq$Freq[freq$Var1=="M"])-(freq$Freq[freq$Var1=="M"]/2)
freq$lab.pos[freq$Var1=="LAR"] = cumsum(freq$Freq[freq$Var1=="LAR"])-(freq$Freq[freq$Var1=="LAR"]/2)
freq$lab.pos[freq$Var1=="BL2"] = cumsum(freq$Freq[freq$Var1=="BL2"])-(freq$Freq[freq$Var1=="BL2"]/2)
freq$lab.pos[freq$Var1=="BL1"] = cumsum(freq$Freq[freq$Var1=="BL1"])-(freq$Freq[freq$Var1=="BL1"]/2)
freq$Label<-paste0(round(freq$Freq,0),"%")
freq$Label[freq$Label=="0%"]<-""
test<-format(fisher.test(conta,simulate.p.value=FALSE, workspace = 2000000)$p.value,scientific = T, digits=2)

g1<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=mycols[c(1,2)]) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="bottom") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "white",size=2.5) +
  annotate('text', 2, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=1, size=2.5) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme(axis.title = element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))

# multivariate

df.temp<-df[!is.na(df$IO_cat_qPCR) & df$TNBCtypes!="ND",]

fit1<-glm(pCR.num~IO_cat_qPCR+TNBCtypes, data=df.temp[df.temp$Arm=="CT/A",], family="binomial")
fit2<-glm(pCR.num~IO_cat_qPCR+TNBCtypes, data=df.temp[df.temp$Arm=="CT",], family="binomial")
df.forest<-data.frame(Variable=rep(names(fit1$coefficients),2),
                      OR=c(fit1$coefficients, fit2$coefficients),
                      CIlow=c(confint(fit1)[,1], confint(fit2)[,1]),
                      CIhi=c(confint(fit1)[,2], confint(fit2)[,2]),
                      Arm=rep(c("CT/A","CT"), each=6))
df.forest<-df.forest[df.forest$Variable!="(Intercept)",]
df.forest$Variable[df.forest$Variable=="IO_cat_qPCRIO+"]<-"IO+"
df.forest$Variable<-gsub("TNBCtypes","",df.forest$Variable)
df.forest$Variable<-factor(df.forest$Variable, levels=rev(unique(df.forest$Variable)))

g2<-ggplot(data=df.forest, aes(x=Variable, y=OR, color=Arm)) +
  geom_pointrange(aes(ymin=CIlow, ymax=CIhi), size=0.3)+
  geom_hline(yintercept =0, linetype=2) +
  facet_wrap(~Arm) +
  xlab("") +
  ylab("log(OR)") +
  coord_flip() +
  theme_pubr(border=T, legend="none") +
  scale_color_manual(values=carto_pal(12,"Bold")[7:8]) +
  theme(axis.title = element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8))


# RNA-Seq IO score

df.temp<-df[df$TNBCtypes!="ND",]
df.temp$TNBCtypes<-droplevels(df.temp$TNBCtypes)
conta<-table(df.temp$IO_cat_NGS,df.temp$TNBCtypes)
freq<-round(t(conta)/colSums(conta)*100,1)
freq<-data.frame(freq)
freq<-freq[order(freq$Var1,freq$Var2,decreasing=T),]
freq$lab.pos[freq$Var1=="ND"] = cumsum(freq$Freq[freq$Var1=="ND"])-(freq$Freq[freq$Var1=="ND"]/2)
freq$lab.pos[freq$Var1=="MSL"] = cumsum(freq$Freq[freq$Var1=="MSL"])-(freq$Freq[freq$Var1=="MSL"]/2)
freq$lab.pos[freq$Var1=="M"] = cumsum(freq$Freq[freq$Var1=="M"])-(freq$Freq[freq$Var1=="M"]/2)
freq$lab.pos[freq$Var1=="LAR"] = cumsum(freq$Freq[freq$Var1=="LAR"])-(freq$Freq[freq$Var1=="LAR"]/2)
freq$lab.pos[freq$Var1=="BL2"] = cumsum(freq$Freq[freq$Var1=="BL2"])-(freq$Freq[freq$Var1=="BL2"]/2)
freq$lab.pos[freq$Var1=="BL1"] = cumsum(freq$Freq[freq$Var1=="BL1"])-(freq$Freq[freq$Var1=="BL1"]/2)
freq$Label<-paste0(round(freq$Freq,0),"%")
freq$Label[freq$Label=="0%"]<-""
test<-format(fisher.test(conta,simulate.p.value=FALSE, workspace = 20000000)$p.value,scientific = T, digits=2)

g3<-ggplot(data=freq,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="", values=mycols[c(1,2)]) +
  xlab("") +
  ylab("% of patients") +
  theme_pubr(legend="bottom") +
  geom_text(aes(x=Var1, y= lab.pos, label = Label), color = "white",size=2.5) +
  annotate('text', 2, 106, 
           label=paste("italic(P)==", test), 
           parse=TRUE, 
           hjust=1, size=2.5) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme(axis.title = element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.position="bottom",
        legend.margin=margin(t=-14,b=0))

# multivariate

df.temp<-df[df$TNBCtypes!="ND",]

fit1<-glm(pCR.num~IO_cat_NGS+TNBCtypes, data=df.temp[df.temp$Arm=="CT/A",], family="binomial")
fit2<-glm(pCR.num~IO_cat_NGS+TNBCtypes, data=df.temp[df.temp$Arm=="CT",], family="binomial")
df.forest<-data.frame(Variable=rep(names(fit1$coefficients),2),
                      OR=c(fit1$coefficients, fit2$coefficients),
                      CIlow=c(confint(fit1)[,1], confint(fit2)[,1]),
                      CIhi=c(confint(fit1)[,2], confint(fit2)[,2]),
                      Arm=rep(c("CT/A","CT"), each=6))
df.forest<-df.forest[df.forest$Variable!="(Intercept)",]
df.forest$Variable[df.forest$Variable=="IO_cat_NGSIO+"]<-"IO+"
df.forest$Variable<-gsub("TNBCtypes","",df.forest$Variable)
df.forest$Variable<-factor(df.forest$Variable, levels=rev(unique(df.forest$Variable)))

g4<-ggplot(data=df.forest, aes(x=Variable, y=OR, color=Arm)) +
  geom_pointrange(aes(ymin=CIlow, ymax=CIhi),size=0.3)+
  geom_hline(yintercept =0, linetype=2) +
  facet_wrap(~Arm) +
  xlab("") +
  ylab("log(OR)") +
  coord_flip() +
  theme_pubr(border=T, legend="none") +
  scale_color_manual(values=carto_pal(12,"Bold")[7:8]) +
  theme(axis.title = element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8))


myplot1 <- arrangeGrob(g1,top = textGrob("A",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot3 <- arrangeGrob(g2,top = textGrob("C",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot2 <- arrangeGrob(g3,top = textGrob("B",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))
myplot4 <- arrangeGrob(g4,top = textGrob("D",
                                         x = unit(0, "npc"),
                                         y = unit(0.2, "npc"),
                                         just=c("left","top"),
                                         gp=gpar(col="black",fontsize=10,fontface="bold")))

pdf(paste0(outputPath, "/Supplementary_Figure_5.pdf"),width=7,height=5)
grid.arrange(myplot1,myplot3, myplot2, myplot4)
dev.off()

#*************************************************************
# Supplementary Table 2 - Comparison RT-qPCR - NGS cohorts
#*************************************************************

suppTab2<-data.frame(Characteristic=rep("", 30), RNASeq_cohort="", RTqPCR_cohort="", pValue="")
tmp<-table(df$Clinical_stage, useNA="ifany")
tmpPerc<-round(table(df$Clinical_stage, useNA="ifany")/nrow(df)*100,1)
suppTab2$Characteristic[1]<-"Clinical stage"
suppTab2$Characteristic[2:3]<-names(tmp)
suppTab2$RNASeq_cohort[2:3]<-paste0(tmp, " (",tmpPerc,"%)")
tmp<-table(df2$Clinical_stage, useNA="ifany")
tmpPerc<-round(table(df2$Clinical_stage, useNA="ifany")/nrow(df2)*100,1)
suppTab2$RTqPCR_cohort[2:3]<-paste0(tmp, " (",tmpPerc,"%)")
suppTab2$pValue[2]<-round(chisq.test(cbind(table(df$Clinical_stage), table(df2$Clinical_stage)))$p.value,3)

suppTab2$Characteristic[4]<-"PD-L1 (IHC)"
tmp<-table(df$PDL1_IHC, useNA="ifany")
tmpPerc<-round(table(df$PDL1_IHC, useNA="ifany")/nrow(df)*100,1)
suppTab2$Characteristic[5:6]<-c("Negative", "Positive")
suppTab2$RNASeq_cohort[5:6]<-paste0(tmp, " (",tmpPerc,"%)")
tmp<-table(df2$PDL1_IHC, useNA="ifany")
tmpPerc<-round(table(df2$PDL1_IHC, useNA="ifany")/nrow(df2)*100,1)
suppTab2$RTqPCR_cohort[5:6]<-paste0(tmp, " (",tmpPerc,"%)")
suppTab2$pValue[5]<-round(chisq.test(cbind(table(df$PDL1_IHC), table(df2$PDL1_IHC)))$p.value,3)

suppTab2$Characteristic[7]<-"sTILs (IHC)"
tmp<-table(df$sTILs_IHC, useNA="ifany")
tmpPerc<-round(table(df$sTILs_IHC, useNA="ifany")/nrow(df)*100,1)
suppTab2$Characteristic[8:10]<-names(tmp)
suppTab2$RNASeq_cohort[8:10]<-paste0(tmp, " (",tmpPerc,"%)")
tmp<-table(df2$sTILs_IHC, useNA="ifany")
tmpPerc<-round(table(df2$sTILs_IHC, useNA="ifany")/nrow(df2)*100,1)
suppTab2$RTqPCR_cohort[8:10]<-paste0(tmp, " (",tmpPerc,"%)")
suppTab2$pValue[8]<-round(chisq.test(cbind(table(df$sTILs_IHC), table(df2$sTILs_IHC)))$p.value,3)

suppTab2$Characteristic[11]<-"Treatment arm"
tmp<-table(df$Arm, useNA="ifany")
tmpPerc<-round(table(df$Arm, useNA="ifany")/nrow(df)*100,1)
suppTab2$Characteristic[12:13]<-names(tmp)
suppTab2$RNASeq_cohort[12:13]<-paste0(tmp, " (",tmpPerc,"%)")
tmp<-table(df2$Arm, useNA="ifany")
tmpPerc<-round(table(df2$Arm, useNA="ifany")/nrow(df2)*100,1)
suppTab2$RTqPCR_cohort[12:13]<-paste0(tmp, " (",tmpPerc,"%)")
suppTab2$pValue[12]<-round(chisq.test(cbind(table(df$Arm), table(df2$Arm)))$p.value,3)

suppTab2$Characteristic[14]<-"Response"
tmp<-table(df$pCR, useNA="ifany")
tmpPerc<-round(table(df$pCR, useNA="ifany")/nrow(df)*100,1)
suppTab2$Characteristic[15:16]<-names(tmp)
suppTab2$RNASeq_cohort[15:16]<-paste0(tmp, " (",tmpPerc,"%)")
tmp<-table(df2$pCR, useNA="ifany")
tmpPerc<-round(table(df2$pCR, useNA="ifany")/nrow(df2)*100,1)
suppTab2$RTqPCR_cohort[15:16]<-paste0(tmp, " (",tmpPerc,"%)")
suppTab2$pValue[15]<-round(chisq.test(cbind(table(df$pCR), table(df2$pCR)))$p.value,3)

suppTab2$Characteristic[17]<-"IO score (NGS)"
tmp<-table(df$IO_cat_NGS, useNA="ifany")
tmpPerc<-round(table(df$IO_cat_NGS, useNA="ifany")/nrow(df)*100,1)
suppTab2$Characteristic[18:19]<-names(tmp)
suppTab2$RNASeq_cohort[18:19]<-paste0(tmp, " (",tmpPerc,"%)")
tmp<-table(df2$IO_cat_NGS, useNA="ifany")
tmpPerc<-round(table(df2$IO_cat_NGS, useNA="ifany")/nrow(df2)*100,1)
suppTab2$RTqPCR_cohort[18:19]<-paste0(tmp, " (",tmpPerc,"%)")
suppTab2$pValue[18]<-round(chisq.test(cbind(table(df$IO_cat_NGS), table(df2$IO_cat_NGS)))$p.value,3)

suppTab2$Characteristic[20]<-"IO score (RT-qPCR)"
tmp<-table(df$IO_cat_qPCR, useNA="ifany")
tmpPerc<-round(table(df$IO_cat_qPCR, useNA="ifany")/nrow(df)*100,1)
suppTab2$Characteristic[21:23]<-names(tmp)
suppTab2$RNASeq_cohort[21:23]<-paste0(tmp, " (",tmpPerc,"%)")
tmp<-table(df2$IO_cat_qPCR, useNA="ifany")
tmpPerc<-round(table(df2$IO_cat_qPCR, useNA="ifany")/nrow(df2)*100,1)
suppTab2$RTqPCR_cohort[21:23]<-c(paste0(tmp, " (",tmpPerc,"%)"),0)
suppTab2$pValue[21]<-round(chisq.test(cbind(table(df$IO_cat_qPCR), table(df2$IO_cat_qPCR)))$p.value,3)

suppTab2$Characteristic[24]<-"TNBC subtypes"
tmp<-table(df$TNBCtypes, useNA="ifany")
tmpPerc<-round(table(df$TNBCtypes, useNA="ifany")/nrow(df)*100,1)
suppTab2$Characteristic[25:30]<-names(tmp)
suppTab2$RNASeq_cohort[25:30]<-paste0(tmp, " (",tmpPerc,"%)")
tmp<-table(df2$TNBCtypes, useNA="ifany")
tmpPerc<-round(table(df2$TNBCtypes, useNA="ifany")/nrow(df2)*100,1)
suppTab2$RTqPCR_cohort[25:30]<-paste0(tmp, " (",tmpPerc,"%)")
suppTab2$pValue[25]<-round(chisq.test(cbind(table(df$TNBCtypes), table(df2$TNBCtypes)))$p.value,3)
write.xlsx(suppTab2, file=paste0(outputPath, "/Supplementary_Table_2.xlsx"), rowNames=F)


#*************************************************************
# Table 1 - IO Score RT-qPCR
#*************************************************************

# Univariate

df2<-df[!is.na(df$IO_cat_qPCR),]
fit.overall<-glm(pCR.num~IO_cat_qPCR, family =binomial, data=df2)
fit.cta<-glm(pCR.num~IO_cat_qPCR, family =binomial, data=df2[df2$Arm=="CT/A",])
fit.ct<-glm(pCR.num~IO_cat_qPCR, family =binomial, data=df2[df2$Arm=="CT",])
fit.interact<-glm(pCR.num~IO_cat_qPCR*Arm, family =binomial, data=df2)

res.overall<-data.frame(Odds_Overall = round(exp(coef(fit.overall))[2],2), CI2.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[2,1],2),
                        CI97.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[2,2],2),P_Overall = round(summary(fit.overall)$coefficients[2,4],4))
rownames(res.overall)<-c("IO+ vs IO-")

res.cta<-data.frame(Odds_CTA = round(exp(coef(fit.cta))[2],2), CI2.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[2,1],2),
                    CI97.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[2,2],2),P_CTA = round(summary(fit.cta)$coefficients[2,4],4))
rownames(res.cta)<-c("IO+ vs IO-")

res.ct<-data.frame(Odds_CT = round(exp(coef(fit.ct))[2],2), CI2.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[2,1],2),
                   CI97.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[2,2],2),P_CT = round(summary(fit.ct)$coefficients[2,4],4))
rownames(res.ct)<-c("IO+ vs IO-")

res.uni<-data.frame(Analysis="Univariate",cbind(res.overall,res.cta,res.ct))
res.uni$P_interaction<-round(summary(fit.interact)$coefficients[4,4],4)
res.uni$Variable<-rownames(res.uni)

# Multivariate

df2$Arm<-factor(df2$Arm,levels=c("CT","CT/A"))

fit.overall<-glm(pCR.num~IO_cat_qPCR+sTILs_IHC+PDL1_IHC, family =binomial, data=df2)
fit.cta<-glm(pCR.num~IO_cat_qPCR+sTILs_IHC+PDL1_IHC, family =binomial, data=df2[df2$Arm=="CT/A",])
fit.ct<-glm(pCR.num~IO_cat_qPCR+sTILs_IHC+PDL1_IHC, family =binomial, data=df2[df2$Arm=="CT",])
fit.interact<-glm(pCR.num~IO_cat_qPCR*Arm+sTILs_IHC+PDL1_IHC, family =binomial, data=df2)

res.overall<-as.data.frame(cbind(round(exp(cbind("Odds ratio" = coef(fit.overall), confint.default(fit.overall, level = 0.95))),2), "P-value" = round(summary(fit.overall)$coefficients[,4],4))[-1,])
colnames(res.overall)<-paste0(c("Odds", "CI2.5","CI97.5","P"),"_Overall")
rownames(res.overall)<-c("IO+ vs IO-", "sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.cta<-as.data.frame(cbind(round(exp(cbind("Odds ratio" = coef(fit.cta), confint.default(fit.cta, level = 0.95))),2), "P-value" = round(summary(fit.cta)$coefficients[,4],4))[-1,])
colnames(res.cta)<-paste0(c("Odds", "CI2.5","CI97.5","P"),"_CTA")
rownames(res.cta)<-c("IO+ vs IO-", "sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.ct<-as.data.frame(cbind(round(exp(cbind("Odds ratio" = coef(fit.ct), confint.default(fit.ct, level = 0.95))),2), "P-value" = round(summary(fit.ct)$coefficients[,4],4))[-1,])
colnames(res.ct)<-paste0(c("Odds", "CI2.5","CI97.5","P"),"_CT")
rownames(res.ct)<-c("IO+ vs IO-", "sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.multi<-data.frame(Analysis="Multivariable",cbind(res.overall,res.cta,res.ct))
res.multi$P_interaction<-c(round(summary(fit.interact)$coefficients[6,4],4),"","")
res.multi$Variable<-rownames(res.multi)
res<-rbind(res.uni,res.multi)
res<-res[,c(1,15,2:14)]
write.xlsx(res,file=paste0(outputPath,"/Table_1.xlsx"), rowNames=F)


#*************************************************************
# Supplementary Table 3 - IO Score NGS
#*************************************************************

# Univariate

fit.overall<-glm(pCR.num~IO_cat_NGS, family =binomial, data=df)
fit.cta<-glm(pCR.num~IO_cat_NGS, family =binomial, data=df[df$Arm=="CT/A",])
fit.ct<-glm(pCR.num~IO_cat_NGS, family =binomial, data=df[df$Arm=="CT",])
fit.interact<-glm(pCR.num~IO_cat_NGS*Arm, family =binomial, data=df)

res.overall<-data.frame(Odds_Overall = round(exp(coef(fit.overall))[2],2), CI2.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[2,1],2),
                        CI97.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[2,2],2),P_Overall = round(summary(fit.overall)$coefficients[2,4],4))
rownames(res.overall)<-c("IO+ vs IO-")

res.cta<-data.frame(Odds_CTA = round(exp(coef(fit.cta))[2],2), CI2.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[2,1],2),
                    CI97.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[2,2],2),P_CTA = round(summary(fit.cta)$coefficients[2,4],4))
rownames(res.cta)<-c("IO+ vs IO-")

res.ct<-data.frame(Odds_CT = round(exp(coef(fit.ct))[2],2), CI2.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[2,1],2),
                   CI97.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[2,2],2),P_CT = round(summary(fit.ct)$coefficients[2,4],4))
rownames(res.ct)<-c("IO+ vs IO-")

res.uni<-data.frame(Analysis="Univariate",cbind(res.overall,res.cta,res.ct))
res.uni$P_interaction<-round(summary(fit.interact)$coefficients[4,4],4)
res.uni$Variable<-rownames(res.uni)

# Multivariable

df$Arm<-factor(df$Arm,levels=c("CT","CT/A"))

fit.overall<-glm(pCR.num~IO_cat_NGS+sTILs_IHC+PDL1_IHC, family =binomial, data=df)
fit.cta<-glm(pCR.num~IO_cat_NGS+sTILs_IHC+PDL1_IHC, family =binomial, data=df[df$Arm=="CT/A",])
fit.ct<-glm(pCR.num~IO_cat_NGS+sTILs_IHC+PDL1_IHC, family =binomial, data=df[df$Arm=="CT",])
fit.interact<-glm(pCR.num~IO_cat_NGS*Arm+sTILs_IHC+PDL1_IHC, family =binomial, data=df)

res.overall<-as.data.frame(cbind(round(exp(cbind("Odds ratio" = coef(fit.overall), confint.default(fit.overall, level = 0.95))),2), "P-value" = round(summary(fit.overall)$coefficients[,4],4))[-1,])
colnames(res.overall)<-paste0(c("Odds", "CI2.5","CI97.5","P"),"_Overall")
rownames(res.overall)<-c("IO+ vs IO-", "sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.cta<-as.data.frame(cbind(round(exp(cbind("Odds ratio" = coef(fit.cta), confint.default(fit.cta, level = 0.95))),2), "P-value" = round(summary(fit.cta)$coefficients[,4],4))[-1,])
colnames(res.cta)<-paste0(c("Odds", "CI2.5","CI97.5","P"),"_CTA")
rownames(res.cta)<-c("IO+ vs IO-", "sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.ct<-as.data.frame(cbind(round(exp(cbind("Odds ratio" = coef(fit.ct), confint.default(fit.ct, level = 0.95))),2), "P-value" = round(summary(fit.ct)$coefficients[,4],4))[-1,])
colnames(res.ct)<-paste0(c("Odds", "CI2.5","CI97.5","P"),"_CT")
rownames(res.ct)<-c("IO+ vs IO-", "sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.multi<-data.frame(Analysis="Multivariate",cbind(res.overall,res.cta,res.ct))
res.multi$P_interaction<-c(round(summary(fit.interact)$coefficients[6,4],4),"","")
res.multi$Variable<-rownames(res.multi)
res<-rbind(res.uni,res.multi)
res<-res[,c(1,15,2:14)]
write.xlsx(res,file=paste0(outputPath,"/Supplementary_Table_3.xlsx"), rowNames=F)

#*************************************************************
# Table 2 - TNBC subtypes
#*************************************************************

df2<-df[df$TNBCtypes!="ND",]

# univariate analysis (LAR as reference)

fit.overall<-glm(pCR.num~TNBCtypes, family =binomial, data=df2)
fit.cta<-glm(pCR.num~TNBCtypes, family =binomial, data=df2[df2$Arm=="CT/A",])
fit.ct<-glm(pCR.num~TNBCtypes, family =binomial, data=df2[df2$Arm=="CT",])
fit.interact<-glm(pCR.num~TNBCtypes*Arm, family =binomial, data=df2)

res.overall<-data.frame(Odds_Overall = round(exp(coef(fit.overall))[-1],2), CI2.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[-1,1],2),
                        CI97.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[-1,2],2),P_Overall = round(summary(fit.overall)$coefficients[-1,4],4))
rownames(res.overall)<-paste(gsub("TNBCtypes","",rownames(res.overall)),"vs LAR")

res.cta<-data.frame(Odds_CTA = round(exp(coef(fit.cta))[-1],2), CI2.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[-1,1],2),
                    CI97.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[-1,2],2),P_CTA = round(summary(fit.cta)$coefficients[-1,4],4))
rownames(res.cta)<-paste(gsub("TNBCtypes","",rownames(res.cta)),"vs LAR")

res.ct<-data.frame(Odds_CT = round(exp(coef(fit.ct))[-1],2), CI2.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[-1,1],2),
                   CI97.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[-1,2],2),P_CT = round(summary(fit.ct)$coefficients[-1,4],4))
rownames(res.ct)<-paste(gsub("TNBCtypes","",rownames(res.ct)),"vs LAR")

res.uni<-data.frame(Analysis="Univariate",cbind(res.overall,res.cta,res.ct))
res.uni$P_interaction<-round(summary(fit.interact)$coefficients[7:10,4],4)
res.uni$Variable<-rownames(res.uni)

# multivariate analysis (LAR as reference)

fit.overall<-glm(pCR.num~TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2)
fit.cta<-glm(pCR.num~TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2[df2$Arm=="CT/A",])
fit.ct<-glm(pCR.num~TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2[df2$Arm=="CT",])
fit.interact<-glm(pCR.num~TNBCtypes*Arm+sTILs_IHC+PDL1_IHC, family =binomial, data=df2)

res.overall<-data.frame(Odds_Overall = round(exp(coef(fit.overall))[-1],2), CI2.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[-1,1],2),
                        CI97.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[-1,2],2),P_Overall = round(summary(fit.overall)$coefficients[-1,4],4))
rownames(res.overall)<-paste(gsub("TNBCtypes","",rownames(res.overall)), "vs LAR")

res.cta<-data.frame(Odds_CTA = round(exp(coef(fit.cta))[-1],2), CI2.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[-1,1],2),
                    CI97.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[-1,2],2),P_CTA = round(summary(fit.cta)$coefficients[-1,4],4))
rownames(res.cta)<-paste(gsub("TNBCtypes","",rownames(res.cta)), "vs LAR")

res.ct<-data.frame(Odds_CT = round(exp(coef(fit.ct))[-1],2), CI2.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[-1,1],2),
                   CI97.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[-1,2],2),P_CT = round(summary(fit.ct)$coefficients[-1,4],4))
rownames(res.ct)<-paste(gsub("TNBCtypes","",rownames(res.ct)), "vs LAR")

res.multi<-data.frame(Analysis="Multivariate",cbind(res.overall,res.cta,res.ct))
res.multi$P_interaction<-c(round(summary(fit.interact)$coefficients[9:12,4],4),"","")
res.multi$Variable<-rownames(res.multi)
res<-rbind(res.uni,res.multi)
res<-res[,c(1,15,2:14)]
rownames(res)[grep("sTILs",rownames(res))]<-"sTILs-High vs sTILs-Low/Mid"
rownames(res)[grep("PDL1",rownames(res))]<-"PD-L1+ vs PD-L1-"
write.xlsx(res,file=paste0(outputPath, "/Table_2.xlsx"), rowNames=F)

#*************************************************************
# Supplementary Table 4 - TNBC subtypes (one vs all)
#*************************************************************

# univariate analysis (one-vs-all)

res.overall<-data.frame(Subtype=unique(df2$TNBCtypes),Odds_Overall = 0, CI2.5_Overall=0,CI97.5_Overall=0,P_Overall = 0)
for(i in unique(df2$TNBCtypes)){
  TNBCtmp<-factor(ifelse(df2$TNBCtypes==i, i,"Other"),levels=c("Other",i))
  fit<-glm(pCR.num~TNBCtmp, family =binomial, data=df2)
  res.overall[res.overall$Subtype==i,-1]<-c(round(exp(coef(fit))[-1],2), round(exp(confint.default(fit, level = 0.95))[-1,1],2),
                                            round(exp(confint.default(fit, level = 0.95))[-1,2],2),round(summary(fit)$coefficients[-1,4],3))
}
res.cta<-data.frame(Subtype=unique(df2$TNBCtypes),Odds_CTA = 0, CI2.5_CTA=0,CI97.5_CTA=0,P_CTA = 0)
for(i in unique(df2$TNBCtypes)){
  df2tmp<-df2[df2$Arm=="CT/A",]
  TNBCtmp<-factor(ifelse(df2tmp$TNBCtypes==i, i,"Other"),levels=c("Other",i))
  fit<-glm(pCR.num~TNBCtmp, family =binomial, data=df2tmp)
  res.cta[res.cta$Subtype==i,-1]<-c(round(exp(coef(fit))[-1],2), round(exp(confint.default(fit, level = 0.95))[-1,1],2),
                                    round(exp(confint.default(fit, level = 0.95))[-1,2],2),round(summary(fit)$coefficients[-1,4],3))
}
res.ct<-data.frame(Subtype=unique(df2$TNBCtypes),Odds_CT = 0, CI2.5_CT=0,CI97.5_CT=0,P_CT = 0)
for(i in unique(df2$TNBCtypes)){
  df2tmp<-df2[df2$Arm=="CT",]
  TNBCtmp<-factor(ifelse(df2tmp$TNBCtypes==i, i,"Other"),levels=c("Other",i))
  fit<-glm(pCR.num~TNBCtmp, family =binomial, data=df2tmp)
  res.ct[res.ct$Subtype==i,-1]<-c(round(exp(coef(fit))[-1],2), round(exp(confint.default(fit, level = 0.95))[-1,1],2),
                                  round(exp(confint.default(fit, level = 0.95))[-1,2],2),round(summary(fit)$coefficients[-1,4],3))
}
res.int<-data.frame(Subtype=unique(df2$TNBCtypes),P_Interaction = 0)
for(i in unique(df2$TNBCtypes)){
  TNBCtmp<-factor(ifelse(df2$TNBCtypes==i, i,"Other"),levels=c("Other",i))
  fit<-glm(pCR.num~TNBCtmp*Arm, family =binomial, data=df2)
  res.int$P_Interaction[res.ct$Subtype==i]<-round(summary(fit)$coefficients[4,4],3)
}

res.uni<-data.frame(Analysis="Univariate",cbind(res.overall,res.cta[,-1],res.ct[,-1]))
res.uni$P_interaction<-res.int$P_Interaction
res.uni<-res.uni[order(res.uni$Subtype),]
res.uni$Subtype<-paste(res.uni$Subtype, "vs others")
write.xlsx(res.uni,file=paste0(outputPath, "/Supplementary_Table_4.xlsx"), rowNames=F)

#*************************************************************
# Table 3 - Full model with all variables (IO score RT-qPCR)
#*************************************************************

df2<-df[df$TNBCtypes!="ND" & !is.na(df$IO_cat_qPCR),]
fit.overall<-glm(pCR.num~IO_cat_qPCR*Arm+TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2)
fit.cta<-glm(pCR.num~IO_cat_qPCR+TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2[df2$Arm=="CT/A",])
fit.ct<-glm(pCR.num~IO_cat_qPCR+TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2[df2$Arm=="CT",])

res.overall<-data.frame(Odds_Overall = round(exp(coef(fit.overall))[-1],2), CI2.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[-1,1],2),
                        CI97.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[-1,2],2),P_Overall = round(summary(fit.overall)$coefficients[-1,4],4))
res.overall$Variable<-c("IO+ vs IO-","CT/A vs CT","BL1 vs LAR","BL2 vs LAR","M vs LAR","MSL vs LAR","sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-","IO+:CT/A")

res.cta<-data.frame(Odds_CTA = round(exp(coef(fit.cta))[-1],2), CI2.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[-1,1],2),
                    CI97.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[-1,2],2),P_CTA = round(summary(fit.cta)$coefficients[-1,4],4))
res.cta$Variable<-c("IO+ vs IO-","BL1 vs LAR","BL2 vs LAR","M vs LAR","MSL vs LAR","sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.ct<-data.frame(Odds_CT = round(exp(coef(fit.ct))[-1],2), CI2.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[-1,1],2),
                   CI97.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[-1,2],2),P_CT = round(summary(fit.ct)$coefficients[-1,4],4))
res.ct$Variable<-c("IO+ vs IO-","BL1 vs LAR","BL2 vs LAR","M vs LAR","MSL vs LAR","sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.multi<-merge(res.overall,res.cta,by="Variable",all=T)
res.multi<-merge(res.multi,res.ct,by="Variable",all=T)
res.multi<-res.multi[c(4,3,1,2,6,7,9,8,5),]
res.multi$P_Interaction<-c(res.multi$P_Overall[nrow(res.multi)], rep("", nrow(res.multi)-1))
res.multi<-res.multi[-nrow(res.multi),]
write.xlsx(res.multi,file=paste0(outputPath, "/Table_3.xlsx"), rowNames=F)

#*************************************************************
# Table 3 - Full model with all variables (IO score NGS)
#*************************************************************

df2<-df[df$TNBCtypes!="ND",]
fit.overall<-glm(pCR.num~IO_cat_NGS*Arm+TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2)
fit.cta<-glm(pCR.num~IO_cat_NGS+TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2[df2$Arm=="CT/A",])
fit.ct<-glm(pCR.num~IO_cat_NGS+TNBCtypes+sTILs_IHC+PDL1_IHC, family =binomial, data=df2[df2$Arm=="CT",])

res.overall<-data.frame(Odds_Overall = round(exp(coef(fit.overall))[-1],2), CI2.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[-1,1],2),
                        CI97.5_Overall=round(exp(confint.default(fit.overall, level = 0.95))[-1,2],2),P_Overall = round(summary(fit.overall)$coefficients[-1,4],4))
res.overall$Variable<-c("IO+ vs IO-","CT/A vs CT","BL1 vs LAR","BL2 vs LAR","M vs LAR","MSL vs LAR","sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-","IO+:CT/A")

res.cta<-data.frame(Odds_CTA = round(exp(coef(fit.cta))[-1],2), CI2.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[-1,1],2),
                    CI97.5_CTA=round(exp(confint.default(fit.cta, level = 0.95))[-1,2],2),P_CTA = round(summary(fit.cta)$coefficients[-1,4],4))
res.cta$Variable<-c("IO+ vs IO-","BL1 vs LAR","BL2 vs LAR","M vs LAR","MSL vs LAR","sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.ct<-data.frame(Odds_CT = round(exp(coef(fit.ct))[-1],2), CI2.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[-1,1],2),
                   CI97.5_CT=round(exp(confint.default(fit.ct, level = 0.95))[-1,2],2),P_CT = round(summary(fit.ct)$coefficients[-1,4],4))
res.ct$Variable<-c("IO+ vs IO-","BL1 vs LAR","BL2 vs LAR","M vs LAR","MSL vs LAR","sTILs-High vs sTILs-Low/Mid","PD-L1+ vs PD-L1-")

res.multi<-merge(res.overall,res.cta,by="Variable",all=T)
res.multi<-merge(res.multi,res.ct,by="Variable",all=T)
res.multi<-res.multi[c(4,3,1,2,6,7,9,8,5),]
res.multi$P_Interaction<-c(res.multi$P_Overall[nrow(res.multi)], rep("", nrow(res.multi)-1))
res.multi<-res.multi[-nrow(res.multi),]
write.xlsx(res.multi,file=paste0(outputPath, "/Supplementary_Table_5.xlsx"), rowNames=F)


