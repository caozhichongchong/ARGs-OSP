###---------------------format MGD mothertable---------------------------------
Metadata=read.delim('Metadata_habitat.txt',header=T) #Table S2
ARG=read.delim('SARGv1.0',header=T)
temp = list.files(pattern="MGD.ARG.blastx.id50.hit50*")
for (m in c(2:(length(temp))))
{
  Unique=read.delim(temp[m],header=F)
  Unique=Unique[rev(order(as.numeric(as.matrix(Unique$V13)))),]
  Unique=subset(Unique, !duplicated(Unique$V2))
  write.table(Unique,paste('MGD.ARG.blastx.id50.hit50.new1.',m,sep = ''),sep='\t',row.names=F,quote=F,col.names=F)
  Unique=Unique[,c(1,3,4,5,12)]
  Unique=as.matrix(Unique)
  Unique[,3]=as.integer(as.numeric(Unique[,3]))
  Unique=Unique[which(as.numeric(Unique[,3])>=50),]
  Unique=Unique[which(as.numeric(Unique[,4])>=17),]
  Unique=data.frame(Unique)
  Unique=merge(Unique,Metadata,by.x='V1',by.y='sampleID2',all.x=T)
  Unique=Unique[,c(1:5,7,8)]
  Unique=merge(Unique,ARG[,1:3],by.x='V3',by.y='seqID1',all.x=T)
  Unique=Unique[,c(1,8,9,2,7,6,3,4,5)]
  Unique=as.matrix(Unique)
  Unique[which(is.na(Unique[,4])),4]='None'
  write.table(Unique,paste('MGD.ARG.blastx.id50.hit50.new2.',m,sep = ''),sep='\t',row.names=F,quote=F,col.names=F)
}


###---------------MGD world map-------------------
library(ggplot2)  # FYI you need v2.0
library(dplyr)    # yes, i could have not done this and just used 'subset' instead of 'filter'
library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
library(ggthemes) # theme_map and tableau colors
dat=read.delim('Metadata_location.txt',header=T)
library(plyr)
dat=count(dat, c("eco_type","lat","long"))
world <- map_data("world")
#plot
gg <- ggplot()
gg <- gg + geom_map(data=world, map=world,
                    aes(x=long, y=lat, map_id=region),
                    color="white", fill="#7f7f7f", size=0.05, alpha=1/4)

gg <- gg + geom_point(data=dat, 
                      aes(x=as.numeric(as.matrix(dat$long)), 
                          y=as.numeric(as.matrix(dat$lat)), 
                          color=eco_type,
                          size=as.numeric(as.matrix(dat$freq))), 
                      alpha=0.8)
gg <- gg + scale_color_brewer(palette="Set1")
gg <- gg + theme(strip.background=element_blank())
gg


###--------------WGD rarefaction curve-----------------------------------------
Mt=read.delim('WGD_ARG.mothertable',header=T)
library(plyr)
library(ggplot2)
Count1=count(Mt, c("ARG","Genome","phylum"))
Count1=merge(Count1,ARG[,c(1,6)],by.x='ARG',by.y='V1',all.x=T)
Count3=count(Ge, c("assembly_accession"))
Count3$Genum=row.names(Count3)
Count1=merge(Count1,Count3[,c(1,3)],by.x='Genome',
             by.y='assembly_accession',all.x=T)
write.table(Count1,'WGD_ARG.mothertable.rarefaction',quote=F,row.names = F,sep='\t')


###------------------WGD taxonomy conserve----------------------------
Mt=read.delim('WGD_ARG.mothertable',header=T)
library(plyr)
Summary=count(data.frame(Mt), c("ARG", "subtype","type","Geneid","Genome","phylum","class","order","family","genus","species","strain",'Pathogen'))
Unique1=subset(data.frame(Summary),!duplicated(Summary[,1]))
Unique1=as.matrix(Unique1)
Unique1[,13:14]=0
for (j in 1:nrow(Summary))
{
  for (i in 1:nrow(Unique1))
    if (Summary[j,1]==Unique1[i,1])
    {Unique1[i,13]=as.numeric(Unique1[i,13])+as.numeric(Summary[j,13])
    Unique1[i,14]=as.numeric(Unique1[i,14])+as.numeric(Summary[j,14])
    for (k in seq(12, 6, by=-1))
    {
      if (Unique1[i,k]!=Summary[j,k])
      {
        Unique1[i,k]='None'
      }
    }
    break
    }
}
for (i in 1:nrow(Unique1))
  for (k in seq(12, 6, by=-1))
  {
    if(Unique1[i,k]!='None' && !grepl('_',Unique1[i,k]))
      break
    if (grepl('_',Unique1[i,k]))
    {
      Unique1[i,k]='None'
    }
  }
Unique=Unique1
write.table(Unique,'WGD_ARG.mothertable.conserve.txt',row.names = F,sep='\t',quote=F)


###--------------MGD cooccurrence of ARG and IntI1------------------------------
Metadata=read.delim('Metadata_habitat.txt',header=T)
Mt=read.delim('MGD_ARG_cellnumber.mothertable.txt',header=T)
Mt2=read.delim('MGD_IntI1_cellnumber.mothertable.txt',header=T)
library(plyr)
Count1=subset(data.frame(Mt$sampleID,Mt$Total_ARGpercell))
Count2=subset(data.frame(Mt2$sampleID,Mt2$Total_Intpercell))
Count1=data.frame(Count1)
Mt=data.frame(Mt)
Count1=Count1[-c(1:3),]
Count2=Count2[-c(1),]
Count=merge(Count1[which(Count1$Mt.Total_ARGpercell!='None'),],
            Count2[which(Count2$Mt2.Total_Intpercell!='None'),],
            by.x='Mt.sampleID',by.y='Mt2.sampleID',all=T)
Count=merge(Count,Metadata[,c(1,2)],by.x='Mt.sampleID',by.y='sampleID',all.X=T)
Count$Mt.Total_ARGpercell[which(is.na( Count$Mt.Total_ARGpercell))]=0
Count=as.matrix(Count)
Count[which(is.na(Count[,3])),3]=0
Count=data.frame(Count)
Count$Color=Count$eco_type
Types=unique(Count$eco_type)
Count$Mt.Total_ARGpercell=as.numeric(as.matrix(Count$Mt.Total_ARGpercell))
Count$Mt2.Total_Intpercell=as.numeric(as.matrix(Count$Mt2.Total_Intpercell))

#plot
library(RColorBrewer)
library("ggplot2")
library('fBasics')
ggplot(Count,aes(x=log(as.numeric(as.matrix(Count$Mt.Total_ARGpercell)),10),
                 y=log(as.numeric(as.matrix(Count$Mt2.Total_Intpercell)),10))) + 
  geom_point(aes(colour=Count$Color),shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_color_brewer(palette="Set1")+
  scale_x_continuous(limits = c(-4, 1))+
  scale_y_continuous(limits = c(-4, 1))
Countnew=Count
Countnew=Countnew[which(Countnew$Mt.Total_ARGpercell >0),]
Countnew=Countnew[which(Countnew$Mt2.Total_Intpercell >0),]

#linear regression for allsamples
df=data.frame(
  x=log(Countnew$Mt.Total_ARGpercell,10),
  y=log(Countnew$Mt2.Total_Intpercell,10))
m <- lm(y ~ x, df)
format(coef(m)[2], digits = 2)
format(summary(m)$r.squared, digits = 3)

Count$Mt.Total_ARGpercell=as.numeric(as.matrix(Count$Mt.Total_ARGpercell))
Count$Mt2.Total_Intpercell=as.numeric(as.matrix(Count$Mt2.Total_Intpercell))

#linear regression for each eco-type
Regression=matrix(0,nrow=7,ncol=3)
row.names(Regression)=Types[1:7]
colnames(Regression)=c("Inter",'k','R2')
for (i in 1:(length(Types)-1))
{
  Countnew=Count[which(Count$eco_type==Types[i]),]
  Countnew=Countnew[which(Countnew$Mt.Total_ARGpercell >0),]
  Countnew=Countnew[which(Countnew$Mt2.Total_Intpercell >0),]
  df=data.frame(
    x=log(Countnew$Mt.Total_ARGpercell,10),
    y=log(Countnew$Mt2.Total_Intpercell,10))
  m <- lm(y ~ x, df)
  
  Regression[i,1]=format(coef(m)[1], digits = 2)
  Regression[i,2]=format(coef(m)[2], digits = 2)
  Regression[i,3]=format(summary(m)$r.squared, digits = 3)
}
write.table(Regression,
            'MGD_ARG_IntI1_regression.txt',quote=F,sep='\t')

#plot
Count=Countcopy
Count$Color=Count$eco_type
Count$ratio=0
Count$radom=0
Types=unique(Count$eco_type)
Count=as.matrix(Count)
Count[,6]=as.numeric(Count[,2])/as.numeric(Count[,3])
for (i in 1:8)
  Count[which(Count[,4]==Types[i]),7]=runif(length(which(Count[,4]==Types[i])),i-0.5,i+0.5)
Count=data.frame(Count)
Count$ratio=as.numeric(as.matrix(Count$ratio))
Count$radom=as.numeric(as.matrix(Count$radom))


ggplot(Count,aes(x=Count$radom,
                 y=log(Count$ratio+1)))+ ##or x2,y2, or x3,y3
  geom_point(aes(colour=Count$Color),shape=16,size=2,
             alpha = 0.6,stroke=0.1)+ 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_color_brewer(palette="Set1")


###------------------------MGD and WGD rarefaction---------------------------
Rare=read.delim('ARG_rarefaction.average.txt',header=F)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggplot2)
ggplot(Rare,aes(x=(as.numeric(as.matrix(Rare$V1))),
                y=(as.numeric(as.matrix(Rare$V2))))) + ##or x2,y2, or x3,y3
  geom_point(#aes(color=Rare$V3,fill=Rare$V3)
    color=cbPalette[8]
    ,shape=16,size=2,
    alpha = 0.2,stroke=0)+ ###shape=16
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_y_continuous(limits = c(0, 3150),breaks = scales::pretty_breaks(n = 6))+
  scale_color_brewer(palette="Set1")+
  geom_smooth(aes(x = (as.numeric(as.matrix(Rare$V1))), 
                  y = (as.numeric(as.matrix(Rare$V2)))), method = "lm",
              formula = y ~ poly(x, 21), se = FALSE,color=cbPalette[6],size=0.5)

