library(ggplot2)
library(reshape2)

data=read.csv(file='~/Documents/Projects/Cirrus_OCT_Python/Data/etdrs_data.csv',
              colClasses=c('integer',
                           'integer',
                           'character',
                           'character',
                           'character',
                           'numeric',
                           'character',
                           'character'),
              na.strings=c('','--'))


data$group[data$group>0]=1
data$Eye<-'OS'
data$Eye[data$id %in% c('sn2828','sn31475','sn31476','sn31489')] <- 'OS'
data$Eye[data$id %in% c('sn2821','sn2833','sn31469','sn31482')] <- 'OD'

# load the foreign data
dat_f<-read.csv(file='Data/foreignEtdrs.csv',
                stringsAsFactors=FALSE)

# extract only LE data
dat_f <- subset(dat_f,Eye=='OS')

# removed 2016-02-14 changing to display LE data from sickkids subjects
# # swap nasal/temporal where eye is OS
# dat_f_copy <- dat_f
# dat_f_copy[dat_f$Region=='InnerNasal' & dat_f$Eye=='OS','Region']<-'InnerTemporal'
# dat_f_copy[dat_f$Region=='InnerTemporal' & dat_f$Eye=='OS','Region']<-'InnerNasal'
# dat_f_copy[dat_f$Region=='OuterNasal' & dat_f$Eye=='OS','Region']<-'OuterTemporal'
# dat_f_copy[dat_f$Region=='OuterTemporal' & dat_f$Eye=='OS','Region']<-'OuterNasal'
# dat_f_copy$Eye <- 'OD'
# dat_f <- dat_f_copy

# calculate the ring values for the foreign data
 # cast data from long to wide
x<-dcast(dat_f,Eye+Layer~Region,value.var="Thickness", mean)
 #calculate the new regions
x$Inner <- mean(c(x$InnerInferior,x$InnerNasal,x$InnerSuperior,x$InnerTemporal))
x$Outer <- mean(c(x$OuterInferior,x$OuterNasal,x$OuterSuperior,x$OuterTemporal))
x$Macular <- mean(c(x$Inner,x$Outer,x$Fovea))
# recast back to long
x<-melt(x,
          id.vars=c("Eye","Layer"),
          measure.vars=c("Fovea",
                         "InnerInferior","InnerNasal","InnerSuperior","InnerTemporal",
                         "OuterInferior","OuterNasal","OuterSuperior","OuterTemporal",
                         "Inner","Outer","Macular"),
          variable.name="Region",
          value.name="Thickness")
x$Region<-levels(x$Region)[x$Region]
dat_f<-x



dat_f$group<-2
dat_f$id<-'foreign'
dat_f$subject<-'foreign'
dat_f$test<-2
dat_f$Thickness<-as.numeric(dat_f$Thickness)
#calculate the inner + outer complexes
for(region in unique(dat_f$Region)){
  subData<-subset(dat_f,Region==region)
  InnerComplex<-sum(subData[subData$Layer %in% c('GCL','IPL','INL','OPL'),'Thickness'])
  OuterComplex<-sum(subData[subData$Layer %in% c('ONL','PhR'),'Thickness'])
  OuterComplex2<-sum(subData[subData$Layer %in% c('ONL','PhR','OPL'),'Thickness'])
  Total<-sum(subData[subData$Layer %in% c('GCL','IPL','INL','ONL','PhR','OPL'),'Thickness'])
  IC1<-sum(subData[subData$Layer %in% c('GCL','IPL'),'Thickness'])
  IC2<-sum(subData[subData$Layer %in% c('GCL','IPL','INL'),'Thickness'])
  newdf<-data.frame(Eye=rep('OS',6),
                    Layer=c('Inner Complex','Outer Retina','Outer Complex + OPL','Total GCL-Bruchs','GCL + IPL','GCL + IPL + INL'),
                    Region=rep(region,6),
                    Thickness=c(InnerComplex,OuterComplex,OuterComplex2,Total,IC1,IC2),
                    group=rep(2,6),
                    id=rep('foreign',6),
                    subject=rep('foreign',6),
                    test=rep(2,6))
  
  dat_f<-rbind(dat_f,newdf)
}



#merge the datasets
ids<-c(data$id,dat_f$id)
subjects<-c(data$subject,dat_f$subject)
eyes<-c(data$Eye,dat_f$Eye)
groups<-c(data$group,dat_f$group)
tests<-c(data$test,dat_f$test)
layers<-c(data$layer,dat_f$Layer)
regions<-c(data$variable,dat_f$Region)
thicknesses<-c(data$value,dat_f$Thickness)

data<-data.frame(id=ids,
                 subject=subjects,
                 eye=eyes,
                 group=groups,
                 test=tests,
                 layer=layers,
                 region=regions,
                 thickness=thicknesses)

# # Create ring averages
# # cast data from long to wide
# x<-dcast(data,group+subject+eye+test+layer~region,value.var="thickness", mean)
# #calculate the new regions
# x$Inner <- x$InnerInferior + x$InnerNasal + x$InnerSuperior + x$InnerTemporal
# x$Outer <- x$OuterInferior + x$OuterNasal + x$OuterSuperior + x$OuterTemporal
# x$Macular <- x$Inner + x$Outer + x$Fovea
# # recast back to long
# data<-melt(x,
#            id.vars=c("subject","layer","eye","test"),
#            measure.vars=c("Fovea","InnerInferior","InnerNasal","InnerSuperior",
#                           "InnerTemporal","OuterInferior","OuterNasal","OuterSuperior",
#                           "Inner","Outer","Macular"),
#            variable.name="region",
#            value.name="thickness")

data$id<-factor(data$id)
data$subject<-factor(data$subject, 
                     levels=c('2337335','2337333','foreign'),
                     labels=c('Family A, II-1','Family A, II-2', 'Family B, II-2'))
data$eye<-factor(data$eye)
data$group<-factor(data$group,levels=c(0,1,2),labels=c('control','sickkids','foreign'))
data$test<-factor(data$test, levels=c(0,1,2,3))
data$layer<-factor(data$layer)
data$region<-factor(data$region)

data<-subset(data,eye=='OS')

data[data$group=='foreign','test']<-3

#Only subset the ring regions
#data<-subset(data,region %in% c('Fovea','Inner','Outer','Macular'))
#data$region <- factor(data$region,levels=c('Fovea','Inner','Outer','Macular'))
data<-subset(data,region %in% c('Inner','Outer','Macular'))
data$region <- factor(data$region,levels=c('Inner','Outer','Macular'))

data$pos_i = as.numeric(data$region)
data$pos_i[data$test==1 & data$group=='sickkids']=data$pos_i[data$test==1 & data$group=='sickkids']-0.2
data$pos_i[data$test==2 & data$group=='sickkids']=data$pos_i[data$test==2 & data$group=='sickkids']+0.2

# calculate summary statistics on the control data
data_c<-subset(data,group=='control')
x <- ddply(data_c,.(layer,region),summarize,ymin=quantile(thickness,0.05),ymax=quantile(thickness,0.95),med=median(thickness))
x$xmin <- as.numeric(x$region) - 0.5
x$xmax <- as.numeric(x$region) + 0.5

default_theme <- theme_bw() +
  theme(axis.text = element_text(size=12),
        title = element_text(size=12,face='bold'),
        axis.title = element_text(size=14,face='bold'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))

for(layer in levels(data$layer)){
  subdata <- data[data$layer==layer,]
  subdata_x <- x[x$layer==layer,]
  p <- ggplot(subdata)
  #p <- p+geom_boxplot(data=subset(subdata,group=='control'),outlier.shape=NA)
#   p <- p+geom_rect(data=subset(subdata_x,region=='Fovea'),
#                    aes(xmin=0.6,xmax=1.4,ymin=ymin,ymax=ymax),
#                    colour="grey20",
#                    alpha=0.2)
#   p <- p+geom_segment(data=subset(subdata_x,region=='Fovea'),
#                       aes(x=0.6,xend=1.4,y=med,yend=med))
  p <- p+geom_rect(data=subset(subdata_x,region=='Inner'),
                   aes(xmin=0.6,xmax=1.4,ymin=ymin,ymax=ymax),
                   colour="grey20",
                   alpha=0.2)
  p <- p+geom_segment(data=subset(subdata_x,region=='Inner'),
                      aes(x=0.6,xend=1.4,y=med,yend=med))
  p <- p+geom_rect(data=subset(subdata_x,region=='Outer'),
                   aes(xmin=1.6,xmax=2.4,ymin=ymin,ymax=ymax),
                   colour="grey20",
                   alpha=0.2)
  p <- p+geom_segment(data=subset(subdata_x,region=='Outer'),
                      aes(x=1.6,xend=2.4,y=med,yend=med))
  p <- p+geom_rect(data=subset(subdata_x,region=='Macular'),
                   aes(xmin=2.6,xmax=3.4,ymin=ymin,ymax=ymax),
                   colour="grey20",
                   alpha=0.2)
  p <- p+geom_segment(data=subset(subdata_x,region=='Macular'),
                      aes(x=2.6,xend=3.4,y=med,yend=med))
  p <- p+geom_point(data=subset(subdata,group!='control'),
                    aes(x=pos_i,
                        y=thickness,
                        shape=test,
                        color=subject,
                        size=4))
  p <- p+ggtitle(layer)
  p <- p+scale_x_continuous(name='Region',
                            breaks=c(1,2,3,4),
                            labels=c('Fovea','Inner','Outer','Macular'))
  p <- p+scale_y_continuous(name=expression(paste("Thickness (",mu,"m)")))
  p <- p + default_theme
  print(p)
  filename = file.path('Output','Etdrs','Rings',paste0(layer,'_test-all.png'))
  ggsave(filename,p,width=4.5,height=5.5,units=c("in"))
}


for(layer in levels(data$layer)){
  subdata <- data[data$layer==layer & data$test==2,]
  subdata_x <- x[x$layer==layer,]
  p <- ggplot(subdata)
  #p <- p+geom_boxplot(data=subset(subdata,group=='control'),outlier.shape=NA)
#   p <- p+geom_rect(data=subset(subdata_x,region=='Fovea'),
#                    aes(xmin=0.6,xmax=1.4,ymin=ymin,ymax=ymax),
#                    colour="grey20",
#                    alpha=0.2)
#   p <- p+geom_segment(data=subset(subdata_x,region=='Fovea'),
#                       aes(x=0.6,xend=1.4,y=med,yend=med))
  p <- p+geom_rect(data=subset(subdata_x,region=='Inner'),
                   aes(xmin=0.6,xmax=1.4,ymin=ymin,ymax=ymax),
                   colour="grey20",
                   alpha=0.2)
  p <- p+geom_segment(data=subset(subdata_x,region=='Inner'),
                      aes(x=0.6,xend=1.4,y=med,yend=med))
  p <- p+geom_rect(data=subset(subdata_x,region=='Outer'),
                   aes(xmin=1.6,xmax=2.4,ymin=ymin,ymax=ymax),
                   colour="grey20",
                   alpha=0.2)
  p <- p+geom_segment(data=subset(subdata_x,region=='Outer'),
                      aes(x=1.6,xend=2.4,y=med,yend=med))
  p <- p+geom_rect(data=subset(subdata_x,region=='Macular'),
                   aes(xmin=2.6,xmax=3.4,ymin=ymin,ymax=ymax),
                   colour="grey20",
                   alpha=0.2)
  p <- p+geom_segment(data=subset(subdata_x,region=='Macular'),
                      aes(x=2.6,xend=3.4,y=med,yend=med))
  p <- p+geom_point(data=subset(subdata,group!='control'),
                    aes(x=as.numeric(region),
                        y=thickness, 
                        color=subject,
                        size=4))
  p <- p+ggtitle(layer)
  p <- p+scale_x_continuous(name='Region',
                            breaks=c(1,2,3,4),
                            labels=c('Fovea','Inner','Outer','Macular'))
  p <- p+scale_y_continuous(name=expression(paste("Thickness (",mu,"m)")))
  p <- p + default_theme
  print(p)
  filename = file.path('Output','Etdrs','Rings',paste0(layer,'_test-2.png'))
  ggsave(filename,p)
}

  
