plot.negative.slope_vl.measurements<-function(vl,missing,re_vl){

  
re_vl$cohort<-as.factor(sapply(sapply(rownames(re_vl),function(x) strsplit(x,split="[[:digit:]]")),"[",1))
re_vl$negative_slope<-as.factor(is.na(re_vl$spvl));levels(re_vl$negative_slope)=paste(c('positive','negative'),'slope')
df<-melt(apply(table(re_vl$cohort,re_vl$negative_slope),1,prop.table))


p1<-ggplot(df,aes(value,fill=Var1,x=Var2))+
  geom_bar(position='stack',stat='identity')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1),legend.position="none")+
  xlab('Cohort')+
  ylab('Fraction')+
  ggtitle('Negative slope spvl_model per cohort')
  
X<-as.data.frame(table(vl$patient_id));colnames(X)[1]<-'patient_id';X$negative_slope<-as.factor(X$patient_id%in%missing)

p2<-ggplot(X,aes(Freq,fill=negative_slope))+
  geom_histogram(position = "stack", binwidth=1)+
  theme_bw()+
  xlab("number of vl measurements")+
  ylab("patients")+
  scale_x_discrete(breaks=c(2,seq(5,range(X$Freq)[2],5)))+
  scale_fill_discrete(labels = paste(c("positive","negative"), "slope"))+
  theme(legend.position=c(0.8,0.7),legend.title=element_blank())+
  ggtitle("Negative slope spvl_model")

multiplot(p1,p2)
}

