library(tidyverse)
library(MuMIn)
library(ggplot2)
library(dplyr)


k <- list.files("./Project 8 Data", full.names = T,pattern = ".csv")
print(k)
k.l <- list()

for(i in k){
  met.dat <- unlist(strsplit(i,"_"))
  sub <- met.dat[2]
  ang <- as.numeric(met.dat[3])
  exp <- gsub("\\..+","",met.dat[4])
  k.l[[i]] <- read_delim(i,delim = " ", col_names = c("Reading","Force","Unit"), id="Experiment",progress = FALSE) %>%select(Force)%>%
    mutate(sub=sub,ang=ang,exp=exp)
}

data <- do.call(rbind, k.l)
data <- data%>%
  group_by(sub,exp,ang)%>%
  summarise(max.force=max(abs(Force), na.rm=TRUE), n=n())

data.joined <- data%>%
  group_by(sub,exp)%>%
  summarize(max.f2 = max(max.force))%>%
  left_join(data)%>%
  mutate(norm.force=max.force/max.f2)

data.con <- filter(data.joined,exp=="control")
data.fat <- filter(data.joined,exp=="fatigue")

ang.con <- data.con$ang
normF.con <- data.con$norm.force

ang.fat <- data.fat$ang
normF.fat <- data.fat$norm.force

#control first
poly.m2.con <- lm(normF.con~poly(ang.con,2)) #second order
poly.m3.con <- lm(normF.con~poly(ang.con,3)) #third order
poly.m4.con <- lm(normF.con~poly(ang.con,4)) #fourth order

AICc(poly.m2.con,poly.m3.con,poly.m4.con) #the fourth order model fits best

x.pred<- seq(45,157.5,length.out = 1000) #define 1000 angles from our range

normF.pred.con <- predict(poly.m4.con,newdata = data.frame(ang.con=x.pred)) #predict the force using 1000 angles

qplot(ang.con,normF.con)+geom_point(aes(x=x.pred,y=normF.pred.con),col="red")+geom_point(aes(x=x.pred[which.max(normF.pred.con)],y=normF.pred.con[which.max(normF.pred.con)]),size=5,col="blue")

x.pred[which.max(normF.pred.con)] #theta_max is 157.5

#fatigue
poly.m2.fat <- lm(normF.fat~poly(ang.fat,2)) #second order
poly.m3.fat <- lm(normF.fat~poly(ang.fat,3)) #third order
poly.m4.fat <- lm(normF.fat~poly(ang.fat,4)) #fourth order

AICc(poly.m2.fat,poly.m3.fat,poly.m4.fat) #the fourth order model fits best


normF.pred.fat <- predict(poly.m4.fat,newdata = data.frame(ang.fat=x.pred)) #predict the force using 1000 angles

qplot(ang.fat,normF.fat)+geom_point(aes(x=x.pred,y=normF.pred.fat),col="red")+geom_point(aes(x=x.pred[which.max(normF.pred.fat)],y=normF.pred.fat[which.max(normF.pred.fat)]),size=5,col="blue")


x.pred[which.max(normF.pred.fat)]-x.pred[which.max(normF.pred.con)] #0 Degree Shift

data.joined%>%
  ggplot(aes(ang,norm.force,col=exp))+geom_point()


AICs <- data.joined%>%
  group_by(sub,exp)%>%
  summarize(
    m2=AICc(lm(norm.force~poly(ang,2))), #second order
    m3=AICc(lm(norm.force~poly(ang,3))), #third order
    m4=AICc(lm(norm.force~poly(ang,4))) #fourth order
  )%>%
  pivot_longer(m2:m4,names_to="model",values_to="AICc")%>%
  print()

fits <- data.joined%>%
  group_by(sub,exp)%>%
  summarize(
    m2=predict(lm(norm.force~poly(ang,2)),newdata=data.frame(ang=x.pred)), #second order
    m3=predict(lm(norm.force~poly(ang,3)),newdata=data.frame(ang=x.pred)), #third order
    m4=predict(lm(norm.force~poly(ang,4)),newdata=data.frame(ang=x.pred)) #fourth order
  )%>%
  pivot_longer(m2:m4,names_to="model")%>%
  group_by(sub,exp,model)%>%
  summarize(theta_max=x.pred[which.max(value)])%>%
  print()

best.models <- fits%>%
  left_join(AICs)%>%
  group_by(sub,exp)%>%
  mutate(best=AICc==min(AICc))%>%
  filter(best==TRUE)%>%
  dplyr::select(-best)%>%
  print()

anova(lm(theta_max~exp,best.models))

best.models%>%
  pivot_wider(id_cols=sub,names_from = exp,values_from=theta_max)%>%
  mutate(shift=fatigue-control)%>%
  ungroup()%>%
  summarize(mean.shift=abs(mean(na.rm=TRUE,shift)),se.shift=sd(na.rm=TRUE,shift)/sqrt(length(shift)))

