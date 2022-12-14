library(tidyverse)
library(Momocs)
library(tibble)
library(ape)

f <- list.files("class_out_data",pattern=".txt|.csv",full.names = TRUE)
Species_Scale_Bar <- read_csv("Scale Bar  - Group 1.csv")

Species_Scale_Bar <- Species_Scale_Bar %>% 
  filter(Scale_Bar=="1")


out <- read_delim(f[1],delim="\t") %>% 
  as.matrix()
#make a large df with vroom
out.df <- vroom::vroom(f, id = "filename")

out.df <- out.df %>% 
  mutate(wing=gsub("XY_.+_(hindwing|forewing)\\..+","\\1",basename(filename))) %>% 
  na.omit()

#make list
outs.l <- sapply(f,function(x) out.df %>% filter(filename==x) %>% select(X,Y) %>% as.matrix)

#extract wing info
wings <- gsub("XY_.+_(hindwing|forewing)\\..+","\\1",basename(names(outs.l)))
outs <-  outs.l %>% 
  Out(fac=list(wing=wings)) %>% 
  coo_flipx()
forewings <- outs %>% 
  filter(wing=="forewing")

hindwings <- outs %>% 
  filter(wing=="hindwing")


fore.min <- forewings %>% 
  coo_nb() %>% 
  min()

forewing.gp <- forewings %>%
  coo_interpolate(fore.min) %>% 
  fgProcrustes() 

hind.min <- hindwings %>% 
  coo_nb() %>% 
  min()

hindwing.gp <- hindwings %>% 
  coo_interpolate(hind.min) %>% 
  coo_slide(id=1) %>% 
  coo_align()  %>%
  fgProcrustes()

hind.area <- hindwing.gp%>%
  coo_trans(10,10) %>% #just moves every outline up and right 10 units so areas are all positive
  coo_area()


hind.w <- hindwing.gp%>%
  coo_length()

hindwing.AR <- hind.w^2/hind.area
hindwing.AR.DF <- data.frame(xy.file=basename(names(hindwing.AR))) %>% 
  mutate(identifier=gsub("XY_|_hindwing|_forewing|.txt","",xy.file)) %>% 
  mutate(hindwing.AR)


fore.area <- forewing.gp%>%
  coo_trans(10,10) %>% #just moves every outline up and right 10 units so areas are all positive
  coo_area()
fore.w <- forewing.gp %>%
  coo_length()
forewing.AR <- fore.w^2/fore.area
forewing.AR.DF <- data.frame(xy.file=basename(names(forewing.AR))) %>% 
  mutate(identifier=gsub("XY_|_hindwing|_forewing|.txt","",xy.file)) %>% 
  mutate(forewing.AR)

#forewing.AR.DF <- forewing.AR.DF%>% 
# left_join(Species_Scale_Bar) %>%
# na.omit()



#PCA Analysis

lep.tree <- ape::read.tree("lep_tree2.tre")
lep.tree <- ladderize(lep.tree)

lep.tree$tip.label <- gsub("_"," ",lep.tree$tip.label)


lep.sp <- read_csv("lep_image_data.csv")

out.data <- tibble(xy.file=basename(names(outs))) %>% 
  mutate(identifier=gsub("XY_|_hindwing|_forewing|.txt","",xy.file)) %>% 
  left_join(lep.sp)

hindwing.AR.DF <- hindwing.AR.DF %>% 
  left_join(lep.sp)

forewing.AR.DF <- forewing.AR.DF %>% 
  left_join(lep.sp)

forewing.pca <- forewings %>%
  coo_interpolate(fore.min) %>%
  coo_align()  %>%
  coo_slide(id=1) %>% 
  fgProcrustes() %>% 
  efourier(norm=FALSE) %>% 
  PCA()

hindwing.pca <-hindwings %>% 
  coo_interpolate(hind.min) %>% 
  coo_align()  %>%
  coo_slide(id=1) %>% 
  fgProcrustes() %>% 
  efourier(norm=FALSE) %>% 
  PCA()

out.data <- tibble(xy.file=basename(names(outs))) %>% 
  mutate(identifier=gsub("XY_|_hindwing|_forewing|.txt","",xy.file)) %>% 
  left_join(lep.sp)

hindwing.pca2 <-  tibble(xy.file=basename(rownames(hindwing.pca$x)),PC1=hindwing.pca$x[,1],PC2=hindwing.pca$x[,2]) %>% 
  left_join(out.data)

forewing.pca2 <-  tibble(xy.file=basename(rownames(forewing.pca$x)),PC1=forewing.pca$x[,1],PC2=forewing.pca$x[,2])%>% 
  left_join(out.data)


hindwing.AR.PCA <- hindwing.AR.DF %>% 
  left_join(hindwing.pca2)%>% 
  na.omit()


forewing.AR.PCA <- forewing.AR.DF %>% 
  left_join(forewing.pca2)%>% 
  na.omit()

hindwing.AR.PCA %>% 
  ggplot(aes(x=hindwing.AR,y=PC1))+geom_point()+geom_smooth(method="lm")

forewing.AR.PCA %>% 
  ggplot(aes(x=forewing.AR,y=PC1))+geom_point()+geom_smooth(method="lm")

hindwing.AR.PCA %>% 
  ggplot(aes(x=hindwing.AR,y=PC2))+geom_point()+geom_smooth(method="lm")

forewing.AR.PCA %>% 
  ggplot(aes(x=forewing.AR,y=PC2))+geom_point()+geom_smooth(method="lm")




#evolutionary data
drops <- lep.tree$tip.label[!lep.tree$tip.label%in%unique(out.data$species)]

lep.tree2 <- drop.tip(lep.tree,drops)

plot(lep.tree2,cex=0.5)

hind.AR <- hindwing.AR.DF %>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(hindwing.AR=mean(hindwing.AR)) %>% 
  pull

names(hind.AR) <-  hindwing.AR.DF%>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(hindwing.AR=mean(hindwing.AR)) %>% 
  pull(species)


fore.AR <- forewing.AR.DF %>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(forewing.AR=mean(forewing.AR)) %>% 
  pull

names(fore.AR) <-  forewing.AR.DF%>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(forewing.AR=mean(forewing.AR)) %>% 
  pull(species)

library(phytools)

foreAR.BM<-brownie.lite(lep.tree2,fore.AR*10)
hindAR.BM<-brownie.lite(lep.tree2,hind.AR*10)

foreAR.BM$sig2.single

library(RRphylo)

hindAR.RR <- RRphylo(tree=lep.tree2,y=hind.AR)
hindAR.RR$rates

hindAR.SS<- search.shift(RR=hindAR.RR,status.type="clade")

hindAR.SS$single.clades

plot(lep.tree2)
nodelabels(node = as.numeric(rownames(hindAR.SS$single.clades)),text = rownames(hindAR.SS$single.clades))

hindAR.plot <- plotShift(RR=hindAR.RR,SS=hindAR.SS)

foreAR.plot <- plotShift(RR=hindAR.RR,SS=hindAR.SS)
hindAR.plot$plotClades()

library(ggtree)
library(wesanderson)

plot_AR <- function(tre=NULL,SS=NULL,tax=NULL){
  
  
  nodes <- as.numeric(rownames(SS$single.clades))
  
  pal <- wes_palette("Zissou1",n=length(nodes))
  sp <- list()
  for(i in nodes){
    sp.i <- extract.clade(tre,i)$tip.label
    
    #print(head(tax))
    sub.names <- lapply(tax,function(x) x[x%in%sp.i]) 
    
    in.clades <- lapply(sub.names,function(x) length(x)>0) 
    all.of.clade <- lapply(sub.names,function(x) all(sapply(sp.i,function(z) z%in%x))) 
    
    high.clade <- names(sub.names)[last(which(all.of.clade==T))]
    all.clades <- names(sub.names)[which(in.clades==T)]
    crown <- ""
    if(high.clade!=last(names(sub.names))) crown <- "crown-"
    
    sub.clades <- NULL
    if(length(grepl("oidea",all.clades))>0) sub.clades <- all.clades[grepl("oidea",all.clades)]
    
    high.clade2 <- paste0(crown,high.clade,": ",paste0(sub.clades,collapse = "+"))
    sp[[paste0(i)]] <- tibble(n=i,species=sp.i,clade=high.clade2)
    
  }
  
  
  d <- do.call(rbind,sp)%>% 
    rename(label=species) 
  
  d2<- d %>% rename(clade_name=clade) 
  
  p <- ggtree(tre)+ scale_y_reverse()
  
  p$data <- p$data %>% left_join(d) %>% left_join(tibble(node=nodes,SS$single.clades) %>% mutate(shift=ifelse(rate.difference>0,"+","-")))
  
  p <-  p+geom_tiplab(aes(col=clade),geom="text",size=1.2)+
    geom_cladelab(data=d2,mapping=aes(node=n,col=clade_name,label=clade_name),offset=1,size=1.5)+
    geom_hilight(data=d2,mapping = aes(node = n,fill=clade_name),alpha = 0.01)+
    scale_fill_manual(values = pal)+
    scale_color_manual(values = pal)+
    theme(legend.position = "none")+geom_nodepoint(mapping=aes(subset = shift =="-"), size=5, shape=25,fill='blue',color='blue',alpha=0.7)+
    geom_nodepoint(mapping=aes(subset = shift =="+"), size=5, shape=24, fill='red',color='red',alpha=0.7)
  p <- p+xlim(NA,6)
  res <- tibble(n=nodes,SS$single.clades) %>% left_join(d %>% select(n,clade) %>% unique)
  
  return(list(plot=p,res=res))
  
}

tax.names <- readRDS("Lep_classification.RDS")

hindAR.res <- plot_SS(lep.tree2,hindAR.SS,tax = tax.names)

hindAR.res$plot

hindAR.res$res

hindAR.pic <- pic(hind.AR,phy = lep.tree2)
foreAR.pic <- pic(fore.AR,phy = lep.tree2)

AR.pic <- tibble(
  hind=hindAR.pic,
  fore=foreAR.pic
)

AR.pic %>% 
  ggplot(aes(x=fore,y=hind))+geom_point()+geom_smooth(method="lm")

summary(lm(hind~fore,AR.pic))
