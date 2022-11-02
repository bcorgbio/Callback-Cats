#Module 6 - adam script

library(tidyverse)
library(Momocs)

f <- list.files("class_out_data",pattern=".txt|.csv",full.names = TRUE)

out <- read_delim(f[1],delim="\t") %>% 
  as.matrix()

#make a large df with vroom
out.df <- vroom::vroom(f, id = "filename")

#add wing info
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

#visualize each group
forewings <- outs %>% 
  filter(wing=="forewing")

hindwings <- outs %>% 
  filter(wing=="hindwing")

forewings %>% 
  stack()
hindwings %>% 
  stack()

#Procrustes: fit samples together
fore.min <- forewings %>% 
  coo_nb() %>% 
  min()

forewings %>%
  coo_interpolate(fore.min) %>% 
  fgProcrustes() %>% 
  stack()

hind.min <- hindwings %>% 
  coo_nb() %>% 
  min()

hindwings %>% 
  coo_interpolate(hind.min) %>% 
  coo_slide(id=1) %>% 
  coo_align()  %>%
  fgProcrustes() %>%
  stack()

#perform EFA
forewings %>%
  coo_interpolate(fore.min) %>% 
  coo_align()  %>%
  fgProcrustes() %>% 
  efourier(norm=FALSE) 

hindwings %>% 
  coo_interpolate(hind.min) %>% 
  coo_align()  %>%
  fgProcrustes() %>% 
  efourier(norm=FALSE) 

#perform PCA
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

hindwings %>% 
  coo_interpolate(hind.min) %>% 
  coo_align()  %>%
  coo_slide(id=1) %>% 
  fgProcrustes() %>% 
  stack
#visualize PCA
forewing.pca %>% 
  plot_PCA(title = "forewings")
hindwing.pca %>% 
  plot_PCA(title = "hindwings")

#COMPARATIVE ANALYSIS
library(ape)

lep.tree <- ape::read.tree("lep_tree2.tre")

lep.tree <- ladderize(lep.tree)
plot(lep.tree,cex=0.1)

lep.tree$tip.label <- gsub("_"," ",lep.tree$tip.label)
basename(names(outs))[1:5]

lep.sp <- read_csv("lep_image_data.csv")
head(lep.sp)
head(lep.sp$identifier)

out.data <- tibble(xy.file=basename(names(outs))) %>% 
  mutate(identifier=gsub("XY_|_hindwing|_forewing|.txt","",xy.file)) %>% 
  left_join(lep.sp)
head(out.data)

head(hindwing.pca$x,1)

#First 2 PCs for each wing stored with species info
hindwing.pca2 <-  tibble(xy.file=basename(rownames(hindwing.pca$x)),PC1=hindwing.pca$x[,1],PC2=hindwing.pca$x[,2]) %>% 
  left_join(out.data)
forewing.pca2 <-  tibble(xy.file=basename(rownames(forewing.pca$x)),PC1=forewing.pca$x[,1],PC2=forewing.pca$x[,2])%>% 
  left_join(out.data)

#EVOLUTIONARY RATES
drops <- lep.tree$tip.label[!lep.tree$tip.label%in%unique(out.data$species)]

lep.tree2 <- drop.tip(lep.tree,drops)

plot(lep.tree2,cex=0.5)

#PC1s
hind.pc1 <- hindwing.pca2 %>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC1=mean(PC1)) %>% 
  pull

names(hind.pc1) <-  hindwing.pca2%>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC1=mean(PC1)) %>% 
  pull(species)

fore.pc1 <- forewing.pca2 %>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC1=mean(PC1)) %>% 
  pull(PC1)

names(fore.pc1) <-  forewing.pca2 %>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC1=mean(PC1)) %>% 
  pull(species)

#PC2s
hind.pc2 <- hindwing.pca2 %>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC2=mean(PC2)) %>% 
  pull(PC2)

names(hind.pc2) <-  hindwing.pca2%>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>%
  summarize(PC2=mean(PC2)) %>% 
  pull(species)

fore.pc2 <- forewing.pca2 %>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC2=mean(PC2)) %>% 
  pull(PC2)

names(fore.pc2) <-  forewing.pca2 %>% 
  filter(species%in% lep.tree2$tip.label) %>% 
  group_by(species) %>% 
  summarize(PC2=mean(PC2)) %>% 
  pull(species)

library(phytools)
forePC1.BM<-brownie.lite(lep.tree2,fore.pc1*10)
hindPC1.BM<-brownie.lite(lep.tree2,hind.pc1*10)

forePC2.BM<-brownie.lite(lep.tree2,fore.pc2*10)
hindPC2.BM<-brownie.lite(lep.tree2,hind.pc2*10)

forePC1.BM$sig2.single

#Shifts in evolutionary rate
library(RRphylo)
hindPC1.RR <- RRphylo(tree=lep.tree2,y=hind.pc1)
hindPC1.RR$rates

hindPC1.SS<- search.shift(RR=hindPC1.RR,status.type="clade")

hindPC1.SS$single.clades

plot(lep.tree2)
nodelabels(node = as.numeric(rownames(hindPC1.SS$single.clades)),text = rownames(hindPC1.SS$single.clades))

hindPC1.plot <- plotShift(RR=hindPC1.RR,SS=hindPC1.SS)

forePC1.plot <- plotShift(RR=hindPC1.RR,SS=hindPC1.SS)
hindPC1.plot$plotClades()

library(ggtree)
library(wesanderson)

plot_SS <- function(tre=NULL,SS=NULL,tax=NULL){
  
  
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
hindPC1.res <- plot_SS(lep.tree2,hindPC1.SS,tax = tax.names)

hindPC1.res$plot



