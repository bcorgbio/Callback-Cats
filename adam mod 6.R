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
