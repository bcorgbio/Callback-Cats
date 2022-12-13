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

hind.w <- hindwings %>%
  coo_length()

hindwing.AR <- hind.w^2/hind.area
hindwing.AR.DF <- data.frame(names(hindwing.AR),unname(hindwing.AR))
colnames(hindwing.AR.DF) <- c("Filename","AR")


fore.area <- forewing.gp%>%
  coo_trans(10,10) %>% #just moves every outline up and right 10 units so areas are all positive
  coo_area()
fore.w <- forewings %>%
  coo_length()
forewing.AR <- fore.w^2/fore.area
forewing.AR.DF <- data.frame(names(forewing.AR),unname(forewing.AR))
colnames(forewing.AR.DF) <- c("Filename","AR")


