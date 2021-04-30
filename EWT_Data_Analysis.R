#packages needed
library("vegan")
library("tidyverse")
library("iNEXT")
library("magrittr")
library("reshape")
library("wiqid")
library("betapart")
library("unmarked")

#import colne data
colne <- read.csv("ewt1_sites.csv")
colnames(colne)

#import blackwater data
blackwater <- read.csv("ewt2_sites.csv")
colnames(blackwater)

#importing site data
site <- read.csv("site_data.csv")
colnames(site)[1] <- c("Site")

#combining datasets
combined <- cbind(blackwater, colne[,4:18], deparse.level = 1)
colnames(combined)

#species accumulation curve - carnivore
carnivora <- subset(combined, Taxon=="Carnivora")
carnivora <- carnivora[,-(1:2)]

#Remove rownames and change to species names
row.names(carnivora) =NULL
carnivora <- column_to_rownames(carnivora, var = "species")

#combine dataset with itself in order to plot curves looking at both rivers and overall
carnivora <- cbind(carnivora, carnivora)
carnivora <- t(carnivora) #transpose data

#transform reads into Presense/absense data
carnivora <- carnivora %>% decostand(., "pa")

#combine with site data
carnivora <- cbind(site$River, carnivora, deparse.level = 1)
colnames(carnivora)[1] <- c("River")
carnivora <- carnivora[-(1:5),] # remove beaver experiment
rownames(carnivora)
carnivora <- as.data.frame(carnivora[-(26:30),])

#split data according to system
richlist.carnivora <- carnivora %>%
  group_split(River) %>%
  setNames(c('Blackwater','Colne', 'Overall')) #set names of each matrix in the list (check order!!!!)

richlist.carnivora <- lapply(richlist.carnivora, t) #transpose -  this will also convert them to matrices
richlist.carnivora <- lapply(richlist.carnivora, function(x) x[-1,]) # remove first row with names of loc
richlist.carnivora[] <- lapply(richlist.carnivora, type.convert, as.is = TRUE)
str(richlist.carnivora)

#creates species accumulation data looking at species richness (q=0) with all curves exprapolated to 50 sample sites
chao.div.carnivora <- iNEXT(richlist.carnivora, q=0, datatype = 'incidence_raw', endpoint = 50)
chao.div.carnivora$AsyEst 

#rarefaction curve - carnivora
sample.carnivora <- ggiNEXT(chao.div.carnivora, type=1)+
  theme_bw()+
  ylab("Species Richness")+
  xlab("Number of sites")+
  scale_color_manual(values = c("orange", "deepskyblue3", "black"))+
  scale_fill_manual(values = c("orange", "skyblue", "grey"))+
  ylim(c(1,15))
sample.carnivora
ggsave("acc_curve_carnivora.png")

#species accumulation curve - Artiodactyla
artiodactyla <- subset(combined, Taxon=="Artiodactyla")
artiodactyla <- artiodactyla[,-(1:2)]
row.names(artiodactyla) =NULL
artiodactyla <- column_to_rownames(artiodactyla, var = "species")
artiodactyla <- cbind(artiodactyla, artiodactyla)
artiodactyla <- t(artiodactyla)
artiodactyla <- artiodactyla %>% decostand(., "pa")
artiodactyla <- cbind(site$River, artiodactyla, deparse.level = 1)
colnames(artiodactyla)[1] <- c("River")
artiodactyla <- artiodactyla[-(1:5),]
artiodactyla <- as.data.frame(artiodactyla[-(26:30),])

richlist.artiodactyla <- artiodactyla %>%
  group_split(River) %>%
  setNames(c('Blackwater','Colne', 'Overall')) 

richlist.artiodactyla <- lapply(richlist.artiodactyla, t) 
richlist.artiodactyla <- lapply(richlist.artiodactyla, function(x) x[-1,]) 
richlist.artiodactyla[] <- lapply(richlist.artiodactyla, type.convert, as.is = TRUE)
str(richlist.artiodactyla)
chao.div.artiodactyla <- iNEXT(richlist.artiodactyla, q=0, datatype = 'incidence_raw', endpoint = 50)
chao.div.artiodactyla$AsyEst
sample.artiodactyla <- ggiNEXT(chao.div.artiodactyla, type=1)+
  theme_bw()+
  ylab("Species Richness")+
  xlab("Number of sites")+
  scale_color_manual(values = c("orange", "deepskyblue3", "black"))+
  scale_fill_manual(values = c("orange","skyblue", "grey"))+
  ylim(c(1,15))
sample.artiodactyla
ggsave("acc_curve_artiodactyla.png")

#species accumulation curve - Eulipotyphla
Eulipotyphla <- subset(combined, Taxon=="Eulipotyphla")
Eulipotyphla <- Eulipotyphla[,-(1:2)]
row.names(Eulipotyphla) =NULL
Eulipotyphla <- column_to_rownames(Eulipotyphla, var = "species")
Eulipotyphla <- cbind(Eulipotyphla, Eulipotyphla)
Eulipotyphla <- t(Eulipotyphla)
Eulipotyphla <- Eulipotyphla %>% decostand(., "pa")
Eulipotyphla <- cbind(site$River, Eulipotyphla, deparse.level = 1)
colnames(Eulipotyphla)[1] <- c("River")
Eulipotyphla <- Eulipotyphla[-(1:5),]
Eulipotyphla <- as.data.frame(Eulipotyphla[-(26:30),])

richlist.Eulipotyphla <- Eulipotyphla %>%
  group_split(River) %>%
  setNames(c('Blackwater','Colne', 'Overall')) 

richlist.Eulipotyphla <- lapply(richlist.Eulipotyphla, t) 
richlist.Eulipotyphla <- lapply(richlist.Eulipotyphla, function(x) x[-1,]) 
richlist.Eulipotyphla[] <- lapply(richlist.Eulipotyphla, type.convert, as.is = TRUE)
str(richlist.Eulipotyphla)
chao.div.Eulipotyphla <- iNEXT(richlist.Eulipotyphla, q=0, datatype = 'incidence_raw', endpoint = 50)
chao.div.Eulipotyphla$AsyEst
sample.Eulipotyphla <- ggiNEXT(chao.div.Eulipotyphla, type=1)+
  theme_bw()+
  ylab("Species Richness")+
  xlab("Number of sites")+
  scale_color_manual(values = c("orange", "deepskyblue3", "black"))+
  scale_fill_manual(values = c("orange","skyblue","grey"))+
  ylim(c(1,15))
sample.Eulipotyphla
ggsave("acc_curve_Eulipotyphla.png")

#species accumulation curve - Lagomorpha
Lagomorpha <- subset(combined, Taxon=="Lagomorpha")
Lagomorpha <- Lagomorpha[,-(1:2)]
row.names(Lagomorpha) =NULL
Lagomorpha <- column_to_rownames(Lagomorpha, var = "species")
Lagomorpha <- cbind(Lagomorpha, Lagomorpha)
Lagomorpha <- t(Lagomorpha)
Lagomorpha <- Lagomorpha %>% decostand(., "pa")
Lagomorpha <- cbind(site$River, Lagomorpha, deparse.level = 1)
colnames(Lagomorpha)[1] <- c("River")
Lagomorpha <- Lagomorpha[-(1:5),]
Lagomorpha <- as.data.frame(Lagomorpha[-(26:30),])

richlist.Lagomorpha <- Lagomorpha %>%
  group_split(River) %>%
  setNames(c('Blackwater','Colne', 'Overall')) 

richlist.Lagomorpha <- lapply(richlist.Lagomorpha, t) 
richlist.Lagomorpha <- lapply(richlist.Lagomorpha, function(x) x[-1,]) 
richlist.Lagomorpha[] <- lapply(richlist.Lagomorpha, type.convert, as.is = TRUE)
str(richlist.Lagomorpha)
chao.div.Lagomorpha <- iNEXT(richlist.Lagomorpha, q=0, datatype = 'incidence_raw', endpoint = 50)
chao.div.Lagomorpha$AsyEst
sample.Lagomorpha <- ggiNEXT(chao.div.Lagomorpha, type=1)+
  theme_bw()+
  ylab("Species Richness")+
  xlab("Number of sites")+
  scale_color_manual(values = c("orange", "deepskyblue3", "black"))+
  scale_fill_manual(values = c("orange","skyblue","grey"))+
  ylim(c(1,15))
sample.Lagomorpha
ggsave("curve_Lagomorpha.png")

#species accumulation curve - Rodentia
Rodentia <- subset(combined, Taxon == "Rodentia")
Rodentia <- Rodentia[,-(1:2)]
row.names(Rodentia) =NULL
Rodentia <- column_to_rownames(Rodentia, var = "species")
Rodentia <- cbind(Rodentia, Rodentia)
Rodentia <- t(Rodentia)
Rodentia <- Rodentia %>% decostand(., "pa")
Rodentia <- cbind(site$River, Rodentia, deparse.level = 1)
colnames(Rodentia)[1] <- c("River")
Rodentia <- Rodentia[-(1:5),]
Rodentia <- as.data.frame(Rodentia[-(26:30),])

richlist.Rodentia <- Rodentia %>%
  group_split(River) %>%
  setNames(c('Blackwater','Colne', 'Overall'))

richlist.Rodentia <- lapply(richlist.Rodentia, t) 
richlist.Rodentia <- lapply(richlist.Rodentia, function(x) x[-1,]) 
richlist.Rodentia[] <- lapply(richlist.Rodentia, type.convert, as.is = TRUE)
str(richlist.Rodentia)
chao.div.Rodentia <- iNEXT(richlist.Rodentia, q=0, datatype = 'incidence_raw', endpoint = 50)
chao.div.Rodentia$AsyEst
sample.Rodentia <- ggiNEXT(chao.div.Rodentia, type=1)+
  theme_bw()+
  ylab("Species Richness")+
  xlab("Number of sites")+
  scale_color_manual(values = c("orange", "deepskyblue3", "black"))+
  scale_fill_manual(values = c("orange","skyblue","grey"))+
  ylim(c(1,15))
sample.Rodentia
ggsave("acc_curve_Rodentia.png")

#accumulation curve - overall diversity 
comb <- combined[,-(1:2)]
comb <- column_to_rownames(comb, var = "species")
comb <- cbind(comb, comb)
comb <- t(comb)
comb <- comb %>% decostand(., "pa")
comb <- cbind(site$River, comb, deparse.level = 1)
colnames(comb)[1] <- c("River")
comb <- comb[-(1:5),]
comb <- as.data.frame(comb[-(26:30),])

richlist <- comb %>%
  group_split(River) %>%
  setNames(c('Blackwater','Colne', 'Overall')) 

richlist <- lapply(richlist, t) 
richlist <- lapply(richlist, function(x) x[-1,]) 
richlist[] <- lapply(richlist, type.convert, as.is = TRUE)
str(richlist)
chao.div <- iNEXT(richlist, q=0, datatype = 'incidence_raw', endpoint = 200, knots = 200)
chao.div$AsyEst
sample <- ggiNEXT(chao.div, type=1)+
  theme_bw()+
  ylab("Species Richness")+
  xlab("Number of sites")+
  scale_color_manual(values = c("orange", "skyblue", "black"))+
  scale_fill_manual(values = c("orange","skyblue","grey"))+
  ylim(c(1,40))
sample
ggsave("acc_curve_combined.png")

#Tranforming Blackwater data for ChaoII
blackwater_sampling <- blackwater[,-(4:8)]
blackwater_sampling <- blackwater_sampling[,-(1:2)]
blackwater_sampling[,2:11] %<>% decostand(., "pa") 
blackwater_sampling.t <- t(blackwater_sampling[,-1]); colnames(blackwater_sampling.t) <- blackwater_sampling[,1]

#melt data according to species
blackwater_sampling.t.melted <- melt(blackwater_sampling.t); colnames(blackwater_sampling.t.melted) <- c("Site", "species", "presence")

#Tranforming Colne data for ChaoII
colne_sampling <- colne[,-(1:2)]
colne_sampling[,2:16] %<>% decostand(., "pa")
colne_sampling.t <- t(colne_sampling[,-1]); colnames(colne_sampling.t) <- colne_sampling[,1]

#melt data according to species
colne_sampling.t.melted <- melt(colne_sampling.t); colnames(colne_sampling.t.melted) <- c("Site", "species", "presence")

#combine datasets 
combined_data <- merge(colne_sampling.t.melted, blackwater_sampling.t.melted , by = c("Site","species"), suffixes = c(".colne",".blackwater"), all = T)

#turns NAs into 0s
combined_data[is.na(combined_data)] <- 0

#merge with site data and add sum column
combined_data <- merge(unique(site[1:25,]), combined_data, by = "Site")
combined_data$presence.blackwater[combined_data$presence.blackwater > 0] <- 2
combined_data$sum <- combined_data$presence.colne + combined_data$presence.blackwater #to have a code for find in both sampling
data.chao.II <- combined_data

#creating chao matrix
result.chao.II <- data.frame(type = c("Total","Colne", "Blackwater"),
                             method = c("eDNA", "eDNA", "eDNA"),
                             Species_richness = NA, 
                             chaoII = NA,
                             chaoII.low = NA,
                             chaoII.up = NA)

#fill in chaoII matrix
result.chao.II[1,3] <- data.chao.II %>% 
  droplevels() %>% 
  cast(data=., Site ~ species, value = "sum", add.missing = T, fill = 0, fun.aggregate = "sum") %>% 
  decostand(x = ., method = "pa") %>% 
  specpool() %>% 
  '[' (1)

result.chao.II[2,3] <- subset(data.chao.II, River == "Colne") %>% 
  droplevels() %>% 
  cast(data=., Site ~ species, value = "sum", add.missing = T, fill = 0, fun.aggregate = "sum") %>% 
  decostand(x = ., method = "pa") %>% 
  specpool() %>% 
  '[' (1)

result.chao.II[3,3] <- subset(data.chao.II, River == "Blackwater") %>% 
  droplevels() %>% 
  cast(data=., Site ~ species, value = "sum", add.missing = T, fill = 0, fun.aggregate = "sum") %>% 
  decostand(x = ., method = "pa") %>% 
  specpool() %>% 
  '[' (1)

result.chao.II[1,4:6] <- data.chao.II %>% 
  droplevels() %>% 
  cast(data=., Site ~ species, value = "sum", add.missing = T, fill = 0, fun.aggregate = "sum") %>% 
  decostand(x = ., method = "pa") %>% 
  t %>%
  richChao2(correct = T) %>% 
  '[' (1:3)

result.chao.II[2,4:6] <- subset(data.chao.II, River == "Colne") %>% 
  droplevels() %>% 
  cast(data=., Site ~ species, value = "sum", add.missing = T, fill = 0, fun.aggregate = "sum") %>% 
  decostand(x = ., method = "pa") %>% 
  t %>%
  richChao2(correct = T) %>% 
  '[' (1:3)

result.chao.II[3,4:6] <- subset(data.chao.II, River == "Blackwater") %>% 
  droplevels() %>% 
  cast(data=., Site ~ species, value = "sum", add.missing = T, fill = 0, fun.aggregate = "sum") %>% 
  decostand(x = ., method = "pa") %>% 
  t %>%
  richChao2(correct = T) %>% 
  '[' (1:3)

#Chao matrix
write.csv(result.chao.II, "chao_matrix_rode.csv")

#NMDs with elipses
species <- column_to_rownames(combined, var = "species") # converts the "Species" column to rownames
head(species)
species <- species[,-(1:2)]
species %<>% decostand(., "pa")
species.t <- as.data.frame(t(species)) # transpose data

#Jaccard distances between sites
species_dist <- vegdist(species.t[6:30,], method = "jaccard") 
species_dist

#remove beaver enclousure
ewt.nmds <- species.t[-(1:5),]

#NMDs function 
set.seed(1000)
ewt.mds <- metaMDS(ewt.nmds, distance = "jaccard", autotransform = F)
ewt.mds

#stressplot to show how well the NMDs data fits the actual data
stressplot(ewt.mds)
plot(ewt.mds) #basic plot
plot(ewt.mds, type = "t") #basic plot with site names
site_data <- site[(6:30),]#edit site data

#transform data
site.scrs <- as.data.frame(scores(ewt.mds, display = "sites")) 
site.scrs <- cbind(site.scrs, River = site_data$River)  
site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) 
head(site.scrs)

#change river data to a factor
site.scrs$River <- as.factor(site.scrs$River)
str(site.scrs)

#nmds plot
nmds.plot.ewt <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2, fill=River))+ 
  geom_point(shape=21, size=5, colour="black", alpha=0.8,aes(NMDS1, NMDS2, colour = factor(site.scrs$River)))+ 
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "River")+ 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))
nmds.plot.ewt<- nmds.plot.ewt+scale_fill_manual(values = c("orange", "deepskyblue3"))
nmds.plot.ewt

#function to create ellipses
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the management factor
df_ell.ewt.river <- data.frame() #sets up a data frame before running the function.
for(g in levels(site.scrs$River)){
  df_ell.ewt.river <- rbind(df_ell.ewt.river, cbind(as.data.frame(with(site.scrs [site.scrs$River==g,],
                                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,River=g))
}

#adding ellipses and title
nmds.plot.ewt+ 
  geom_path(data = df_ell.ewt.river, aes(x = NMDS1, y = NMDS2, group = River, colour=River))+ #this is the ellipse, separate ones by Site. 
  labs(title = "Ordination plot with ellipses")+
  scale_color_manual(values = c("orange", "deepskyblue3"))
ggsave("nmds_plot.png")

#create data frame for PARMANOVA
nmds.spec <- as.data.frame(scores(ewt.mds)) %>% 
  rownames_to_column(var = "Site") %>% 
  mutate(River = ifelse(Site %in% c("C01", "C02", "C03", "C04","C05","C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14", "C15"), "Colne", 
                        ifelse(Site %in% c("B06", "B07", "B08", "B09", "B10", "B11", "B12", "B13", "B14", "B15"), "Blackwater", NA))) 
#PERMANOVA
set.seed(1000)
adonis(species.t[6:30,]~River, data = nmds.spec, method = "jaccard", permutations = 1000)

#JITTER_PLOT data
combined_jit <- combined[,-(1:2)]

#species to rownames
jit_data <- column_to_rownames(combined_jit, var = "species")
str(jit_data)
jit_data <- t(jit_data) # transpose data
jit_data <- jit_data[-(1:5),] #remove beaver enclosure 
jit_data <- jit_data %<>% decostand(., "pa") #PA <- reads

#row sums 
species_rich <- as.data.frame(rowSums(jit_data))
colnames(species_rich)[1] <- c("Species_Richness")

#combined with site data
species_richness <- cbind(site_data, species_rich)

#jitter boxplot
ewt_box <- ggplot(species_richness, aes(x=River, y=Species_Richness, colour=River))+
  geom_boxplot(size=1, colour="black")+
  geom_jitter(size=4, alpha=0.6)+
  theme_classic()+
  ylab("Species Richness")+
  theme(legend.position = "none", text = element_text(size=15))+
  scale_color_manual(values = c("orange","deepskyblue3"))
ewt_box
ggsave("box_jitter.png")

#proportional bubble charts
prop <- as.data.frame(t(combined[,4:33]))

#calculating proportions
sum <- as.data.frame(rowSums(prop, na.rm = FALSE))
colnames(sum) <- c("sum")
prop <- cbind(prop, sum)
cols <- 1:25
prop <- cbind(prop, prop[cols]/prop$sum)
prop <- prop[,27:51]
prop <- as.data.frame(t(prop))
prop <- cbind(colne[,2],prop)
colnames(prop)[1] <- c("species")

#melt data acording to species
species_melt <- melt.data.frame(prop); colnames(species_melt) <- c("species","site", "reads")

#colne data
colne_melt <- species_melt[376:750,]
colne_melt <- subset(colne_melt, reads>0) #remove rows with 0 reads
c <- as.data.frame(round_any(colne_melt$reads, 0.25, f=ceiling)) # round proportion to nearest 25%
colnames(c) <- c("Reads")
colne_melt <- cbind(colne_melt, c)

#blackwater data
blackwater_melt <- species_melt[126:375,]
blackwater_melt <- subset(blackwater_melt, reads>0)
b <- as.data.frame(round_any(blackwater_melt$reads, 0.25, f=ceiling))
colnames(b) <- c("Reads")
blackwater_melt <- cbind(blackwater_melt, b)

#beaver data
beaver_melt <- species_melt[1:125,]
beaver_melt <- subset(beaver_melt, reads>0)
b <- as.data.frame(round_any(beaver_melt$reads, 0.25, f=ceiling))
colnames(b) <- c("Reads")
beaver_melt <- cbind(beaver_melt, b)

#colne proportional bubble chart
bubble_colne <- ggplot(colne_melt, aes(x=site, y=species, size=Reads))+
  geom_point(alpha = 0.5, shape=21, fill="deepskyblue3")+
  ylab('Species')+
  xlab("Location - River Colne")+
  theme_bw()+
  scale_size_continuous(breaks = c(0.25, 0.5, 0.75,1), 
                        name = 'Proportion of Reads', 
                        range = c(5,14))+
  theme(axis.text.y = element_text(face = "italic"))
bubble_colne
ggsave("colne_prop_bubble_chart.png")

#blackwater proportional bubble chart
bubble_blackwater <- ggplot(blackwater_melt, aes(x=site, y=species, size=Reads))+
  geom_point(alpha = 0.5, shape=21, fill="deepskyblue3")+
  ylab('Species')+
  xlab("Location - River Blackwater")+
  theme_bw()+
  scale_size_continuous(breaks = c(0.25, 0.5, 0.75,1), 
                        name = 'Proportion of Reads', 
                        range = c(5,14))+
  theme(axis.text.y = element_text(face = "italic"))
bubble_blackwater
ggsave("blackwater_prop_bubble_chart.png")

#beaver proportional bubble chart
bubble_beaver <- ggplot(beaver_melt, aes(x=site, y=species, size=Reads))+
  geom_point(alpha = 0.5, shape=21, fill="deepskyblue3")+
  ylab('Species')+
  xlab("Location - River Blackwater (Beaver Enclosure)")+
  theme_bw()+
  scale_size_continuous(breaks = c(0.25,0.5,0.75,1), 
                        name = 'Proportion of Reads', 
                        range = c(5,14))+
  theme(axis.text.y = element_text(face = "italic"))
bubble_beaver
ggsave("beaver_prop_bubble_chart.png")

#Set Working files
colne_beta <- colne[,-(1:2)]
colne_beta <- column_to_rownames(colne_beta, var = "species")
colne_beta <- t(colne_beta)


blackwater_beta <- blackwater[,-(1:2)]
blackwater_beta <- column_to_rownames(blackwater_beta, var = "species")
blackwater_beta <- t(blackwater_beta)
blackwater_beta <- blackwater_beta[-(1:5),]

combined_beta <- combined[,-(1:2)]
combined_beta <- column_to_rownames(combined_beta, var = "species")
combined_beta <- t(combined_beta)
combined_beta <- combined_beta[-(1:5),]
all_b <- blackwater_beta %<>% decostand(., "pa")
all_b <- colSums(all_b)
all_b <- as.data.frame(all_b)
all_c <- colne_beta %<>% decostand(., "pa")
all_c <- colSums(all_c)
all_c <- as.data.frame(all_c)
all <- cbind(all_b, all_c)
all <- t(all)
rownames(all)=c("blackwater","colne")

#Get incidence based dataframes

### Blackwater
dim(blackwater_beta)
colnames(blackwater_beta)
dataBlack<-(blackwater_beta>0)*1      # remove the species presents less than 1 time and create a presence/absence table
dataBlack<-dataBlack[,apply(dataBlack,2,function(dataBlack) !all(dataBlack==0))] #remove species not detected - all zeros
rich_black=colSums(dataBlack) # summarize the number of OTUs in each location
rich_black

### Colne

dim(colne_beta)
colnames(colne_beta)
dataColne<-(colne_beta>0)*1      # remove the species presents less than 1 time and create a presence/absence table
dataColne<-dataColne[,apply(dataColne,2,function(dataColne) !all(dataColne==0))] #remove species not detected - all zeros
rich_colne=colSums(dataColne) # summarize the number of OTUs in each location
rich_colne

### Combined Dataset

dim(combined_beta)
colnames(combined_beta)
datacomb<-(combined_beta>0)*1      # remove the species presents less than 1 time and create a presence/absence table
datacomb<-datacomb[,apply(datacomb,2,function(datacomb) !all(datacomb==0))] #remove species not detected - all zeros
rich_comb=colSums(datacomb) # summarize the number of OTUs in each location
rich_comb


### All (comparing datasets by richness)

dim(all)
colnames(all)
dataall<-(all>0)*1      # remove the species presents less than 1 time and create a presence/absence table
dataall<-dataall[,apply(dataall,2,function(dataall) !all(dataall==0))] #remove species not detected - all zeros
rich_all=colSums(dataall) # summarize the number of OTUs in each location
rich_all


## betapart analysis

#### Multiple site measures 

dataBlack.multi <- beta.multi(dataBlack, index.family = "sor") # Jaccard index, for Sorensen change for sor
dataColne.multi <- beta.multi(dataColne,index.family = "sor")
datacomb.multi <- beta.multi(datacomb, index.family = "sor")
dataall.multi <- beta.multi(dataall,index.family = "sor")

# sampling across equal sites 
dataBlack.samp <- beta.sample(dataBlack, sites=10, samples=100) 
dataColne.samp <- beta.sample(dataColne, sites=10, samples=100)
datacomb.samp <- beta.sample(datacomb, sites=10, samples=100)


# plotting the distributions of components 
dist.b <- dataBlack.samp$sampled.values 
dist.c <- dataColne.samp$sampled.values
plot(density(dist.b$beta.SOR), xlim=c(0,0.9), ylim=c(0, 32), xlab='Beta diversity', main='', col='darkorange2',lwd=3) 
lines(density(dist.b$beta.SNE), col='darkorange2',lty=5, lwd=)
lines(density(dist.b$beta.SIM),col='darkorange2',lty=3,lwd=2) 
lines(density(dist.c$beta.SOR),col='cadetblue3', lty=1,lwd=3) 
lines(density(dist.c$beta.SNE), col='cadetblue3', lty=5, lwd=3) 
lines(density(dist.c$beta.SIM), col='cadetblue3', lty=3, lwd=2) #
legend("topleft",
       c("Blackwater","Colne"),
       fill=c("darkorange2","cadetblue3")
)

dist.comb <- datacomb.samp$sampled.values 
plot(density(dist.comb$beta.SOR), xlim=c(0,0.9), ylim=c(0, 20), xlab='Beta diversity', main='', col='darkorange2',lwd=3) 
lines(density(dist.comb$beta.SNE), col='darkorange2',lty=5, lwd=)
lines(density(dist.comb$beta.SIM),col='darkorange2',lty=3,lwd=2) 



## pairwise 
pair.b <- beta.pair(dataBlack) 
pair.c <- beta.pair(dataColne) 
pair.comb <- beta.pair(datacomb) 
pair.all <- beta.pair(dataall)

#blackwater
write.table(as.matrix(pair.b$beta.sim),"blackwater-sim.txt")
write.table(as.matrix(pair.b$beta.sne),"blackwater-sne.txt")
write.table(as.matrix(pair.b$beta.sor),"blackwater-sor.txt")

#colne
write.table(as.matrix(pair.c$beta.sim),"colne-sim.txt")
write.table(as.matrix(pair.c$beta.sne),"colne-sne.txt")
write.table(as.matrix(pair.c$beta.sor),"colne-sor.txt")

#combined dataset
write.table(as.matrix(pair.comb$beta.sim),"comb-sim.txt")
write.table(as.matrix(pair.comb$beta.sne),"comb-sne.txt")
write.table(as.matrix(pair.comb$beta.sor),"comb-sor.txt")

#combined dataset_all
write.table(as.matrix(pair.all$beta.sim),"all-sim.txt")
write.table(as.matrix(pair.all$beta.sne),"all-sne.txt")
write.table(as.matrix(pair.all$beta.sor),"all-sor.txt")

# plotting clusters 
dist.b <- dataBlack.samp$sampled.values
plot(hclust(pair.b$beta.sim, method="average"), hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sim]), line=0.3)
dist.b <- dataBlack.samp$sampled.values
plot(hclust(pair.b$beta.sne, method="average"), hang=-1, main='', sub='', xlab='')
title(xlab=expression(beta[sne]), line=0.3) 

dist.c <- dataColne.samp$sampled.values
plot(hclust(pair.c$beta.sim, method="average"), hang=-1, main='', sub='', xlab='') 
title(xlab=expression(beta[sim]), line=0.3) 
dist.c <- dataColne.samp$sampled.values
plot(hclust(pair.c$beta.sne, method="average"), hang=-1, main='', sub='', xlab='') 
title(xlab=expression(beta[sne]), line=0.3) 

dist.comb <- datacomb.samp$sampled.values
plot(hclust(pair.comb$beta.sim, method="average"), hang=-1, main='', sub='', xlab='') 
title(xlab=expression(beta[sim]), line=0.3)
dist.comb <- datacomb.samp$sampled.values
plot(hclust(pair.comb$beta.sne, method="average"), hang=-1, main='', sub='', xlab='') 
title(xlab=expression(beta[sne]), line=0.3) 

#occupancy+detection probability
occupancy_detection <- read.csv("occupancy_detection.csv")

#selection species - DO FOR ALL SPECIES!!!!
occupancy_detection2 <- subset(occupancy_detection, Species=="Sorex minutus")

#remove first colum of site names
dets <- occupancy_detection2[,-1]
dets <- t(dets)

#Build unmarkedFrame object from detection data
dets2 <- unmarkedFrameOccu(y = dets)
summary(dets)

#model of constant detection with no covariates
m1 <- occu(~1 ~1, data = dets2)

#Back transform estimate of psi (occupancy)
m1.psi <- backTransform(m1, type = "state")

#95% CI of psi
m1CI.psi <- confint(m1.psi)

#Back transform estimate of p (detection probability)
m1.p <- backTransform(m1, type = "det")

#95% CI of p
m1CI.p <- confint(m1.p)

#River Colne only
dets_col <- dets[1:15,]
dets_col <- dets_col[,-1]
dets2_col <- unmarkedFrameOccu(y = dets_col)
summary(dets2_col)
m1_col <- occu(~1 ~1, data = dets2_col)
m1.psi_col <- backTransform(m1_col, type = "state")
m1CI.psi_col <- confint(m1.psi_col)
m1.p_col <- backTransform(m1_col, type = "det")
m1CI.p_col <- confint(m1.p_col)

#River Blackwater only
dets_bla <- dets[16:25,]
dets_bla <- dets_bla[,-1]
dets2_bla <- unmarkedFrameOccu(y = dets_bla)
summary(dets2_bla)
m1_bla <- occu(~1 ~1, data = dets2_bla)
m1.psi_bla <- backTransform(m1_bla, type = "state")
m1CI.psi_bla <- confint(m1.psi_bla)
m1.p_bla <- backTransform(m1_bla, type = "det")
m1CI.p_bla <- confint(m1.p_bla)

#create a data frame containing occupancy and detection probability with 95% confidence intervals
occ_det <- data.frame(system = c("Combined", "Colne", "Blackwater"),
                      occupancy = NA,
                      occupancy.low = NA, 
                      occupancy.up = NA,
                      detection_prob = NA,
                      detection_prob.low = NA,
                      detection_prob.up=NA)
occ_det[1,2] <- round(m1.psi@estimate,2)
occ_det[1,3] <- round(m1CI.psi[1,1],2)
occ_det[1,4] <- round(m1CI.psi[1,2],2)
occ_det[1,5] <- round(m1.p@estimate,2)
occ_det[1,6] <- round(m1CI.p[1,1],2)
occ_det[1,7] <- round(m1CI.p[1,2],2)  
occ_det[2,2] <- round(m1.psi_col@estimate,2)
occ_det[2,3] <- round(m1CI.psi_col[1,1],2)
occ_det[2,4] <- round(m1CI.psi_col[1,2],2)
occ_det[2,5] <- round(m1.p_col@estimate,2)
occ_det[2,6] <- round(m1CI.p_col[1,1],2)
occ_det[2,7] <- round(m1CI.p_col[1,2],2)  
occ_det[3,2] <- round(m1.psi_bla@estimate,2)
occ_det[3,3] <- round(m1CI.psi_bla[1,1],2)
occ_det[3,4] <- round(m1CI.psi_bla[1,2],2)
occ_det[3,5] <- round(m1.p_bla@estimate,2)
occ_det[3,6] <- round(m1CI.p_bla[1,1],2)
occ_det[3,7] <- round(m1CI.p_bla[1,2],2) 
view(occ_det)
