#EWT_Library1
data <- read.table("EWT1_start.csv",sep=",",head=T)
colnames(data)

#subset of data containing just Brachyteles hypoxanthus
pc_1 <- subset(data,scientific_name=="Brachyteles hypoxanthus")

#find row with highest PC reads
pc_1 <- pc_1[1,]
y <- sum(pc_1[,16:103])
z <- sum(pc_1[,104:107])
pc <- y/z

# Select different groups of columns

data_taxo <- data[,1:15] 
data_samples <- data[,16:89]
data_blank <- data[,91:103] 
data_PC <- data[,104:107] 
data_seqs <- data[,108:109]


# Select amount of reads to remove based on tagjumping percentage (using positive control samples)

calculate_reads <- (data_samples*pc) #calculate the amount of reads to remove for each MOTU based on the previously defined percentage 
calculate_reads <- round(calculate_reads)
data_samples <- data_samples-calculate_reads

# Reconstruct the matrix

data_seqs$total_reads <- rowSums(data_samples,na.rm = T)
ewt1_PC <- data.frame(data_taxo,data_blank,data_samples,data_seqs)
sum(ewt1_PC$total_reads)

# Write the new output

write.table(ewt1_PC,"EWT1_PC_Corrected.csv",row.names=F,sep=";")

ewt1_pc <- subset(ewt1_PC, family_name!="Atelidae") # removing human reads 

sum(ewt1_pc$total_reads)

data <- ewt1_pc
colnames(data)

# Select different groups of columns

data_taxo <- data[,1:15] 
data_blank <- data[,16:28] 
data_samples <- data[,29:102]
data_seqs <- data[,103:104]


# Select blank reads to remove
data_blank$max <- do.call(pmax, data_blank)
blank_remove <- data_blank$max


# remove blank reads from samples for every row
for (i in 1:nrow(data_samples)){
  data_samples[i,] <- data_samples[i,] - data_blank$max[i]
  data_samples[i,][data_samples[i,]<2] <- 0
}

# Reconstruct the matrix
data_seqs$total_reads <- rowSums(data_samples,na.rm = T)

ewt1_blanks <- data.frame(data_taxo,data_samples,data_seqs)

# Save the new matrix

write.table(ewt1_blanks,"EWT1_Blanks_Corrected.csv",row.names=F,sep=";")

#Additional filtering steps to remove reads belonging to humans and domestic animals
sum(ewt1_blanks$total_reads) #obtain the total of reads after removing reads from the blanks 
ewt1_mammals <- subset(ewt1_blanks, class_name=="Mammalia") #select only mammals
sum(ewt1_mammals$total_reads)
write.csv(ewt1_mammals, "EWT1_Mammal_MOTUs.csv")
ewt1_mammals_nonhuman <- subset(ewt1_mammals, family_name!="Hominidae") # removing human reads 
sum(ewt1_mammals_nonhuman$total_reads)
write.csv(ewt1_mammals_nonhuman, "EWT1_Human_Corrected.csv")
ewt1_mammals_wild <- subset(ewt1_mammals_nonhuman, genus_name!= "Ovis")# removing domestic animals 
ewt1_mammals_wild1 <- subset(ewt1_mammals_wild, genus_name!= "Bos")
ewt1_mammals_wild2 <- subset(ewt1_mammals_wild1, genus_name!= "Sus")
ewt1_mammals_wild3 <- subset(ewt1_mammals_wild2, genus_name!= "Equus")
ewt1_mammals_wild4 <- subset(ewt1_mammals_wild3, genus_name!= "Canis")
ewt1_mammals_wild5 <- subset(ewt1_mammals_wild4, genus_name!= "Felis")
sum(ewt1_mammals_wild5$total_reads)
write.csv(ewt1_mammals_wild5, "EWT1_Wild_Mammals.csv")

#filter out MOTUs with a best identity <98% and fewer than 5 reads in total
ewt1_mammals_fin <- subset(ewt1_mammals_wild5, best_identity >= 0.98)
ewt1_mammals_final <- subset(ewt1_mammals_fin, total_reads > 5)
sum(ewt1_mammals_final$total_reads)

#final filtered csv file
write.table(ewt1_mammals_final,"EWT1_Final.csv",row.names=F,sep=";")

#EWT_Library2
data <- read.table("EWT2_start.csv",sep=",",head=T)
colnames(data)

#subset of data containing Brachyteles hypoxanthus only
pc_2 <- subset(data,scientific_name=="Brachyteles hypoxanthus")

#find row with highest PC reads
pc_2 <- pc_2[1,]
y <- sum(pc_2[,16:103])
z <- sum(pc_2[,104:107])
pc <- y/z

# Select different groups of columns
data_taxo <- data[,1:15] 
data_samples <- data[,16:90]
data_blank <- data[,91:103] 
data_PC <- data[,104:107] 
data_seqs <- data[,108:109]

# Select amount of reads to remove based on tagjumping percentage (using positive control samples)
calculate_reads <- (data_samples*pc) #calculate the amount of reads to remove for each MOTU based on the previously defined percentage 
calculate_reads <- round(calculate_reads)
data_samples <- data_samples-calculate_reads

# Reconstruct the matrix
data_seqs$total_reads <- rowSums(data_samples,na.rm = T)
ewt2_PC <- data.frame(data_taxo,data_blank,data_samples,data_seqs)
sum(ewt2_PC$total_reads)

# Write the new output
write.table(ewt2_PC,"EWT2_PC_Corrected.csv",row.names=F,sep=";")
ewt2_pc <- subset(ewt2_PC, family_name!="Atelidae") # removing human reads 
sum(ewt2_pc$total_reads)
data <- ewt2_pc
colnames(data)

# Select different groups of columns
data_taxo <- data[,1:15] 
data_blank <- data[,16:28] 
data_samples <- data[,29:103]
data_seqs <- data[,104:105]

# Select blank reads to remove
data_blank$max <- do.call(pmax, data_blank)
blank_remove <- data_blank$max
blank_remove

# remove blank reads from samples for every row
for (i in 1:nrow(data_samples)){
  data_samples[i,] <- data_samples[i,] - data_blank$max[i]
  data_samples[i,][data_samples[i,]<2] <- 0
}

# Reconstruct the matrix
data_seqs$total_reads <- rowSums(data_samples,na.rm = T)
ewt2_blanks <- data.frame(data_taxo,data_samples,data_seqs)

# Save the new matrix
write.table(ewt2_blanks,"EWT2_blankscorrected.csv",row.names=F,sep=";")

#Additional filtering steps to remove reads belonging to humans and domestic animals
sum(ewt2_blanks$total_reads) #obtain the total of reads after removing reads from the blanks 
ewt2_mammals <- subset(ewt2_blanks, class_name=="Mammalia") #select only mammals
sum(ewt2_mammals$total_reads)
write.table(ewt2_mammals, "EWT2_Mammal_MOTUs.csv", row.names = F, sep=";")
ewt2_mammals_nonhuman <- subset(ewt2_mammals, family_name!="Hominidae") # removing human reads 
sum(ewt2_mammals_nonhuman$total_reads)
write.table(ewt2_mammals_nonhuman, "EWT2_Human_Corrected.csv", row.names = F, sep=";")
ewt2_mammals_wild <- subset(ewt2_mammals_nonhuman, genus_name!= "Ovis")# removing domestic animals 
ewt2_mammals_wild1 <- subset(ewt2_mammals_wild, genus_name!= "Bos")
ewt2_mammals_wild2 <- subset(ewt2_mammals_wild1, genus_name!= "Sus")
ewt2_mammals_wild3 <- subset(ewt2_mammals_wild2, genus_name!= "Equus")
ewt2_mammals_wild4 <- subset(ewt2_mammals_wild3, genus_name!= "Canis")
ewt2_mammals_wild5 <- subset(ewt2_mammals_wild4, genus_name!= "Felis")
sum(ewt2_mammals_wild5$total_reads)
write.csv(ewt2_mammals_wild5, "EWT2_Wild_Mammals.csv")

#filtering to exclude MOTUs with a best identity <98% and fewer than 5 reads in total
ewt2_mammals_fin <- subset(ewt2_mammals_wild5, best_identity >= 0.98)
ewt2_mammals_final <- subset(ewt2_mammals_fin, total_reads > 5)
sum(ewt2_mammals_final$total_reads)

#creating fina filtered csv file
write.table(ewt2_mammals_final,"EWT2_Final.csv",row.names=F,sep=";")