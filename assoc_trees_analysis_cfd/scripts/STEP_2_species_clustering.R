#Catherine dell'Olio

#SPECIES CLUSTERING
#Emma suggested I run an NMDS to get a better idea if clusters of species are present and thus might weigh together instead of separately

#you will need the following packages:
library(dendextend)
library(ggdendro)
library(ggplot2)
library(vegan)

# Source - https://stackoverflow.com/a/55778423
# Posted by Moritz Schmid
# Retrieved 2026-04-15, License - CC BY-SA 4.0
# calculate the dissimilarity matrix with vegdist so you can use the sorenson/bray 
#method
distBray <- vegdist(df[1:nrow(df), 2:ncol(df)], method = "bray") 

# calculate the clusters with ward.D2
clust1 <- hclust(distBray, method = "ward.D2")
clust1
clust1[["labels"]] <-abundance_DA$focal_tree

#plot the cluster dendrogram with dendextend
dend <- clust1 %>% as.dendrogram %>%
  set("branches_k_color", k = 5) %>% set("branches_lwd", 0.5)  %>%  set("clear_leaves") %>% set("labels_colors", k = 5)  %>% set("leaves_cex", 0.5) %>%
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dend)
ggplot(ggd1, horiz = TRUE)

#group 1 2 and 3 are more related than group 4 and 5. group 1 and 2 are more related to each other than to 3.
group1 <- c("109", "22", "104", "125", "231", "96", "91", "182", "238", "236", "23", "234", "169", "68", "93", "196", "71", "64", "100", "101", "67a", "587", "126", "230", "127", "175")
group2 <- c("92", "78", "124", "72", "123")
group3 <- c("50", "45", "68", "27", "60", "6", "7", "228", "130", "61", "203", "57", "207")
group4 <- c("226", "44", "40", "223")
group5 <- c("199", "53", "582", "39", "34", "216", "157", "36")

group1df <- df %>% filter(focal_tree %in% group1)
group2df <- df %>% filter(focal_tree %in% group2)
group3df <- df %>% filter(focal_tree %in% group3)
group4df <- df %>% filter(focal_tree %in% group4)
group5df <- df %>% filter(focal_tree %in% group5)

#listing these groups together with a column to show the grouping.
#first I will start with two treatments: groups 1 2 and 3 (trt 0) versus groups 4 and 5 (trt 1)

group1df <- group1df %>% mutate(Trt = "0")
group2df <- group2df %>% mutate(Trt = "0")
group3df <- group3df %>% mutate(Trt = "0")
group4df <- group4df %>% mutate(Trt = "1")
group5df <- group5df %>% mutate(Trt = "1")
NMDS_abundances <- rbind(group1df, group2df, group3df, group4df, group5df)

data <- NMDS_abundances[2:(ncol(df)-1)]
trt <- NMDS_abundances[ncol(NMDS_abundances)]

dim1<-metaMDS(data, k = 1, trymax = 100)
s1<-dim1[[22]]
s1
#s1 = 0.42
dim2<-metaMDS(data,k=2,trymax = 100) #ran 20 times
s2<-dim2[[22]]
s2
#s2 = 0.25
#this is borderline but really too high for good interpretation
dim3<-metaMDS(data,k=3,trymax = 100) #20 times
s3<-dim3[[22]]
s3
#s3 = 0.18

#NMDS plot
ord<-metaMDS(data,k=2) # ran 20 times

ano <-  anosim(data, trt$Trt, distance = "bray", permutations = 9999)
plot(ano)

data.scores <- as.data.frame(scores(ord)$sites)
data.scores$trt <- trt$Trt

## this next section is how to make connected centroids (helps visualize comparison)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ trt, data = data.scores, FUN = mean)

segs <- merge(data.scores, setNames(cent, c('Site','oNMDS1','oNMDS2')),
              by = trt$trt, sort = FALSE)

all.trt.nmds <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes( shape = factor(trt), colour = factor(trt))) + 
  geom_segment(data = segs,
               mapping = aes(xend = 5, yend = 5)) + 
  geom_point(data = cent, size = 5, colour=c("goldenrod","darkblue"), shape=c(0,1)) +
  scale_shape_manual(values=c(0,1,2))+
  xlim(-3,3)+
  ylim(-3,3)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "none", 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2))+
  labs(x = "NMDS1", y = "NMDS2")  + 
  scale_colour_manual(values = c("goldenrod", "darkblue"))
plot(all.trt.nmds)

simper_data <- simper(data, trt)
summary(simper_data)