#############################################################################################################
## Plot and modify genetic map
##
## 
## 
## Last Modified Mat 2, 2023
#################################################################################################
#Import libraries
library(tidyverse)

gen_map <- read_csv("data/genomic_data/CE10_v2_genome_map.csv") %>% 
  filter(chr2.0 %in% c("CE10_chr1",
                       "CE10_chr2",
                       "CE10_chr3",
                       "CE10_chr4",
                       "CE10_chr5",
                       "CE10_chr6",
                       "CE10_chr7",
                       "CE10_chr8"))

#################################################################################################

#Make initial plot
ggplot(gen_map, aes(x=bp2.0, y=cM)) + 
  geom_point(aes(), size =1)+
#  geom_text(hjust=-.15, vjust=-.2)+
#  scale_y_continuous(name="Lambda Slope")+
#  scale_x_continuous(name="BF>30 Offset")+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))+
facet_wrap(~chr2.0,ncol=4)

#################################################################################################
##Fix Inverted Chromosomes

#Separate by chromosome
gen_map_chr1 <- gen_map %>% filter(chr2.0 == "CE10_chr1")
gen_map_chr2 <- gen_map %>% filter(chr2.0 == "CE10_chr2")
gen_map_chr3 <- gen_map %>% filter(chr2.0 == "CE10_chr3")
gen_map_chr4 <- gen_map %>% filter(chr2.0 == "CE10_chr4")
gen_map_chr5 <- gen_map %>% filter(chr2.0 == "CE10_chr5")
gen_map_chr6 <- gen_map %>% filter(chr2.0 == "CE10_chr6")
gen_map_chr7 <- gen_map %>% filter(chr2.0 == "CE10_chr7")
gen_map_chr8 <- gen_map %>% filter(chr2.0 == "CE10_chr8")

#gen_map_chr1 <- gen_map_chr1[order(gen_map_chr1$bp2.0),]
#gen_map_chr2 <- gen_map_chr2[order(gen_map_chr1$bp2.0),]
#gen_map_chr3 <- gen_map_chr3[order(gen_map_chr1$bp2.0),]
#gen_map_chr4 <- gen_map_chr4[order(gen_map_chr1$bp2.0),]
#gen_map_chr5 <- gen_map_chr5[order(gen_map_chr1$bp2.0),]
#gen_map_chr6 <- gen_map_chr6[order(gen_map_chr1$bp2.0),]
#gen_map_chr7 <- gen_map_chr7[order(gen_map_chr1$bp2.0),]
#gen_map_chr8 <- gen_map_chr8[order(gen_map_chr1$bp2.0),]



#Flip inverted chromosomes
gen_map_chr2 <- gen_map_chr2 %>% mutate(bp2.0=max(gen_map_chr2$bp2.0)-bp2.0)
gen_map_chr6 <- gen_map_chr6 %>% mutate(bp2.0=max(gen_map_chr6$bp2.0)-bp2.0)
gen_map_chr7 <- gen_map_chr7 %>% mutate(bp2.0=max(gen_map_chr7$bp2.0)-bp2.0)
gen_map_chr8 <- gen_map_chr8 %>% mutate(bp2.0=max(gen_map_chr8$bp2.0)-bp2.0)

#Plot afterwards
gen_map_update <- rbind(gen_map_chr1,
                       gen_map_chr2,
                       gen_map_chr3,
                       gen_map_chr4,
                       gen_map_chr5,
                       gen_map_chr6,
                       gen_map_chr7,
                       gen_map_chr8)


#Make initial plot
ggplot(gen_map_chr2, aes(x=bp2.0, y=cM)) + 
  geom_point(aes(), size =1)+
  #  geom_text(hjust=-.15, vjust=-.2)+
  #  scale_y_continuous(name="Lambda Slope")+
  #  scale_x_continuous(name="BF>30 Offset")+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))+
  facet_wrap(~chr2.0,ncol=4)

#################################################################################################
## Export and manually curate
#write_csv(gen_map_update, "data/genomic_data/gen_map_update.csv")

#Import after manual curation

gen_map_order <- read_csv("data/genomic_data/gen_map_order.csv")

#Make new plot
ggplot(gen_map_order, aes(x=bp2.0, y=cM)) + 
  geom_point(aes(), size =1)+
  #  geom_text(hjust=-.15, vjust=-.2)+
  #  scale_y_continuous(name="Lambda Slope")+
  #  scale_x_continuous(name="BF>30 Offset")+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    #legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.position = c(0.85, 0.25),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))+
  facet_wrap(~chr2.0,ncol=4)








