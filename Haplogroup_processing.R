setwd("/home/jnguy254/projects/def-ccastel/SharedResources/CLSA_private/results/genomics")

library(data.table)
suppressMessages(library(ggpubr))

hp <- fread("haplogroups.txt")

# add in subgroupings for haplogroups
hp$subgroup1 <- substr(hp$Haplogroup, start = 1, stop = 1)
hp$subgroup2 <- substr(hp$Haplogroup, start = 1, stop = 2)
write.table(hp, file="haplogroups.txt", row.names = FALSE)

# load in epi_responses
responses <- fread("/home/jnguy254/projects/def-ccastel/SharedResources/CLSA_private/data/baseline/2203013_WesternU_CCastellani_Baseline_CoPv7.csv")
epi_responses <- responses[-which(is.na(responses$ADM_EPIGEN2_COM)),]
epi_responses <- epi_responses[-which(is.na(epi_responses$ADM_GWAS3_COM)),]


hp_epi <- hp[which(hp$ADM_GWAS_COM %in% epi_responses$ADM_GWAS3_COM),]

write.table(hp_epi, file="haplogroups_epi.txt", row.names = FALSE)

# density of haplogroup quality scores
p1 <- ggplot(hp, aes(x=Quality)) + 
  geom_density(alpha = 0.25, fill = "#BCE5DD") + 
  geom_rug(linewidth = 1.5, col = "#1D4B42") +
  geom_vline(xintercept = 0.53, linetype = "longdash", color = "#1D4B42") +
  annotate("text", x=0.52, y = 3.7, label = "q = 0.53 | N = 20976", angle = 90) +
  geom_vline(xintercept = 0.625, linetype = "longdash", color = "#1D4B42") +
  annotate("text", x=0.615, y = 3.7, label = "q = 0.625 | N = 13598", angle = 90) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size=12, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "N = 26626")
p2 <- ggplot(hp_epi, aes(x=Quality)) + 
  geom_density(alpha = 0.25, fill = "#BCE5DD") + 
  geom_rug(linewidth = 1.5, col = "#1D4B42") +
  geom_vline(xintercept = 0.54, linetype = "longdash", color = "#1D4B42") +
  annotate("text", x=0.53, y = 2.4, label = "q = 0.54 | N = 1109", angle = 90) +
  geom_vline(xintercept = 0.635, linetype = "longdash", color = "#1D4B42") +
  annotate("text", x=0.625, y = 2.4, label = "q = 0.635 | N = 747", angle = 90) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size=12, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "N = 1445")
png("HP_QualityScores.png", width=350, height=175, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, ncol = 2)
dev.off()