library('tidyverse')
library('gridExtra')
library(gghighlight)

# work with dates in radiocarbon space

# read in dates
dates_orig <- read.delim("data/processed/master_dates_file.txt",
                    sep="\t", stringsAsFactors = F)

# add 'processed' names ----
# read in mammal taxonomy file
mammalTaxa <- read.delim("data/raw/TaxonomyMatchingFile.txt", sep="\t")
matchedTaxa<- mammalTaxa[match(dates_orig$prelim_taxon_name, mammalTaxa$OriginalName),
                         c('RevisedName', 'Genus')]
dates <- cbind(dates_orig, matchedTaxa)

# filter out "no" dates and Box 999 dates
dates <- dates[-which(dates$UseSample=="N"),]
dates <- dates[-which(dates$box=="999"),]

# Add column indicating a "too old" flag
GreaterThan <- rep(0, nrow(dates))
GreaterThan[grep(">", dates$X14C_age_BP)] <- 1
dates <- cbind(dates, GreaterThan)
dates$X14C_age_BP <- gsub(">", "", dates$X14C_age_BP) #remove the > from the old dates and convert to numeric from factor
dates$X14C_age_BP <- as.numeric(dates$X14C_age_BP)

# Rename 14C
colnames(dates)[which(colnames(dates)=="X14C_age_BP")] <- "C14_age_BP"
colnames(dates)[which(colnames(dates)=="X14C_age_error")] <- "C14_age_error"

# re-order factor levels for plotting
dates$box <- factor(dates$box, levels = rev(c("14", "7b", "13", "1", "HC", "999", "4", "10")))
dates$RevisedName <- factor(dates$RevisedName, levels = rev(c("Sylvilagus sp", "Sylvilagus bachmani", "Sylvilagus audubonii", "Lepus sp", "Otospermophilus beecheyi", "Neotoma sp", "Thomomys sp", "Mustela frenata", "Canis latrans")))

# select only the columns I need
dates <- select(dates, "UCIAMS_Number", "C14_age_BP", "C14_age_error", "box", "Canister", "RevisedName", "Genus", "GreaterThan")
highlight_df <- dates[which(dates$Canister == "misc"),]
#highlight_df <- dates[c(which(is.na(dates$Canister)), which(dates$Canister == "misc")),]
  
########## PLOTTING ###########
## BY BOX ##
pdf(file="output/boxplot-alldates.pdf", height=6, width=12)
base <- ggplot(dates) +
  geom_boxplot(aes(x = box, y = C14_age_BP, fill = box), alpha = 0.7, show.legend = FALSE) +
  geom_point(aes(x = box, y = C14_age_BP, fill = box), 
             size = 3, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +  #
#  geom_point(data=highlight_df, aes(x = box, y = C14_age_BP), fill="gray", size = 3, shape = 21, show.legend = FALSE) +
  xlab("Box") + 
  ylab("C14 years BP")
all <- base + 
  theme_light()+
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=28,face="bold"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
all + coord_flip() + scale_y_reverse()
dev.off()

## BY TAXON ##
pdf(file="output/boxplot-alltaxa.pdf", height=6, width=12)
base <- ggplot(dates) +
  geom_boxplot(aes(x = RevisedName, y = C14_age_BP, fill = RevisedName), 
               alpha = 0.7, show.legend = FALSE) +
  geom_point(aes(x = RevisedName, y = C14_age_BP, fill = RevisedName), 
             size = 3, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +
  xlab("Taxon") + 
  ylab("C14 years BP") 
all <- base + 
  theme_light()+
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=28,face="bold"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
  
all + scale_y_reverse() + coord_flip()
dev.off()

## start working with individual taxa. For this, focus only on the P23 boxes
P23_dates_taxa <- dates[c(which(dates$box == "1"), 
                     which(dates$box == "7b"), 
                     which(dates$box == "13"), 
                     which(dates$box == "14")),]
P23_dates_taxa <- P23_dates_taxa[c(which(P23_dates_taxa$RevisedName == 
                                      "Sylvilagus audubonii"), 
                              which(P23_dates_taxa$RevisedName==
                                      "Sylvilagus bachmani"), 
                              which(P23_dates_taxa$RevisedName==
                                      "Sylvilagus sp"), 
                              which(P23_dates_taxa$RevisedName==
                                      "Otospermophilus beecheyi")),]
P23_dates_taxa <- droplevels(P23_dates_taxa)
#levels(P23_dates_taxa$box)<- c("14", "7b", "13", "1")

## NEED to remove legend OR remove labels from X axis and change legend title (better option)
base <- ggplot(P23_dates_taxa, aes(y = RevisedName, x=C14_age_BP))
plot <- base + 
  geom_boxplot(aes(y = RevisedName, x = C14_age_BP, fill = RevisedName), 
               alpha = 0.7) +
  geom_point(aes(y = RevisedName, x = C14_age_BP, fill = RevisedName), 
             size = 3, shape = 21, position = position_jitterdodge()) +
  facet_grid(box ~ .) +
  theme_light() +
  ylab("Taxon") + 
  xlab("C14 years BP")
plot + scale_x_reverse() 

# Plot at genus level to group Sylvilagus together
base <- ggplot(P23_dates_taxa, aes(y = Genus, x=C14_age_BP))
plot <- base + 
  geom_boxplot(aes(y = Genus, x = C14_age_BP, fill = Genus), 
               alpha = 0.7) +
  geom_point(aes(y = Genus, x = C14_age_BP, fill = Genus), 
             size = 3, shape = 21, position = position_jitterdodge()) +
  facet_grid(fct_rev(box) ~ .) +
  theme_light() +
  ylab("Taxon") + 
  xlab("C14 years BP")
pdf(file="output/boxplot-Sylv_Oto.pdf", height=6, width=12)
plot + scale_x_reverse()
dev.off()

# Plot only the P23 boxes
P23_dates <- dates[c(which(dates$box == "1"), 
which(dates$box == "7b"), 
which(dates$box == "13"), 
which(dates$box == "14")),]
P23_dates$box <- factor(P23_dates$box,
levels = c('1','13', '14', '7b'),ordered = TRUE)

base <- ggplot(P23_dates) +
  geom_boxplot(aes(x = box, y = C14_age_BP, fill = box), alpha = 0.7, show.legend = FALSE) +
  geom_point(aes(x = box, y = C14_age_BP, fill = box), 
             size = 3, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +  #
  #  geom_point(data=highlight_df, aes(x = box, y = C14_age_BP), fill="gray", size = 3, shape = 21, show.legend = FALSE) +
  xlab("Box") + 
  ylab("C14 years BP")
all <- base + 
  theme_light()+
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=28,face="bold"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

pdf(file="output/boxplot-P23_dates.pdf", height=6, width=12)
all + coord_flip() + scale_y_reverse()
dev.off()


anova <- aov(C14_age_BP~RevisedName, data=P23_dates_taxa[which(P23_dates_taxa$box=="1"),])
anova <- aov(C14_age_BP~RevisedName, data=P23_dates_taxa[which(P23_dates_taxa$box=="14"),])
anova <- aov(C14_age_BP~RevisedName, data=P23_dates_taxa[which(P23_dates_taxa$box=="13"),])
anova <- aov(C14_age_BP~RevisedName, data=P23_dates_taxa[which(P23_dates_taxa$box=="7b"),])

# Plot each box separately and run stats
P23_boxes <- c("14", "7b", "13", "1")
Focal_taxa <- c("Sylvilagus audubonii", "Sylvilagus bachmani", "Otospermophilus beecheyi")
for (i in 1:length(P23_boxes)){
  
  # create box/taxon dataframes
Saudu <- P23_dates %>% filter(RevisedName == "Sylvilagus audubonii", box==P23_boxes[i])
Sbach14 <- P23_dates %>% filter(RevisedName == "Sylvilagus bachmani", box==P23_boxes[i])
Obeech24 <- P23_dates %>% filter(RevisedName == "Otospermophilus beecheyi", box==P23_boxes[i])


base <- ggplot(P23_dates) 
base +
  geom_boxplot(aes(x = box, y = C14_age_BP, fill = box), alpha = 0.7) +
  geom_point(aes(x = box, y = C14_age_BP, fill = box), 
             size = 3, shape = 21, position = position_jitterdodge()) +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right", 
        legend.text=element_text(size=12)) +
  theme_light() +
  xlab("Box") + 
  ylab("C14 years BP")


}

Sylv <- P23_dates %>% filter(RevisedName == "Sylvilagus sp")

# Box 14
P23_dates_taxa[which(is.na(P23_dates_taxa$Canister)),]





### Plot only boxes 1, 7b, and 13

# set colors
scale_fill_manual() for box plot, bar plot, violin plot, etc
scale_color_manual() for lines and points
# Box plot
bp + scale_fill_manual(values=c("yellow", "teal", "purple"))
# Scatter plot
sp + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# Plot only the P23 boxes
P23_dates_sub <- dates[c(which(dates$box == "1"), 
                     which(dates$box == "7b"), 
                     which(dates$box == "13")),]
P23_dates_sub$box <- factor(P23_dates_sub$box,
                        levels = c('1','13', '7b'),ordered = TRUE)

base <- ggplot(P23_dates_sub) +
  geom_boxplot(aes(x = box, y = C14_age_BP, fill = box), alpha = 0.7, show.legend = FALSE) +
  geom_point(aes(x = box, y = C14_age_BP, fill = box), 
             size = 3, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +  
  xlab("Box") + 
  ylab("C14 years BP")
all <- base + 
  scale_fill_brewer(palette="Set1")+
  theme_light()+
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

pdf(file="output/boxplot-P23_dates_no14.pdf", height=6, width=12)
all + coord_flip() + scale_y_reverse()
dev.off()
