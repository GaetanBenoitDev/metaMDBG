
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(gcookbook)
library(dplyr)

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

dataFilename = paste0(script.basename, "/data_performances.csv")
outputFilename = paste0(script.basename, "/fig_runtime.png")
print(outputFilename)

png(outputFilename, 800, 800)
#pdf(getwd() + "/test.pdf")

#df <- data.frame(dose=c("D0.5", "D1", "D2", "D0.5", "D1", "D2"),
#                len=c(4.2, 10, 29.5, 8, 11, 12),
#                supp=c("OO", "OO", "OO", "JJ", "JJ", "JJ"))

df <- read.csv(dataFilename, header=TRUE)
df$Software = factor(df$Software, levels=unique(df$Software))
df$Dataset = factor(df$Dataset, levels=unique(df$Dataset))
#df$Type = factor(df$Type, levels=rev(unique(df$Type)))

head(df)
#lala.lala

#ce <- cabbage_exp %>% arrange(Date, rev(Cultivar))

# Calculate y position, placing it in the middle
#ce <- ce %>%
#  group_by(df$Software) %>%
#  mutate(label_y = cumsum(df$Count) - 0.5 * df$Count)

ggplot(df, aes(x=Software, y=Runtime)) +
    geom_bar(stat="identity", width = 0.8, fill="tomato1") + #, color="black" , color="brown4"
    scale_y_continuous(expand = c(0, 0), limits=c(0, 650)) + #
    scale_x_discrete(expand = c(0.3, 0)) +
    theme(axis.title.x=element_blank()) +
    theme(text = element_text(size = 30)) + 
    #theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + #Hide y axis
    theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
    theme(legend.title=element_blank()) + 
    #theme(panel.border = element_rect(color = "black", size = 1)) +
    ylab("Runtime (h)")+
    facet_wrap(~Dataset, nrow = 1, strip.position="bottom") +
    theme(panel.spacing = unit(2, "lines"))  + 
    theme(axis.line.y = element_line(color="black", size = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.5)) +
    theme(strip.placement = 'outide', strip.background = element_blank()) +
    #scale_fill_manual(values = c("#99CCFF", "#33FF33")) +
    #scale_fill_manual(values = c("darkolivegreen4", "darkolivegreen3")) +#+ theme(strip.placement = 'bottom') + #+ theme_classic() #+ theme_classic() +  
    geom_text(aes(label = round(Runtime, digits = 0)), vjust=-1, size=6)  +
    theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #Remove facet background
    theme(panel.background = element_blank()) #+ #Remove grey background
    #theme(legend.position = c(0.25,0.9)) + 
    #guides(fill = guide_legend(reverse=TRUE))

    #geom_text(aes(y = label_y, label = Count), vjust = 1.5, colour = "white")
        
  
  
  
  #+ theme(panel.spacing = unit(0, "lines"))

#ggplot(data=df, aes(x=dose, y=Count, fill=supp)) +
#  geom_bar(stat="identity")#+
  #geom_text(aes(y=len, label=len), vjust=1.6, color="white", size=3.5)+
  #scale_fill_brewer(palette="Paired")
  #+ theme_minimal()


dev.off()