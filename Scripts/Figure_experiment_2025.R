# Date: June 22nd 2025
# Purpose: Create graphs comparing differences in transmission experiment with wounding
# Seagrass in Palau

library(tidyverse)
library(ggplot2)
library(readxl)
library(DHARMa)
library(ggpattern)

exp_data <- read_xlsx("Data/Experiment_2025.xlsx") #%>% 
  #filter(Questionable == 0)


##### Statistics ###############################################################
exp_data 

library(glmmTMB)

# Do the analysis but separate out lesioned and wounded as two treatments in model
m2 <- glmmTMB(Disease ~ lesioned*wounded, data = exp_data, family = binomial)
sim_m2 <- simulateResiduals(m2, plot = F)
plot(sim_m2)
summary(m2)

# 
exp(1.4663) # Odds
exp(1.4663)/(1 + exp(1.4663)) # Probability

anova_test <- car::Anova(m2, type = "III")
car::Anova(m2, type = "II")

##### Plotting #################################################################

plot2 <- exp_data %>% 
  group_by(Type) %>% 
  summarize(proportion_disease = sum(Disease)/n()) %>% 
  ungroup() %>% 
  mutate(health = c("Healthy", "Healthy", "Lesioned", "Lesioned"),
         wounded = c("Not Wounded", "Wounded", "Not Wounded", "Wounded"),
         Treatment = c("C", "D", "A", "B"),
         Type = str_replace_all(Type, "_", " / ")) %>%
  arrange(-proportion_disease) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  ggplot(aes(fill=wounded, y=proportion_disease, x=reorder(health, proportion_disease, decreasing = TRUE))) + 
    geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("#2C5F2D", "bisque")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.9),
        axis.line.y = element_line(color = "grey80"),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.title = element_text(color = "grey30")) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(y = 0,label = Treatment), vjust = 1.5, position=position_dodge(width = 1), stat="identity", size = 3.3, color = "#4D4D4DFF") +
  geom_hline(yintercept = 0, color = "grey80") +
  labs(x = NULL, y = "Percent of blades with lesions on Day 6")
plot2 

plot3 <- exp_data %>% 
  group_by(Type) %>% 
  summarize(proportion_disease = sum(Disease)/n()) %>% 
  ungroup() %>% 
  mutate(health = c("Healthy", "Healthy", "Lesioned", "Lesioned"),
         wounded = c("Not Wounded", "Wounded", "Not Wounded", "Wounded"),
         Treatment = c("C", "D", "A", "B"),
         Type = str_replace_all(Type, "_", " / ")) %>%
  arrange(-proportion_disease) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  ggplot(aes(fill=wounded, y=proportion_disease, x=reorder(health, proportion_disease, decreasing = TRUE))) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("#2C5F2D", "bisque")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line.y = element_line(color = "grey80"),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.7),
        axis.text.x = element_text(color = "black")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.6), n.breaks = 7) +
  geom_text(aes(y = 0,label = Treatment), vjust = 1.5, position=position_dodge(width = 1), stat="identity", size = 3.3, color = "#4D4D4DFF") +
  labs(x = "Treatment Group", y = "Percent of blades with lesions on Day 6") +
  geom_hline(yintercept = 0, color = "grey80") 
  

library(ggpubr)
library(rstatix)

stat.test <- broom::tidy(anova_test) %>%
  rename(group1 = term) %>%
  # Example: compare everything to the first level or specific groups
  mutate(group2 = "Healthy", 
         p = round(p.value, 3), 
         y.position = 0.55) %>% 
  filter(group1 == "lesioned") %>% 
  mutate(group1 = "Lesioned")

plot_w_sig <- plot3 + stat_pvalue_manual(
  stat.test,  
  label = "p = {p}", 
  tip.length = 0) 

plot_w_sig
#### Export #################################################################### 

ggsave(plot_w_sig, 
       file = "Figures/experiment_2025_v2.png", 
       width = 4, height = 4, units = "in") 

