library(tidyverse)
library(lme4)
library(ggeffects) 
library(MuMIn)
library(DHARMa)
library(ggpubr)
library(wesanderson)
library(extrafont)
library(svglite)

# update working dir
setwd("")

loadfonts(device = "postscript") 

#### data cleaning ####

# read in dataframe
raw_data <- read_csv("data_raw.csv")
                     
data_clean <- raw_data %>%
  replace(is.na(.), 0)  %>%
  select(-sm, -sens) %>%
  rename(replicate = rep)
  
# convert to binomial 

binom_data <- pmap_dfr(data_clean,
                      function(timepoint, phage, glucose, replicate, crispr, 
                               moi, cell_density, phage_density) {
                        data.frame(timepoint = timepoint,
                                   phage = phage,
                                   glucose = glucose,
                                   replicate = replicate,
                                   moi=moi,
                                   cell_density=cell_density,
                                   phage_density=phage_density,
                                   crispr = c(rep(1, crispr),
                                             rep(0, 24 - crispr)))})

model_df <- binom_data %>%
  unite(treatment, phage, glucose, replicate, sep = "_", remove = FALSE) %>%
  mutate(treatment = as.character(treatment)) %>%
  mutate(product = cell_density * phage_density) %>%
  mutate(log_foi = log10(product)) %>%
  mutate(log_cell = log10(cell_density)) %>%
  mutate(log_phage = log10(phage_density))

######### model for each timepoint #########

# make dfs #

tp1_df <- model_df %>%
  filter(timepoint == "1")

tp2_df <- model_df %>%
  filter(timepoint == "2")

tp3_df <- model_df %>%
  filter(timepoint == "3")

### timepoint 1 ###

glm_tp1_max <- glmer(crispr ~ phage + glucose + log_cell + log_phage + (1|treatment), 
                 data=tp1_df, na.action = 'na.fail', family="binomial")

# model summary
summary(glm_tp1_max)

# model selection
dredge1 <- dredge(glm_tp1_max)

# final model = simplest within 2 delta AICs

glm_tp1 <- glmer(crispr ~ phage + log_phage + (1|treatment), 
                 data=tp1_df, na.action = 'na.fail', family="binomial")

# model summary
summary(glm_tp1)

# test dispersion of final model
testDispersion(glm_tp1)
simulationOutput <- simulateResiduals(fittedModel = glm_tp1)
plot(simulationOutput)

### timepoint 2 ###

glm_tp2_max <- glmer(crispr ~ phage + glucose + log_cell + log_phage + (1|treatment), 
                     data=tp2_df, na.action = 'na.fail', family="binomial")

# model summary
summary(glm_tp2_max)

# model selection
dredge2 <- dredge(glm_tp2_max)

# final model = simplest within 2 delta AICs

glm_tp2 <- glmer(crispr ~ log_cell + log_phage + (1|treatment), 
                     data=tp2_df, na.action = 'na.fail', family="binomial")

# model summary
summary(glm_tp2)

# test dispersion of final model
testDispersion(glm_tp2)
simulationOutput <- simulateResiduals(fittedModel = glm_tp2)
plot(simulationOutput)

### timepoint 3 ###

glm_tp3_max <- glmer(crispr ~ phage + glucose + log_cell + log_phage + (1|treatment), 
                     data=tp3_df, na.action = 'na.fail', family="binomial")

# model summary
summary(glm_tp3_max)

# model selection
dredge3 <-dredge(glm_tp3_max)

# final model = simplest within 2 delta AICs

glm_tp3 <- glmer(crispr ~ phage + log_cell + log_phage + (1|treatment), 
                 data=tp3_df, na.action = 'na.fail', family="binomial")

summary(glm_tp3)

# test dispersion of final model
testDispersion(glm_tp3)
simulationOutput <- simulateResiduals(fittedModel = glm_tp3)
plot(simulationOutput)

######### prediction plots #########

### timepoint 1 ###

pred_tp1 <- ggpredict(glm_tp1, terms = c("log_phage", "phage [all]"), 
                                type='fixed', ci.lvl = 0.95) 


plot_tp1 <- ggplot(pred_tp1, aes(x = x, y = predicted, colour=group, fill=group)) +
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.5, colour = NA) +
  labs(x = expression("log"[10]*" phage density"), y = "CRISPR-Cas prob.", colour=expression(atop("Initial phage dose"," (pfu ml"^-1*")")),
       fill=expression(atop("Initial phage dose"," (pfu ml"^-1*")"))) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 6)) +
  scale_fill_manual(values = wes_palette("Cavalcanti1")) +
  scale_colour_manual(values = wes_palette("Cavalcanti1"))

### timepoint 2 ###

pred_tp2 <- ggpredict(glm_tp2, terms = c("log_phage", "log_cell"), 
                      type='fixed', ci.lvl = 0.95) 

plot_tp2 <- ggplot(pred_tp2, aes(x = x, y = predicted, colour=group, fill=group)) +
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),  alpha=0.5, colour = NA) +
  labs(x = expression("log"[10]*" phage density"), y = "CRISPR-Cas prob.", colour=expression("log"[10]*" cell density"),
       fill=expression("log"[10]*" cell density")) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 6)) +
  scale_fill_manual(values = wes_palette("Cavalcanti1")) +
  scale_colour_manual(values = wes_palette("Cavalcanti1"))

### timepoint 3 ###

pred_tp3 <- ggpredict(glm_tp3, terms = c("log_phage", "log_cell", "phage"), 
                      type='fixed', ci.lvl = 0.95) 

plot_tp3 <- ggplot(pred_tp3, aes(x = x, y = predicted, colour=group, fill=group)) +
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),  alpha=0.5, colour = NA) +
  labs(x = expression("log"[10]*" phage density"), y = "CRISPR-Cas prob.", colour=expression("log"[10]*" cell density"),
       fill=expression("log"[10]*" cell density")) +
  facet_grid(~facet) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 6)) +
  scale_fill_manual(values = wes_palette("Cavalcanti1")) +
  scale_colour_manual(values = wes_palette("Cavalcanti1"))


max_pred <- ggarrange(plot_tp1, plot_tp2, plot_tp3,
                            ncol=1,
                            nrow=3)

ggsave(plot=max_pred, "all_timepoints_selected_models.svg", width = 20, height = 30, units="cm")

## AIC tables ##
write.csv(dredge1, "AIC_tp1.csv")
write.csv(dredge2, "AIC_tp2.csv")
write.csv(dredge3, "AIC_tp3.csv")

## model effects tables ##

tp1_summary <- summary(glm_tp1)
tp1_table <- tp1_summary$coefficients
write.csv(tp1_table, file = "model1_tp1_effects.csv")

tp2_summary <- summary(glm_tp2)
tp2_table <- tp2_summary$coefficients
write.csv(tp2_table, file = "model1_tp2_effects.csv")

tp3_summary <- summary(glm_tp3)
tp3_table <- tp3_summary$coefficients
write.csv(tp3_table, file = "model1_tp3_effects.csv")

