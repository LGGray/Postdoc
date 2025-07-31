library(readxl)
library(dplyr)
library(tidyverse)

embryo <- read_excel('~/Downloads/43587_2025_856_MOESM2_ESM.xlsx', sheet = 7)
young <- read_excel('~/Downloads/43587_2025_856_MOESM2_ESM.xlsx', sheet = 9)
adult <- read_excel('~/Downloads/43587_2025_856_MOESM2_ESM.xlsx', sheet = 3)
aged <- read_excel('~/Downloads/43587_2025_856_MOESM2_ESM.xlsx', sheet = 11)

age_list <- list(embryo, young, adult, aged)

Xist_list <- list()
for(df in 1:length(age_list)) {
    Xist_TPM <- subset(age_list[[df]], name == 'Xist' & sex != 'male')
    organ <- split(Xist_TPM, Xist_TPM$organ)
    Xist_list[[df]] <- organ
}
names(Xist_list) <- c('embryo', 'young', 'adult', 'aged')

# Extract TPM values for each organ and age group
brain <- lapply(Xist_list, function(x){
    x[['Br']]$TPM
})
brain <- bind_rows(brain) %>% 
    mutate(organ = 'brain')

heart <- lapply(Xist_list, function(x){
    x[['He']]$TPM
})
heart <- bind_rows(heart) %>% 
    mutate(organ = 'heart')

kidney <- lapply(Xist_list, function(x){
    x[['Ki']]$TPM
})
kidney <- bind_rows(kidney) %>% 
    mutate(organ = 'kidney')

liver <- lapply(Xist_list, function(x){
    x[['Li']]$TPM
})
liver <- bind_rows(liver) %>% 
    mutate(organ = 'liver')

lung <- lapply(Xist_list, function(x){
    x[['Lu']]$TPM
})
lung <- bind_rows(lung) %>% 
    mutate(organ = 'lung')

spleen <- lapply(Xist_list, function(x){
    x[['Sp']]$TPM
})
spleen <- bind_rows(spleen) %>% 
    mutate(organ = 'spleen')
spleen$embryo <- 0

muscle <- lapply(Xist_list, function(x){
    x[['Mu']]$TPM
})
muscle <- bind_rows(muscle) %>% 
    mutate(organ = 'muscle')
muscle$embryo <- 0

plot_data <- rbind(brain, heart, kidney, liver, lung, spleen, muscle)

# Reshape data
long_data <- plot_data %>%
  pivot_longer(cols = c(embryo, young, adult, aged), 
               names_to = "age", 
               values_to = "value")

# Remove TPM = 0
long_data <- subset(long_data, value > 0)

# Set factors
long_data$age <- factor(long_data$age, levels = c("embryo", "young", "adult", "aged"))
long_data$organ <- factor(long_data$organ, levels = c("brain", "heart", "kidney", "liver", "lung", "spleen", "muscle"))

organ_colors <- c("#DCB465", "#8B1812","#244C51", "#C97A41", "#A77A76", "#555463", "#97A092")

pdf('../Xist_TPM.pdf')
ggplot(long_data, aes(x = age, y = value, fill = organ)) +
  geom_boxplot() +
  scale_fill_manual(values = organ_colors) +
  # show points
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ organ, scales = "free_y") +
  labs(title = "Xist expression across development", x = "Age", y = "Xist TPM") +
  theme_minimal() +
  theme(legend.position = 'none')
dev.off()

brain_model <- lm(log2(value + 1) ~ age, data = subset(long_data, organ == "brain"))
summary(brain_model)
heart_model <- lm(log2(value + 1) ~ age, data = subset(long_data, organ == "heart"))
summary(heart_model)
kidney_model <- lm(log2(value + 1) ~ age, data = subset(long_data, organ == "kidney"))
summary(kidney_model)
liver_model <- lm(log2(value + 1) ~ age, data = subset(long_data, organ == "liver"))
summary(liver_model)
lung_model <- lm(log2(value + 1) ~ age, data = subset(long_data, organ == "lung"))
summary(lung_model)
spleen_model <- lm(log2(value + 1) ~ age, data = subset(long_data, organ == "spleen"))
summary(spleen_model)
muscle_model <- lm(log2(value + 1) ~ age, data = subset(long_data, organ == "muscle"))
summary(muscle_model)


#### cardiac cell types ####
cardiac_celltypes <- read_excel('~/Downloads/43587_2025_856_MOESM2_ESM.xlsx', sheet = 5)
Xist_TPM_cardiac <- subset(cardiac_celltypes, name == 'Xist' & sex == 'female')

pdf('../Xist_TPM_cardiac.pdf')
ggplot(Xist_TPM_cardiac, aes(x = celltype, y = TPM, fill = celltype)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Xist expression in cardiac cell types", x = "Cell Type", y = "Xist TPM") +
  theme_minimal() +
  theme(legend.position = 'none')
dev.off()

cardiac_model <- lm(log2(TPM + 1) ~ celltype, data = Xist_TPM_cardiac)
summary(cardiac_model)
cardiac_model_anova <- anova(cardiac_model)
print(cardiac_model_anova)