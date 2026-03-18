
library(tidyverse)
library(printr)
library(vegan)
library(broom)
library(gglmannotate)
comm=read.csv("vOTU_lytic.csv",header = T,row.names = 1)
distances <- vegdist(comm, method = "canberra", upper = TRUE) 
canberra <- distances %>%
  tidy() %>%
  dplyr::rename(query_sample = item1, reference_sample = item2) %>%
  mutate_at(c("query_sample", "reference_sample"), as.character)
write_tsv(canberra, "canberra_distances_lytic_vOTU.tsv")

metadata <- read.csv("../metadata1.csv")
references <- metadata %>%
  filter(type == "skin") %>%
  select(-type) %>%
  add_count(site, timepoint, name = "pool_size") %>%
  filter(pool_size >= 4) %>%
  select(-pool_size) %>%
  rename_all(~ str_c("reference_", .))

matches <- canberra %>%
  left_join(metadata %>% rename_with(~ str_c("query_", .))) %>%
  filter(query_type %in% c("household surface", "public surface")) %>%
  right_join(references) %>%                         # keep only ref combos
  group_by(query_sample, reference_site, reference_timepoint) %>%
  slice_min(distance, n = 1, with_ties = FALSE) %>%  # single winner
  ungroup() %>%
  mutate(is_accurate = query_location == reference_location)

write_tsv(matches, "watanabe_matches_lytic_vOTU.tsv")

plot_data=read.csv("watanabe_plot_data_summary.csv",header = T)

plot_data$query_site=factor(plot_data$query_site, levels = c("Door knob","Bed headboard","Park/campus handrail", "Subway exit handrail"))

ggplot(plot_data, aes(x = delay*12, y = accuracy,color=organism)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "loess") +scale_color_manual(values = c("#5DA143","#F25C5E","#1B4275","#B31455"))+
  facet_wrap(~ query_site, ncol = 4) +
  labs(
    title = "Watanabe matching accuracy by sampling delay",
    subtitle = "Broken down by query site",
    caption = "Sampling delay is query time point minus reference time point",
    x = "Sampling delay",
    y = "Matching accuracy"
  )+theme_bw()




# overall accuracy rate:
matches %>%
  dplyr::summarise(accuracy = sum(is_accurate) / dplyr::n())

matches %>%
  group_by(query_type) %>%
  dplyr::summarise(accuracy = sum(is_accurate) / dplyr::n())

matches %>%
  group_by(query_site) %>%
  dplyr::summarise(accuracy = sum(is_accurate) / dplyr::n())

matches %>%
  group_by(reference_time, query_type, query_time) %>%
  dplyr::summarise(accuracy = 100 * (sum(is_accurate) / dplyr::n()),.groups = "drop") %>%  # <-- FIX
  ungroup() %>%
  mutate(accuracy = signif(accuracy, 2)) %>%
  mutate_at(c("query_type", "reference_time", "query_time"), str_to_sentence) %>%
  mutate(query_type = fct_relabel(query_type, ~ str_wrap(.x, width = 8))) %>%
  mutate_at(c("reference_time"), ~ glue("{.} reference")) %>%
  mutate_at(c("query_time"), ~ glue("{.} query")) %>%
  mutate_at(c("query_time", "reference_time"), ~ str_wrap(.x, width = 8)) %>%
  mutate_at(c("query_time", "reference_time"), fct_rev) %>%
  ggplot(aes(x = query_type, y = accuracy, label = glue("{accuracy}%"))) +
  geom_col(position = "dodge") +
  facet_grid(reference_time ~ query_time) +
  geom_bar_text(position = "dodge", min.size = 0) +
  labs(
    x = "Query site",
    y = "Accuracy (%)"
  ) +
  theme_classic()

ggsave("matching_diurnality_lytic_vOTU_type.pdf", width = 3.5, height = 3, dpi = 300)

matches %>%
  filter(query_type == "household surface") %>%
  dplyr::count(query_time, is_accurate) %>%
  spread(is_accurate, n) %>%
  as.data.frame() %>%
  column_to_rownames("query_time") %>%
  as.matrix() %>%
  chisq.test() %>%
  tidy()

matches %>%
  filter(query_type == "household surface") %>%
  dplyr::count(query_time, is_accurate) %>%
  spread(is_accurate, n) %>%
  mutate(accuracy = 100 * `TRUE` / (`TRUE` + `FALSE`))


matches %>%
  filter(query_type == "public surface") %>%
  dplyr::count(query_time, is_accurate) %>%
  spread(is_accurate, n) %>%
  as.data.frame() %>%
  column_to_rownames("query_time") %>%
  as.matrix() %>%
  chisq.test() %>%
  tidy()

matches %>%
  filter(query_type == "public surface") %>%
  dplyr::count(query_time, is_accurate) %>%
  spread(is_accurate, n) %>%
  mutate(accuracy = 100 * `TRUE` / (`TRUE` + `FALSE`))


matches %>%
  filter(query_type == "household surface") %>%
  dplyr::count(reference_time, is_accurate) %>%
  spread(is_accurate, n) %>%
  as.data.frame() %>%
  column_to_rownames("reference_time") %>%
  as.matrix() %>%
  chisq.test() %>%
  tidy()

matches %>%
  filter(query_type == "household surface") %>%
  dplyr::count(reference_time, is_accurate) %>%
  spread(is_accurate, n) %>%
  mutate(accuracy = 100 * `TRUE` / (`TRUE` + `FALSE`))

matches %>%
  filter(query_type == "public surface") %>%
  dplyr::count(reference_time, is_accurate) %>%
  spread(is_accurate, n) %>%
  as.data.frame() %>%
  column_to_rownames("reference_time") %>%
  as.matrix() %>%
  chisq.test() %>%
  tidy()

matches %>%
  filter(query_type == "public surface") %>%
  dplyr::count(reference_time, is_accurate) %>%
  spread(is_accurate, n) %>%
  mutate(accuracy = 100 * `TRUE` / (`TRUE` + `FALSE`))

#################
bOTU <- read_tsv("../watanabe_matches_rMAG_new.tsv") %>%
  mutate(method = "bOTU") %>%
  select(method, query_sample, query_site, query_type, query_location,
         query_day, query_time, query_timepoint, reference_site,
         reference_location, reference_day, reference_time, reference_timepoint,
         is_accurate)

vOTU <- read_tsv("../watanabe_matches_vOTU.tsv") %>%
  mutate(method = "vOTU") %>%
  select(method, query_sample, query_site, query_type, query_location,
         query_day, query_time, query_timepoint, reference_site,
         reference_location, reference_day, reference_time, reference_timepoint,
         is_accurate)
lysogenic <- read_tsv("watanabe_matches_lysogenic_vOTU.tsv.tsv") %>%
  mutate(method = "lysogenic") %>%
  select(method, query_sample, query_site, query_type, query_location,
         query_day, query_time, query_timepoint, reference_site,
         reference_location, reference_day, reference_time, reference_timepoint,
         is_accurate)
lytic <- read_tsv("watanabe_matches_lytic_vOTU.tsv") %>%
  mutate(method = "lytic") %>%
  select(method, query_sample, query_site, query_type, query_location,
         query_day, query_time, query_timepoint, reference_site,
         reference_location, reference_day, reference_time, reference_timepoint,
         is_accurate)
matching <- bind_rows(vOTU,bOTU, lysogenic, lytic) %>%
  mutate(query_hours = (24 * (query_day - 1)) + if_else(query_time == "morning", 0, 12)) %>%
  mutate(reference_hours = (24 * (reference_day - 1)) + if_else(reference_time == "morning", 0, 12)) %>%
  mutate(delay = query_hours - reference_hours)

plot_data<-matching %>%
  group_by(method, query_site, delay) %>%
  summarise(accuracy = sum(is_accurate) / n()) %>%
  ungroup() %>%
  mutate_at("accuracy", ~ . * 100) %>%
  mutate_at("query_site", str_to_sentence)
  
plot_data$query_site=factor(plot_data$query_site, levels = c("Door knob","Bed headboard","Park/campus handrail", "Subway exit handrail"))
  ggplot(plot_data,aes(x = delay, y = accuracy, colour = method)) +
  geom_point() +scale_color_manual(values = c("#1B4275","#5DA143","#F25C5E","#B31455"))+
  geom_smooth(method = "loess") +
  facet_wrap(~ query_site, nrow = 1) +
  labs(x = "Sampling delay (hours)",
       y = "Accuracy (%)",
       colour = "Method") +
  theme_classic()




