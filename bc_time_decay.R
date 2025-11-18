library(vegan)
comm=read.csv("residential_342_rMAG_normed_coverage.csv",header = T,row.names = 1)
distances <- vegdist(comm, method = "bray", upper = TRUE) 
bc <- distances %>%
       tidy() %>%
       rename(query_sample = item1, reference_sample = item2) %>%
       mutate_at(c("query_sample", "reference_sample"), as.character)
write_tsv(bc, "bc_distances_rMAG.tsv")

matches <- bc %>%
       left_join(
             metadata %>% rename_all(~ str_c("query_", .)))
matches1 <- matches %>%
       left_join(
             metadata %>% rename_all(~ str_c("reference_", .)))
matches2 <- matches1 %>%
       mutate(delay = query_timepoint - reference_timepoint)

matches_filtered <- matches2 %>%
       filter(query_location == reference_location &
                                 query_site == reference_site)
matches_filtered$reference_site=factor(matches_filtered$reference_site,levels=c("left palm","right palm","door knob","bed headboard","park/campus handrail","subway exit handrail"))
ggplot(matches_filtered,aes(x=delay*24,y=1-distance))+geom_point()+facet_wrap(~reference_site)+geom_smooth()+theme_bw()
