##
## Filter type (10L, Whole water), sterivex vs membrane
# Prokaryotes
# OTUs
prokaryotes_membrane_vs_sterivex_10L_whole_water <-
  prok_diversity_metadata %>% 
  filter(method=="filter",effected_volume == 10) %>% 
  mutate(device = ifelse(device == "membrane", "Membrane", "Sterivex"))

prok_wilc_dev_10L_ww <-
  prokaryotes_membrane_vs_sterivex_10L_whole_water %>% 
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ device) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x = "device") %>% 
  mutate(Test = "Mann-Whitney",
         Taxonomic_group = "Prokaryotic")

## Sterivex vs membrane, protists
# OTUS
protist_membrane_vs_sterivex_10L_whole_water <-
  prot_diversity_metadata %>% 
  filter(method == "filter", effected_volume == 10) %>% 
  mutate(device = ifelse(device == "membrane", "Membrane", "Sterivex"))

prot_wilc_m_vs_st_10L_ww <-
  protist_membrane_vs_sterivex_10L_whole_water %>% 
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ device) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x = "device")%>% 
  mutate(Test = "Mann-Whitney",
         Taxonomic_group = "Protist") 

# Plot (Sterivex vs membrane, 10L, WW, protist)

## Whole water vs size fractions, 10L
# OTUs
prokaryotes_ww_vs_sf_10L_membrane <-
  prok_diversity_metadata %>% 
  filter(effected_volume == 10,device == "membrane")

# Kruskal test
prokaryotes_kruskal_test_ww_vs_sf_10L_membrane <-
  prokaryotes_ww_vs_sf_10L_membrane %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  kruskal_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance()%>% 
  mutate(Sequencing_strategy = "MetaB16SV4V5",
         Test = "Kruskal-Wallis",
         Taxonomic_group = "Prokaryotes")

prokaryotes_ww_vs_sf_10L_membrane %>% 
  group_by(Sequencing_strategy) %>% 
  kruskal_effsize(Species_richness ~ size_fraction)

# post-hoc test (Dunn)
prokaryotes_dunn_ww_vs_sf_10L_membrane <- 
  prokaryotes_ww_vs_sf_10L_membrane %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  dunn_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x="size_fraction") %>% 
  mutate(Sequencing_strategy = "MetaB16SV4V5",
         Test = "Dunn",
         Taxonomic_group = "Prokaryotes")


## Wilcoxon test for metagenomes
prokaryotes_wilcox_test_ww_vs_sf_10L_membrane <-
  prokaryotes_ww_vs_sf_10L_membrane %>% 
  filter(Sequencing_strategy != "MetaB16SV4V5") %>% 
  wilcox_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(Sequencing_strategy = "Metagenome",
         Test = "Mann-Whitney",
         Taxonomic_group = "Prokaryotes")

# Protists, whole water vs size fractions, 10L
protist_ww_sf_10L_membrane <- 
  prot_diversity_metadata %>% 
  filter(device == "membrane",effected_volume == 10)

# Kruskall wallis test, protist, ww vs sf, 10L
protist_ww_sf_10L_membrane_kruskall_test <- 
  protist_ww_sf_10L_membrane %>% 
  group_by(Sequencing_strategy) %>%
  kruskal_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(Test = "Kruskal-Wallis",
         Taxonomic_group = "Protist")

# Dunn test not necessary, because Kruskall wallis was not significant.  

# quick overview, protists, 10L, ww vs sf
protist_ww_sf_10L_membrane %>% 
  group_by(Sequencing_strategy,size_fraction) %>% 
  summarise(min(Species_richness),
            max(Species_richness))

## Size fractions, 100L, membrane, prokaryotes
prokaryotes_size_fractions_100L_membrane <- 
  prok_diversity_metadata %>% 
  filter(device == "membrane",effected_volume == 100)

## 16S (3 groups to compare)
# Kruskal test, size fractions, membrane prokaryotes, (16S only)
prokaryotes_size_fractions_100L_membrane_kruskal_test_16S <-
  prokaryotes_size_fractions_100L_membrane %>%
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  kruskal_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(Sequencing_strategy = "MetaB16SV4V5",
         Test = "Kruskal-Wallis",
         Taxonomic_group = "Prokaryotic")

# post-hoc test (Dunn) size fractions 100L prokaryotes, (16S only)
prokaryotes_size_fractions_100L_membrane_dunn_test_16S <- 
  prokaryotes_size_fractions_100L_membrane %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  dunn_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x="size_fraction") %>%
  mutate(Sequencing_strategy = "MetaB16SV4V5",
         Test = "Dunn",
         Taxonomic_group = "Prokaryotic")

prokaryotes_size_fractions_100L_membrane_dunn_test_16S$xmin <- c(1,1,2)
prokaryotes_size_fractions_100L_membrane_dunn_test_16S$xmax <- c(2,3,3)

# Size fractions, 100L, membrane, 16S only

## Some quick metrics
prokaryotes_size_fractions_100L_membrane %>% 
  group_by(Sequencing_strategy,size_fraction) %>% 
  summarize(median(Species_richness))

## metagenomes (2 groups to compare)
# Size fractions, 100L, membrane, Metagenomes
prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes <-
  prokaryotes_size_fractions_100L_membrane %>% 
  filter(Sequencing_strategy == "MetaG") %>%
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x = "size_fraction") %>% 
  mutate(Sequencing_strategy = "Metagenome",
         Test = "Mann-Whitney",
         Taxonomic_group = "Prokaryotic")

prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes$xmin <- 1
prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes$xmax <- 2

## Size fractions, 100L, membrane, protists
protist_size_fractions_100L_membrane <- 
  prot_diversity_metadata %>% 
  filter(device == "membrane",effected_volume == 100)

## Kruskall wallis test
protist_size_fractions_100L_membrane_kruskall <- 
  protist_size_fractions_100L_membrane %>% 
  group_by(Sequencing_strategy) %>% 
  kruskal_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(test = "Kruskal-Wallis",
         Taxonomic_group = "Protist")

# post-hoc test for 18S (Dune test, protist, size fractions, 100L, membrane)
protist_size_fractions_100L_membrane_dunn <-
  protist_size_fractions_100L_membrane %>% 
  group_by(Sequencing_strategy) %>% 
  dunn_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(group = "size_fraction") %>% 
  mutate(test = "Dunn",
         Taxonomic_group = "Protist")

## 2.5L vs 10L, WW, sterivex
# Prokaryotes
prok_2.5L_vs_10L_wilcoxon <- 
  prok_diversity_metadata %>% 
  filter(device == "sterivex", planned_volume != "1L") %>%
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ effected_volume) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  mutate(Test = "Mann-Whitney",
         Taxonomic_group = "Prokaryotic")

#quick overview
prok_diversity_metadata %>% 
  filter(device == "sterivex", planned_volume != "1L") %>% 
  group_by(Sequencing_strategy,effected_volume) %>% 
  summarise(min(Species_richness),
            max(Species_richness),
            median(Species_richness),
            IQR(Species_richness),
            n())

# Protist 2.5L vs 4*2.5L, sterivex, WW
protist_2.5_vs_10_st_ww <- 
  prot_diversity_metadata %>% 
  filter(device == "sterivex", planned_volume != "1L")

# Wilcoxon test
prot_2.5_vs_10_st_ww_wilcox <- 
  protist_2.5_vs_10_st_ww %>% 
  group_by(Sequencing_strategy) %>%
  wilcox_test(Species_richness ~ planned_volume) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  mutate(Test = "Mann-Whitney",
         Taxonomic_group = "Protist")



