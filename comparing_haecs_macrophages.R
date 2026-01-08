######## comparing LPS treated macrophages to HAECs with IL1B
library(tidyverse)
library(dplyr)
library(ggplot2)
set.seed(12347)

## read in macropahge results:
macrophages <- read.delim("/Volumes/data3/anna/splicing_analysis_2024/GSE147310/human_monocyte_derived_macrophages_leafcutter_summary.txt", as.is = T, stringsAsFactors = F)

# read in HAEC result
haecs <- read.delim("/Volumes/data3/anna/splicing_analysis_2024/20240321_annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth.txt", as.is = T, stringsAsFactors = F)

## join dataframes
macrophages %>%
  mutate(chrjunction = unlist(vapply(strsplit(macrophages$intron, ":"), function(x) paste(x[c(1,2,3)], collapse = "-"), character(1L)))) -> macrophages
macrophages%>% 
  select(intron, chrjunction, genes, deltapsi, p.adjust, spliceType) %>%
  full_join(., haecs[,c("intron", "chrjunction", "deltapsi", "p.adjust", "spliceType")], by = "chrjunction", suffix = c("_macro", "_haec")) -> together

## summarise data
together %>%
  mutate(macrosig = ifelse(!is.na(p.adjust_macro), ifelse(p.adjust_macro < 0.05 & abs(deltapsi_macro) > 0.05, "sig", "nonsig"), "nonsig"),
         haecsig = ifelse(!is.na(p.adjust_haec), ifelse(p.adjust_haec < 0.05 & abs(deltapsi_haec) > 0.05, "sig", "nonsig"), "nonsig"), 
         category = ifelse(macrosig == "sig" & haecsig == "sig", "both significant", 
                           ifelse(haecsig == "sig" & !macrosig == "sig", "haec only", 
                                  ifelse(!haecsig == "sig" & macrosig == "sig", "macro only", "not sig"))),
         spliceType = ifelse(category == "both significant", spliceType_haec, 
                             ifelse(category == "haec only", spliceType_haec, spliceType_macro)),
         spliceType = factor(ifelse(spliceType %in% c("cryptic_fiveprime", "cryptic_threeprime", " cryptic_unanchored", "novel annotated pair", "unknown_strand", NA), "cryptic", spliceType), levels = c("A5", "A3", "AF", "AL", 'MX', "RI", "SE", "cryptic"))
         ) -> together

palette <- c("#76E1B2", "#E56E63", "#BB46E0", "#DCC7D8", "#D0DEAF", "#DDAB69", "#83A1DA", "#96DADD", "#9272D5", "#85E359", "#DE77B6", "#7E7D73", "#D8E068")

ggplot(together, aes(x = ifelse(is.na(deltapsi_haec), 0, deltapsi_haec), y = ifelse(is.na(deltapsi_macro), 0, deltapsi_macro), color = spliceType)) + geom_point() + theme_bw() + theme(text = element_text(size = 12))  + xlab("HAEC control vs IL1B deltaPSI") + ylab("Macrophage control vs LPS deltaPSI") + ggtitle("comparing inflammation-induced differential \nsplicing in endothelial and macrophage cells") + scale_color_manual(values = palette)

together %>% select(!spliceType_macro) %>% select(!spliceType_haec)  %>%
  pivot_longer(cols = !c(chrjunction, intron_haec, intron_macro, spliceType, genes, macrosig, haecsig, category), 
  names_to = c('.value', 'grp'), 
  names_sep = "_") -> together_long

ggplot(together_long[!together_long$spliceType == "cryptic" & !is.na(together_long$spliceType) & !together_long$category == "not sig",], aes(y = spliceType, fill = category)) + geom_bar(stat = "count", position = "fill") + scale_fill_manual(values = c("#76E1B2","#83A1DA","#DE77B6")) + theme_bw() + xlab("Percent of DSGs") + theme(text = element_text(size = 12)) + ggtitle("Comparing DSGs between \nmacrophages and HAECs")

library(ggvenn)
df <- list(
  HAECs = haecs$chrjunction[haecs$p.adjust < 0.05 & abs(haecs$deltapsi) > 0.05],
  Macrophages = macrophages$chrjunction[macrophages$p.adjust < 0.05 & abs(macrophages$deltapsi) > 0.05]
)

ggvenn(df, fill_color = c("#DE77B6","#83A1DA"),
       stroke_size = 0.5, set_name_size = 4)
