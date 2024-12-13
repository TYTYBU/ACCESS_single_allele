setwd("/Users/tianyu/Downloads")
library(tidyverse)
library(tidyquant)
library(patchwork)
library(ggbeeswarm)
library(ggpubr)
library(ggseqlogo)

ls.plt = list()

### occupancy schematic plots
# HepG2 USF1
schematicDir = '/Users/tianyu/Downloads/ACCESS_proj/schemetics_data/'
dir.create(paste0(schematicDir, "/plots/"), recursive = T, showWarnings = F)

# motif feature plot
df.edit_frac_unlabeled = read_csv(paste0(schematicDir, 'HepG2_USF1_ENCSR000BGM_GTCACGTGRCS.edit_frac_features.csv'))
df.edit_frac_unlabeled = df.edit_frac_unlabeled %>% mutate(edit_frac_ma = zoo::rollmean(edit_frac, k=3, align = "center", na.pad = T))

df.motif_features = read_csv(paste0(schematicDir, 'HepG2_USF1_ENCSR000BGM_GTCACGTGRCS.motif_features.csv'))
r_peak_pos_list = as.numeric(unlist(str_split(gsub(pattern = "\\[|\\]|,", replacement = "", df.motif_features$r_peak_pos_list), pattern=" ")))
l_peak_pos_list = as.numeric(unlist(str_split(gsub(pattern = "\\[|\\]|,", replacement = "", df.motif_features$l_peak_pos_list), pattern=" ")))
center_pos_list = as.numeric(unlist(str_split(gsub(pattern = "\\[|\\]|,", replacement = "", df.motif_features$center_pos_list), pattern=" ")))
left_peak = df.motif_features$left_peak
footprint = df.motif_features$footprint
right_peak = df.motif_features$right_peak
df.motif_features_expanded = data.frame(relative_pos = c(l_peak_pos_list, center_pos_list, r_peak_pos_list),
                                        feature = c(rep("left_peak", length(l_peak_pos_list)), rep("footprint", length(center_pos_list)), rep("right_peak", length(r_peak_pos_list))),
                                        mean_editFrac = c(rep(left_peak, length(l_peak_pos_list)), rep(footprint, length(center_pos_list)), rep(right_peak, length(r_peak_pos_list))))

plt = ggplot() + 
  geom_rect(data=df.motif_features_expanded, aes(xmin=relative_pos-0.5, xmax=relative_pos+0.5, fill=feature), ymin=-Inf, ymax=Inf, alpha=0.5) + 
  geom_segment(data=df.motif_features_expanded, aes(x=relative_pos-0.5, xend=relative_pos+0.5, y=mean_editFrac, yend=mean_editFrac, color=feature), linewidth=1) + 
  geom_line(data=df.edit_frac_unlabeled, aes(x=relative_pos, y=edit_frac_ma), color='black', linewidth=0.5) + 
  geom_line(data=df.edit_frac_unlabeled, aes(x=relative_pos, y=edit_frac), color='black', linewidth=0.25) + 
  scale_color_manual(values=c('#d8b365', '#5ab4ac', '#5ab4ac')) +
  scale_fill_manual(values=c('#d8b365', '#5ab4ac', '#5ab4ac')) +
  theme_void() + theme(legend.position = "none")
ggsave(filename = paste0(schematicDir, "/plots/HepG2_USF1_motif_features_simplified.pdf"), plot = plt, height = 2, width = 4)

plt = ggplot() + 
  geom_rect(data=df.motif_features_expanded, aes(xmin=relative_pos-0.5, xmax=relative_pos+0.5, fill=feature), ymin=-Inf, ymax=Inf, alpha=0.5) + 
  geom_segment(data=df.motif_features_expanded, aes(x=relative_pos-0.5, xend=relative_pos+0.5, y=mean_editFrac, yend=mean_editFrac, color=feature), linewidth=1) + 
  geom_line(data=df.edit_frac_unlabeled, aes(x=relative_pos, y=edit_frac_ma), color='black', linewidth=0.5) + 
  geom_line(data=df.edit_frac_unlabeled, aes(x=relative_pos, y=edit_frac), color='black', linewidth=0.25) + 
  scale_color_manual(values=c('#d8b365', '#5ab4ac', '#5ab4ac')) +
  scale_fill_manual(values=c('#d8b365', '#5ab4ac', '#5ab4ac')) +
  xlim(-80, 80) + 
  theme_void() + theme(legend.position = "none")
ggsave(filename = paste0(schematicDir, "/plots/HepG2_USF1_motif_features_simplified_cropped.pdf"), plot = plt, height = 2, width = 4)
ls.plt[['schematic_feature']] = plt

# pre-defined labels
df.edit_frac_labeled = read_csv(paste0(schematicDir, 'HepG2_USF1_ENCSR000BGM_GTCACGTGRCS.predefined_labels_editFrac.csv'))
df.edit_frac_labeled = df.edit_frac_labeled %>% mutate(edit_frac_ma = zoo::rollmean(edit_frac, k=3, align = "center", na.pad = T))

plt = ggplot() + 
  # geom_rect(data=df.motif_features_expanded, aes(xmin=relative_pos-0.5, xmax=relative_pos+0.5, fill=feature), ymin=-Inf, ymax=Inf, alpha=0.5) + 
  geom_line(data=df.edit_frac_labeled, aes(x=relative_pos, y=edit_frac_ma, color=read_type), linewidth=0.5) + 
  geom_line(data=df.edit_frac_labeled, aes(x=relative_pos, y=edit_frac, color=read_type), linewidth=0.25) + 
  # scale_color_manual(values=c('red', 'green', 'blue')) +
  # scale_fill_manual(values=c('#e78ac3', '#a6d854', '#ffd92f')) +
  theme_void() + theme(legend.position = "none")
# ggsave(filename = paste0(schematicDir, "/plots/HepG2_USF1_prelabeled_simplified.pdf"), plot = plt, height = 2, width = 4)
ls.plt[['schematic_predefined']] = plt



### bound probability - ChIPseq score correlation overview
df.stats = read_csv('/Users/tianyu/Downloads/TF_summary.features_filtered.csv')
eval_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/evaluation_3stateModel_v2/'

# deep-learning rho violin plot
ind = which(df.stats$cell_type == "HepG2" & df.stats$TF == "USF1")
ind2 = which(df.stats$cell_type == "K562" & df.stats$TF == "STAT5A")
df.stats$shape_size = 0.2
df.stats$shape_size[c(ind, ind2)] = 0.8
df.stats2 = df.stats[c(ind, ind2),]

plt = ggplot() + 
  geom_violin(data=df.stats, aes(x=model_rho, y=cell_type, color=cell_type, fill=cell_type), width = 1, bw=0.02, alpha=0.5) + 
  # geom_quasirandom(orientation='y', size=0.5, varwidth = T) +
  geom_beeswarm(data=df.stats, aes(x=model_rho, y=cell_type, color=cell_type, size=shape_size), cex=2.5) +
  geom_text(data=df.stats2, aes(x = model_rho, y=cell_type, label = TF, color=cell_type), hjust = -0.3, vjust=0.5, size=8/3.5) + 
  scale_size_identity() + 
  scale_color_manual(values = c("blue", "red")) + 
  scale_fill_manual(values = c("blue", "red")) + 
  labs(x='Spearman rho of \nOccuPIE mean allelic bound probability vs. ChIP-seq signal', y='Cell type') + 
  xlim(-0.1, 0.8) + 
  theme_classic() +
  theme(legend.position = "none")
# ggsave(filename = paste0(eval_statsDir, "plots/deepLearn_rho_violin.pdf"), plot = plt, height = 3, width = 6)
ls.plt[['cor_violin']] = plt

### bound probability - ChIPseq score correlation
### composite edit fraction for three labels
eval_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/evaluation_3stateModel_v2/'
df.selected_TFs = read_csv(paste0(eval_statsDir, 'selected_TF_list.csv'))
suffix_str = '.edit_frac.csv'

fnames = list.files(eval_statsDir, pattern = suffix_str)
temp = gsub(pattern = suffix_str, replacement = '', fnames)
df.eval_files = str_split(temp, pattern = '_', simplify = T, n = 3)
colnames(df.eval_files) = c('cell_type', 'TF_name', 'motif_str')
df.eval_files = data.frame(df.eval_files)
df.eval_files$edit_frac_fnames = fnames

suffix_str2 = '.read_type_frac.csv'
fnames2 = gsub(pattern = suffix_str, replacement = suffix_str2, fnames)
df.eval_files$read_tyoe_frac_fnames = fnames2

ind = which(df.selected_TFs$cell_type == "HepG2" & df.selected_TFs$TF_name == "USF1")
ind2 = which(df.selected_TFs$cell_type == "K562" & df.selected_TFs$TF_name == "STAT5A")
for (i in c(ind, ind2)){
  TF_name = df.selected_TFs$TF_name[i]
  cell_type = df.selected_TFs$cell_type[i]
  ind = which(df.eval_files$cell_type == cell_type & df.eval_files$TF_name == TF_name)
  motif_str = df.eval_files$motif_str[ind]
  
  # edit fraction plot
  df.edit_frac = read_csv(paste0(eval_statsDir, df.eval_files$edit_frac_fnames[ind]))
  df.edit_frac = df.edit_frac %>% mutate(edit_frac_ma = zoo::rollmean(edit_frac, k=3, align = "center", na.pad = T))
  df.edit_frac = df.edit_frac %>% filter(label == 'all_channel')
  df.edit_frac$read_type[df.edit_frac$read_type == 'recently bound'] = 'Unbound\nAccessible'
  df.edit_frac$read_type[df.edit_frac$read_type == 'unbound'] = 'Unbound\nInaccessible'
  df.edit_frac$read_type[df.edit_frac$read_type == 'bound'] = 'Bound'
  x_max = max(df.edit_frac$relative_pos)
  y_max = max(df.edit_frac$edit_frac)
  
  plt1 = ggplot() +
    # geom_line(data=df.edit_frac, aes(x=relative_pos, y=edit_frac_ma, color = read_type), linewidth=0.2) + 
    geom_line(data=df.edit_frac, aes(x=relative_pos, y=edit_frac, color = read_type), linewidth=0.2) + 
    # geom_point(size = 0.2, alpha = 0.25) + 
    # geom_ma(ma_fun = SMA, n = 5, linetype = 1, size=0.5) + 
    labs(x = "Relative position", y = "Edit fraction", color = NULL, title = TF_name) + 
    ylim(0, 1) + 
    theme_classic() + 
    theme(
      plot.title = element_text(color = ifelse(TF_name == "USF1", "blue", "red"), hjust=0.5),
      legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1),
      legend.background = element_rect(fill = NA),  # Transparent background
      legend.key.height = unit(0.8, "cm"), legend.margin = margin(t=0)
    )
  ls.plt[[paste0(c(cell_type, TF_name, "editFrac"), collapse = "_")]] = plt1
  
  # correlation plot
  df.read_type_frac = read_csv(paste0(eval_statsDir, df.eval_files$read_tyoe_frac_fnames[ind]))
  df.read_type_frac = df.read_type_frac %>% 
    filter(read_type == 'bound' & chipseq_norm_signalVal > 0) %>%
    mutate(chipseq_norm_signalVal_log10 = log10(chipseq_norm_signalVal))
  res_spearman = cor.test(df.read_type_frac$read_frac_predicted_prob, df.read_type_frac$chipseq_norm_signalVal_log10, method = 'spearman')
  res_pearson = cor.test(df.read_type_frac$read_frac_predicted_prob, df.read_type_frac$chipseq_norm_signalVal_log10, method = 'pearson')
  # cor_str = paste0('rho=', signif(res_spearman$estimate, digits = 2), ', p=', signif(res_spearman$p.value, digits = 2), '\nr=', signif(res_pearson$estimate, digits = 2), ', p=', signif(res_pearson$p.value, digits = 2))
  cor_str = paste0('rho=', signif(res_spearman$estimate, digits = 2), '\np=', signif(res_spearman$p.value, digits = 2))
  x_max = max(df.read_type_frac$read_frac_predicted_prob)
  y_min = min(df.read_type_frac$chipseq_norm_signalVal_log10)
  
  plt2 = ggplot() + 
    geom_point(data=df.read_type_frac, aes(x=read_frac_predicted_prob, y=chipseq_norm_signalVal_log10), size=0.25, alpha=0.25) + 
    labs(x="OccuPIE mean allelic bound probability", y="Log10 norm. ChIP-seq signal", title = TF_name) + 
    geom_label(data=df.read_type_frac, x=x_max, y=y_min, label=cor_str, hjust=1, vjust=0, size=8/3.5, fill = "white", label.size=0, alpha=0.5) + 
    geom_smooth(data=df.read_type_frac, aes(x=read_frac_predicted_prob, y=chipseq_norm_signalVal_log10), method="lm", color="black", linewidth=0.5, se=FALSE) +
    theme_classic() + 
    theme(plot.title = element_text(color = ifelse(TF_name == "USF1", "blue", "red"), hjust=0.5))
  ls.plt[[paste0(c(cell_type, TF_name, "corr"), collapse = "_")]] = plt2
  
  # plt = (plt1 | plt2)
  # ggsave(filename = paste0(eval_statsDir, "plots/", cell_type, "/", cell_type, "_", TF_name, "_", motif_str, ".comboPlt.pdf"), height = 6, width = 3)
}

plt_edit_frac = (ls.plt$K562_STAT5A_editFrac | ls.plt$HepG2_USF1_editFrac) + 
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
  
plt_cor_dotplot = (ls.plt$K562_STAT5A_corr | ls.plt$HepG2_USF1_corr) + 
  plot_layout(axis_titles = "collect")
  
plt_right_column = ls.plt$cor_violin + plt_edit_frac + plt_cor_dotplot +
  plot_layout(heights = c(2,2,2)) &
  theme(
    text = element_text(size = 8), 
    axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(linewidth=0.3),
    legend.text = element_text(size = 8), legend.title = element_text(size = 8),
    plot.title = element_text(size = 8),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5)
  )

ggsave(filename = paste0(eval_statsDir, "plots/eval_comboPlt.pdf"), plot=plt_right_column, height = 6, width = 4, units = "in")




### supplementary figure 6
dir.create(paste0(eval_statsDir, "suplementary_fig6"), showWarnings = F, recursive = T)

# pre-defined rho vs. deep-learning rho 
df.stats = read_csv('/Users/tianyu/Downloads/TF_summary.features_filtered.csv')
eval_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/evaluation_3stateModel_v2/'

xy_max = max(df.stats$rho, df.stats$model_rho)
xy_min = min(df.stats$rho, df.stats$model_rho)
plt = df.stats %>% ggplot(aes(x=rho, y=model_rho, color=cell_type)) + 
  geom_point(size=0.25, alpha=1) + 
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray", linewidth=0.5) + 
  labs(x='Pre-defined spearman rho', y='Deep learning spearman rho', color=NULL) + 
  xlim(xy_min, xy_max) + ylim(xy_min, xy_max) + 
  theme_classic() + 
  theme(
    legend.position = "inside", 
    legend.position.inside = c(1, 0), 
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = NA),  # Transparent background
    legend.key.height = unit(0.5, "cm")
  )
ggsave(filename = paste0(eval_statsDir, "plots/suplementary_fig6/preDefined_vs_deepLearn.jpg"), height = 3, width = 3)

# bound probability - ChIPseq score correlation
# composite edit fraction for three labels
eval_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/evaluation_3stateModel_v2/'
df.selected_TFs = read_csv(paste0(eval_statsDir, 'selected_TF_list.csv'))
suffix_str = '.edit_frac.csv'

fnames = list.files(eval_statsDir, pattern = suffix_str)
temp = gsub(pattern = suffix_str, replacement = '', fnames)
df.eval_files = str_split(temp, pattern = '_', simplify = T, n = 3)
colnames(df.eval_files) = c('cell_type', 'TF_name', 'motif_str')
df.eval_files = data.frame(df.eval_files)
df.eval_files$edit_frac_fnames = fnames

suffix_str2 = '.read_type_frac.csv'
fnames2 = gsub(pattern = suffix_str, replacement = suffix_str2, fnames)
df.eval_files$read_tyoe_frac_fnames = fnames2

ind = which(df.selected_TFs$cell_type == "K562" & df.selected_TFs$TF_name == "MITF")
ind2 = which(df.selected_TFs$cell_type == "HepG2" & df.selected_TFs$TF_name == "CREB1")
for (i in c(ind, ind2)){
  TF_name = df.selected_TFs$TF_name[i]
  cell_type = df.selected_TFs$cell_type[i]
  ind = which(df.eval_files$cell_type == cell_type & df.eval_files$TF_name == TF_name)
  motif_str = df.eval_files$motif_str[ind]
  
  # edit fraction plot
  df.edit_frac = read_csv(paste0(eval_statsDir, df.eval_files$edit_frac_fnames[ind]))
  df.edit_frac = df.edit_frac %>% mutate(edit_frac_ma = zoo::rollmean(edit_frac, k=3, align = "center", na.pad = T))
  df.edit_frac = df.edit_frac %>% filter(label == 'all_channel')
  df.edit_frac$read_type[df.edit_frac$read_type == 'recently bound'] = 'Unbound\nAccessible'
  df.edit_frac$read_type[df.edit_frac$read_type == 'unbound'] = 'Unbound\nInaccessible'
  df.edit_frac$read_type[df.edit_frac$read_type == 'bound'] = 'Bound'
  x_max = max(df.edit_frac$relative_pos)
  y_max = max(df.edit_frac$edit_frac)
  
  plt1 = ggplot() +
    geom_line(data=df.edit_frac, aes(x=relative_pos, y=edit_frac, color = read_type), linewidth=0.2) + 
    labs(x = "Relative position", y = "Edit fraction", color = NULL, title = paste0(cell_type, " ", TF_name)) + 
    ylim(0, 1) + 
    theme_classic() + 
    theme(
      legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1),
      legend.background = element_rect(fill = NA),  # Transparent background
      legend.key.height = unit(0.8, "cm"), legend.margin = margin(t=0)
    )
  
  # correlation plot
  df.read_type_frac = read_csv(paste0(eval_statsDir, df.eval_files$read_tyoe_frac_fnames[ind]))
  df.read_type_frac = df.read_type_frac %>% 
    filter(read_type == 'bound' & chipseq_norm_signalVal > 0) %>%
    mutate(chipseq_norm_signalVal_log10 = log10(chipseq_norm_signalVal))
  res_spearman = cor.test(df.read_type_frac$read_frac_predicted_prob, df.read_type_frac$chipseq_norm_signalVal_log10, method = 'spearman')
  res_pearson = cor.test(df.read_type_frac$read_frac_predicted_prob, df.read_type_frac$chipseq_norm_signalVal_log10, method = 'pearson')
  # cor_str = paste0('rho=', signif(res_spearman$estimate, digits = 2), ', p=', signif(res_spearman$p.value, digits = 2), '\nr=', signif(res_pearson$estimate, digits = 2), ', p=', signif(res_pearson$p.value, digits = 2))
  cor_str = paste0('rho=', signif(res_spearman$estimate, digits = 2), '\np=', signif(res_spearman$p.value, digits = 2))
  x_max = max(df.read_type_frac$read_frac_predicted_prob)
  y_min = min(df.read_type_frac$chipseq_norm_signalVal_log10)
  
  plt2 = ggplot() + 
    geom_point(data=df.read_type_frac, aes(x=read_frac_predicted_prob, y=chipseq_norm_signalVal_log10), size=0.25, alpha=0.25) + 
    labs(x="OccuPIE mean allelic bound probability", y="Log10 norm. ChIP-seq signal", title = paste0(cell_type, " ", TF_name)) + 
    geom_label(data=df.read_type_frac, x=x_max, y=y_min, label=cor_str, hjust=1, vjust=0, size=8/3.5, fill = "white", label.size=0, alpha=0.5) + 
    geom_smooth(data=df.read_type_frac, aes(x=read_frac_predicted_prob, y=chipseq_norm_signalVal_log10), method="lm", color="black", linewidth=0.5, se=FALSE) +
    theme_classic()
  
  plt = (plt1 | plt2) 
  ggsave(filename = paste0(eval_statsDir, "plots/suplementary_fig6/", cell_type, "_", TF_name, ".comboPlt.jpg"), height = 3.5, width = 7)
}






# plt = plt1 / plt2 + plot_layout(heights = c(5, 10))
# ggsave(filename = paste0(cobind_statsDir, 'plots/K562_STAT5A_center_dist_21_shared_read_edit_map.png'), plot = plt, width = 6, height = 8)
# 









# plt_logo = ggplot() + 
#   # geom_rect(data = df.motif_pair_seqs, xmin=2.5, xmax=11.5, ymin=-Inf, ymax=Inf, fill='lightgrey', alpha=0.25) + 
#   # geom_rect(data = df.motif_pair_seqs, xmin=23.5, xmax=32.5, ymin=-Inf, ymax=Inf, fill='lightgrey', alpha=0.25) + 
#   add_full_rect(xmin = 2.5, xmax = 11.5, fill = "lightgrey", alpha = 0.5) +
#   add_full_rect(xmin = 23.5, xmax = 32.5, fill = "lightgrey", alpha = 0.5) +
#   geom_logo(data = df.motif_pair_seqs$sub_seq) + theme_logo(base_size = 11) +
#   scale_x_continuous(expand = c(0, 0), labels = seq(-10, selected_center_dist+10, 5), breaks = seq(1, 10+1+selected_center_dist+10, 5)) + 
#   scale_y_continuous(expand = c(0, 0)) + 
#   theme(panel.grid = element_blank())
# plt_logo
#   
# ggsave(filename = paste0(out_dir, '/seqLogo_plot.png'), plot = plt_logo, width = 5.6, height = 2)
# 

