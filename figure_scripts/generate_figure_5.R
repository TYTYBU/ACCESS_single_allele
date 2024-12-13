setwd("/Users/tianyu/Downloads")
library(tidyverse)
library(tidyquant)
library(patchwork)
library(ggbeeswarm)
library(ggpubr)
library(ggseqlogo)


### median co-occupancy vs. significance dotplot
order_TF_pairs_by_name <- function(TF1, TF2){
  TF_pair = c(TF1, TF2)
  sorted_TF_pair = sort(TF_pair)
  return(paste0(sorted_TF_pair, collapse = "-"))
}

cobind_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/cobinding_v2_read_thres50_rho0/'
inf_y_pos = 350

cell_type = "K562"
df.cobinding_bound_median_all = read_csv(paste0(cobind_statsDir, 'cobinding_bound_median_all_', cell_type, '.csv'))
df.cobinding_bound_median_all$TF_pairs = apply(df.cobinding_bound_median_all, 1, function(row) order_TF_pairs_by_name(row[['TF1']], row[['TF2']]))
df.cobinding_bound_median_all = df.cobinding_bound_median_all %>% group_by(TF_pairs) %>% slice_head()
df.cobinding_bound_median_all = df.cobinding_bound_median_all %>% mutate(log10_corrected_p = -log10(delta_oe_pval_corrected))
df.cobinding_bound_median_all$log10_corrected_p[is.infinite(df.cobinding_bound_median_all$log10_corrected_p)] = inf_y_pos

# label specific TF groups
df.cobinding_bound_median_all$label = "Other"
TF_list = c("bHLHE40", "MITF", "USF2", "MAX", "USF1")
ind = which((df.cobinding_bound_median_all$TF1 %in% TF_list) | (df.cobinding_bound_median_all$TF2 %in% TF_list))
df.cobinding_bound_median_all$label[ind] = "bHLH-LZ"
TF_list = c("NEUROD1", "ETS1", "MAFF", "REST", "GABPA", "ATF3")
ind = which((df.cobinding_bound_median_all$TF1 %in% TF_list) | (df.cobinding_bound_median_all$TF2 %in% TF_list))
df.cobinding_bound_median_all$label[ind] = "Repressor"
ind = which((df.cobinding_bound_median_all$TF_pairs == "STAT5A-STAT5A"))
df.cobinding_bound_median_all$label[ind] = "STAT5A-STAT5A"
# ind = which(df.cobinding_bound_median_all$delta_oe_median>0.01 & df.cobinding_bound_median_all$delta_oe_median<0.02 & (-log10(df.cobinding_bound_median_all$delta_oe_pval_corrected)>100))
# df.cobinding_bound_median_all$label[ind] = "selected"
df.cobinding_bound_median_all$label = factor(df.cobinding_bound_median_all$label, levels = c("Other", "bHLH-LZ", "Repressor", "STAT5A-STAT5A"))

plt = ggplot() + 
  geom_point(data=df.cobinding_bound_median_all, aes(x=delta_oe_median, y=log10_corrected_p, color=label), size=0.8, alpha=0.75, stroke=0) + 
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.25) + 
  scale_color_manual(values = c("lightgray", "#377eb8", "#e41a1c", "#4daf4a")) + 
  scale_y_continuous(expand = c(0.01,0), breaks=c(seq(0, 300, 100), inf_y_pos), labels = c(seq(0, 300, 100), "Inf.")) + 
  labs(x='Median Obs. - Exp.', y='-log10 P-value', color=NULL) + 
  theme_classic(base_size = 8) + 
  theme(
    legend.position = "inside", legend.position.inside = c(1, 0.01), legend.justification = c(1, 0),
    legend.background = element_rect(fill = NA),  # Transparent background
    legend.key.height = unit(0.3, "cm"), legend.margin = margin(t=0,l=0,b=0.5,r=0), legend.spacing.x = unit(0.02, "cm"), legend.key.width= unit(0.25, "cm"),
    plot.title = element_text(hjust = 0.5, margin = margin(b = 10)) # Adjust top margin
  )
ggsave(filename = paste0(cobind_statsDir, "plots/", cell_type, "_median_coocc_dotplot.pdf"), plot = plt, height = 1.8, width = 3.5, units = "in")




### Observed vs. Expected dotplot
order_TF_pairs_by_start <- function(start1, start2){
  if (start1 <= start2){
    return(paste0(start1, "-", start2))
  } else {
    return(paste0(start2, "-", start1))
  }
}

# plot the TF pairs with top and bottom 5 observed-expected
cobind_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/cobinding_v2_read_thres50_rho0/'
suffix_str = '.cobinding_stats.csv.gz'
cell_type = "HepG2"
TF1_name = "USF1"
fnames = list.files(paste0(cobind_statsDir, cell_type, "_footprint_mean"), pattern = suffix_str)
temp = str_split(sub(pattern=suffix_str, replacement="", fnames), pattern = "_", n = 3, simplify = T)
df.fnames = data.frame(cell_type=temp[,1], TF_name=temp[,2], motif_str=temp[,3], fnames)

df.cobinding_bound_median_all = read_csv(paste0(cobind_statsDir, 'cobinding_bound_median_all_', cell_type, '.csv'))
df.cobinding_bound_median_TF1 = df.cobinding_bound_median_all %>% filter(TF1 == TF1_name) %>% arrange(delta_oe_pval_corrected)
ind = which(df.fnames$cell_type == cell_type & df.fnames$TF_name == TF1_name)
df.cobinding_bound = read_csv(paste0(cobind_statsDir, cell_type, '/', df.fnames$fnames[ind]))
df.cobinding_bound = df.cobinding_bound %>% filter(read_type_pair == "bound--bound")

dir.create(paste0(cobind_statsDir, "plots/"), recursive = T, showWarnings = F)

TF2_name = "BHLHE40"
i = which(df.cobinding_bound_median_TF1$TF1 == TF1_name & df.cobinding_bound_median_TF1$TF2 == TF2_name)
median_coocc = df.cobinding_bound_median_TF1$delta_oe_median[i]
pval_corrected = df.cobinding_bound_median_TF1$delta_oe_pval_corrected[i]

df.cobinding_bound_TF1_TF2 = df.cobinding_bound %>% filter(TF1 == TF1_name & TF2 == TF2_name)
if (TF1_name == TF2_name){
  df.cobinding_bound_TF1_TF2$TF_pairs = apply(df.cobinding_bound_TF1_TF2, 1, function(row) order_TF_pairs_by_start(row[['start1']], row[['start2']]))
  df.cobinding_bound_TF1_TF2 = df.cobinding_bound_TF1_TF2 %>% group_by(TF_pairs) %>% slice_head()
}
xy_max = ceiling(max(c(df.cobinding_bound_TF1_TF2$expected_probability, df.cobinding_bound_TF1_TF2$observed_probability)))

plt = ggplot() + 
  geom_point(data=df.cobinding_bound_TF1_TF2, aes(x=expected_probability, y=observed_probability), alpha=0.25, size=0.75, color='black', stroke=NA) + 
  geom_abline(slope=1, intercept=0, linetype="dashed", color="black", linewidth=0.25) + 
  xlim(0, xy_max) + ylim(0, xy_max) + 
  annotate("text", x=xy_max, y=0, label=paste0("Median Obs.-Exp.\nCo-occ. = ", signif(median_coocc, digits=2), "\nP = ", signif(pval_corrected, digits=2)), hjust=1, vjust=0, size=2) + 
  # scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  # labs(x=expression(P[expected]), y=expression(P[observed]), title=paste0(cell_type, " ", TF1_name, "-", TF2_name, " Co-occ.")) + 
  labs(x=expression(P[expected]), y=expression(P[observed]), title=NULL) + 
  theme_classic(base_size = 8) + 
  theme(plot.title = element_text(size = 8))
ggsave(filename = paste0(cobind_statsDir, "plots/", TF1_name, "_", TF2_name, ".coocc_dotplot.pdf"), plot = plt, height = 1.5, width = 1.8)




### Periodicity vs. significance
cobind_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/cobinding_v2_read_thres50_rho0/'
inf_y_pos = 3.5

df.fft_stats_all = read_csv(paste0(cobind_statsDir, "/cobind_data_new/TF1-all.fft_stats_all.csv"))
df.fft_stats_all = df.fft_stats_all %>% mutate(log10_pval_corrected = -log10(pval_bonferroni))
df.fft_stats_all$log10_pval_corrected[is.infinite(df.fft_stats_all$log10_pval_corrected)] = inf_y_pos

plt = ggplot() +
  geom_point(data=df.fft_stats_all, aes(x=signal_noise_ratio, y=log10_pval_corrected, color=cell_type), alpha=0.75, size=0.8, stroke=0) + 
  scale_color_manual(values = c("#377eb8", "#e41a1c")) + 
  scale_y_continuous(expand = c(0.01,0), breaks=c(seq(0, 3, 1), inf_y_pos), labels = c(seq(0, 3, 1), "Inf.")) + 
  geom_hline(yintercept = 2, linetype='dashed', linewidth = 0.25) + 
  geom_vline(xintercept = 1, linetype='dashed', linewidth = 0.25) + 
  labs(x="Periodicity strength", y="-log10 P-value", color=NULL) + 
  theme_classic(base_size = 8) + 
  theme(
    legend.position = "inside", legend.position.inside = c(1, 0.01), legend.justification = c(1, 0),
    legend.background = element_rect(fill = NA),  # Transparent background
    legend.key.height = unit(0.3, "cm"), legend.margin = margin(t=0,l=0,b=0,r=0), legend.spacing.x = unit(0.01, "cm"), legend.key.width= unit(0.25, "cm")
  )
ggsave(filename = paste0(cobind_statsDir, "plots/periodicity_dotplot.pdf"), plot = plt, height = 1.6, width = 1.8, units = "in")



### TF1-all spacial-resolved co-occupancy plots
get_period_shift <- function(df.fft_stats, df.cobinding_bound_by_pos, period_steps=100, peak=F){
  min_period_shift_list = c()
  df.trough_pos = c()
  df.temp_all = c()
  
  for (i in 1:nrow(df.fft_stats)){
    side = df.fft_stats$side[i]
    period = df.fft_stats$period[i]
    
    if (is.na(period)){
      min_period_shift_list = c(min_period_shift_list, NA)
      next
    } else {
      period_shift_list = seq(0, period, length.out = period_steps+1)
      period_shift_list = period_shift_list[-length(period_shift_list)]
    }
    
    if (side == 'Negative'){
      df.cobinding_bound_by_pos_side = df.cobinding_bound_by_pos %>% filter(center_dist < 0)
    } else {
      df.cobinding_bound_by_pos_side = df.cobinding_bound_by_pos %>% filter(center_dist > 0)
    }
    pos_min = min(df.cobinding_bound_by_pos_side$center_dist)
    pos_max = max(df.cobinding_bound_by_pos_side$center_dist)
    
    median_delta_oe_by_pos_list = c()
    for (period_shift in period_shift_list){
      period_shift_pos_list = round(seq(pos_min+period_shift, pos_max, period))
      temp = df.cobinding_bound_by_pos_side %>% filter(center_dist %in% period_shift_pos_list)
      median_delta_oe_by_pos_list = c(median_delta_oe_by_pos_list, median(temp$delta_oe_smooth, na.rm = T))
    }
    df.temp = data.frame(period_shift_list, median_delta_oe_by_pos_list, side)
    df.temp_all = rbind(df.temp_all, df.temp)
    
    min_period_shift = period_shift_list[ifelse(peak, which.max(median_delta_oe_by_pos_list), which.min(median_delta_oe_by_pos_list))]
    min_period_shift_list = c(min_period_shift_list, min_period_shift)
    df_trough_pos_side = data.frame(center_dist = round(seq(pos_min+min_period_shift, pos_max, period)), side)
    df.trough_pos = rbind(df.trough_pos,df_trough_pos_side )
  }
  
  df.fft_stats$period_shift = min_period_shift_list
  return(list(df.fft_stats, df.trough_pos, df.temp_all))
}

cobind_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/cobinding_v2_read_thres50_rho0/'
df.cobinding_bound_by_pos_all = read_csv(paste0(cobind_statsDir, "/cobind_data_new/TF1-all.cobinding_bound_by_pos.csv.gz"))
df.fft_stats_all = read_csv(paste0(cobind_statsDir, "/cobind_data_new/TF1-all.fft_stats_all.csv"))
window_size = 7

# generate plot for selected TFs
df.selected_TFs = read_csv(paste0(cobind_statsDir, 'selected_TF_list_2.csv'))
dir.create(paste0(cobind_statsDir, 'plots/final_selected_TFs_2/'), recursive = T, showWarnings = F)

ls.plt = list()
for (i in 1:nrow(df.selected_TFs)){
  cell_type = df.selected_TFs$cell_type[i]
  TF_name = df.selected_TFs$TF_name[i]
  
  df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_all %>%
    filter(cell_type == {{cell_type}} & TF_name == {{TF_name}} & bound_filter == '[0.01, 0.99]' & motif_pair_orientation == 'either strand') %>%
    drop_na(delta_oe_by_pos_norm)
  
  df.fft_stats_TF = df.fft_stats_all %>%
    filter(cell_type == {{cell_type}} & TF_name == {{TF_name}})
  
  # calculate moving average
  df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_TF %>% mutate(delta_oe_smooth = signal::sgolayfilt(delta_oe_by_pos_norm, p=2, n=window_size))
  df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_TF %>% mutate(side = ifelse(center_dist>0, "Positive", "Negative"))
  
  # get period shift
  temp = get_period_shift(df.fft_stats_TF, df.cobinding_bound_by_pos_TF, peak=F)
  df.fft_stats_TF = temp[[1]]
  df.fft_stats_TF$side = factor(df.fft_stats_TF$side, levels = c('Positive', 'Negative'))
  df.trough_pos = temp[[2]]
  
  # specific y-axis range
  y_min = df.selected_TFs$y_min[i]
  y_max = df.selected_TFs$y_max[i]
  y_diff = y_max - y_min
  
  # prepare periodicity labels
  x_pos = df.trough_pos %>% filter(side == "Positive") %>% arrange(center_dist) %>% slice_head(n=2) %>% pull(center_dist)
  x_neg = df.trough_pos %>% filter(side == "Negative") %>% arrange(desc(center_dist)) %>% slice_head(n=2) %>% pull(center_dist)
  df.fft_stats_TF = df.fft_stats_TF %>%
    mutate(T_str = paste0("T[Per] == ", round(period*2)/2), # round period to nearest 0.5
           P_str = paste0(ifelse(pval_bonferroni>0, "P[Per] < ", "P[Per] < 0.001"), round(pval_bonferroni, digits = 3)),
           x_label = ifelse(side == "Positive", x_pos[1]*0.95, x_neg[1]*0.95),
           x_arrow1 = ifelse(side == "Positive", x_pos[1], x_neg[1]), x_arrow2 = ifelse(side == "Positive", x_pos[2], x_neg[2]),
           y_label1 = y_min+y_diff*0.2, y_label2 = y_min, y_arrow=y_min+y_diff*0.2)
  
  # make plot
  plt =  ggplot() + 
    geom_point(data = df.cobinding_bound_by_pos_TF, aes(x=center_dist, y=delta_oe_by_pos_norm, color=side), size=0.2, alpha=0.5) + 
    geom_line(data = df.cobinding_bound_by_pos_TF, aes(x=center_dist, y=delta_oe_smooth, color=side), linewidth=0.3) + 
    labs(x='Adjacent motif displacement', y='Norm. Median (Obs. - Exp.)', title = paste0(cell_type, " ", TF_name)) + 
    geom_vline(data = df.trough_pos, aes(xintercept=center_dist, color=side), linetype='dashed', linewidth=0.25) +
    scale_x_continuous(breaks = seq(-100, 100, 10), labels = seq(-100, 100, 10), expand=c(0, 0)) + 
    coord_cartesian(ylim = c(y_min, y_max)) +
    theme_light(base_size = 8) +
    geom_segment(data = df.fft_stats_TF, aes(x=x_arrow1, y=y_arrow, xend=x_arrow2, yend=y_arrow, color=side), arrow=arrow(length = unit(0.1, "cm"), type="closed", ends="both"), linewidth=0.2) +
    geom_text(data = df.fft_stats_TF, aes(x=x_label, y=y_label1, label=T_str, color=side, hjust=ifelse(side == "Positive", 1, 0)), vjust=0, parse=TRUE, size=2) + 
    geom_text(data = df.fft_stats_TF, aes(x=x_label, y=y_label2, label=P_str, color=side, hjust=ifelse(side == "Positive", 1, 0)), vjust=0, parse=TRUE, size=2) + 
    scale_color_manual(values = c("#f1a340", "#998ec3")) +
    theme(legend.position ='none', panel.grid = element_blank(), plot.title = element_text(size = 8))
  
  # ggsave(filename = paste0(cobind_statsDir, 'plots/final_selected_TFs_2/', cell_type, "_", TF_name, ".spatial_coocc_slim_trough_lines.pdf"), width = 10, height = 1.5)
  ls.plt[[paste(cell_type, TF_name, collapse = " ")]] = plt
}

plt_all = wrap_plots(ls.plt, ncol = 2) + plot_layout(axis_titles = "collect")
ggsave(filename = paste0(cobind_statsDir, "plots/4TFs.spatial_coocc_slim.pdf"), plot = plt_all, width = 11, height = 2, units="in")




# # composite periodic plots
# cobind_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/cobinding_v2_read_thres50_rho0'
# df.cobinding_bound_by_pos_all = read_csv(paste0(cobind_statsDir, "/composite_periodicity_test/composite_1/TF1-all.cobinding_bound_by_pos.csv.gz"))
# df.cobinding_bound_by_pos_all = df.cobinding_bound_by_pos_all %>% mutate(center_dist = center_dist_phase_corrected)
# df.fft_stats_all = read_csv(paste0(cobind_statsDir, "/composite_periodicity_test/composite_1/TF1-all.fft_stats_all.csv"))
# window_size = 7
# 
# for (T_str in c("T_1.8nt", "T_0.9nt")){
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_all[which(df.cobinding_bound_by_pos_all$T_str == T_str),]
#   df.fft_stats_TF = df.fft_stats_all[which(df.fft_stats_all$T_str == T_str),]
#   
#   # calculate moving average
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_TF %>% mutate(delta_oe_smooth = signal::sgolayfilt(delta_oe_by_pos_norm, p=2, n=window_size))
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_TF %>% mutate(side = ifelse(center_dist>0, "Positive", "Negative"))
#   
#   # get period shift
#   temp = get_period_shift(df.fft_stats_TF, df.cobinding_bound_by_pos_TF, peak = F)
#   df.fft_stats_TF = temp[[1]]
#   df.fft_stats_TF$side = factor(df.fft_stats_TF$side, levels = c('Positive', 'Negative'))
#   df.trough_pos = temp[[2]]
#   
#   # specific y-axis range
#   y_min = min(df.cobinding_bound_by_pos_TF$delta_oe_by_pos_norm, na.rm = T)
#   y_max = max(df.cobinding_bound_by_pos_TF$delta_oe_by_pos_norm, na.rm = T)
#   y_diff = y_max - y_min
#   
#   # prepare periodicity labels
#   x_pos = df.trough_pos %>% filter(side == "Positive") %>% arrange(center_dist) %>% slice_head(n=2) %>% pull(center_dist)
#   x_neg = df.trough_pos %>% filter(side == "Negative") %>% arrange(desc(center_dist)) %>% slice_head(n=2) %>% pull(center_dist)
#   df.fft_stats_TF = df.fft_stats_TF %>%
#     mutate(T_str = paste0("T[Per] == ", round(period*2)/2), # round period to nearest 0.5
#            P_str = paste0(ifelse(pval>0, "P[Per] < ", "P[Per] < 0.001"), signif(pval, digits = 2)),
#            x_label = ifelse(side == "Positive", x_pos[1]*0.95, x_neg[1]*0.95),
#            x_arrow1 = ifelse(side == "Positive", x_pos[1], x_neg[1]), x_arrow2 = ifelse(side == "Positive", x_pos[2], x_neg[2]),
#            y_label1 = y_min+y_diff*0.2, y_label2 = y_min, y_arrow=y_min+y_diff*0.2)
#   
#   # make plot
#   plt =  ggplot() + 
#     geom_point(data = df.cobinding_bound_by_pos_TF, aes(x=center_dist, y=delta_oe_by_pos_norm, color=side), size=0.2) + 
#     geom_line(data = df.cobinding_bound_by_pos_TF, aes(x=center_dist, y=delta_oe_smooth, color=side)) + 
#     labs(x='Adjacent motif displacement', y='Norm. Median (Obs. - Exp.)') + 
#     geom_vline(data = df.trough_pos, aes(xintercept=center_dist, color=side), linetype='dashed', size=0.25) +
#     geom_segment(data = df.fft_stats_TF, aes(x=x_arrow1, y=y_arrow, xend=x_arrow2, yend=y_arrow, color=side), arrow=arrow(length = unit(0.15, "cm"), type="closed", ends="both"), size=0.5) +
#     geom_text(data = df.fft_stats_TF, aes(x=x_label, y=y_label1, label=T_str, color=side, hjust=ifelse(side == "Positive", 1, 0)), vjust=0, parse=TRUE, size=3) + 
#     geom_text(data = df.fft_stats_TF, aes(x=x_label, y=y_label2, label=P_str, color=side, hjust=ifelse(side == "Positive", 1, 0)), vjust=0, parse=TRUE, size=3) + 
#     scale_x_continuous(breaks = seq(-100, 100, 10), labels = seq(-100, 100, 10)) + 
#     coord_cartesian(ylim = c(y_min, y_max)) +
#     theme_light() + theme(legend.position ='none', panel.grid = element_blank())
#   
#   ggsave(filename = paste0(cobind_statsDir, "/composite_periodicity_test/composite_1/", T_str, "_composite.spatial_coocc_slim.pdf"), plot = plt, width = 10, height = 2)
#   
# }
# 
# 
# # composite periodic plots (median delta_oe level)
# cobind_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/cobinding_v2_read_thres50_rho0'
# df.cobinding_bound_by_pos_all = read_csv(paste0(cobind_statsDir, "/composite_periodicity_test/composite_2/TF1-all.cobinding_bound_by_pos.csv.gz"))
# df.cobinding_bound_by_pos_all = df.cobinding_bound_by_pos_all %>% mutate(center_dist = center_dist_phase_corrected) %>% drop_na()
# df.fft_stats_all = read_csv(paste0(cobind_statsDir, "/composite_periodicity_test/composite_2/TF1-all.fft_stats_all.csv"))
# window_size = 7
# 
# for (T_str in c("T_1.8nt", "T_0.9nt")){
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_all[which(df.cobinding_bound_by_pos_all$T_str == T_str),]
#   df.fft_stats_TF = df.fft_stats_all[which(df.fft_stats_all$T_str == T_str),]
#   
#   # calculate moving average
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_TF %>% mutate(delta_oe_smooth = signal::sgolayfilt(delta_oe_by_pos_norm_by_pos, p=2, n=window_size))
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_TF %>% mutate(side = ifelse(center_dist>0, "Positive", "Negative"))
#   
#   # get period shift
#   temp = get_period_shift(df.fft_stats_TF, df.cobinding_bound_by_pos_TF, peak = F)
#   df.fft_stats_TF = temp[[1]]
#   df.fft_stats_TF$side = factor(df.fft_stats_TF$side, levels = c('Positive', 'Negative'))
#   df.trough_pos = temp[[2]]
#   
#   # specific y-axis range
#   y_min = min(df.cobinding_bound_by_pos_TF$delta_oe_by_pos_norm_by_pos, na.rm = T)
#   y_max = max(df.cobinding_bound_by_pos_TF$delta_oe_by_pos_norm_by_pos, na.rm = T)
#   y_diff = y_max - y_min
#   
#   # make plot
#   plt =  ggplot() + 
#     geom_point(data = df.cobinding_bound_by_pos_TF, aes(x=center_dist, y=delta_oe_by_pos_norm_by_pos, color=side), size=0.2) + 
#     geom_line(data = df.cobinding_bound_by_pos_TF, aes(x=center_dist, y=delta_oe_smooth, color=side)) + 
#     labs(x='Adjacent motif displacement', y='Norm. Median (Obs. - Exp.)') + 
#     scale_x_continuous(breaks = seq(-100, 100, 10), labels = seq(-100, 100, 10)) + 
#     coord_cartesian(ylim = c(y_min, y_max)) +
#     theme_light() + theme(legend.position ='none', panel.grid = element_blank())
#   
#   ggsave(filename = paste0(cobind_statsDir, "/composite_periodicity_test/composite_2/", T_str, "_composite.spatial_coocc_slim.pdf"), plot = plt, width = 10, height = 2)
# }


# 
# # generate plot for all TFs
# cobind_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/cobinding_v2_read_thres50_rho0/'
# df.cobinding_bound_by_pos_all = read_csv(paste0(cobind_statsDir, "/cobind_data_new/TF1-all.cobinding_bound_by_pos.csv.gz"))
# df.fft_stats_all = read_csv(paste0(cobind_statsDir, "/cobind_data_new/TF1-all.fft_stats_all.csv"))
# window_size = 7
# 
# df.fft_stats_max = df.fft_stats_all %>% group_by(cell_type, TF_name) %>% 
#   slice_min(pval_bonferroni) %>% slice_max(signal_noise_ratio) %>% ungroup() %>%
#   arrange(cell_type, pval_bonferroni, desc(signal_noise_ratio))
# 
# for (i in 1:nrow(df.fft_stats_max)){
#   cell_type = df.fft_stats_max$cell_type[i]
#   TF_name = df.fft_stats_max$TF_name[i]
#   
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_all %>%
#     filter(cell_type == {{cell_type}} & TF_name == {{TF_name}} & bound_filter == '[0.01, 0.99]' & motif_pair_orientation == 'either strand') %>%
#     drop_na(delta_oe_by_pos_norm)
#   
#   df.fft_stats_TF = df.fft_stats_all %>%
#     filter(cell_type == {{cell_type}} & TF_name == {{TF_name}})
#   
#   # get period shift
#   temp = get_period_shift(df.fft_stats_TF, df.cobinding_bound_by_pos_TF)
#   df.fft_stats_TF = temp[[1]]
#   df.fft_stats_TF$side = factor(df.fft_stats_TF$side, levels = c('Positive', 'Negative'))
#   df.trough_pos = temp[[2]]
#   
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_TF %>% mutate(delta_oe_smooth = signal::sgolayfilt(delta_oe_by_pos_norm, p=2, n=window_size))
#   df.cobinding_bound_by_pos_TF = df.cobinding_bound_by_pos_TF %>% mutate(side = ifelse(center_dist>0, "Positive", "Negative"))
#   
#   df.fft_stats_pos = df.fft_stats_TF %>% filter(side == 'Positive') %>%
#     mutate(foot_str = paste0("Dominant T: ", round(period, digits=2), ", phase: +", round(period_shift, digits=2)), 
#            head_str = paste0("Strength: ", round(signal_noise_ratio, digits=3), ", corrected P: ", round(pval_bonferroni, digits=2)))
#   df.fft_stats_neg = df.fft_stats_TF %>% filter(side == 'Negative') %>%
#     mutate(foot_str = paste0("Dominant T: ", round(period, digits=2), ", phase: +", round(period_shift, digits=2)), 
#            head_str = paste0("Strength: ", round(signal_noise_ratio, digits=3), ", corrected P: ", round(pval_bonferroni, digits=2)))
#   
#   y_min = min(df.cobinding_bound_by_pos_TF$delta_oe_by_pos_norm)
#   y_max = max(df.cobinding_bound_by_pos_TF$delta_oe_by_pos_norm)
#   y_diff = y_max - y_min
#   x_min = min(df.cobinding_bound_by_pos_TF$center_dist)
#   x_max = max(df.cobinding_bound_by_pos_TF$center_dist)
#   
#   plt =  ggplot() + 
#     geom_point(data = df.cobinding_bound_by_pos_TF, aes(x=center_dist, y=delta_oe_by_pos_norm, color=side), size=0.2) + 
#     geom_line(data = df.cobinding_bound_by_pos_TF, aes(x=center_dist, y=delta_oe_smooth, color=side)) + 
#     labs(x='Adjacent motif displacement', y='Norm. Median (Obs. - Exp.)', title = paste0(cell_type, " ", TF_name)) + 
#     theme_light() + theme(legend.position ='none', panel.grid = element_blank())
#   
#   plt = plt + 
#     geom_vline(data=df.trough_pos, aes(xintercept=center_dist, color=side), linetype='dashed', size=0.25) +
#     annotate('label', x=x_min, y=y_min, label=df.fft_stats_neg$foot_str, hjust=0, vjust=0) +
#     annotate('label', x=x_min, y=y_max, label=df.fft_stats_neg$head_str, hjust=0, vjust=1) +
#     annotate('label', x=x_max, y=y_min, label=df.fft_stats_pos$foot_str, hjust=1, vjust=0) +
#     annotate('label', x=x_max, y=y_max, label=df.fft_stats_pos$head_str, hjust=1, vjust=1)
#   
#   dir.create(paste0(cobind_statsDir, 'cobind_data_new/plots/'), recursive = T, showWarnings = F)
#   ggsave(filename = paste0(cobind_statsDir, 'cobind_data_new/plots/', cell_type, "_", TF_name, ".cobinding.pdf"), width = 8, height = 3)
# }






### STAT5A-STAT5A read edit map
plt1_width = 7.2
plt1_height = 2.8
selected_center_dist = 21
cobind_statsDir = '/Users/tianyu/Downloads/ACCESS_proj/cobinding_v2_read_thres50_rho0/'
out_dir = paste0(cobind_statsDir, 'plots/K562_STAT5A_center_dist_21_shared_read_edit_map/')
dir.create(out_dir, recursive = T, showWarnings = F)

# raster map
df.edits_long = read_csv(paste0(cobind_statsDir, "K562_footprint_mean/K562_STAT5A_center_dist_21_shared_read_edits_long.csv.gz"))
df.edits_long = df.edits_long %>% filter(relative_pos >= -100 & relative_pos <= (100+selected_center_dist))
df.edits_long$read_edit_type = gsub(pattern="_", replacement = "-", df.edits_long$read_edit_type)
df.side = df.edits_long %>% group_by(read_edit_type, index) %>% summarise(delta_oe_prob = mean(delta_oe_prob))
df.side$x_pos = 0

plt_raster = df.edits_long %>% ggplot(aes(x=relative_pos, y=index, fill=edit)) +
  # geom_tile(color = NA) + 
  geom_raster() + 
  # scale_fill_gradient(low = "white", high = "red") +
  scale_fill_viridis_b() +
  labs(x="Position relative to 5' STAT5A motif", y="Read index") +
  # geom_vline(xintercept = c(0, selected_center_dist), linetype='dashed', linewidth=0.2, color='white') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  facet_wrap("read_edit_type", nrow = 2, scales = 'free_y', strip.position = "right") + 
  theme_light(base_size = 8) + 
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        strip.placement = 'outside', strip.clip="off",
        plot.margin = margin(0,0,0,2))

plt_sideBar = df.side %>% ggplot(aes(x=x_pos, y=index, fill=delta_oe_prob)) +
  # geom_tile(color = NA) +
  geom_raster() + 
  scale_fill_gradient2(
    name = NULL,
    low = "blue", mid = "white", high = "red", midpoint = mean(df.side$delta_oe_prob), 
    guide = guide_colorbar(barwidth = 0.5)) +
  labs(x="", y='Allele Obs.-Exp. Co-occ.') +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), position = "right") +
  facet_wrap("read_edit_type", nrow = 2, scales = 'free_y', strip.position = "left") +
  theme_light(base_size = 8) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(0,2,0,2))

plt1 = (plt_raster | plt_sideBar) + plot_layout(widths = c(24, 1))
# ggsave(filename = paste0(out_dir, '/raster_plot.png'), plot = plt1, width = 8, height = 4)
ggsave(filename = paste0(out_dir, '/raster_plot.pdf'), plot = plt1, width = plt1_width, height = plt1_height, units="in")


# edit fraction plot
df.edit_frac = read_csv(paste0(cobind_statsDir, "K562_footprint_mean/K562_STAT5A_center_dist_21_shared_read_edit_frac.csv"))
df.edit_frac = df.edit_frac %>% mutate(edit_frac_ma = zoo::rollmean(edit_frac, k=5, align = "center", na.pad = T)) %>% 
  filter(relative_pos >= -100 & relative_pos <= (100+selected_center_dist))

highlight_range1 = c(-8, 0)
highlight_range2 = c(13, 21)
plt_editFrac = ggplot() + 
  annotate("rect", xmin = highlight_range1[1], xmax = highlight_range1[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
  annotate("rect", xmin = highlight_range2[1], xmax = highlight_range2[2], ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_point(data=df.edit_frac, aes(x=relative_pos, y=edit_frac), size=0.2, alpha=0.5) +
  geom_line(data=df.edit_frac, aes(x=relative_pos, y=edit_frac_ma), linewidth=0.2) + 
  # geom_ma(data=df.edit_frac, aes(x=relative_pos, y=edit_frac), ma_fun = SMA, n = 5, linetype=1, size=0.1) + 
  labs(y='Edit fraction', x=NULL) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme_light(base_size = 8) + 
  theme(plot.margin = margin(0,2,0,2))

# seqLogo plot
df.motif_pair_seqs = read_csv(paste0(cobind_statsDir, "K562_footprint_mean/K562_STAT5A_center_dist_21_motif_pair_seqs.csv"))

flank_len = 150
span_from_center = 10
x_label_step = 5
center1 = 1+flank_len
center2 = 1+flank_len+selected_center_dist
df.motif_pair_seqs$sub_seq = substr(df.motif_pair_seqs$motif_pair_seqs, start = center1-span_from_center, stop = center2+span_from_center)

library(grid)
add_full_rect <- function(xmin, xmax, fill, alpha) {
  annotation_custom(
    grob = rectGrob(
      x = unit(0.5, "npc"), y = unit(0.5, "npc"),
      width = unit(1, "npc"), height = unit(1, "npc"),
      gp = gpar(fill = fill, alpha = alpha, col = NA)
    ),
    xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf
  )
}

plt_logo = ggplot() + 
  annotate("rect", xmin = highlight_range1[1]+span_from_center+1-0.5, xmax = highlight_range1[2]+span_from_center+1+0.5, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
  annotate("rect", xmin = highlight_range2[1]+span_from_center+1-0.5, xmax = highlight_range2[2]+span_from_center+1+0.5, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_logo(data = df.motif_pair_seqs$sub_seq) + 
  # theme_logo(base_size = 11) +
  theme_light(base_size = 8) + 
  scale_x_continuous(expand = c(0, 0), 
                     labels = seq(-span_from_center, selected_center_dist+span_from_center, x_label_step), 
                     breaks = seq(1, (span_from_center+1)+selected_center_dist+span_from_center, x_label_step)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(panel.grid = element_blank(),
        plot.margin = margin(2,2,0,2))

plt2 = (plt_logo / plt_editFrac) + plot_layout(heights = c(1, 2))
ggsave(filename = paste0(out_dir, '/seqLogo_editFrac_plot.pdf'), plot = plt2, width = round(plt1_width*0.825, digits = 2), height = plt1_height/2*1.5, units="in")

