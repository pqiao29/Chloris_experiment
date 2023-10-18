res_compare <- readRDS("result/compare.rds")

I_Chloris_RDR <- unlist(lapply(res_compare, function(x) x[1, 1]))
I_Chloris_both10 <- unlist(lapply(res_compare, function(x) x[2, 1]))
I_inferCNV <- unlist(lapply(res_compare, function(x) x[3, 1]))
I_HB_RDR <- unlist(lapply(res_compare, function(x) x[4, 1]))
I_HB_both10 <- unlist(lapply(res_compare, function(x) x[5, 1]))
I_HB_both50 <- unlist(lapply(res_compare, function(x) x[6, 1]))
I_HB_both100 <- unlist(lapply(res_compare, function(x) x[7, 1]))


H_Chloris_RDR <- unlist(lapply(res_compare, function(x) x[1, 2]))
H_Chloris_both10 <- unlist(lapply(res_compare, function(x) x[2, 2]))
H_inferCNV <- unlist(lapply(res_compare, function(x) x[3, 2]))
H_HB_RDR <- unlist(lapply(res_compare, function(x) x[4, 2]))
H_HB_both10 <- unlist(lapply(res_compare, function(x) x[5, 2]))
H_HB_both50 <- unlist(lapply(res_compare, function(x) x[6, 2]))
H_HB_both100 <- unlist(lapply(res_compare, function(x) x[7, 2]))

K_Chloris_RDR <- unlist(lapply(res_compare, function(x) x[1, 3]))
K_Chloris_both10 <- unlist(lapply(res_compare, function(x) x[2, 3]))
K_inferCNV <- unlist(lapply(res_compare, function(x) x[3, 3]))


### Clustering
N <- length(res_compare)
df_I <- data.frame("accuracy" = c(I_HB_RDR, I_HB_both10, I_HB_both100, I_inferCNV, I_Chloris_RDR, I_Chloris_both10), 
                   "label" = rep(c("HB_RDR", "HB_both10", "HB_both100", "inferCNV", "RDR", "full"), each = N), 
                   "xlabel" = rep(c(0, 0.3, 0.6, 1.3, 2.2, 2.8), each = N),
                   "width" = c(rep(1.5, 3*N), rep(2, N), rep(0.8, 2*N)))

method_label_order <- c("HB_RDR", "HB_both10", "HB_both100", "inferCNV", "RDR", "full")
df_I$label <- factor(df_I$label, levels = method_label_order, labels = method_label_order)
y_ticks <- c(median(I_HB_RDR), median(I_HB_both10), median(I_HB_both100), median(I_inferCNV), median(I_Chloris_RDR))
y_ticks <- sort(unique(y_ticks), decreasing = TRUE)

gg_cluster <- ggplot(df_I, aes(x = xlabel, y = accuracy, fill = factor(label), color = label, width = width))  +
  geom_boxplot(aes(weight = width), outlier.size = 1, varwidth = TRUE, alpha = 0.4) +
  geom_violin(alpha = 1, position = position_dodge(c(0.8))) +
  scale_fill_manual(values = c("#F7FCB9", "#ADDD8E",  "#41AB5D", "#DFC27D", "#F4A582",  "#4393C3")) +
  scale_color_manual(values = c(rep("#006837", 3),  "#8C510A", "#B2182B", "#2166AC")) +
  geom_segment(aes(x = -1, y = y_ticks[1], xend = 3.5, yend = y_ticks[1]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  geom_segment(aes(x = -1, y = y_ticks[2], xend = 1.8, yend = y_ticks[2]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  geom_segment(aes(x = -1, y = y_ticks[3], xend = 1, yend = y_ticks[3]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  geom_segment(aes(x = -1, y = y_ticks[4], xend = 0.2, yend = y_ticks[4]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  geom_segment(aes(x = -1, y = y_ticks[5], xend = 0.6, yend = y_ticks[5]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  theme(legend.position = "none", 
        legend.text = element_blank(),
        legend.title = element_blank(), 
        legend.key.size = unit(1, 'cm'),
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(),  
        panel.background = element_blank())

ggsave(gg_cluster, width = 9, height = 8, file = "result/clustering.png")


### States
df_H <- data.frame("accuracy" = c(H_HB_RDR, H_HB_both10, H_HB_both100, H_inferCNV, H_Chloris_RDR, H_Chloris_both10), 
                   "label" = rep(c("HB_RDR", "HB_both10", "HB_both100", "inferCNV", "RDR", "full"), each = N), 
                   "xlabel" = rep(c(0, 0.7, 1.2, 2.3, 3.4, 4.2), each = N),
                   "width" = c(rep(1, 2*N), rep(0.7, N), rep(1.2, 3*N)))

df_H$label <- factor(df_H$label, levels = method_label_order, labels = method_label_order)

y_ticks <- c(median(H_HB_RDR), median(H_inferCNV), median(H_Chloris_RDR),  median(H_Chloris_both10))
y_ticks <- sort(unique(y_ticks), decreasing = TRUE)

gg_state <- ggplot(df_H, aes(x = xlabel, y = accuracy, fill = factor(label), color = label, width = width))  +
  geom_violin(alpha = 0.4, position = position_dodge(c(0.8))) +
  geom_boxplot(aes(weight = width), outlier.size = 1, varwidth = TRUE) +
  scale_fill_manual(values = c("#F7FCB9", "#ADDD8E",  "#41AB5D", "#DFC27D", "#F4A582",  "#4393C3")) +
  scale_color_manual(values = c(rep("#006837", 3),  "#8C510A", "#B2182B", "#2166AC")) +
  geom_segment(aes(x = 3.8, y = y_ticks[1], xend = 5, yend = y_ticks[1]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  geom_segment(aes(x = 3, y = y_ticks[2], xend = 5, yend = y_ticks[2]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  geom_segment(aes(x = 1, y = y_ticks[3], xend = 5, yend = y_ticks[3]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  geom_segment(aes(x = -1, y = y_ticks[4], xend = 5, yend = y_ticks[4]), color = "gray30", linetype="dotted", linewidth=0.6) + 
  theme(legend.position = "none", 
        legend.text = element_blank(),
        legend.title = element_blank(), 
        legend.key.size = unit(1, 'cm'),
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(),  
        panel.background = element_blank())

ggsave(gg_state, width = 9, height = 8, file = "result/CNstate.png")

#### K
## =============================== load data ==============================================
load("~/Projects/HMM4/sims/results/inferCNV_r1/infercnv_K5_outlier5_a.Rdata")
res_fin <- res_fin[unlist(lapply(res_fin, function(x) length(x) > 1))]
n_inferCNV <- length(res_fin)
K_infercnv_RDR <- unlist(lapply(res_fin, function(x) x[1, 4]))
K_Chloris_RDR <- unlist(lapply(res_fin, function(x) x[2, 4]))



## =============================== reformat ==============================================
possible_Ks <- c(3, 7, 4, 6, 5)

K_cnt_infercnv_RDR <- NULL
K_cnt_Chloris_RDR <- NULL
K_cnt_Chloris_both <- NULL
for(k in possible_Ks){
  K_cnt_infercnv_RDR <- rbind(K_cnt_infercnv_RDR, c(k, sum(K_inferCNV == k)))
  K_cnt_Chloris_RDR <- rbind(K_cnt_Chloris_RDR, c(k*10, sum(K_Chloris_RDR == k)))
  K_cnt_Chloris_both <- rbind(K_cnt_Chloris_both, c(k*100, sum(K_Chloris_both10 == k)))
}

K_cnt_infercnv_RDR[, 2] <- K_cnt_infercnv_RDR[, 2]/N
K_cnt_Chloris_RDR[, 2] <- K_cnt_Chloris_RDR[, 2]/N
K_cnt_Chloris_both[, 2] <- K_cnt_Chloris_both[, 2]/N

K_cnt_all <- rbind(K_cnt_infercnv_RDR, K_cnt_Chloris_RDR, K_cnt_Chloris_both)
methods <- c("inferCNV", "Chloris (RDR)", "Chloris (both)")
method_label <- rep(methods, each = length(possible_Ks))
K_true_cnt <- c(K_cnt_infercnv_RDR[length(possible_Ks), 2], K_cnt_Chloris_RDR[length(possible_Ks), 2], K_cnt_Chloris_both[length(possible_Ks), 2])


## =============================== plot ==============================================
df <- data.frame("K" = K_cnt_all[, 1], "count" = K_cnt_all[, 2], "method" = method_label)
df$method <- factor(df$method, levels = methods, labels = 1:length(methods))
df$method <- as.numeric(df$method)
K_label <- rep(possible_Ks, length(methods)) * rep(c(1, 10, 100), each = length(possible_Ks))
df$K <- factor(df$K, levels = K_label, labels = K_label)

bar_colors <- c("#faf4ee", "#faf4ee", "#fae5ac", "#fae5ac", "#f6bd6a", 
                "#F4A582", "#F4A582", "#F4A582", "#F4A582", "#B2182B", 
                "#D1E5F0", "#D1E5F0", "#92b5d8", "#92b5d8", "#1b3461")

gg_K <- ggplot(data=df, aes(x = method, y = count, fill = K, width = 1)) +
  geom_bar(stat="identity") + 
  #scale_x_continuous(breaks = c(1:length(methods)), labels = c("inferCNV", "Chloris(RDR)", "Chloris(both)")) + 
  geom_segment(aes(x = 0.5, y = 0, xend = 0.5, yend = 1), color = "white", size = 0.75) + 
  geom_segment(aes(x = 1.5, y = 0, xend = 1.5, yend = 1), color = "white", size = 1.5) + 
  geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 1), color = "white", size = 1.5) + 
  geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = 1), color = "white", size = 0.75) + 
  annotate("text", x = c(1:length(K_true_cnt)), y = K_true_cnt - 0.06, label= paste0(round(K_true_cnt*100), "%"), color = "white",  fontface = "bold", size = 8) +  
  annotate("text", x = c(1:length(K_true_cnt)), y = 0.1, label = c("inferCNV\n(RDR)", "Our model\n(RDR)", "Our model\n(RDR + BAF)"), color = "white", size = 8) +  
  scale_fill_manual(values = bar_colors) +
  theme(legend.position="none",
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.background = element_blank())


ggsave(gg_K, width = 6, height = 5, file = "result/K.png")


