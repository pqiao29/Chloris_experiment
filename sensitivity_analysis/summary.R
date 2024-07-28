
t = 10 ### Adjust parameter t /in {10, 20, 50. 100}
ws = c(1:5, 10, 20, 30, 40, 50:55)
count <- method <- K <- NULL
for(w in ws){
    res <- readRDS(paste0("result/w", w, "T", t, ".rds"))
    res_length <- unlist(lapply(res, length))
    res <- res[res_length > 1]
    
    tmp_RDR <- table(unlist(lapply(res, function(x) min(x[1, 1], 8))))
    tmp_both <- table(unlist(lapply(res, function(x) min(x[2, 1], 8))))
    
    count <- c(count, as.matrix(tmp_RDR)/length(res), as.matrix(tmp_both)/length(res))
    method <- c(method, rep(paste0("RDR", w), length(as.matrix(tmp_RDR))), rep(paste0("both", w), length(as.matrix(tmp_both))))
    K <- c(K, names(tmp_RDR), as.numeric(names(tmp_both))*10)
}

K_color_order <- c(8:5, 2:4)
df <- data.frame("method" = method, "count" = count)
K_label <- rep(K_color_order, 2) * rep(c(1, 10), each = 7)
df$K <- factor(as.numeric(K), levels = K_label, labels = K_label)
labels <- paste0(rep(c("RDR", "both"), length(ws)), rep(ws, each = 2))
df$method <- factor(df$method, levels = labels, labels = 1:length(labels))
df$method <- as.numeric(df$method)
df$count[df$K %in% c(2:4, 20, 30, 40)] <- -df$count[df$K %in% c(2:4, 20, 30, 40)] 

bar_colors <- c("#B4C3CE", "#95A5AF", "#597790", "#264C6B", "#B4C3CE", "#95A5AF", "#597790",
                "#FDDBC7", "#eda880", "#d97e4a", "#CC4C02", "#FDDBC7", "#eda880", "#d97e4a")[c(K_color_order, K_color_order*10) %in% as.numeric(K)]
gg <- ggplot(data=df, aes(x = method, y = count, fill = K, width = 1)) +
  geom_bar(stat="identity")  +  scale_fill_manual(values = bar_colors) + 
  scale_x_continuous(breaks = seq(1.5, length(labels), length.out = length(ws)), labels = ws, expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(-1, 1, 0.1), labels = paste0(seq(-1, 1, 0.1)*100, "%")) + 
  geom_segment(aes(x = 0.5, y = -1, xend = 0.5, yend = 1), color = "white", linewidth = 0.75) +
  geom_segment(aes(x = 2.5, y = -1, xend = 2.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 4.5, y = -1, xend = 4.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 6.5, y = -1, xend = 6.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 8.5, y = -1, xend = 8.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 10.5, y = -1, xend = 10.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 12.5, y = -1, xend = 12.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 14.5, y = -1, xend = 14.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 16.5, y = -1, xend = 16.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 18.5, y = -1, xend = 18.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 20.5, y = -1, xend = 20.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 22.5, y = -1, xend = 22.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 24.5, y = -1, xend = 24.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 26.5, y = -1, xend = 26.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 28.5, y = -1, xend = 28.5, yend = 1), color = "white", linewidth = 1.5) + 
  geom_segment(aes(x = 30.5, y = -1, xend = 30.5, yend = 1), color = "white", linewidth = 0.75) + 
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.y = element_blank(), 
    panel.background = element_blank()) 

ggsave(gg, width = 10, height = 6, file = paste0("T", t, ".png"))

