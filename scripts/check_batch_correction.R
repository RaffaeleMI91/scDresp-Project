pre_corrected_extended <- read.csv("/Volumes/iorio/Raffaele/SCDRESP_data/data/pre_corrected_psclones.csv")
pre_corrected_extended$X <- NULL
post_corrected_extended <- read.csv("/Volumes/iorio/Raffaele/SCDRESP_data/data/post_corrected_psclones.csv")
post_corrected_extended$X <- NULL
polyclonal_CCLS <- read.csv("/Volumes/iorio/Raffaele/SCDRESP_data/data/filtered_clustering_res.csv")
clones_IDs <- read.csv("/Volumes/iorio/Raffaele/SCDRESP_data/data/CLONE_IDs.csv")
clones_IDs <- clones_IDs$Clone_ID
pre_corrected_extended$tissue <- NULL

clones_to_remove <- setdiff(polyclonal_CCLS$CellLine, pre_corr_ext_bulk$SangerModelID)

ccls <- polyclonal_CCLS$CellLine
ccls <- ccls[!ccls %in% c("SIDM00128", "SIDM00690")]

pre_corr_ext_bulk <- pre_corrected_extended %>% filter(type=="bulk", SangerModelID %in% polyclonal_CCLS$CellLine)
pre_corr_ext_psclones <- pre_corrected_extended %>% filter(type=="sc pseudoclones", SangerModelID %in% polyclonal_CCLS$CellLine)
rownames(pre_corr_ext_psclones) <- clones_IDs
rownames(pre_corr_ext_bulk) <- ccls[!ccls %in% c("SIDM00128", "SIDM00690")]


# Gower dist
install.packages("StatMatch")
G_all_pre <- StatMatch::gower.dist(rbind.data.frame(
  apply(subset(pre_corr_ext_bulk, select=-c(tissue2,SangerModelID,type)),2,as.numeric),
  apply(subset(pre_corr_ext_psclones, select=-c(tissue2,SangerModelID,type)),2,as.numeric)
))

rownames(G_all_pre) <- c(rownames(pre_corr_ext_bulk),rownames(pre_corr_ext_psclones))
colnames(G_all_pre) <- c(rownames(pre_corr_ext_bulk),rownames(pre_corr_ext_psclones))

G_all_pre <- G_all_pre[
  !(rownames(G_all_pre) %in% c("SIDM00128_clone_0", "SIDM00128_clone_1", "SIDM00690_clone_0", "SIDM00690_clone_1")),
  !(colnames(G_all_pre) %in% c("SIDM00128_clone_0", "SIDM00128_clone_1", "SIDM00690_clone_0", "SIDM00690_clone_1"))
]


post_corr_ext_bulk <- post_corrected_extended %>% filter(type=="bulk", SangerModelID %in% polyclonal_CCLS$CellLine)
post_corr_ext_psclones <- post_corrected_extended %>% filter(type=="sc pseudoclones", SangerModelID %in% polyclonal_CCLS$CellLine)

rownames(post_corr_ext_psclones) <- clones_IDs
rownames(post_corr_ext_bulk) <- ccls[!ccls %in% c("SIDM00128", "SIDM00690")]

G_all_post <- StatMatch::gower.dist(rbind.data.frame(
  apply(subset(post_corr_ext_bulk, select=-c(tissue,SangerModelID,type)),2,as.numeric),
  apply(subset(post_corr_ext_psclones, select=-c(tissue,SangerModelID,type)),2,as.numeric)
))

rownames(G_all_post) <- c(rownames(post_corr_ext_bulk),rownames(post_corr_ext_psclones))
colnames(G_all_post) <- c(rownames(post_corr_ext_bulk),rownames(post_corr_ext_psclones))

G_all_post <- G_all_post[
  !(rownames(G_all_post) %in% c("SIDM00128_clone_0", "SIDM00128_clone_1", "SIDM00690_clone_0", "SIDM00690_clone_1")),
  !(colnames(G_all_post) %in% c("SIDM00128_clone_0", "SIDM00128_clone_1", "SIDM00690_clone_0", "SIDM00690_clone_1"))
]



distances_pre <- list()
distances_post <- list()
for (cl in ccls){
  g_pre <- G_all_pre[grep(cl, rownames(G_all_pre)), grep(cl, rownames(G_all_pre))]
  g_post <-G_all_post[grep(cl, rownames(G_all_post)), grep(cl, rownames(G_all_post))]
  g_pre_upper <- g_pre[upper.tri(g_pre)]
  g_post_upper <- g_post[upper.tri(g_post)]
  distances_pre[[cl]] <- g_pre_upper
  distances_post[[cl]] <- g_post_upper
}

wilcox.test(as.numeric(unlist(distances_pre)), as.numeric(unlist(distances_post)))
t.test(as.numeric(unlist(distances_pre)),
       as.numeric(unlist(distances_post)))

boxplot(as.numeric(unlist(distances_pre)),
        as.numeric(unlist(distances_post)),
        )

distances_pre_numeric <- as.numeric(unlist(distances_pre))
distances_post_numeric <- as.numeric(unlist(distances_post))

distances_df <- data.frame(
  distance = c(distances_pre_numeric, distances_post_numeric),
  group = rep(c("Pre", "Post"), c(length(distances_pre_numeric), length(distances_post_numeric)))
)
distances_df$group <- factor(distances_df$group, levels = c("Pre", "Post"))
# Plot a boxplot with ggplot
ggplot(distances_df, aes(x = group, y = distance, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Gower distance") +
  ggtitle("Intra-similarity ditributions pre- and post-correction") + theme_classic() 
  #scale_fill_manual(values = c("Pre" = "grey", "Post" = "")) 
