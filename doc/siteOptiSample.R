## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,   # wider images
  fig.height = 6,  # taller images
  fig.retina = 2   # higher resolution
)

## ----setup--------------------------------------------------------------------
library(siteOptiSample)

library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(scico)
library(GGally)
library(patchwork)

## -----------------------------------------------------------------------------

# Import structural metrics raster
data(san_lorenzo_rast)

# Import polygons of existing plots
data(san_lorenzo_plots)

# Prepare data objects
r <- unwrap(san_lorenzo_rast)
p <- san_lorenzo_plots

## -----------------------------------------------------------------------------
# Visualize correlation among structural layers
ggcorr(as.data.frame(r), label = TRUE)

## -----------------------------------------------------------------------------
old_ext <- extractPlotMetrics(r, p, fun = mean)

## -----------------------------------------------------------------------------
old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)

## -----------------------------------------------------------------------------
old_dist_euclidean <- pcaDist(old_pca$r_pca$x, old_pca$p_pca, 
  n_pca = 3, k = 1, method = "euclidean")

old_dist_mahalanobis <- pcaDist(old_pca$r_pca$x, old_pca$p_pca, 
  n_pca = 3, k = 1, method = "mahalanobis")

old_dist_list <- list(
  "euclidean" = old_dist_euclidean,
  "mahalanobis" = old_dist_mahalanobis)

## -----------------------------------------------------------------------------
data.frame(
  euclidean = c(old_dist_euclidean), 
  mahalanobis = c(old_dist_mahalanobis)) %>% 
  ggplot(., aes(x = euclidean, y = mahalanobis)) + 
  geom_abline(linetype = 2, colour = "red") + 
  geom_point(alpha = 0.5) + 
  geom_smooth(colour = "blue") + 
  geom_smooth(method = "lm", colour = "green") + 
  theme_classic() +
  labs(
    x = "Euclidean distance",
    y = "Mahalanobis distance") + 
  ggtitle("Distance in structural space from each pixel to nearest plot")

## -----------------------------------------------------------------------------
data.frame(
  euclidean = c(old_dist_euclidean), 
  mahalanobis = c(old_dist_mahalanobis)) %>% 
  pivot_longer(everything()) %>% 
  ggplot(., aes(x = value, fill = name)) + 
    geom_histogram(position = "identity", alpha = 0.5) + 
    scale_fill_discrete(name = "Distance") + 
    theme_classic() + 
    labs(
      x = "Distance",
      y = "Number of pixels")

## -----------------------------------------------------------------------------
# Prepare matrices with raw variables 
vars_raw <- names(r)[c(1,3,5)]
r_vars_raw <- values(r)[,vars_raw]
old_ext_vars_raw <- old_ext[,vars_raw]

# Run distance function with raw variables
old_dist_raw <- pcaDist(r_vars_raw, old_ext_vars_raw,
  n_pca = 3, k = 1, method = "euclidean")

# Run distance function with weighted raw variables
old_dist_raw_weights <- pcaDist(r_vars_raw, old_ext_vars_raw, w = c(1, 5, 10), 
  n_pca = 3, k = 1, method = "euclidean")

# Compare weighted and unweighted distances
data.frame(
  raw = c(old_dist_raw), 
  raw_weights = c(old_dist_raw_weights)) %>% 
  ggplot(., aes(x = raw, y = raw_weights)) + 
  geom_abline(linetype = 2, colour = "red") + 
  geom_point(alpha = 0.5) + 
  geom_smooth(colour = "blue") + 
  geom_smooth(method = "lm", colour = "green") + 
  theme_classic() +
  labs(
    x = "Distance with unweighted raw data",
    y = "Distance with weighted raw data") + 
  ggtitle("Distance in structural space from each pixel to nearest plot")

## -----------------------------------------------------------------------------
# Prepare rasters with distances
old_dist_r_list <- lapply(old_dist_list, function(x) { 
  out <- r[[1]]
  values(out)[!is.na(values(out))] <- x
  return(out)
})
names(old_dist_r_list) <- names(old_dist_list)

vis_map_list <- lapply(names(old_dist_r_list), function(x) { 
  ggplot() +
    geom_spatraster(data = old_dist_r_list[[x]]) + 
    scale_fill_scico(name = paste0("Relative distance to nearest plot"), 
      palette = "bamako",
      limits = c(0, 8),
      oob = scales::squish) +
    geom_sf(data = p, 
      colour = "#E74B5E", fill = NA) + 
    guides(
      fill = guide_colourbar(
        title.position = "top",
        barwidth = 20,
        barheight = 1)) + 
    theme_bw() +
    ggtitle(x) +
    labs(
      x = NULL, 
      y = NULL) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12))
})
names(vis_map_list) <- names(old_dist_r_list)

wrap_plots(vis_map_list) + 
  plot_layout(
    nrow = 1, 
    guides = "collect") & 
  theme(legend.position = "bottom")

## -----------------------------------------------------------------------------
# Extract PCA values from r_pca
pca_pt_r <- as.data.frame(old_pca$r_pca$x)
pca_pt_r$type <- "Landscape value" 

# Extract PCA values from p_pca
pca_pt_p <- as.data.frame(old_pca$p_pca)
pca_pt_p$type <- "Existing plot" 

# Combine PCA values
pca_pt_all <- rbind(pca_pt_r, pca_pt_p)
pca_pt_all$type <- factor(pca_pt_all$type, 
  levels = c("Landscape value", "Existing plot"))

ggplot() +
  geom_point(data = pca_pt_all, 
    aes(x = PC1, y = PC2, 
      size = type, colour = type, alpha = type),
      shape = "circle") + 
  scale_colour_discrete(name = NULL) + 
  scale_size_manual(name = NULL, 
    values = c("Landscape value" = 1, "Existing plot" = 3)) + 
  scale_alpha_manual(name = NULL, 
    values = c("Landscape value" = 1, "Existing plot" = 0.5)) + 
  labs(x = "PC1", y = "PC2") + 
  theme_bw() +
  ggtitle("Structural space coverage") +
  theme(plot.title = element_text(hjust = 0.5))

## -----------------------------------------------------------------------------
# Run pixel classification
old_pca_classif <- classifLandscape(
  r_pca = old_pca$r_pca$x, 
  p = old_pca$p_pca, 
  n_pca = 3, 
  ci = 0.95
)

# Prepare raster with pixel classifications
r_classif <- r[[1]]
values(r_classif)[!is.na(values(r_classif))] <- old_pca_classif$group
r_classif_cls <- data.frame(
  id = c(1, 2), 
  classif = c("Poorly represented", "Well-represented by existing plots"))
levels(r_classif) <- r_classif_cls

# Create map of pixels and their classification
ggplot() +
  geom_spatraster(data = r_classif) + 
  scale_fill_discrete(name = NULL, guide = guide_legend(nrow = 3),
    na.value = NA, na.translate = FALSE) +
  theme_bw() +
  ggtitle("Landscape representativeness") +
  labs(
    x = NULL, 
    y = NULL) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12))

## -----------------------------------------------------------------------------
old_pca_classif_df <- as.data.frame(old_pca$r_pca$x)
old_pca_classif_df$group <- old_pca_classif$group

ggplot(old_pca_classif_df, aes(PC1, PC2, colour = group)) +
  geom_point(alpha = 0.8) +
  scale_colour_discrete(name = NULL) +
  theme_bw() + 
  labs(x = "PC1", y = "PC2") + 
  theme(legend.position = "bottom")

## -----------------------------------------------------------------------------
# Test with existing plots
test1 <- plotSelect(
  r = r,
  p = p,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE,
  n_pca = 3,
  coord = NULL,
  method = meanminSelect
)

test1$type <- "With existing plots"

# Test without existing plots
test2 <- plotSelect(
  r = r,
  p = NULL,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE, 
  n_pca = 3,
  method = meanminSelect)

test2$type <- "Without existing plots"

# Test without plots PCA 
test3 <- plotSelect(
  r = r,
  p = NULL,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = FALSE, 
  n_pca = NULL,
  method = meanminSelect)
test3$type <- "Using raw raster values"

# Test using all columns
test4 <- plotSelect(
  r = r,
  p = NULL,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE, 
  n_pca = 12,
  method = meanminSelect)
test4$type <- "Using all PCA axes"

# Combine tests
test_1_2_3_4 <- list(
  test1,
  test2, 
  test3,
  test4)

# Create a raster to use as a background
r_bg <- r[[1]]
r_bg[!is.na(values(r_bg))] <- 1
r_bg <- as.factor(r_bg)

# Compare outputs
wrap_plots(lapply(test_1_2_3_4, function(x) {
  ggplot() + 
    geom_spatraster(data = r_bg, show.legend = FALSE) +
    scale_fill_manual(values = "grey", na.value = NA, na.translate = FALSE) + 
    geom_sf(data = x, fill = NA, colour = "red") + 
    theme_void() + 
    ggtitle(unique(x$type))
}))

## -----------------------------------------------------------------------------
# Test with mask layer
r_mask <- r[[1]]
r_mask[values(r_mask) < 35] <- NA_real_
r_mask[!is.na(values(r_mask))] <- 1
r_mask <- as.factor(r_mask)

test5 <- plotSelect(
  r = r,
  p = NULL,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = r_mask,
  pca = TRUE, 
  n_pca = 3,
  method = meanminSelect)
test5$type <- "Using a mask layer"

# Should always overlap mask layer   
ggplot() + 
  geom_spatraster(data = r_mask, show.legend = FALSE) + 
  scale_fill_manual(values = "grey", na.value = NA, na.translate = FALSE) + 
  geom_sf(data = test5, colour = "red", fill = NA) +
  theme_classic()

## -----------------------------------------------------------------------------
# Test with minimax algorithm
test6 <- plotSelect(
  r = r,
  p = p,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE, 
  n_pca = 3,
  method = minimaxSelect)
test6$type <- "minimax"

test_1_6 <- list(
  test1, 
  test6)
test_1_6[[1]]$type <- "meanmin"

wrap_plots(lapply(test_1_6, function(x) {
  ggplot() + 
    geom_spatraster(data = r_bg, show.legend = FALSE) +
    scale_fill_manual(values = "grey", na.value = NA, na.translate = FALSE) + 
    geom_sf(data = x, fill = NA, colour = "red") + 
    theme_void() + 
    ggtitle(unique(x$type))
}))

## -----------------------------------------------------------------------------
# Compare PCA positions using minimax and meanmin algorithms

# Extract metrics from raster for each existing plot
test1_ext <- extractPlotMetrics(r, test1, fun = NULL)
test6_ext <- extractPlotMetrics(r, test6, fun = NULL)
# Using FUN = NULL returns the values of all values within each pixel

# Put pixels and existing plots in the same PCA space
test1_pca <- as.data.frame(PCALandscape(r, test1_ext, center = TRUE, scale. = TRUE)$p_pca)
test1_pca$type <- "meanmin"
test6_pca <- as.data.frame(PCALandscape(r, test6_ext, center = TRUE, scale. = TRUE)$p_pca)
test6_pca$type <- "minimax"

pca_pt_test_1_6 <- rbind(pca_pt_r, test1_pca, test6_pca, pca_pt_p)
pca_pt_test_1_6$type <- factor(pca_pt_test_1_6$type,
  levels = c("Landscape value", "meanmin", "minimax", "Existing plot"))

ggplot() +
  geom_point(data = pca_pt_test_1_6, 
    aes(x = PC1, y = PC2, 
      size = type, colour = type, alpha = type),
      shape = "circle") + 
  scale_colour_discrete(name = NULL) + 
  scale_size_manual(name = NULL, 
    values = c("Existing plot" = 3, "meanmin" = 3, "minimax" = 3, "Landscape value" = 1)) + 
  scale_alpha_manual(name = NULL, 
    values = c("Existing plot" = 0.5, "meanmin" = 0.5, "minimax" = 0.5, "Landscape value" = 1)) + 
  labs(x = "PC1", y = "PC2") + 
  theme_bw() +
  ggtitle("Structural space coverage") +
  theme(plot.title = element_text(hjust = 0.5))

## -----------------------------------------------------------------------------
# Calculate distance from each pixel to nearest plot in PCA space
# Using meanmin and minimax algorithms
test1_dist_euclidean <- pcaDist(old_pca$r_pca$x, 
  as.matrix(test1_pca[,grepl("PC", colnames(test1_pca))]), 
  n_pca = 3, k = 1, method = "euclidean")

test6_dist_euclidean <- pcaDist(old_pca$r_pca$x, 
  as.matrix(test6_pca[,grepl("PC", colnames(test6_pca))]), 
  n_pca = 3, k = 1, method = "euclidean")
 
# Create histograms to compare resulting distances
data.frame(
  existing = c(old_dist_euclidean), 
  meanmin = c(test1_dist_euclidean),
  minimax = c(test6_dist_euclidean)) %>% 
  pivot_longer(everything()) %>% 
  mutate(name = factor(name, 
    levels = c("existing", "meanmin", "minimax"),
    labels = c("Existing plot", "meanmin", "minimax"))) %>% 
  ggplot(., aes(x = value, fill = name)) + 
    geom_histogram(position = "identity", colour = "black") + 
    facet_wrap(~name, scales = "fixed", ncol = 1) + 
    theme_classic() + 
    theme(legend.position = "none") + 
    labs(
      x = "Distance",
      y = "Number of pixels")

## -----------------------------------------------------------------------------
# Extract structural metrics as a non-spatial matrix
r_df <- values(r)
valid_rows <- complete.cases(r_df)
r_df <- r_df[valid_rows,]
r_coord <- crds(r)
r_df_all <- cbind(r_df, r_coord)

# Extract existing plots as row indices in r_df
p_vect <- vect(p)
p_vect <- project(p_vect, crs(r))
r_cov <- cells(r, p_vect)[,"cell"]
p_ind <- match(r_cov, which(valid_rows))
p_ind <- p_ind[!is.na(p_ind)]

## -----------------------------------------------------------------------------
all(PCALandscape(r_df, old_ext, center = TRUE, scale. = TRUE)$r_pca$x == old_pca$r_pca$x)

## -----------------------------------------------------------------------------
test_df1 <- plotSelect(
  r = r_df_all,
  p = p_ind,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE,
  n_pca = 3,
  coord = c("x", "y"),
  method = meanminSelect
)
test_df1_coord <- as.data.frame(r_coord[test_df1,])

test_r1 <- plotSelect(
  r = r,
  p = p,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE,
  n_pca = 3,
  method = meanminSelect
)

test_df2 <- plotSelect(
  r = r_df_all,
  p = p_ind,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE,
  n_pca = 3,
  coord = c("x", "y"),
  method = minimaxSelect
)
test_df2_coord <- as.data.frame(r_coord[test_df2,])

test_r2 <- plotSelect(
  r = r,
  p = p,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE,
  n_pca = 3,
  method = minimaxSelect
)

p_df_r1 <- ggplot() +
  geom_spatraster(data = r[[1]]) + 
  scale_fill_scico(name = "Canopy height (zmax)", palette = "bamako", 
    na.value = NA) +
  geom_sf(data = p, fill = NA, colour = "red") + 
  geom_sf(data = test_r1, fill = NA, colour = "blue", linewidth = 1) + 
  geom_point(data = test_df1_coord, aes(x = x, y = y ), shape = 21, colour = "black", fill = "cyan", size = 3) +
  theme_bw() +
  ggtitle("meanminSelect") + 
  labs(
    x = NULL, 
    y = NULL) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12))
  
p_df_r2 <- ggplot() +
  geom_spatraster(data = r[[1]]) + 
  scale_fill_scico(name = "Canopy height (zmax)", palette = "bamako", 
    na.value = NA) +
  geom_sf(data = p, fill = NA, colour = "red") + 
  geom_sf(data = test_r2, fill = NA, colour = "blue", linewidth = 1) + 
  geom_point(data = test_df2_coord, aes(x = x, y = y ), shape = 21, colour = "black", fill = "cyan", size = 3) +
  theme_bw() +
  ggtitle("minimaxSelect") + 
  labs(
    x = NULL, 
    y = NULL) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12))

wrap_plots(p_df_r1, p_df_r2) + 
  plot_layout(
    nrow = 1,
    guides = "collect")

## -----------------------------------------------------------------------------
het_q_list <- lapply(seq(0.2, 1, 0.2), function(x) { 
  out <- plotSelect(
    r = r,
    p = p,
    n_plots = 10,
    p_new_dim = c(100, 100),
    r_mask = NULL,
    pca = TRUE,
    n_pca = 3,
    het_q = x,
    method = minimaxSelect
  )
  out$het_q <- x
  out
})

het_q_poly <- bind_rows(het_q_list)

ggplot() +
  geom_spatraster(data = r[[1]]) + 
  scale_fill_scico(name = "Canopy height (zmax)", palette = "bamako", 
    na.value = NA) +
  geom_sf(data = p, fill = NA, colour = "red") + 
  geom_sf(data = het_q_poly, aes(colour = het_q), fill = NA, linewidth = 1) + 
  scale_colour_scico(name = "Heterogeneity threshold", palette = "roma") + 
  theme_bw() +
  labs(
    x = NULL, 
    y = NULL) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12))

