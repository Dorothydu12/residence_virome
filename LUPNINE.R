library(tidyverse)

# 1. Extract metadata and species abundance matrix
abund_wide <- read.csv("all_combin.csv", row.names = 1, header = TRUE)
meta  <- abund_wide %>% select(Date, sample_type)
species_mat <- abund_wide %>% select(-Date, -sample_type)

# 2. Assign a unique ID to each sample
meta <- meta %>% mutate(sample_id = row_number())

# 3. Convert to long format: sample_id | species | abundance
long <- species_mat %>%
  dplyr::mutate(sample_id = row_number()) %>%
  pivot_longer(cols = -sample_id,
               names_to  = "species",
               values_to = "abundance") %>%
  left_join(meta, by = "sample_id")

# Construct time point string (can be changed to any desired unit)
long <- long %>%
  mutate(time_point = paste(time_point))

# Check: Each sample must appear only once!
stopifnot(anyDuplicated(long$sample_id) == 0)

# Convert back to wide format: rows = samples, columns = species
wide <- long %>%
  dplyr::select(sample_id, species, abundance) %>%
  pivot_wider(names_from = species, values_from = abundance)

# Split list by time points
lst <- long %>%
  dplyr::select(sample_id, time_point) %>%
  distinct() %>%
  split(.$time_point) %>%
  purrr::map(~ pull(., sample_id))    # Explicitly use purrr::map

# Slice species matrix sequentially by time point
arrays <- purrr::map(lst, ~ wide %>% filter(sample_id %in% .x) %>% dplyr::select(-sample_id) %>% as.matrix())

# Combine into a 3D array: samples × species × time
abund_3d <- abind::abind(arrays, along = 3)
dimnames(abund_3d) <- list(
  NULL,                      # Sample names can be omitted
  colnames(arrays[[1]]),     # Species names
  names(arrays)              # Time point names
)

## 1. Loop through time points: Retain only species with variance > 0 and not all NA across all time points
good_sp <- colnames(abund_3d)  # Assume all are good initially

for (t in 1:dim(abund_3d)[3]) {
  mat_t <- abund_3d[,,t]
  # Remove columns that are all NA or have zero variance
  keep <- colnames(mat_t)[apply(mat_t, 2, function(x)
    var(x, na.rm = TRUE) > 0 & !all(is.na(x)))]
  good_sp <- intersect(good_sp, keep)
}

## 2. Filter out low-quality/invalid species
abund_3d_clean <- abund_3d[, good_sp, ]

## 3. Check dimensions
dim(abund_3d_clean)

## 4. Run LUPINE analysis
library(LUPINE)
res <- LUPINE(abund_3d, is.transformed = TRUE, ncomp = 1, cutoff = 0.001)


### Impute missing data
library(splines)       # Handles ns() or smooth.spline()
library(compositions)  

## 1. Load data -------------------------------------------------------------
## Assume your_data is an n x p matrix/data.frame, where columns represent components
## Assume there is a variable t (length = n) related to sample order, such as time, age, or sampling sequence
data(your_data)         # Your own composition matrix
t <- your_data$t         # Predictor variable (e.g., sampling sequence)
X  <- as.matrix(your_data[, -match("t", colnames(your_data))])  # Exclude column t, leaving only components

## 2. CLR transformation --------------------------------------------------------------
## zCompositions::clr() automatically replaces zeros with a small value and performs CLR
clr.x <- clr(X)          # Returns a matrix of the same dimension; missing positions remain NA

## 3. Loop through columns (components) for spline interpolation ---------------------------
imputed_clr <- clr.x   
for(j in 1:ncol(clr.x)){
  y <- clr.x[, j]
  idx  <- which(!is.na(y))
  nobs <- length(idx)
  if(nobs < 4){             # Use linear interpolation if too few observations
    imputed_clr[is.na(y), j] <- approx(t[idx], y[idx],
                                       xout = t[is.na(y)])$y
    next
  }
  # Adaptive degrees of freedom: at most nobs-1, capped at 5
  df <- min(5, nobs - 1)
  fit <- smooth.spline(t[idx], y[idx], df = df)
  imputed_clr[is.na(y), j] <- predict(fit, x = t[is.na(y)])$y
}


ivi.list <- vector("list", 17)

for(i in 1:17){
  g <- graph_from_adjacency_matrix(res[[i]],
                                   mode     = "undirected",
                                   weighted = NULL)
  
  ivi.df <- data.frame(ivi(graph = g))
  
  ivi.df$subgraph_id <- i+1
  
  ivi.list[[i]] <- ivi.df
}

IVI.all <- do.call(cbind, ivi.list)
IVI.all <- IVI.all[ !grepl("subgraph_id", names(IVI.all), ignore.case = TRUE) ]
ivi_cols <- grep("^ivi", names(IVI.all))
names(IVI.all)[ivi_cols] <- paste0("d", 2:18)
write.csv(IVI.all, "IVI_all_subgraphs.csv", row.names = FALSE)

head(IVI.all)


g_tp2 <- graph_from_adjacency_matrix(res[[1]], mode ="undirected", weighted = NULL)
write.graph(g_tp1, "/Users/shicong/Desktop/Residence-3#/LUPINE/all_combine/tp18.gml", format = "gml")

# Calculate influential nodes
IVI_tp2 <- data.frame(ivi(graph = g_tp2))

IVI_night <- cbind(IVI_tp2, IVI_tp3, IVI_tp4, IVI_tp5, IVI_tp6, IVI_tp7, IVI_tp8, IVI_tp9, IVI_tp10, IVI_tp11, IVI_tp12, IVI_tp13, IVI_tp14, IVI_tp15, IVI_tp16, IVI_tp17, IVI_tp18) 

pca.res <- prcomp(t(IVI_night), scale. = TRUE)
explained_ratio <- pca.res$sdev^2 / sum(pca.res$sdev^2)
pca_sample <- data.frame(pca.res$x[ ,1:2])

ggplot(data = pca_sample, aes(x = PC1, y = PC2, color = group$group)) + 
  geom_point(size = 3) + 
  geom_text(label = row.names(pca_sample)) + 
  labs(x = paste('PCA1:', 38.9, '%'), y = paste('PCA2:', 13.7, '%'), color = '')
