library(tidyverse)
# 1. 拆出元数据和物种
abund_wide <- read.csv("all_combin.csv",row.names = 1,header = T)
meta  <- abund_wide %>% select(Date,sample_type)
species_mat <- abund_wide %>% select(-Date, -sample_type)

# 2. 给每个样本一个唯一 ID
meta <- meta %>% mutate(sample_id = row_number())

# 3. 变长表：sample_id | species | abundance
long <- species_mat %>%
  dplyr::mutate(sample_id = row_number()) %>%
  pivot_longer(cols = -sample_id,
               names_to  = "species",
               values_to = "abundance") %>%
  left_join(meta, by = "sample_id")

# 构造时间点字符串（可改成你要的任何单位）
long <- long %>%
  mutate(time_point = paste(time_point))

# 检查：每个样本只能出现一次！
stopifnot(anyDuplicated(long$sample_id) == 0)

# 先变宽：行=样本，列=物种
wide <- long %>%
  dplyr::select(sample_id, species, abundance) %>%
  pivot_wider(names_from = species, values_from = abundance)

# 再按时间点拆列表
lst <- long %>%
  dplyr::select(sample_id, time_point) %>%
  distinct() %>%
  split(.$time_point) %>%
  purrr::map(~ pull(., sample_id))   # 显式用 purrr::map

# 按顺序把物种矩阵切片
arrays <- purrr::map(lst, ~ wide %>% filter(sample_id %in% .x) %>% dplyr::select(-sample_id) %>% as.matrix())

# 合并成 3D 数组：样本 × 物种 × 时间
abund_3d <- abind::abind(arrays, along = 3)
dimnames(abund_3d) <- list(
  NULL,                                    # 样本名可不要
  colnames(arrays[[1]]),                   # 物种名
  names(arrays)                            # 时间点名
)

## 1. 按时间点循环：只保留在所有时间点里方差都>0 且 非全NA 的物种
good_sp <- colnames(abund_3d)  # 先假设全部是好的

for (t in 1:dim(abund_3d)[3]) {
  mat_t <- abund_3d[,,t]
  # 去掉全NA或方差=0的列
  keep <- colnames(mat_t)[apply(mat_t, 2, function(x)
    var(x, na.rm = TRUE) > 0 & !all(is.na(x)))]
  good_sp <- intersect(good_sp, keep)
}

## 2. 剔除坏物种
abund_3d_clean <- abund_3d[, good_sp, ]

## 3. 再看维度
dim(abund_3d_clean)

## 4. 重新跑
library(LUPINE)
res <- LUPINE(abund_3d, is.transformed = TRUE, ncomp = 1,cutoff = 0.001)

###impute data
library(splines)         # 负责 ns() 或 smooth.spline()
library(compositions)  
## 1. 读入数据 -------------------------------------------------------------
## 假设 your_data 是一个 n 行×p 列矩阵/数据框，列名就是成分
## 再假设有一个与样本顺序有关的变量 t（长度 = n），如时间、年龄、采样序号
data(your_data)          # 你自己的成分矩阵
t <- your_data$t         # 预测变量（如采样序号）
X  <- as.matrix(your_data[, -match("t", colnames(your_data))])  # 去掉 t 列，只剩成分

## 2. clr 变换 --------------------------------------------------------------
## 用 zCompositions::clr() 会自动把 0 替换为一个极小值并完成 clr
clr.x <- clr(X)          # 返回同样维度的矩阵，缺失位置仍是 NA

## 3. 按列（成分）循环做样条插值 -------------------------------------------
imputed_clr <- clr.x   
for(j in 1:ncol(clr.x)){
  y <- clr.x[, j]
  idx  <- which(!is.na(y))
  nobs <- length(idx)
  if(nobs < 4){              # 太少就线性插值
    imputed_clr[is.na(y), j] <- approx(t[idx], y[idx],
                                       xout = t[is.na(y)])$y
    next
  }
  # 自适应自由度：最多用 nobs-1，且不超过 5
  df <- min(5, nobs - 1)
  fit <- smooth.spline(t[idx], y[idx], df = df)
  imputed_clr[is.na(y), j] <- predict(fit, x = t[is.na(y)])$y
}


ivi.list <- vector("list", 17)

for(i in 1:17){

  g <- graph_from_adjacency_matrix(res[[i]],
                                   mode      = "undirected",
                                   weighted  = NULL)
  

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




g_tp2 <- graph_from_adjacency_matrix(res[[1]], mode ="undirected", weighted =NULL)
write.graph(g_tp1,"/Users/shicong/Desktop/Residence-3#/LUPINE/all_combine/tp18.gml", format="gml")

#calculate infulutial nodes
IVI_tp2 <- data.frame(ivi(graph = g_tp2))

IVI_night <- cbind(IVI_tp2,IVI_tp3,IVI_tp4,IVI_tp5,IVI_tp6,IVI_tp7,IVI_tp8,IVI_tp9,IVI_tp10,IVI_tp11,IVI_tp12,IVI_tp13,IVI_tp14,IVI_tp15,IVI_tp16,IVI_tp17,IVI_tp18) 
pca.res <- prcomp(t(IVI_night), scale. = TRUE)
explained_ratio <- pca.res$sdev^2 / sum(pca.res$sdev^2)
pca_sample <- data.frame(pca.res$x[ ,1:2])
ggplot(data = pca_sample, aes(x =PC1, y = PC2,color=group$group)) + geom_point(size = 3)+ geom_text(label=row.names(pca_sample))+labs(x =  paste('PCA1:', 38.9, '%'), y = paste('PCA2:', 13.7, '%'), color = '')



