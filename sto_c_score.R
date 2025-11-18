
library(dplyr)
###STO
abund=read.csv("/Users/shicong/Desktop/Residence-3#/ecology/SAD/vOTU.csv",header = T,row.names = 1)
pres <- as.data.frame((abund > 0) + 0)
sample <- rownames(pres)
abund <- cbind(sample = sample, pres)
#colnames(abund)[1] <- "sample"
abund=merge(abund,metadata,"sample")
abund=abund[,-1]

test <- abund %>%
  filter(location == "residence 1", site == "left palm") %>%
  tibble::column_to_rownames(var = "timepoint") %>%  # 把 timepoint 变成行名
  dplyr::select(1:342)  
mat=t(test)
st_tbl <- mat %>%
  as.data.frame() %>%
  mutate(species = rownames(.)) %>%
  pivot_longer(-species, names_to = "tp", values_to = "pres") %>%
  mutate(tp = as.integer(tp)) %>%        # 确保时间顺序
  arrange(tp) %>%
  group_by(species) %>%
  mutate(
    cum.pres = cumsum(pres),
    occup    = cum.pres / row_number()   # Oᵢ(t)
  ) %>%
  group_by(tp) %>%
  summarise(mean.occup = mean(occup), .groups = "drop")
write.csv(st_tbl,"../STO/r1_park_campus_handrail.csv")

# 回归得 STO 斜率
sto_fit  <- lm(mean.occup ~ tp, data = st_tbl)
sto_slope <- coef(sto_fit)[["tp"]]

ggplot(st_tbl,aes(y=mean.occup,x= tp)) +
  geom_point(alpha = .15, size = 5) +
  stat_poly_eq(formula = y ~ x,               # 简单线性
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.., sep = "*\", \"*")),
               coef.digits = 2,                # 保留 3 位小数
               parse = TRUE,
               size=3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE)+theme_bw()

library(vegan)
library(EcoSimR)
###C-socre
abund_sub <- abund %>% 
  filter(site == "right palm") %>% 
  dplyr::select(1:718)          # 指定包前缀
##342
# 再去掉列和为 0 的列
test <- abund_sub %>% 
  dplyr::select(where(~ sum(.x, na.rm = TRUE) != 0)) 
test_num <- as.matrix(as.data.frame(lapply(test, as.numeric)))

###null_model
set.seed(888)
null.out <- cooc_null_model(test_num,
                            algo   = "sim9",   # 固定-固定
                            metric = "c_score",
                            nReps  = 10000)
summary(null.out) 
