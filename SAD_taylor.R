library(sads)
abund=read.csv("vOTU.csv",header = T,row.names = 1)##col=species
abund=clr(abund)
sample <- rownames(abund)
abund <- cbind(sample = sample, abund)
abund <- merge(abund,metadata,"sample")
row.names(abund) <- abund$sample
abund1=abund[,-1]
###vOTU:718
###rMAG:342
regional.abund <- abund %>%
  filter(location == "residence 1",site == "leftpalm") %>%dplyr::select(1:718) %>%
  as.matrix() %>%                       # 先变成矩阵
  t() %>%                               # 转置：物种→行
  .[rowSums(.) > 0, ] %>%               # 去掉行和为 0 的物种
  t()   
regional.abund <- colSums(regional.abund)
fit.ln  <- fitsad(regional.abund, sad = "lnorm")
fit.bs  <- fitsad(regional.abund, sad = "bs")
fit.weibull <- fitsad(regional.abund, sad = "weibull")

AICtab(fit.ln, fit.bs,fit.weibull,nobs = sum(regional.abund),weights=TRUE)
best <- fit.weibull
obs.rad=rad(regional.abund)
## Predicted frequencies from a fitted model
fit.weibull.r <- radpred(fit.weibull)
rad.df <- data.frame(
  rank   = obs.rad$rank,
  obs    = obs.rad$abund,
  pred   = fit.weibull.r$abund,
  residual = obs.rad$abund - fit.weibull.r$abund
)
write.csv(rad.df,"model/bacteria_park_campushandrail_handrail.csv")
best

## Rank-abundance plot with observed and predicted frequencies
plot(rad(regional.abund))
lines(fit.ln.r,col="red")
###model R2
df=read.csv("r1_leftpalm.csv",header = T)
ss_res=sum((log10(df$obs)-log10(df$pred))^2)
ss_tot =sum((log10(df$obs)-log10(mean(df$pred)))^2)
r_squared = 1 - ss_res / ss_tot



###Taylor’s power law. log(σ^2)=log(a)+b⋅log(μ)
library(tidyverse)
library(broom)
library(reshape)
rMAG=read.csv("/Users/shicong/Desktop/Residence-3#/SAD/residential_342_rMAG_normed_coverage.csv",header = T)
dat=melt(rMAG)
colnames(dat)=c("sample","species","abund")
vOTU=read.csv("/Users/shicong/Desktop/Residence-3#/SAD/vOTU.csv",header = T)
dat=melt(vOTU)
colnames(dat)=c("sample","species","abund")
metadata=read.csv("metadata.csv",header = T)
dat1=merge(dat,metadata,"sample")
####基于residence location和site
tl_fit <- dat1 %>%
  group_by(location, site, species) %>%
  summarise(
    mu  = mean(abund),   # 该 location-site 内该物种的平均丰度
    var = var(abund),    # 该 location-site 内该物种的方差
    .groups = "drop"
  ) %>%
  filter(var > 0, mu > 0) %>%
  mutate(
    log_mu  = log10(mu),
    log_var = log10(var)
  )

write.csv(tl_fit,"site/taylor_rMAG.csv")

# 线性回归估计 log(a) 和 b
lm_taylor <- lm(log_var ~ log_mu, data = tl_fit)
##b ≈ 1	随机分布（泊松，方差≈均值）
##1 < b < 2	聚集分布（常见）
##b ≈ 2	强烈聚集（方差随均值平方增长

ggplot(tl_fit, aes(log_mu, log_var)) +
  geom_point(alpha = .7) +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  labs(x = "log₁₀(mean abundance)",
       y = "log₁₀(variance)",
       title = sprintf("Taylor’s power law: b = %.2f", coef(lm_taylor)[2]))