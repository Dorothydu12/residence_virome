library(sads)
abund = read.csv("vOTU.csv", header = TRUE, row.names = 1) ## col = species
abund = clr(abund)
sample <- rownames(abund)
abund <- cbind(sample = sample, abund)
abund <- merge(abund, metadata, "sample")
row.names(abund) <- abund$sample
abund1 = abund[, -1]

### vOTU: 718
### rMAG: 342
regional.abund <- abund %>%
  filter(location == "residence 1", site == "leftpalm") %>% 
  dplyr::select(1:718) %>%
  as.matrix() %>%                      # Convert to matrix first
  t() %>%                              # Transpose: species -> rows
  .[rowSums(.) > 0, ] %>%              # Remove species with row sum of 0
  t()   

regional.abund <- colSums(regional.abund)
fit.ln  <- fitsad(regional.abund, sad = "lnorm")
fit.bs  <- fitsad(regional.abund, sad = "bs")
fit.weibull <- fitsad(regional.abund, sad = "weibull")

AICtab(fit.ln, fit.bs, fit.weibull, nobs = sum(regional.abund), weights = TRUE)
best <- fit.weibull
obs.rad = rad(regional.abund)

## Predicted frequencies from a fitted model
fit.weibull.r <- radpred(fit.weibull)
rad.df <- data.frame(
  rank     = obs.rad$rank,
  obs      = obs.rad$abund,
  pred     = fit.weibull.r$abund,
  residual = obs.rad$abund - fit.weibull.r$abund
)
write.csv(rad.df, "model/bacteria_park_campushandrail_handrail.csv")
best

## Rank-abundance plot with observed and predicted frequencies
plot(rad(regional.abund))
lines(fit.ln.r, col = "red")

### Model R-squared
df = read.csv("r1_leftpalm.csv", header = TRUE)
ss_res = sum((log10(df$obs) - log10(df$pred))^2)
ss_tot = sum((log10(df$obs) - log10(mean(df$pred)))^2)
r_squared = 1 - ss_res / ss_tot


### Taylor's power law: log(σ^2) = log(a) + b · log(μ)
library(tidyverse)
library(broom)
library(reshape)

rMAG = read.csv("/Users/shicong/Desktop/Residence-3#/SAD/residential_342_rMAG_normed_coverage.csv", header = TRUE)
dat = melt(rMAG)
colnames(dat) = c("sample", "species", "abund")

vOTU = read.csv("/Users/shicong/Desktop/Residence-3#/SAD/vOTU.csv", header = TRUE)
dat = melt(vOTU)
colnames(dat) = c("sample", "species", "abund")

metadata = read.csv("metadata.csv", header = TRUE)
dat1 = merge(dat, metadata, "sample")

#### Based on residence location and site
tl_fit <- dat1 %>%
  group_by(location, site, species) %>%
  summarise(
    mu  = mean(abund),   # Mean abundance of the species within this location-site
    var = var(abund),    # Variance of the species within this location-site
    .groups = "drop"
  ) %>%
  filter(var > 0, mu > 0) %>%
  mutate(
    log_mu  = log10(mu),
    log_var = log10(var)
  )

write.csv(tl_fit, "site/taylor_rMAG.csv")

# Linear regression to estimate log(a) and b
lm_taylor <- lm(log_var ~ log_mu, data = tl_fit)
## b ≈ 1     Random distribution (Poisson, variance ≈ mean)
## 1 < b < 2 Clumped/aggregated distribution (common)
## b ≈ 2     Strongly aggregated (variance scales with the square of the mean)

ggplot(tl_fit, aes(log_mu, log_var)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  labs(x = "log₁₀(mean abundance)",
       y = "log₁₀(variance)",
       title = sprintf("Taylor’s power law: b = %.2f", coef(lm_taylor)[2]))
