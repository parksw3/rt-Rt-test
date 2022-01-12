library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(egg)
source("renewal_det.R")

## some arbitrary changes in R
Rfun <- function(t) {
  if (t < 20) {
    1.5
  } else if (t >= 20 & t < 40) {
    1.5 - 1/20 * (t-20)
  } else if (t >= 40 & t < 60) {
    0.5
  } else if (t >= 60 & t < 80) {
    0.5 + 2/20 * (t-60)
  } else if (t >= 80) {
    2.5 - 1/20 * (t-80)
  }
}

gorig <- 3.2
sdorig <- 2.1
cvorig <- sdorig/gorig

gvec <- c(2.4, 3.2, 4)
sdvec <- c(1.6, 2.1, 2.6)
thetavec <- c(0.6, 1, 1.4)

pardata <- expand.grid(gvec, sdvec, thetavec)

reslist <- vector('list', nrow(pardata))

for (i in 1:nrow(pardata)) {
  pp <- pardata[i,]
  cv <- pp[[2]]/pp[[1]]
  
  rr <- renewal_det(Rfun=Rfun,
                    genfun1=function(x) dgamma(x, shape=1/cvorig^2, rate=1/cvorig^2/gorig),
                    genfun2=function(x) dgamma(x, shape=1/cv^2, rate=1/cv^2/pp[[1]]),
                    theta=pp[[3]],
                    tmax=100)
  
  approxR1 <- (1 + cvorig^2 * rr$r1 * gorig)^(1/cvorig^2)
  approxR2 <- (1 + cvorig^2 * rr$r1 * gorig)^(1/cvorig^2) * pp[[3]]
  
  reslist[[i]] <- data.frame(
    time=rr$tvec,
    truer1=rr$r1,
    truer2=rr$r2,
    trueR1=rr$Rt1,
    trueR2=rr$Rt2,
    approxR1=approxR1,
    approxR2=approxR2,
    approxr2=(approxR2^(cv^2) - 1)/(cv^2 * pp[[1]]),
    mean=pp[[1]],
    sd=pp[[2]],
    theta=pp[[3]]
  )
}

resdata <- reslist %>%
  lapply(function(x) tail(x, -1)) %>%
  bind_rows

example <- resdata %>%
  filter(mean==3.2, sd==2.1, theta==1.4)
  
g1 <- ggplot(example) +
  geom_line(aes(time, trueR1, col="Variant 1")) +
  geom_line(aes(time, trueR2, col="Variant 2")) +
  scale_x_continuous("Time") +
  scale_y_continuous("Reproduction number") +
  ggtitle("A") +
  theme(
    legend.title = element_blank()
  )

g2 <- ggplot(example) +
  geom_line(aes(time, truer1, col="Variant 1")) +
  geom_line(aes(time, truer2, col="Variant 2")) +
  scale_x_continuous("Time") +
  scale_y_continuous("Growth rate") +
  ggtitle("B") +
  theme(
    legend.title = element_blank()
  )

g3 <- ggplot(example) +
  geom_line(aes(time, trueR1, col="Variant 1", lty="True")) +
  geom_line(aes(time, trueR2, col="Variant 2", lty="True")) +
  geom_line(aes(time, approxR1, col="Variant 1", lty="Approx")) +
  geom_line(aes(time, approxR2, col="Variant 2", lty="Approx")) +
  scale_x_continuous("Time") +
  scale_y_continuous("Reproduction number") +
  scale_linetype_manual(values=c(2, 1)) +
  ggtitle("C") +
  theme(
    legend.title = element_blank()
  )

g4 <- ggplot(example) +
  geom_line(aes(time, truer1, col="Variant 1", lty="True")) +
  geom_line(aes(time, truer2, col="Variant 2", lty="True")) +
  geom_line(aes(time, approxr2, col="Variant 2", lty="Approx")) +
  scale_x_continuous("Time") +
  scale_y_continuous("Growth rate") +
  scale_linetype_manual(values=c(2, 1)) +
  ggtitle("D") +
  theme(
    legend.title = element_blank()
  )

gtot <- ggarrange(g1, g2, g3, g4)

ggsave("simulate_example.png", gtot, width=12, height=8)

g5 <- ggplot(resdata) +
  geom_line(aes(time, truer1, col=factor(theta)), col="black", lwd=1) +
  geom_line(aes(time, truer2, col=factor(theta))) +
  geom_line(aes(time, approxr2, col=factor(theta)), lty=2) +
  scale_x_continuous("Time") +
  scale_y_continuous("Growth rate") +
  scale_color_viridis_d("Transmission advantage", end=0.8) +
  facet_grid(mean~sd)  +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  )

ggsave("simulate_example2.png", g5, width=12, height=8)

g6 <- ggplot(resdata) +
  # geom_line(aes(time, trueR1, col=factor(theta)), col="black", lwd=1) +
  geom_line(aes(time, trueR2, col=factor(theta))) +
  # geom_line(aes(time, approxR1, col=factor(theta)), lty=2) +
  geom_line(aes(time, approxR2, col=factor(theta)), lty=2) +
  scale_x_continuous("Time") +
  scale_y_continuous("Variant reproduction number") +
  scale_color_viridis_d("Transmission advantage", end=0.8) +
  facet_grid(mean~sd)  +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  )

ggsave("simulate_example3.png", g6, width=12, height=8)
