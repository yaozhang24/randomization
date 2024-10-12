library(readxl)
library(data.table)

## load data

data <- read_excel("../data/australia/S2 Data.xlsx") ## Can be downloaded from https://doi.org/10.1371/journal.pmed.1002412.s002

data <- data.table(data)

#trial 1 or 2
data <- subset(data,hospital=="Dandenong")

trial1 <- subset(data, study1 == 1)
## trial_2 <- subset(data,study1==0)

# What's the realized treatment order?
table(trial1[, c("index_ward", "sw_step", "no_we_exposure")])
cross_over_realized <- c(3, 5, 7, 6, 2, 4)

#
trial1[, treatment_status := as.numeric(sw_step >= cross_over_realized[index_ward])]
trial1[, log_acute_los := log(acute_los)]

library(ggplot2)
library(ggrepel)
## ggplot(trial1) + aes(x = factor(treatment_status), y = log_acute_los) + geom_boxplot() + facet_grid(~ index_ward)
trial1_mean <- trial1[,
                       .(log_acute_los = mean(log(acute_los)),
                         treatment_status = mean(treatment_status)),
                      by = c("index_ward", "sw_step")]

trial1_mean$ward_by_cross_over <- factor(trial1_mean$index_ward)
levels(trial1_mean$ward_by_cross_over) <- c("B", "D", "F", "E", "A", "C")

pl_mean <- ggplot(trial1_mean) + aes(x = sw_step, y = log_acute_los, group = ward_by_cross_over) + geom_line(aes(color = factor(treatment_status))) + geom_label_repel(aes(label = ward_by_cross_over, fill = factor(treatment_status), size = factor(treatment_status)), color = "white", force = 0.001, point.padding = NA, alpha = 0.9) + theme_classic(base_size = 15) + theme(legend.position = "top", legend.title = element_blank()) + xlab("Month") + ylab("Log length of stay") + scale_fill_discrete(labels = c("Regular weekend health service", "No weekend health service")) + scale_color_discrete(labels = c("Regular weekend health service", "No weekend health service")) + scale_size_discrete(range = c(3, 5), labels = c("Regular weekend health service", "No weekend health service"))

ggsave("ward-mean.pdf", pl_mean, width = 10, height = 7)

get_statistic <- function(index_ward,
                          sw_step,
                          log_acute_los,
                          cross_over) {

    treatment_status <- sw_step >= cross_over[index_ward]
    c(lm(log_acute_los ~ treatment_status)$coef[2],
      lm(log_acute_los ~ treatment_status + as.factor(index_ward))$coef[2],
      lm(log_acute_los ~ treatment_status + as.factor(index_ward) +
             sw_step)$coef[2])

}

## Permuting crossover times: individual data, exact p-value
T_obs <- get_statistic(trial1$index_ward,
                       trial1$sw_step,
                       trial1$log_acute_los,
                       cross_over_realized)

## p2 <- sapply(combinat::permn(2:7),
##              get_statistic,
##              index_ward = trial1$index_ward,
##              sw_step = trial1$sw_step,
##              log_acute_los = trial1$log_acute_los,
##              beta = 0)
## p2 <- t(p2)

library(pbapply)
nsim <- 720 * 5
ward_perm <- pbreplicate(nsim, sample(trial1$index_ward), simplify = FALSE, cl = 6)
time_perm <- pbreplicate(nsim, sample(trial1$sw_step), simplify = FALSE, cl = 6)
cross_over_perm <- pbreplicate(nsim / 720, combinat::permn(2:7))

settings <- expand.grid(ward = c(FALSE, TRUE),
                        time = c(FALSE, TRUE),
                        cross_over = c(FALSE, TRUE),
                        beta = seq(-0.2, 0.8, 0.01),
                        sim = 1:nsim)

get_statistic_perm <- function(setting) {

    if (setting$ward) {
        index_ward <- ward_perm[[setting$sim]]
    } else {
        index_ward <- trial1$index_ward
    }

    if (setting$time) {
        sw_step <- time_perm[[setting$sim]]
    } else {
        sw_step <- trial1$sw_step
    }

    if (setting$cross_over) {
        cross_over <- cross_over_perm[[setting$sim]]
    } else {
        cross_over <- cross_over_realized
    }

    get_statistic(index_ward, sw_step, trial1$log_acute_los - setting$beta * trial1$treatment_status, cross_over)

}

library(pbapply)
test_stat <- pbsapply(1:nrow(settings), function(i) get_statistic_perm(settings[i, ]), cl = 6)

test_stat <- t(test_stat)
colnames(test_stat) <- c("T1", "T2", "T3")
sim_results <- cbind(settings, test_stat)

save(sim_results, file = "sim_results3.rda")

################################################################################

load("sim_results2.rda") # 720 * 5
tmp <- sim_results
load("sim_results3.rda") # 720 * 5
sim_results <- rbind(tmp, sim_results)


sim_results$ward <- factor(sim_results$ward)
levels(sim_results$ward) <- c("", "Permute ward")

sim_results$time <- factor(sim_results$time)
levels(sim_results$time) <- c("", "Permute time")

sim_results$cross_over <- factor(sim_results$cross_over)
levels(sim_results$cross_over) <- c("Do not permute crossover time", "Permute crossover time")

sim_results$quasi <- factor(paste(sim_results$ward, sim_results$time))
levels(sim_results$quasi) <- c("", "Permute time", "Permute ward", "Permute time and ward")

sim_results$permute <- factor(paste(sim_results$cross_over, sim_results$ward, sim_results$time))
levels(sim_results$permute) <- c("Nothing", "Time", "Ward", "Time & ward",
                                 "Crossover (randomization)", "Crossover & time", "Crossover & ward", "Crossover, time & ward")
sim_results$permute <- factor(sim_results$permute,
                              c("Nothing", "Crossover (randomization)", "Time", "Ward", "Time & ward",
                                "Crossover & time", "Crossover & ward", "Crossover, time & ward"))

## Quasi-Randomization distributions

library(stringr)

sim_results <- data.table(sim_results)
sim_results <- sim_results[order(permute)]
sim_results[, pretty_permute := str_wrap(permute, width = 10)]
sim_results$pretty_permute <- factor(sim_results$pretty_permute, unique(sim_results$pretty_permute))

pl1 <- ggplot(subset(sim_results, beta == 0 & permute != "Nothing")) + aes(T1) + geom_freqpoly(binwidth = 0.01) + facet_grid(pretty_permute ~ ., switch = "y") + geom_vline(xintercept = T_obs[1], color = "red", linetype = "dashed") + theme_classic(base_size = 15) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), strip.text.y.left = element_text(angle = 0), strip.background.y = element_blank()) + ylab("") + xlab(expression(T[1] ~ "(adjust for nothing)"))
pl1

pl2 <- ggplot(subset(sim_results, beta == 0 & permute != "Nothing")) + aes(T2) + geom_freqpoly(binwidth = 0.01) + facet_grid(pretty_permute ~ .) + geom_vline(xintercept = T_obs[2], color = "red", linetype = "dashed") + theme_classic(base_size = 15) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), strip.text.y.right = element_blank(), strip.background.y = element_blank()) + ylab("") + xlab(expression(T[2] ~ "(adjust for ward)"))
pl2

pl3 <- ggplot(subset(sim_results, beta == 0 & permute != "Nothing")) + aes(T3) + geom_freqpoly(binwidth = 0.01) + facet_grid(pretty_permute ~ .) + geom_vline(xintercept = T_obs[3], color = "red", linetype = "dashed") + theme_classic(base_size = 15) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), strip.text.y.right = element_blank(), strip.background.y = element_blank()) + ylab("") + xlab(expression(T[3] ~ "(adjust for time & ward)"))
pl3


library(cowplot)
pdf("random-dist-new.pdf", width = 10, height = 7)
plot_grid(pl1, pl2, pl3, nrow = 1, rel_widths = c(1.4, 1, 1))
dev.off()


## p-values and confidence intervals
sim_results <- data.table(sim_results)
sim_results <- sim_results[permute != "Nothing"]

pval_pos <- sim_results[,
                    .(p1 = mean(T1 >= T_obs[1] - beta),
                      p2 = mean(T2 >= T_obs[2] - beta),
                      p3 = mean(T3 >= T_obs[3] - beta)),
                    c("beta", "permute")]

pval_neg <- sim_results[,
                    .(p1 = mean(T1 <= T_obs[1] - beta),
                      p2 = mean(T2 <= T_obs[2] - beta),
                      p3 = mean(T3 <= T_obs[3] - beta)),
                    c("beta", "permute")]

cl <- pval_pos[, .(cl1 = max(beta[p1 <= 0.05]),
                   cl2 = max(beta[p2 <= 0.05]),
                   cl3 = max(beta[p3 <= 0.05])
                   ), c("permute")]

cu <- pval_neg[, .(cu1 = min(beta[p1 <= 0.05]),
                   cu2 = min(beta[p2 <= 0.05]),
                   cu3 = min(beta[p3 <= 0.05])
                   ), c("permute")]

results_summary <- merge(merge(pval_pos[beta == 0], cl), cu)

results_summary$permute <- factor(results_summary$permute)

results_summary[, ci1 := paste0("[", round(cl1, 2), ",", round(cu1, 2), "]")]
results_summary[, ci2 := paste0("[", round(cl2, 2), ",", round(cu2, 2), "]")]
results_summary[, ci3 := paste0("[", round(cl3, 2), ",", round(cu3, 2), "]")]
results_summary[, p1r := round(p1, 4)]
results_summary[, p2r := round(p2, 4)]
results_summary[, p3r := round(p3, 4)]

results_summary[, pretty_permute := str_wrap(permute, width = 10)]

library(tables)
booktabs()
toLatex(tabular(Factor(pretty_permute) ~ (p1r + ci1 + p2r + ci2 + p3r + ci3) * Heading() * identity, data = results_summary))

