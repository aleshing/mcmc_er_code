library(multilink)
library(tidyverse)
library(salso)
library(posterior)
library(bayesplot)
library(xtable)
library(ggpubr)

#### Load in helper functions and data ####
source("helper_functions.R")
data(no_dup_data)

#### Set number of iterations for MCMC and warmup #### 
n_iter <- 5000
warmup <- 1000

#### Create comparison data #### 
comparison_list <- create_comparison_data(no_dup_data$records,
                                          types = c("bi", "lv", "lv", "lv",
                                                    "lv", "bi", "bi"),
                                          breaks = list(NA,
                                                        c(0, 0.25, 0.5),
                                                        c(0, 0.25, 0.5),
                                                        c(0, 0.25, 0.5),
                                                        c(0, 0.25, 0.5),
                                                        NA, NA),
                                          file_sizes = no_dup_data$file_sizes,
                                          duplicates = c(0, 0, 0))

#### Specify prior #### 
prior_list <- specify_prior(comparison_list, mus = NA, nus = NA, flat = 0,
                            alphas = rep(1, 7), dup_upper_bound = c(1, 1, 1),
                            dup_count_prior_family = NA,
                            dup_count_prior_pars = NA,
                            n_prior_family = "uniform", n_prior_pars = NA)

#### Initialization #### 
# Specify indexing scheme for initialization
pairs_to_keep_init <-
  (comparison_list$comparisons[, "gname_DL_3"] != TRUE) &
  (comparison_list$comparisons[, "fname_DL_3"] != TRUE)
# Generate 4 partitions for initializing 4 MCMC chains
Z_init <- uniform_index_sample(comparison_list = comparison_list,
                               pairs_to_keep = pairs_to_keep_init,
                               num_samples = 4,
                               duplicates = c(0, 0, 0),
                               seed = 44)

#### Run MCMC #### 
results_extra1 <- vector(mode = "list", length = 4)
results_extra2 <- vector(mode = "list", length = 4)
for(i in 1:4){
  print(i)
  # Gibbs update
  # Twice the number of iterations
  results_extra1[[i]] <- gibbs_sampler(comparison_list, prior_list,
                                       n_iter = n_iter * 2,
                                       Z_init = Z_init[, i],
                                       seed = 43, chaperones_info = NA)
  
  # Chaperones update, nonuniform Chaperones distribution (exact)
  # 5000 updates per Gibbs sampler iteration, 5 restricted Gibbs steps per update
  # Twice the number of iterations
  chaperones_info <- list(chap_type = 1,
                          num_chap_iter = 5000,
                          nonuniform_chap_type = 0,
                          extra_gibbs = FALSE,
                          num_restrict = 5)
  results_extra2[[i]] <- gibbs_sampler(comparison_list, prior_list,
                                       n_iter = n_iter * 2,
                                       Z_init = Z_init[, i],
                                       seed = 43,
                                       chaperones_info = chaperones_info)
}
save(results_extra1, results_extra2,
     file = "illustrative_example_extra_samples.RData")

#### Find reference partitions #### 
# Choose a reference partition: Bayes estimate with respect to Binder's loss
r <- sum(no_dup_data$file_sizes)
for(i in 1:4){
  print(i)
  results_extra1[[i]]$binder <- salso::salso(t(results_extra1[[i]]$partitions)[(warmup + 1):(n_iter * 2), ],
                                             loss = "binder",
                                             maxNClusters = r,
                                             nRuns = 4,
                                             maxZealousAttempts = 0,
                                             nCores = 2,
                                             seconds = 180)
  
  results_extra2[[i]]$binder <- salso::salso(t(results_extra2[[i]]$partitions)[(warmup + 1):(n_iter * 2), ],
                                             loss = "binder",
                                             maxNClusters = r,
                                             nRuns = 4,
                                             maxZealousAttempts = 0,
                                             nCores = 2,
                                             seconds = 180)
}
save(results_extra1, results_extra2,
     file = "illustrative_example_extra_samples.RData")

#### Postprocess MCMC chains ####
summary_chains <- data.frame(.iteration = 1:(n_iter * 2),
                             num_clust = colSums(results_extra1[[1]]$contingency_tables),
                             t(results_extra1[[1]]$contingency_tables),
                             binder_1 =
                               salso::binder(results_extra1[[1]]$binder,
                                             t(results_extra1[[1]]$partitions)),
                             binder_2 =
                               salso::binder(results_extra1[[2]]$binder,
                                             t(results_extra1[[1]]$partitions)),
                             binder_3 =
                               salso::binder(results_extra1[[3]]$binder,
                                             t(results_extra1[[1]]$partitions)),
                             binder_4 =
                               salso::binder(results_extra1[[4]]$binder,
                                             t(results_extra1[[1]]$partitions)),
                             .chain = 1,
                             mcmc = "Gibbs longer")
for(i in 1:4){
  if(i > 1){
    summary_chains <- rbind(summary_chains,
                            data.frame(.iteration = 1:(n_iter * 2),
                                       num_clust = colSums(results_extra1[[1]]$contingency_tables),
                                       t(results_extra1[[1]]$contingency_tables),
                                       binder_1 =
                                         salso::binder(results_extra1[[1]]$binder,
                                                       t(results_extra1[[i]]$partitions)),
                                       binder_2 =
                                         salso::binder(results_extra1[[2]]$binder,
                                                       t(results_extra1[[i]]$partitions)),
                                       binder_3 =
                                         salso::binder(results_extra1[[3]]$binder,
                                                       t(results_extra1[[i]]$partitions)),
                                       binder_4 =
                                         salso::binder(results_extra1[[4]]$binder,
                                                       t(results_extra1[[i]]$partitions)),
                                       .chain = i,
                                       mcmc = "Gibbs longer"))
  }
  
  summary_chains <- rbind(summary_chains,
                          data.frame(.iteration = 1:(n_iter * 2),
                                     num_clust = colSums(results_extra2[[1]]$contingency_tables),
                                     t(results_extra2[[1]]$contingency_tables),
                                     binder_1 =
                                       salso::binder(results_extra2[[1]]$binder,
                                                     t(results_extra2[[i]]$partitions)),
                                     binder_2 =
                                       salso::binder(results_extra2[[2]]$binder,
                                                     t(results_extra2[[i]]$partitions)),
                                     binder_3 =
                                       salso::binder(results_extra2[[3]]$binder,
                                                     t(results_extra2[[i]]$partitions)),
                                     binder_4 =
                                       salso::binder(results_extra2[[4]]$binder,
                                                     t(results_extra2[[i]]$partitions)),
                                     .chain = i,
                                     mcmc = "Nonuniform Chaperones longer"))
}



summary_chains_temp <- as_draws_df(summary_chains %>%
                                     filter(mcmc == "Gibbs longer"))
mcmc_schemes <- unique(summary_chains$mcmc)[-1]
for(i in 1:length(mcmc_schemes)){
  summary_chains_temp <- rbind(summary_chains_temp,
                               as_draws_df(summary_chains %>%
                                             filter(mcmc == mcmc_schemes[i])))
}
summary_chains <- summary_chains_temp

summary_output <- summarize_draws(summary_chains %>%
                                    filter(mcmc == "Gibbs longer") %>%
                                    select(-mcmc)) %>% cbind(mcmc = "Gibbs longer")
mcmc_schemes <- unique(summary_chains$mcmc)[-1]
for(i in 1:length(mcmc_schemes)){
  summary_output <- rbind(summary_output,
                          summarize_draws(summary_chains %>%
                                            filter(mcmc == mcmc_schemes[i]) %>%
                                            select(-mcmc)) %>% cbind(mcmc = mcmc_schemes[i]))
}
save(summary_chains, summary_output,
     file = "illustrative_example_extra_summary.RData")

#### Calculate summaries post-warmup ####
summary_chains_burn <- summary_chains %>% filter(.iteration > warmup)
summary_output_burn <- summarize_draws(summary_chains_burn %>%
                                         filter(mcmc == "Gibbs longer") %>%
                                         select(-mcmc)) %>% cbind(mcmc = "Gibbs longer")
mcmc_schemes <- unique(summary_chains$mcmc)[-1]
for(i in 1:length(mcmc_schemes)){
  summary_output_burn <- rbind(summary_output_burn,
                               summarize_draws(summary_chains_burn %>%
                                                 filter(mcmc == mcmc_schemes[i]) %>%
                                                 select(-mcmc)) %>% cbind(mcmc = mcmc_schemes[i]))
}

#### Gibbs: Get table of results ####
summary_output_burn %>%
  filter(mcmc %in% c("Gibbs longer") &
           !(variable %in% c("X000"))) %>%
  mutate(mcmc  = recode(mcmc,
                        `Gibbs longer` = "Gibbs, 10,000 Iterations"),
         variable  = recode(variable,
                            num_clust = "# of Clusters",
                            binder_1 = "Binder 1",
                            binder_2 = "Binder 2",
                            binder_3 = "Binder 3",
                            binder_4 = "Binder 4")) %>%
  select(variable, mean, q5, q95, rhat, ess_bulk, mcmc) %>%
  rename(Summary = variable, Mean = mean, Lower = q5, Upper = q95,
         Rhat = rhat, ESS = ess_bulk, Sampler = mcmc) %>%
  xtable() %>% print(include.rownames = FALSE)

summary_output_burn %>%
  filter(mcmc %in% c("Gibbs longer") &
           !(variable %in% c("X000"))) %>%
  mutate(mcmc  = recode(mcmc,
                        `Gibbs longer` = "Gibbs, 10,000 Iterations"),
         variable  = recode(variable,
                            num_clust = "# of Clusters",
                            binder_1 = "Binder 1",
                            binder_2 = "Binder 2",
                            binder_3 = "Binder 3",
                            binder_4 = "Binder 4")) %>%
  select(variable, mean, q5, q95, rhat, ess_bulk, mcmc) %>%
  rename(Summary = variable, Mean = mean, Lower = q5, Upper = q95,
         Rhat = rhat, ESS = ess_bulk, Sampler = mcmc) %>%
  xtable(digits = -2) %>% print(include.rownames = FALSE)

#### Chaperones: Get table of results ####
summary_output_burn %>%
  filter(mcmc %in% c("Nonuniform Chaperones longer") &
           !(variable %in% c("X000"))) %>%
  mutate(mcmc  = recode(mcmc,
                        `Nonuniform Chaperones longer` = "Chaperones, 10,000 Iterations"),
         variable  = recode(variable,
                            num_clust = "# of Clusters",
                            binder_1 = "Binder 1",
                            binder_2 = "Binder 2",
                            binder_3 = "Binder 3",
                            binder_4 = "Binder 4")) %>%
  select(variable, mean, q5, q95, rhat, ess_bulk, mcmc) %>%
  rename(Summary = variable, Mean = mean, Lower = q5, Upper = q95,
         Rhat = rhat, ESS = ess_bulk, Sampler = mcmc) %>%
  xtable() %>% print(include.rownames = FALSE)

summary_output_burn %>%
  filter(mcmc %in% c("Nonuniform Chaperones longer") &
           !(variable %in% c("X000"))) %>%
  mutate(mcmc  = recode(mcmc,
                        `Nonuniform Chaperones longer` = "Chaperones, 10,000 Iterations"),
         variable  = recode(variable,
                            num_clust = "# of Clusters",
                            binder_1 = "Binder 1",
                            binder_2 = "Binder 2",
                            binder_3 = "Binder 3",
                            binder_4 = "Binder 4")) %>%
  select(variable, mean, q5, q95, rhat, ess_bulk, mcmc) %>%
  rename(Summary = variable, Mean = mean, Lower = q5, Upper = q95,
         Rhat = rhat, ESS = ess_bulk, Sampler = mcmc) %>%
  xtable(digits = -2) %>% print(include.rownames = FALSE)

#### Chaperones: ESS over time plots ####
ess_plot <- get_ess(summary_chains_burn %>%
                      filter(mcmc == "Gibbs longer") %>%
                      select(-mcmc)) %>% cbind(mcmc = "Gibbs longer")
mcmc_schemes <- unique(summary_chains$mcmc)[-1]
for(i in 1:length(mcmc_schemes)){
  ess_plot <- rbind(ess_plot,
                    get_ess(summary_chains_burn %>%
                              filter(mcmc == mcmc_schemes[i]) %>%
                              select(-mcmc)) %>% cbind(mcmc = mcmc_schemes[i]))
}
ess_plot <- ess_plot %>% mutate(time = time - warmup,
                                time = time * 4, 
                                ess_time = ess / time)

chap_ess1 <- ess_plot %>% filter(mcmc == "Nonuniform Chaperones longer") %>% 
  filter(variable %in% c("num_clust", "binder_1")) %>%
  mutate(variable  = recode(variable,
                            num_clust = "# of Clusters",
                            binder_1 = "Binder 1",
                            binder_2 = "Binder 2",
                            binder_3 = "Binder 3",
                            binder_4 = "Binder 4"),
         Summary = variable) %>%
  ggplot(aes(x = time, y = ess, col = Summary)) + geom_line(size = 2) + 
  geom_hline(yintercept = 400, linetype = "dashed", size = 2) +
  xlab("Iteration") + ylab("ESS") + 
  theme_bw() + theme(text = element_text(size = 25),
                     axis.text.x = element_text(size = 12),
                     axis.text.y = element_text(size = 20))

chap_ess2 <- ess_plot %>% filter(mcmc == "Nonuniform Chaperones longer") %>% 
  filter(variable %in% c("num_clust", "binder_1")) %>%
  mutate(variable  = recode(variable,
                            num_clust = "# of Clusters",
                            binder_1 = "Reference 1",
                            binder_2 = "Reference 2",
                            binder_3 = "Reference 3",
                            binder_4 = "Reference 4"),
         Summary = variable) %>%
  ggplot(aes(x = time, y = ess_time, col = Summary)) + geom_line(size = 2) +
  xlab("Iteration") + ylab("ESS / Iteration") + 
  theme_bw() + theme(text = element_text(size = 25),
                     axis.text.x = element_text(size = 12),
                     axis.text.y = element_text(size = 20))

ggarrange(chap_ess1, chap_ess2, ncol = 2, common.legend = T)
ggsave("plots/chap_ess_extra.pdf", width = 12, height = 6)
