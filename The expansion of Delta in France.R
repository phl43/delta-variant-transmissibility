library(tidyverse)
library(zoo)
library(lubridate)
library(modelsummary)
library(gt)
library(purrr)
library(rstan)
library(brms)

#########################################################################################################
#                                          FUNCTIONS                                                    #
#########################################################################################################

add_var_student <- custom_family(
  "add_var_student", dpars = c("mu", "sigma", "nu", "alpha"),
  links = c("log", "identity", "identity", "identity"),
  lb = c(NA, 0, 1, NA),
  type = "real",
  vars = "vreal1[n]"
)

stan_funs <- "
  real add_var_student_lpdf(real y, real mu, real sigma, real nu, real alpha, real f) {
    real combined_mu = (1 + alpha * f) * mu;
    return student_t_lpdf(y | nu, combined_mu, sigma);
  }
"
stanvars <- stanvar(block = "functions", scode = stan_funs)

fit_model <- function(formula, data, ...) {
  brm(
    formula = formula,
    data = data,
    family = add_var_student,
    stanvars = stanvars,
    chains = 4,
    warmup = 1000,
    iter = 3000,
    cores = parallel::detectCores(),
    control = list(adapt_delta = 0.95),
    seed = 371,
    ...
  )
}

fit_models <- function(data) {
  dynamic_data <- data
  static_data <- map(
    unique(data$period),
    function(p) {
      data %>%
        filter(period == p)
    }
  )
  
  options(mc.cores = parallel::detectCores())
  
  priors <- c(
    prior(gamma(2, 0.1), class = nu),
    prior(normal(0, 1), class = alpha),
    prior(student_t(3, 0, 0.5), class = sigma)
  )
  
  dynamic_results <- list(
    list(
      model = "Base model",
      fit = fit_model(
        R_total | vreal(prevalence_L452R) ~ 1,
        dynamic_data,
        prior = priors
      )
    ),
    list(
      model = "Model with department-level intercept",
      fit = fit_model(
        R_total | vreal(prevalence_L452R) ~ (1 | dep),
        dynamic_data,
        prior = priors
      )
    ),
    list(
      model = "Model with national time-varying component",
      fit = fit_model(
        R_total | vreal(prevalence_L452R) ~ period,
        dynamic_data,
        prior = c(
          priors,
          prior(student_t(3, 0, 0.5), class = "b")
        )
      )
    ),
    list(
      model = "Model with department-level intercept and\nnational time-varying component",
      fit = fit_model(
        R_total | vreal(prevalence_L452R) ~ (1 | dep) + period,
        dynamic_data,
        prior = c(
          priors,
          prior(student_t(3, 0, 0.5), class = "b")
        )
      )
    )
  )
  
  static_results <- map(
    unique(data$period),
    function(period) {
      list(
        period = period,
        fit = fit_model(
          R_total | vreal(prevalence_L452R) ~ 1,
          static_data[[period]],
          prior = priors
        )
      )
    }
  )
  
  list(
    dynamic = dynamic_results,
    static = static_results
  )
}

summarize_model_fit <- function(model_fit) {
  posterior_samples(model_fit) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    group_by(parameter) %>%
    summarize(
      mean = mean(value),
      lower = quantile(value, 0.05),
      upper = quantile(value, 0.95)
    )
}

summarize_results_dynamic <- function(model_fits) {
  map(
    model_fits,
    function(m) {
      summarize_model_fit(m$fit) %>%
        mutate(model = m$model)
    }
  ) %>%
    bind_rows()
}

summarize_results_static <- function(model_fits) {
  map(
    model_fits,
    function(m) {
      summarize_model_fit(m$fit) %>%
        mutate(period = m$period)
    }
  ) %>%
    bind_rows()
}

#########################################################################################################
#                                       DATA PREPARATION                                                #
#########################################################################################################

# I only include data from metropolitan France and exclude overseas territories
overseas_departments <- c(
  "971",
  "972",
  "973",
  "974",
  "975",
  "976",
  "977",
  "978"
)
overseas_regions <- c(
  "01",
  "02",
  "03",
  "04"
)

# source: https://www.data.gouv.fr/en/datasets/donnees-relatives-aux-resultats-des-tests-virologiques-covid-19/
url_department_test_data <- "https://www.data.gouv.fr/en/datasets/r/406c6a23-e283-4300-9484-54e78c8ae675"

incidence_by_department <- read_delim(url_department_test_data, delim = ";") %>%
  filter(cl_age90 == "0" & !(dep %in% overseas_departments)) %>%
  mutate(
    date = ymd(jour),
    cases = P
  ) %>%
  group_by(dep, date) %>%
  mutate(
    smoothed_cases = rollmean(cases, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  ungroup() %>%
  select(
    dep,
    date,
    cases,
    smoothed_cases
  )

national_incidence <- incidence_by_department %>%
  group_by(date) %>%
  summarize(
    cases = sum(cases)
  ) %>%
  mutate(
    smoothed_cases = rollmean(cases, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  select(
    date,
    cases,
    smoothed_cases
  )

# source: https://www.data.gouv.fr/en/datasets/donnees-de-laboratoires-pour-le-depistage-indicateurs-sur-les-mutations/
url_department_mutations_data <- "https://www.data.gouv.fr/en/datasets/r/4d3e5a8b-9649-4c41-86ec-5420eb6b530c"

prevalence_mutations_by_department <- read_delim(url_department_mutations_data, delim = ";") %>%
  filter(!(dep %in% overseas_departments)) %>%
  mutate(
    # the prevalence of B.1.1.7 in the dataset is the average for d to d + 7, so
    # I assign the estimate for a given week to d + 4
    date = ymd(str_extract(semaine, "^[:digit:]{4}-[:digit:]{2}-[:digit:]{2}")) + 3,
    nb_crib_E484K = nb_A0 + nb_A1,
    prevalence_E484K = ifelse(nb_crib_E484K > 0, nb_A1 / nb_crib_E484K, NA),
    prevalence_non_E484K = 1 - prevalence_E484K,
    nb_crib_E484Q = nb_B0 + nb_B1,
    prevalence_E484Q = ifelse(nb_crib_E484Q > 0, nb_B1 / nb_crib_E484Q, NA),
    prevalence_non_E484Q = 1 - prevalence_E484Q,
    nb_crib_L452R = nb_C0 + nb_C1,
    prevalence_L452R = ifelse(nb_crib_L452R > 0, nb_C1 / nb_crib_L452R, NA),
    prevalence_non_L452R = 1 - prevalence_L452R
  ) %>%
  select(
    dep,
    date,
    prevalence_E484K,
    prevalence_non_E484K,
    prevalence_E484Q,
    prevalence_non_E484Q,
    prevalence_L452R,
    prevalence_non_L452R,
    nb_crib_E484K,
    nb_crib_E484Q,
    nb_crib_L452R
  )

national_prevalence_mutations <- prevalence_mutations_by_department %>%
  group_by(date) %>%
  summarize(
    prevalence_E484K = sum(prevalence_E484K * nb_crib_E484K, na.rm = TRUE) / sum(nb_crib_E484K, na.rm = TRUE),
    prevalence_non_E484K = sum(prevalence_non_E484K * nb_crib_E484K, na.rm = TRUE) / sum(nb_crib_E484K, na.rm = TRUE),
    prevalence_E484Q = sum(prevalence_E484Q * nb_crib_E484Q, na.rm = TRUE) / sum(nb_crib_E484Q, na.rm = TRUE),
    prevalence_non_E484Q = sum(prevalence_non_E484Q * nb_crib_E484Q, na.rm = TRUE) / sum(nb_crib_E484Q, na.rm = TRUE),
    prevalence_L452R = sum(prevalence_L452R * nb_crib_L452R, na.rm = TRUE) / sum(nb_crib_L452R, na.rm = TRUE),
    prevalence_non_L452R = sum(prevalence_non_L452R * nb_crib_L452R, na.rm = TRUE) / sum(nb_crib_L452R, na.rm = TRUE),
    nb_crib_L452R = sum(nb_crib_L452R, na.rm = TRUE)
  ) %>%
  select(
    date,
    prevalence_E484K,
    prevalence_non_E484K,
    prevalence_E484Q,
    prevalence_non_E484Q,
    prevalence_L452R,
    prevalence_non_L452R,
    nb_crib_L452R
  )

mean_gt <- 4.8

department_combined_data <- incidence_by_department %>%
  group_by(dep) %>%
  mutate(
    weekly_cases_total = rollsum(cases, 7, fill = rep(0, 6), align = "left"),
    weekly_growth_factor_total = weekly_cases_total / lag(weekly_cases_total, 7),
    R_total = weekly_growth_factor_total ^ (mean_gt / 7)
  ) %>%
  inner_join(prevalence_mutations_by_department %>% mutate(date = date - 3), by = c("dep", "date")) %>%
  inner_join(department_vaccination, by = c("dep", "date")) %>%
  inner_join(department_google_mobility, by = c("dep", "date")) %>%
  mutate(
    weekly_cases_E484K = weekly_cases_total * prevalence_E484K,
    weekly_cases_E484Q = weekly_cases_total * prevalence_E484Q,
    weekly_cases_L452R = weekly_cases_total * prevalence_L452R,
    weekly_cases_non_E484K = weekly_cases_total * prevalence_non_E484K,
    weekly_cases_non_E484Q = weekly_cases_total * prevalence_non_E484Q,
    weekly_cases_non_L452R = weekly_cases_total * prevalence_non_L452R,
    weekly_growth_factor_E484K = weekly_cases_E484K / lag(weekly_cases_E484K, 7),
    weekly_growth_factor_E484Q = weekly_cases_E484Q / lag(weekly_cases_E484Q, 7),
    weekly_growth_factor_L452R = weekly_cases_L452R / lag(weekly_cases_L452R, 7),
    weekly_growth_factor_non_E484K = weekly_cases_non_E484K / lag(weekly_cases_non_E484K, 7),
    weekly_growth_factor_non_E484Q = weekly_cases_non_E484Q / lag(weekly_cases_non_E484Q, 7),
    weekly_growth_factor_non_L452R = weekly_cases_non_L452R / lag(weekly_cases_non_L452R, 7),
    R_E484K = weekly_growth_factor_E484K ^ (mean_gt / 7),
    R_E484Q = weekly_growth_factor_E484Q ^ (mean_gt / 7),
    R_L452R = weekly_growth_factor_L452R ^ (mean_gt / 7),
    R_non_E484K = weekly_growth_factor_non_E484K ^ (mean_gt / 7),
    R_non_E484Q = weekly_growth_factor_non_E484Q ^ (mean_gt / 7),
    R_non_L452R = weekly_growth_factor_non_L452R ^ (mean_gt / 7),
    advantage_E484K = R_E484K / R_non_E484K - 1,
    advantage_E484Q = R_E484Q / R_non_E484Q - 1,
    advantage_L452R = R_L452R / R_non_L452R - 1,
    prevalence_E484K = stats::lag(prevalence_E484K, 7),
    prevalence_E484Q = stats::lag(prevalence_E484Q, 7),
    prevalence_L452R = lag(prevalence_L452R, 7)
  ) %>%
  ungroup() %>%
  filter(
    date %in% seq(ymd("2021-05-31"), ymd("2021-09-27"), by = "7 day") &
    nb_crib_L452R >= 30
  ) %>%
  mutate(
    period = paste0("Week ", as.integer(date - ymd("2021-01-04")) / 7, " to week ", as.integer(date - ymd("2021-01-04")) / 7 + 1),
    t = as.integer(date - ymd("2021-01-04")) / 7 - 21
  ) %>%
  na.omit() %>%
  arrange(period) %>%
  select(
    t,
    period,
    dep,
    weekly_cases_L452R,
    weekly_cases_non_L452R,
    weekly_growth_factor_L452R,
    weekly_growth_factor_non_L452R,
    R_total,
    R_L452R,
    R_non_L452R,
    advantage_L452R,
    prevalence_L452R,
    partial_vaccination_rate,
    complete_vaccination_rate,
    retail_and_recreation,
    transit_stations,
    workplaces,
    residential
  )

department_combined_data$period <- fct_relevel(
  department_combined_data$period,
  c(
    "Week 22 to week 23",
    "Week 23 to week 24",
    "Week 24 to week 25",
    "Week 25 to week 26",
    "Week 26 to week 27",
    "Week 27 to week 28",
    "Week 28 to week 29",
    "Week 29 to week 30",
    "Week 30 to week 31",
    "Week 31 to week 32",
    "Week 32 to week 33",
    "Week 33 to week 34",
    "Week 34 to week 35",
    "Week 35 to week 36",
    "Week 36 to week 37",
    "Week 37 to week 38",
    "Week 38 to week 39"
  )
)

national_combined_data <- national_incidence %>%
  mutate(
    weekly_cases_total = rollsum(cases, 7, fill = rep(0, 6), align = "left"),
    weekly_growth_factor_total = weekly_cases_total / lag(weekly_cases_total, 7),
    R_total = weekly_growth_factor_total ^ (mean_gt / 7)
  ) %>%
  inner_join(national_prevalence_mutations %>% mutate(date = date - 3), by = "date") %>%
  mutate(
    weekly_cases_E484K = weekly_cases_total * prevalence_E484K,
    weekly_cases_E484Q = weekly_cases_total * prevalence_E484Q,
    weekly_cases_L452R = weekly_cases_total * prevalence_L452R,
    weekly_cases_non_E484K = weekly_cases_total * prevalence_non_E484K,
    weekly_cases_non_E484Q = weekly_cases_total * prevalence_non_E484Q,
    weekly_cases_non_L452R = weekly_cases_total * prevalence_non_L452R,
    weekly_growth_factor_E484K = weekly_cases_E484K / lag(weekly_cases_E484K, 7),
    weekly_growth_factor_E484Q = weekly_cases_E484Q / lag(weekly_cases_E484Q, 7),
    weekly_growth_factor_L452R = weekly_cases_L452R / lag(weekly_cases_L452R, 7),
    weekly_growth_factor_non_E484K = weekly_cases_non_E484K / lag(weekly_cases_non_E484K, 7),
    weekly_growth_factor_non_E484Q = weekly_cases_non_E484Q / lag(weekly_cases_non_E484Q, 7),
    weekly_growth_factor_non_L452R = weekly_cases_non_L452R / lag(weekly_cases_non_L452R, 7),
    R_E484K = weekly_growth_factor_E484K ^ (mean_gt / 7),
    R_E484Q = weekly_growth_factor_E484Q ^ (mean_gt / 7),
    R_L452R = weekly_growth_factor_L452R ^ (mean_gt / 7),
    R_non_E484K = weekly_growth_factor_non_E484K ^ (mean_gt / 7),
    R_non_E484Q = weekly_growth_factor_non_E484Q ^ (mean_gt / 7),
    R_non_L452R = weekly_growth_factor_non_L452R ^ (mean_gt / 7),
    advantage_E484K = R_E484K / R_non_E484K - 1,
    advantage_E484Q = R_E484Q / R_non_E484Q - 1,
    advantage_L452R = R_L452R / R_non_L452R - 1,
    prevalence_E484K = lag(prevalence_E484K, 7),
    prevalence_E484Q = lag(prevalence_E484Q, 7),
    prevalence_L452R = lag(prevalence_L452R, 7)
  ) %>%
  filter(
    date %in% seq(ymd("2021-05-31"), ymd("2021-09-27"), by = "7 day")
  ) %>%
  mutate(
    period = paste0("Week ", as.integer(date - 7 - ymd("2021-05-31")) / 7 + 22, " to week ", as.integer(date - ymd("2021-05-31")) / 7 + 22)
  ) %>%
  na.omit() %>%
  arrange(period) %>%
  select(
    period,
    weekly_cases_L452R,
    weekly_cases_non_L452R,
    weekly_growth_factor_L452R,
    weekly_growth_factor_non_L452R,
    R_total,
    R_L452R,
    R_L452R,
    advantage_L452R,
    prevalence_L452R,
    nb_crib_L452R
  )

national_combined_data$period <- fct_relevel(
  national_combined_data$period,
  c(
    "Week 22 to week 23",
    "Week 23 to week 24",
    "Week 24 to week 25",
    "Week 25 to week 26",
    "Week 26 to week 27",
    "Week 27 to week 28",
    "Week 28 to week 29",
    "Week 29 to week 30",
    "Week 30 to week 31",
    "Week 31 to week 32",
    "Week 32 to week 33",
    "Week 33 to week 34",
    "Week 34 to week 35",
    "Week 35 to week 36",
    "Week 36 to week 37",
    "Week 37 to week 38",
    "Week 38 to week 39"
  )
)

#########################################################################################################
#                                     DATA ANALYSIS AND PLOTS                                           #
#########################################################################################################

dir.create("Figures")

ggplot(national_combined_data, aes(x = prevalence_L452R, y = R_total)) +
  geom_point(aes(group = period, color = period)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "Effective reproduction number vs. prevalence of Delta at the national level in France",
    subtitle = "(the prevalence of Delta is inferred from the prevalence of the L452R mutation and weekly growth rates are converted into effective reproduction numbers\nby assuming the generation time distribution has a mean of 4.8 days)"
  ) +
  xlab("Prevalence of Delta") +
  ylab("Delta's effective reproduction number") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")

ggsave("Figures/Effective reproduction number vs. prevalence of Delta at the national level in France.png", width = 12, height = 6)

ggplot(department_combined_data, aes(x = prevalence_L452R, y = R_total)) +
  geom_point(aes(group = period, color = period)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  ylim(min(department_combined_data$R_total), 2.5) +
  theme_bw() +
  ggtitle(
    "Effective reproduction number vs. prevalence of Delta at the department level in France",
    subtitle = "(the prevalence of Delta is inferred from the prevalence of the L452R mutation and weekly growth rates are converted into effective reproduction numbers\nby assuming the generation time distribution has a mean of 4.8 days)"
  ) +
  xlab("Prevalence of Delta") +
  ylab("Delta's effective reproduction number") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")

ggsave("Figures/Effective reproduction number vs. prevalence of Delta at the department level in France.png", width = 12, height = 6)

ggplot(department_combined_data, aes(x = prevalence_L452R, y = R_total)) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE) +
  facet_wrap(~ period, ncol = 2) +
  scale_x_continuous(labels = scales::percent) +
  ylim(min(department_combined_data$R_total), 2.5) +
  theme_bw() +
  ggtitle(
    "Effective reproduction number vs. prevalence of Delta at the department level in France",
    subtitle = "(the prevalence of Delta is inferred from the prevalence of the L452R mutation and weekly growth rates are converted into effective reproduction numbers\nby assuming the generation time distribution has a mean of 4.8 days)"
  ) +
  xlab("Prevalence of Delta") +
  ylab("Effective reproduction number") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")

ggsave("Figures/Effective reproduction number vs. prevalence of Delta at the department level by week.png", width = 12, height = 12)

results <- fit_models(department_combined_data)

dynamic_model_summary <- summarize_results_dynamic(results$dynamic) %>%
  filter(parameter == "alpha") %>%
  select(-parameter)

dynamic_model_summary$model <- fct_relevel(
  dynamic_model_summary$model,
  c(
    "Base model",
    "Model with department-level intercept",
    "Model with national time-varying component",
    "Model with department-level intercept and\nnational time-varying component"
  )
)

static_model_summary <- summarize_results_static(results$static) %>%
  filter(parameter == "alpha") %>%
  select(-parameter)

ggplot(dynamic_model_summary, aes(x = model, y = mean)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .1)) +
  scale_y_continuous(labels = scales::percent) +
  ggtitle(
    "Summary of the results of the analysis based on fitting various models similar to that used in Abbott and Funk (2021) to French data at the department level",
    subtitle = "(error bars are 90% credible intervals)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  xlab("") +
  ylab("Estimate of Delta's transmissibility advantage")

ggsave(
  "Figures/Relationship between the prevalence of Delta and the epidemic's reproduction number (dynamic analysis).png",
  width = 18,
  height = 6
)

ggplot(static_model_summary, aes(x = period, y = mean)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .1)) +
  scale_y_continuous(labels = scales::percent) +
  ggtitle(
    "Summary of the results of the analysis based on fitting various models similar to that used in Abbott and Funk (2021) to French data at the department level (time-sliced analysis)",
    subtitle = "(error bars are 90% credible intervals)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  xlab("") +
  ylab("Estimate of Delta's transmissibility advantage")

ggsave(
  "Figures/Relationship between the prevalence of Delta and the epidemic's reproduction number (static analysis).png",
  width = 18,
  height = 6
)

ggplot(national_combined_data, aes(x = prevalence_L452R, y = advantage_L452R)) +
  geom_point(aes(group = period, color = period)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "Delta's transmission advantage relative to other strains vs. prevalence of Delta at the national level in France",
    subtitle = "(the prevalence of Delta is inferred from the prevalence of the L452R mutation and weekly growth rates are converted into effective reproduction numbers\nby assuming the generation time distribution has a mean of 4.8 days)"
  ) +
  xlab("Prevalence of Delta") +
  ylab("Delta's transmission advantage") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")

ggsave("Figures/Delta's transmission advantage vs. prevalence of Delta at the national level in France.png", width = 12, height = 6)

ggplot(department_combined_data, aes(x = prevalence_L452R, y = R_total)) +
  geom_point(aes(group = period, color = period)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  ylim(min(department_combined_data$R_total), 2.5) +
  theme_bw() +
  ggtitle(
    "Effective reproduction number vs. prevalence of Delta at the department level in France",
    subtitle = "(the prevalence of Delta is inferred from the prevalence of the L452R mutation and weekly growth rates are converted into effective reproduction numbers\nby assuming the generation time distribution has a mean of 4.8 days)"
  ) +
  xlab("Prevalence of Delta") +
  ylab("Delta's effective reproduction number") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")

ggsave("Figures/Effective reproduction number vs. prevalence of Delta at the department level in France.png", width = 12, height = 6)

ggplot(department_combined_data, aes(x = prevalence_L452R, y = advantage_L452R)) +
  geom_point(aes(group = period, color = period)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.5, 3)) +
  theme_bw() +
  ggtitle(
    "Delta's transmission advantage vs. prevalence of Delta at the department level in France",
    subtitle = "(the prevalence of Delta is inferred from the prevalence of the L452R mutation and weekly growth rates are converted into effective reproduction numbers\nby assuming the generation time distribution has a mean of 4.8 days)"
  ) +
  xlab("Prevalence of Delta") +
  ylab("Delta's transmission advantage") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")

ggsave("Figures/Delta's transmission advantage vs. prevalence of Delta at the department level in France.png", width = 12, height = 6)
