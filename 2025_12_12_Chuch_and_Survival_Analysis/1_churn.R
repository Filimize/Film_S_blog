# install.packages(c("lubridate", "ggsurvfit", "gtsummary", "tidycmprsk"))
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
# devtools::install_github("zabore/condsurv")
library(condsurv)
library(readr)
library(survival)


data <- read_csv("synthetic_US_churn.csv")
data$sex <- as.factor(data$sex)
data$plan <- factor(data$plan,
                    levels = c("basic", "standard", "premium"), 
                    ordered = TRUE)
# data$churn <- as.factor(data$churn)

data <- data |> mutate(first_signup_date = ymd(first_signup_date),
                       last_active_date = ymd(last_active_date)
                       )
data <- data |> mutate( duration = as.duration(first_signup_date %--% last_active_date)/ ddays(1)
                        )

Surv(data$duration, data$churn)[1:10]

survfit2(Surv(duration, churn) ~ 1, data = data)|> 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall not churn probability"
  ) + 
  add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(duration, churn) ~ 1, data = data), times = 365.25)

survfit(Surv(duration, churn) ~ 1, data = data) |> 
  tbl_survfit(
    times = 365.25,
    label_header = "**1-year not churn (95% CI)**"
  )

survfit(Surv(duration, churn) ~ 1, data = data) |> 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median not churn (95% CI)**"
  )

survdiff(Surv(duration, churn) ~ sex, data = data)
survdiff(Surv(duration, churn) ~ plan, data = data)

coxph(Surv(duration, churn) ~ sex + plan, data = data)
coxph(Surv(duration, churn) ~ plan, data = data)


# Load packages
library(survival)
library(survminer)

# # Cox model (plan only) — Best model
# cox_plan <- coxph(Surv(duration, churn) ~ plan, data = data)
# summary(cox_plan)

survfit2(Surv(duration, churn) ~ plan, data = data)|> 
  ggsurvfit() +
  labs(
    x = "Subscription Days",
    y = "Probability of not churning",
    title = "Survival Probability by Subscription Plan",
    subtitle = "Kaplan–Meier estimates with 95% confidence intervals",
  ) + 
  add_confidence_interval()


summary(survfit(Surv(duration, churn) ~ plan, data = data), times = 365.25)

survfit(Surv(duration, churn) ~ plan, data = data) |> 
  tbl_survfit(
    times = 365.25,
    label_header = "**1-year not churning (95% CI)**"
  )

survfit(Surv(duration, churn) ~ plan, data = data) |> 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median time of not churning (95% CI)**"
  )


######################################################################
# Model Comparison
cox_plan_sex <- coxph(Surv(duration, churn) ~ plan + sex, data = data)

AIC(cox_plan, cox_plan_sex)
BIC(cox_plan, cox_plan_sex)

# Print comparison
cat("AIC comparison:\n")
print(AIC(cox_plan, cox_plan_sex))

