
data("ckbtools_baseline")
data("ckbtools_endpoints")
data("ckbtools_ascertainments")

baseline <- ckbtools_baseline
endpoints <- ckbtools_endpoints
ascertainments <- ckbtools_ascertainments

baseline$sex <- factor(baseline$is_female,
                       levels = c(0, 1),
                       labels = c("Male", "Female"))

baseline$age_at_study_date <- baseline$age_at_study_date_x100 / 100

baseline$age_at_baseline_group <- cut(baseline$age_at_study_date,
                                      breaks = seq(30, 80, by = 10),
                                      right = FALSE)

baseline$high_school_or_higher <- baseline$highest_education %in% c(3, 4, 5)

baseline$education <- factor(baseline$highest_education,
                             levels = c(0, 1, 2, 3, 4, 5),
                             labels = c("No formal school",
                                        "Primary School",
                                        "Middle School",
                                        "High School",
                                        "Technical School / college",
                                        "University"))

baseline$education_4_groups <- factor(baseline$highest_education,
                                      levels = c(0, 1, 2, 3, 4, 5),
                                      labels = c("No formal school",
                                                 "Primary School",
                                                 "Middle School",
                                                 "High School or higher",
                                                 "High School or higher",
                                                 "High School or higher"))

baseline$smoking <- factor(baseline$smoking_category,
                           levels = c(1, 2, 3, 4),
                           labels = c("Never smoker",
                                      "Occasional smoker",
                                      "Ex regular smoker",
                                      "Smoker"))

baseline$smoking_3_groups <- factor(baseline$smoking_category,
                                    levels = c(1, 2, 3, 4),
                                    labels = c("Never / Occasional smoker",
                                               "Never / Occasional smoker",
                                               "Ex regular smoker",
                                               "Smoker"))

baseline$current_regular_smoker <-  factor(baseline$smoking_category,
                                           levels = c(1, 2, 3, 4),
                                           labels = c("No",
                                                      "No",
                                                      "No",
                                                      "Yes"))

baseline$alcohol <- factor(baseline$alcohol_category,
                           levels = c(1, 2, 3, 4, 5, 6),
                           labels = c("Never regular",
                                      "Ex-regular",
                                      "Occasional",
                                      "Monthly",
                                      "Reduced intake",
                                      "Weekly"))

baseline$alcohol_3_groups <- factor(baseline$alcohol_category,
                                    levels = c(1, 2, 3, 4, 5, 6),
                                    labels = c("Never regular drinker",
                                               "Weekly / Reduced / Ex",
                                               "Occasional / Monthly",
                                               "Occasional / Monthly",
                                               "Weekly / Reduced / Ex",
                                               "Weekly / Reduced / Ex"))

baseline$bmi_grp <- cut(baseline$bmi_calc,
                        c(0, 18.5, 25, 30, Inf),
                        labels = c("<18.5", "18.5-25", "25-30", "30+"),
                        right = FALSE)

baseline$hours_since_last_ate <- baseline$hours_since_last_ate_x10 / 10

baseline$region <- factor(baseline$region_code,
                          levels = c(12, 16, 26, 36, 46, 52, 58, 68, 78, 88),
                          labels = c("Qingdao", "Harbin", "Haikou", "Suzhou", "Liuzhou",
                                     "Sichuan", "Gansu", "Henan", "Zhejiang", "Hunan"))

endpoints$ihd <- endpoints$CKB0003_ep
endpoints$ihd_date <- endpoints$CKB0003_datedeveloped

endpoints$respiratory_diseases <- endpoints$CKB0032_ep
endpoints$respiratory_diseases_date <- endpoints$CKB0032_datedeveloped

endpoints <- endpoints[, c("csid", "ihd", "ihd_date", "respiratory_diseases", "respiratory_diseases_date")]

baseline$all_cause_mortality_date <- baseline$censoring_date
baseline$all_cause_mortality <- as.numeric(grepl("Dead", baseline$censoring_reason))

ckbtools_participant_data <- merge(baseline, endpoints, by = "csid")

csid_in_subcohort <- ascertainments$csid[ascertainments$olinkexpexpan_chd_b1_subcohort == 1]
ckbtools_participant_data$proteomics_subcohort <- as.numeric(ckbtools_participant_data$csid %in% csid_in_subcohort)


usethis::use_data(ckbtools_participant_data, overwrite = TRUE)
