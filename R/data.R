#' AVC: Atrioventricular Canal Repair
#'
#' Survival data for 310 patients who underwent repair of atrioventricular
#' septal defects (congenital heart disease) at the Cleveland Clinic between
#' 1977 and 1993. Exhibits two identifiable hazard phases: an early
#' post-operative risk and a constant late phase.
#'
#' @format A data frame with 310 rows and 11 variables:
#' \describe{
#'   \item{study}{Patient identifier}
#'   \item{status}{NYHA functional class (1--4)}
#'   \item{inc_surg}{Surgical grade of AV valve incompetence}
#'   \item{opmos}{Date of operation (months since January 1967)}
#'   \item{age}{Age at repair (months)}
#'   \item{mal}{Malalignment indicator (0/1)}
#'   \item{com_iv}{Interventricular communication indicator (0/1)}
#'   \item{orifice}{Associated cardiac anomaly indicator (0/1)}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{op_age}{Interaction term: opmos x age}
#' }
#'
#' @source Blackstone, Naftel, and Turner (1986)
#'   \doi{10.1080/01621459.1986.10478314}. Cleveland Clinic Foundation.
#'
#' @examples
#' data(avc)
#' avc <- na.omit(avc)
#'
#' # Kaplan-Meier survival
#' km <- survival::survfit(survival::Surv(int_dead, dead) ~ 1, data = avc)
#' plot(km, xlab = "Months after AVC repair", ylab = "Survival",
#'      main = "AVC: Kaplan-Meier survival estimate")
#'
#' \donttest{
#' # Two-phase hazard fit (early CDF + constant -- what AVC supports)
#' fit <- hazard(
#'   survival::Surv(int_dead, dead) ~ 1, data = avc,
#'   dist = "multiphase",
#'   phases = list(
#'     early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1),
#'     constant = hzr_phase("constant")
#'   ),
#'   fit = TRUE, control = list(n_starts = 5, maxit = 1000)
#' )
#' summary(fit)
#' }
#'
#' @seealso \code{vignette("fitting-hazard-models")},
#'   \code{vignette("prediction-visualization")}
#' @family datasets
"avc"

#' CABGKUL: Primary Isolated Coronary Artery Bypass Grafting (KU Leuven)
#'
#' Survival data for 5,880 patients who underwent primary isolated CABG
#' at KU Leuven, Belgium, between 1971 and July 1987. The simplest dataset
#' structure (intercept-only, right-censored) with large sample size
#' exercising all three temporal hazard phases.
#'
#' @format A data frame with 5880 rows and 2 variables:
#' \describe{
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#' }
#'
#' @source KU Leuven cardiac surgery registry. Primary benchmark dataset for
#'   C binary parity testing.
#'
#' @examples
#' data(cabgkul)
#'
#' # Kaplan-Meier survival
#' km <- survival::survfit(survival::Surv(int_dead, dead) ~ 1, data = cabgkul)
#' plot(km, xlab = "Months after CABG", ylab = "Survival",
#'      main = "CABGKUL: Kaplan-Meier survival (n = 5,880)")
#'
#' \donttest{
#' # Single-phase Weibull fit with parametric overlay
#' fit <- hazard(survival::Surv(int_dead, dead) ~ 1, data = cabgkul,
#'               dist = "weibull", theta = c(mu = 0.10, nu = 1.0), fit = TRUE)
#' t_grid <- seq(0.01, max(cabgkul$int_dead) * 0.9, length.out = 200)
#' surv   <- predict(fit, newdata = data.frame(time = t_grid),
#'                   type = "survival")
#' plot(km, xlab = "Months after CABG", ylab = "Survival",
#'      main = "CABGKUL: Weibull vs. Kaplan-Meier")
#' lines(t_grid, surv, col = "blue", lwd = 2)
#' legend("bottomleft", c("KM", "Weibull"), col = c("black", "blue"),
#'        lty = 1, lwd = c(1, 2))
#' }
#'
#' @seealso \code{vignette("fitting-hazard-models")}
#' @family datasets
"cabgkul"

# `cabgkul` is referenced by name (lazy-loaded dataset) in the package's
# own parity tests; declare it as a known global so static checkers do not
# flag a missing binding.
utils::globalVariables("cabgkul")

#' OMC: Open Mitral Commissurotomy
#'
#' Data for 339 patients who underwent open mitral commissurotomy at the
#' University of Alabama Birmingham. Contains repeated thromboembolic events
#' (up to 3 per patient) with left censoring, exercising the interval
#' censoring likelihood.
#'
#' @format A data frame with 339 rows and 7 variables:
#' \describe{
#'   \item{study}{Patient identifier}
#'   \item{te1}{Indicator for first thromboembolic event}
#'   \item{te2}{Indicator for second thromboembolic event}
#'   \item{te3}{Indicator for third thromboembolic event}
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#'   \item{opdjul}{Operation date (Julian)}
#' }
#'
#' @source University of Alabama Birmingham cardiac surgery registry.
#' @family datasets
"omc"

#' TGA: Transposition of the Great Arteries
#'
#' Survival data for 470 patients who underwent the arterial switch operation
#' for transposition of the great arteries at Boston Children's Hospital and
#' Children's Hospital of Philadelphia. Used for sensitivity analysis and
#' internal validation examples.
#'
#' @format A data frame with 470 rows and 14 variables:
#' \describe{
#'   \item{study}{Patient identifier}
#'   \item{simple}{Simple TGA indicator (0/1)}
#'   \item{dextroin}{D-looped transposition indicator (0/1)}
#'   \item{ca_1rl2c}{Coronary artery pattern indicator}
#'   \item{hyaaproc}{Hybrid approach procedure indicator (0/1)}
#'   \item{no_tca}{No total circulatory arrest indicator (0/1)}
#'   \item{tca_time}{Total circulatory arrest time (minutes)}
#'   \item{age_days}{Age at operation (days)}
#'   \item{arciopyr}{Aortic cross-clamp time per year}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{source}{Source institution (BCH or CHOP)}
#'   \item{ca1_2_l}{Coronary artery configuration (1/2/L)}
#'   \item{opyear}{Year of operation}
#' }
#'
#' @source Boston Children's Hospital and Children's Hospital of Philadelphia.
#' @family datasets
"tga"

#' Valves: Primary Heart Valve Replacement
#'
#' Data for 1,533 patients who underwent primary heart valve replacement.
#' The largest multivariable example dataset with multiple endpoints
#' including death, prosthetic valve endocarditis (PVE), bioprosthesis
#' degeneration, and reoperation.
#'
#' @format A data frame with 1533 rows and 19 variables:
#' \describe{
#'   \item{age_cop}{Age at operation (years)}
#'   \item{nyha}{NYHA functional class (1--4)}
#'   \item{mitral}{Mitral valve position indicator (0/1)}
#'   \item{double_}{Double valve replacement indicator (0/1)}
#'   \item{ao_pinc}{Aortic position, incompetence (0/1)}
#'   \item{black}{Black race indicator (0/1)}
#'   \item{i_path}{Ischemic pathology indicator (0/1)}
#'   \item{nve}{Native valve endocarditis indicator (0/1)}
#'   \item{mechvalv}{Mechanical valve indicator (0/1)}
#'   \item{male}{Male sex indicator (0/1)}
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#'   \item{int_pve}{Follow-up interval to PVE or last contact (months)}
#'   \item{pve}{PVE indicator (1 = PVE, 0 = censored)}
#'   \item{bio}{Bioprosthesis indicator (0/1)}
#'   \item{int_rdg}{Follow-up interval to degeneration or last contact (months)}
#'   \item{reop_dg}{Reoperation for degeneration indicator (0/1)}
#'   \item{int_reop}{Follow-up interval to reoperation or last contact (months)}
#'   \item{reop}{Reoperation indicator (0/1)}
#' }
#'
#' @source Cleveland Clinic Foundation heart valve replacement registry.
#'
#' @examples
#' data(valves)
#' valves_cc <- na.omit(valves)
#'
#' # Kaplan-Meier for two endpoints
#' km_death <- survival::survfit(
#'   survival::Surv(int_dead, dead) ~ 1, data = valves_cc)
#' km_pve <- survival::survfit(
#'   survival::Surv(int_pve, pve) ~ 1, data = valves_cc)
#'
#' plot(km_death, xlab = "Months after valve replacement", ylab = "Survival",
#'      main = "Valves: Death and PVE endpoints")
#' lines(km_pve, col = "red")
#' legend("bottomleft", c("Death", "PVE"), col = c("black", "red"), lty = 1)
#'
#' @seealso \code{vignette("fitting-hazard-models")},
#'   \code{vignette("prediction-visualization")}
#' @family datasets
"valves"

#' US Life Table 2023: All-Interval-Censored SAS Parity Anchor
#'
#' The published NCHS United States life table for 2023, expressed on a
#' synthetic radix of 100,000, as a `PROC HAZARD` job consumes it: one
#' interval-censored row per year of age, weighted by the number of deaths
#' falling in that year. Aggregate published counts only -- no patient-level
#' data and no PHI.
#'
#' This is the reference fixture for `objective = "sas"`. It is the cleanest
#' available anchor for that form for two reasons. Every row is
#' interval-censored, so nothing else dilutes the signal; and every interval
#' is exactly one year wide, so \eqn{\log(u - l) = 0} and the width term
#' switches off, isolating the \eqn{S(u)\,\Delta\Lambda} core. It also settles
#' by itself the rival hypothesis that SAS bakes in a constant width divisor:
#' that reading is off by 593,146 log-likelihood units here.
#'
#' Two rows of the source table -- ages 119--120 and 124--125 -- carry
#' `d_all == 0` and are dropped by the SAS job's own
#' `IF D_ALL=0 THEN DELETE`, leaving 124 of 126. The age grid is therefore not
#' contiguous, which is why the fixture carries explicit `age_l`/`age_u`
#' columns rather than an implied one-year step.
#'
#' @format A data frame with 124 rows and 3 variables:
#' \describe{
#'   \item{age_l}{Lower bound of the age interval, in years (integer)}
#'   \item{age_u}{Upper bound of the age interval, in years (integer);
#'     `age_u - age_l` is exactly 1 on every row}
#'   \item{d_all}{Deaths in the interval on a 100,000 radix. This is the
#'     ICENSOR variable, and it is a *weight*, not an indicator -- a fact the
#'     weighted life table forces and 335 occurrences of
#'     `icensor icens_wt=il_dead` across the SAS corpus corroborate.
#'     Sums to 100000.0125; ranges from 0.2352 to 3620.335}
#' }
#'
#' @source National Center for Health Statistics, United States Life Tables,
#'   2023. Derived from `/studies/general/uslife/table2023` via
#'   `data-raw/make_data.R`; `inst/extdata/uslife2023.csv` is the durable
#'   source, reproducible without a SAS license.
#'
#' @references
#' The SAS reference fit is `distributions/hz.icall.lst` in that study, whose
#' printed log-likelihood is -410414 with `Number of events conserved` =
#' 100000.
#'
#' @examples
#' data(uslife2023)
#'
#' # Every interval is one year wide -- the property that makes this the anchor.
#' stopifnot(all(uslife2023$age_u - uslife2023$age_l == 1))
#'
#' # Deaths by age, on the published 100,000 radix.
#' plot(uslife2023$age_l, uslife2023$d_all, type = "h",
#'      xlab = "Age (years)", ylab = "Deaths per 100,000",
#'      main = "US Life Table 2023: deaths by year of age")
#'
#' @seealso \code{\link{hazard}} for the `objective` argument this fixture
#'   anchors.
#' @family datasets
"uslife2023"
