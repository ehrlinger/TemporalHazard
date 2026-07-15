# TemporalHazard-package.R -- Package-level overview help page (?TemporalHazard)
#
# Provides the front-door documentation: the core multiphase model, the phase
# vocabulary, the SAS/C HAZARD bridge, and a map of the main entry points.
# Modelled on the overview pages of randomForestSRC and ggRandomForests.

#' TemporalHazard: Temporal Parametric Hazard Modeling
#'
#' Native R implementation of the multiphase parametric hazard model of
#' Blackstone, Naftel, and Turner (1986), with a focus on behavioral parity,
#' transparent numerics, and reproducible validation against the original
#' \sQuote{C}/\sQuote{SAS} HAZARD program.  The package fits time-varying hazards
#' as an additive sum of parametric *phases* --- early, constant, and late risk
#' streams --- each carrying its own covariate effects.
#'
#' @section The multiphase model:
#'
#' The total cumulative hazard decomposes additively across \eqn{J} phases:
#'
#' \deqn{H(t \mid \mathbf{x}) = \sum_{j=1}^{J} \mu_j(\mathbf{x}) \,
#'       \Phi_j(t)}
#'
#' where \eqn{\mu_j(\mathbf{x}) = \exp(\alpha_j + \mathbf{x}_j^\top
#' \beta_j)} is the phase-specific log-linear scale (intercept plus
#' covariate effects) and \eqn{\Phi_j(t)} is the temporal shape contributed by
#' phase \eqn{j}, with its own parameters depending on the phase `type` (see the
#' phase vocabulary below).  The instantaneous hazard and survival follow
#' directly:
#'
#' \deqn{h(t \mid \mathbf{x}) = \sum_{j=1}^{J} \mu_j(\mathbf{x}) \,
#'       \varphi_j(t), \qquad
#'       S(t \mid \mathbf{x}) = \exp\!\bigl(-H(t \mid \mathbf{x})\bigr)}
#'
#' with \eqn{\varphi_j = d\Phi_j/dt}.  Each phase is specified with
#' [hzr_phase()] and the model is fit by maximum likelihood with [hazard()].
#'
#' @section Phase vocabulary:
#'
#' \tabular{llll}{
#'   **Phase type** \tab **\eqn{\Phi_j(t)}** \tab **Domain** \tab **Use** \cr
#'   `"cdf"` \tab \eqn{G(t)} \tab \eqn{[0, 1]} \tab Early risk that resolves over time \cr
#'   `"hazard"` \tab \eqn{-\log(1 - G(t))} \tab \eqn{[0, \infty)} \tab Late or aging risk that accumulates \cr
#'   `"g3"` \tab \eqn{G_3(t)} \tab \eqn{[0, \infty)} \tab Unbounded late risk (original C/SAS late phase) \cr
#'   `"constant"` \tab \eqn{t} \tab \eqn{[0, \infty)} \tab Flat background rate (no shape parameters) \cr
#' }
#'
#' Here \eqn{G(t)} is the generalized temporal decomposition CDF computed by
#' [hzr_decompos()], and \eqn{G_3(t)} is the unbounded late-phase intensity from
#' [hzr_decompos_g3()].  See `vignette("mf-mathematical-foundations")` for the
#' full derivation.
#'
#' @section SAS/C HAZARD bridge:
#'
#' The classic three-phase HAZARD model maps directly onto `hzr_phase()` calls:
#'
#' \tabular{lll}{
#'   **HAZARD phase** \tab **Role** \tab **R equivalent** \cr
#'   G1 (early)    \tab Early resolving risk \tab `hzr_phase("cdf", ...)` \cr
#'   G2 (constant) \tab Flat background rate \tab `hzr_phase("constant")` \cr
#'   G3 (late)     \tab Rising late risk     \tab `hzr_phase("g3", ...)` \cr
#' }
#'
#' TemporalHazard generalizes the fixed three-phase structure to \eqn{N} phases
#' of any type.  [hzr_argument_mapping()] gives the full parameter translation
#' table between the SAS/C parameterization and the R arguments.
#'
#' @section Main entry points:
#'
#' \describe{
#'   \item{Model fitting}{[hazard()] --- build and fit single- or multiphase
#'     models; [hzr_phase()] --- specify one phase; [hzr_stepwise()] --- forward,
#'     backward, or bidirectional covariate selection.}
#'   \item{Prediction}{`predict()` on a fitted `hazard` object --- survival,
#'     cumulative hazard, and per-phase decomposed hazard; `summary()` ---
#'     coefficient tables with Wald inference.}
#'   \item{Parametric family}{[hzr_decompos()] --- the early-phase (G1)
#'     decomposition \eqn{G(t)}, \eqn{g(t)}, \eqn{h(t)}; [hzr_decompos_g3()] ---
#'     the late-phase (G3) intensity; [hzr_phase_cumhaz()] and
#'     [hzr_phase_hazard()] --- per-phase \eqn{\Phi(t)} and \eqn{\varphi(t)}.}
#'   \item{Diagnostics}{[hzr_kaplan()], [hzr_nelson()] --- nonparametric
#'     references; [hzr_gof()] --- goodness of fit; [hzr_calibrate()],
#'     [hzr_deciles()] --- calibration; [hzr_bootstrap()] --- resampling CIs;
#'     [hzr_competing_risks()] --- cumulative incidence.}
#' }
#'
#' @section Vignettes:
#'
#' \describe{
#'   \item{`vignette("getting-started")`}{First fit, end to end.}
#'   \item{`vignette("mf-mathematical-foundations")`}{The decomposition family
#'     and multiphase model, with derivations.}
#'   \item{`vignette("fitting-hazard-models")`}{Single-phase through multiphase
#'     fitting.}
#'   \item{`vignette("prediction-visualization")`}{Prediction types and
#'     decomposed-hazard plots.}
#'   \item{`vignette("inference-diagnostics")`}{Bootstrap CIs and diagnostics.}
#'   \item{`vignette("sas-to-r-migration")`}{Translating SAS/C HAZARD code.}
#' }
#'
#' @references
#' Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of time-varying
#' hazard into phases, each incorporating a separate stream of concomitant
#' information. *J Am Stat Assoc.* 1986;81(395):615--624.
#' \doi{10.1080/01621459.1986.10478314}
#'
#' Rajeswaran J, Blackstone EH, Ehrlinger J, Li L, Ishwaran H, Parides MK.
#' Probability of atrial fibrillation after ablation: Using a parametric
#' nonlinear temporal decomposition mixed effects model. *Stat Methods Med Res.*
#' 2018;27(1):126--141. \doi{10.1177/0962280215623583}
#'
#' @keywords internal
"_PACKAGE"
