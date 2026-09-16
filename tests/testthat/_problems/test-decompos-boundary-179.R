# Extracted from test-decompos-boundary.R:179

# test -------------------------------------------------------------------------
tau <- exp(2.6238)
gamma <- 220.08
eta <- 0.0090875
t <- c(0.1, 0.3, 0.5, 0.55)
x <- gamma * log(t / tau)
expect_true(all(x < -700))
alpha <- 1.2174
d <- hzr_decompos_g3(t, tau, gamma, alpha, eta)
ln_scale <- log(gamma) + log(eta) - log(tau) - log(alpha)
expect_equal(log(d$G3), eta * (x - log(alpha)), tolerance = 1e-12)
expect_equal(
    log(d$g3),
    ln_scale + (eta - 1) * (x - log(alpha)) + (gamma - 1) * log(t / tau),
    tolerance = 1e-12
  )
expect_true(all(diff(d$G3) > 0))
t_near <- 10 * exp(-1 / 220)
d_big <- hzr_decompos_g3(t_near, 10, 220, 1e308, 0.5)
expect_equal(log(d_big$G3),
               0.5 * (log(log1p(exp(-1))) - log(1e308)),
               tolerance = 1e-12)
