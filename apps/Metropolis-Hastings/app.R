library(shiny)

log_sum_exp <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

target_spec <- function(name) {
  if (name == "Standard normal") {
    list(
      dim = 1,
      log_density = function(z) dnorm(z[1], log = TRUE),
      density_1d = function(x) dnorm(x),
      limits = c(-4.5, 4.5),
      default_start = c(3, 0),
      default_scale = 1,
      proposal_name = "Symmetric Gaussian random walk",
      proposal_uses_scale = TRUE,
      propose = function(current, scale) {
        candidate <- current + rnorm(1, sd = scale)
        list(candidate = candidate, log_q_forward = 0, log_q_reverse = 0)
      },
      description = "A smooth unimodal target. This is the basic tuning example: proposals that are too small move slowly, while proposals that are too large are often rejected."
    )
  } else if (name == "Bimodal normal mixture") {
    list(
      dim = 1,
      log_density = function(z) {
        a <- log(0.5) + dnorm(z[1], -4, 1, log = TRUE)
        b <- log(0.5) + dnorm(z[1], 4, 1, log = TRUE)
        log_sum_exp(a, b)
      },
      density_1d = function(x) 0.5 * dnorm(x, -4, 1) + 0.5 * dnorm(x, 4, 1),
      limits = c(-8, 8),
      default_start = c(-4, 0),
      default_scale = 1,
      proposal_name = "Symmetric Gaussian random walk",
      proposal_uses_scale = TRUE,
      propose = function(current, scale) {
        candidate <- current + rnorm(1, sd = scale)
        list(candidate = candidate, log_q_forward = 0, log_q_reverse = 0)
      },
      description = "Two well-separated modes expose a major limitation of local random-walk proposals: a chain may explore one mode well but rarely cross to the other."
    )
  } else if (name == "Correlated bivariate normal") {
    rho <- 0.97
    list(
      dim = 2,
      log_density = function(z) {
        -(z[1]^2 - 2 * rho * z[1] * z[2] + z[2]^2) / (2 * (1 - rho^2))
      },
      limits = c(-4, 4),
      default_start = c(3, -3),
      default_scale = 0.5,
      proposal_name = "Symmetric isotropic Gaussian random walk",
      proposal_uses_scale = TRUE,
      propose = function(current, scale) {
        candidate <- current + rnorm(2, sd = scale)
        list(candidate = candidate, log_q_forward = 0, log_q_reverse = 0)
      },
      contour = function(x, y) {
        exp(-(x^2 - 2 * rho * x * y + y^2) / (2 * (1 - rho^2)))
      },
      description = "The target lies along a narrow diagonal ridge. An isotropic proposal does not match that geometry, producing high autocorrelation even when the trace appears active."
    )
  } else if (name == "Banana-shaped target") {
    bend <- 0.18
    list(
      dim = 2,
      log_density = function(z) {
        transformed_y <- z[2] - bend * (z[1]^2 - 4)
        dnorm(z[1], 0, 2, log = TRUE) + dnorm(transformed_y, 0, 0.7, log = TRUE)
      },
      limits = c(-6, 6),
      default_start = c(-5, 4),
      default_scale = 0.8,
      proposal_name = "Symmetric isotropic Gaussian random walk",
      proposal_uses_scale = TRUE,
      propose = function(current, scale) {
        candidate <- current + rnorm(2, sd = scale)
        list(candidate = candidate, log_q_forward = 0, log_q_reverse = 0)
      },
      contour = function(x, y) {
        dnorm(x, 0, 2) * dnorm(y - bend * (x^2 - 4), 0, 0.7)
      },
      description = "A curved target shows why a single global random-walk scale can struggle when the posterior geometry changes with position."
    )
  } else if (name == "Gamma target: asymmetric log-normal walk") {
    list(
      dim = 1,
      log_density = function(z) {
        if (z[1] <= 0) return(-Inf)
        dgamma(z[1], shape = 3, rate = 1, log = TRUE)
      },
      density_1d = function(x) dgamma(x, shape = 3, rate = 1),
      limits = c(0, 12),
      default_start = c(0.4, 0),
      default_scale = 0.55,
      proposal_name = "Asymmetric multiplicative log-normal random walk",
      proposal_uses_scale = TRUE,
      propose = function(current, scale) {
        candidate <- rlnorm(1, meanlog = log(current[1]), sdlog = scale)
        list(
          candidate = candidate,
          log_q_forward = dlnorm(candidate, meanlog = log(current[1]), sdlog = scale, log = TRUE),
          log_q_reverse = dlnorm(current[1], meanlog = log(candidate), sdlog = scale, log = TRUE)
        )
      },
      description = "A multiplicative proposal stays positive but is asymmetric: q(y|x) is not q(x|y). Omitting the Hastings correction gives the wrong stationary distribution."
    )
  } else {
    list(
      dim = 1,
      log_density = function(z) {
        if (z[1] <= 0 || z[1] >= 1) return(-Inf)
        dbeta(z[1], 5, 2, log = TRUE)
      },
      density_1d = function(x) dbeta(x, 5, 2),
      limits = c(0, 1),
      default_start = c(0.5, 0),
      default_scale = 1,
      proposal_name = "Independence proposal Beta(2,2)",
      proposal_uses_scale = FALSE,
      propose = function(current, scale) {
        candidate <- rbeta(1, 2, 2)
        list(
          candidate = candidate,
          log_q_forward = dbeta(candidate, 2, 2, log = TRUE),
          log_q_reverse = dbeta(current[1], 2, 2, log = TRUE)
        )
      },
      description = "The candidate is proposed independently of the current state. Because the proposal is not uniform, its density appears explicitly in the Hastings ratio."
    )
  }
}

run_mh <- function(spec, initial, proposal_scale, iterations) {
  d <- length(initial)
  chain <- matrix(NA_real_, nrow = iterations, ncol = d)
  accepted <- logical(iterations)
  current <- initial
  current_log_density <- spec$log_density(current)
  chain[1, ] <- current

  if (iterations >= 2) {
    for (i in 2:iterations) {
      proposal <- spec$propose(current, proposal_scale)
      candidate <- proposal$candidate
      candidate_log_density <- spec$log_density(candidate)
      log_alpha <- min(0, candidate_log_density - current_log_density +
                         proposal$log_q_reverse - proposal$log_q_forward)
      if (log(runif(1)) < log_alpha) {
        current <- candidate
        current_log_density <- candidate_log_density
        accepted[i] <- TRUE
      }
      chain[i, ] <- current
    }
  }

  colnames(chain) <- paste0("X", seq_len(d))
  list(chain = chain, accepted = accepted)
}

effective_sample_size <- function(x) {
  n <- length(x)
  if (n < 4 || isTRUE(all.equal(var(x), 0))) return(1)
  max_lag <- min(n - 1, max(20, floor(10 * log10(n))))
  rho <- as.numeric(acf(x, lag.max = max_lag, plot = FALSE, demean = TRUE)$acf)[-1]
  if (!length(rho)) return(n)

  pair_sums <- numeric()
  if (length(rho) >= 2) {
    pair_sums <- rho[seq(1, length(rho) - 1, by = 2)] +
      rho[seq(2, length(rho), by = 2)]
    first_nonpositive <- which(pair_sums <= 0)[1]
    if (!is.na(first_nonpositive)) pair_sums <- pair_sums[seq_len(first_nonpositive - 1)]
  }
  tau <- if (length(pair_sums)) 1 + 2 * sum(pair_sums) else 1
  min(n, max(1, n / max(tau, 1)))
}

ui <- fluidPage(
  tags$head(tags$style(HTML("body{background:#f5f7fb;color:#172033}.hero{background:linear-gradient(120deg,#083344,#0f766e);color:white;border-radius:14px;padding:22px 28px;margin:18px 0}.well{background:white;border:0;box-shadow:0 2px 12px #dfe3ec}.metric{font-size:27px;font-weight:750;color:#0f766e}.note{color:#5f6b7a}.warning{color:#b42318;font-weight:650}.good{color:#137333;font-weight:650}.btn-primary{background:#0f766e;border-color:#0f766e}.tab-content{padding-top:18px}"))),
  div(class = "hero",
      h2("Metropolis–Hastings sampler"),
      p("Proposal, accept or reject, and learn from a dependent sample")),
  sidebarLayout(
    sidebarPanel(width = 4,
      selectInput("target", "Target distribution",
                  c("Standard normal", "Bimodal normal mixture",
                    "Correlated bivariate normal", "Banana-shaped target",
                    "Gamma target: asymmetric log-normal walk",
                    "Beta target: independence sampler")),
      sliderInput("iterations", "Number of iterations", 500, 20000, 5000, step = 500),
      sliderInput("burn", "Burn-in", 0, 5000, 500, step = 100),
      conditionalPanel(
        condition = "input.target != 'Beta target: independence sampler'",
        sliderInput("proposal_scale", "Proposal scale", 0.05, 8, 1, step = 0.05)
      ),
      numericInput("start_x", "Initial X₁", 3, step = 0.5),
      conditionalPanel(
        condition = "input.target == 'Correlated bivariate normal' || input.target == 'Banana-shaped target'",
        numericInput("start_y", "Initial X₂", 0, step = 0.5)
      ),
      actionButton("run", "Run chain", class = "btn-primary", width = "100%"),
      hr(),
      p(class = "note", "Change the proposal scale and rerun. A high acceptance rate is not automatically good: many tiny accepted moves can still yield a very small ESS."),
      wellPanel(h5("Proposal"), textOutput("proposal_name")),
      uiOutput("target_note")
    ),
    mainPanel(width = 8,
      fluidRow(
        column(4, wellPanel(h4("Acceptance rate"), div(class = "metric", textOutput("acceptance")))),
        column(4, wellPanel(h4("Post-burn-in draws"), div(class = "metric", textOutput("retained")))),
        column(4, wellPanel(h4("ESS"), div(class = "metric", textOutput("ess_metric"))))
      ),
      tabsetPanel(
        tabPanel("Trace",
          plotOutput("trace_plot", height = 430),
          uiOutput("mixing_message")
        ),
        tabPanel("Autocorrelation",
          selectInput("acf_coordinate", "Coordinate", c("X1", "X2"), width = "180px"),
          plotOutput("acf_plot", height = 410),
          p(class = "note", "Slowly decaying autocorrelation means that nearby iterations carry similar information. ESS converts this dependence into an equivalent number of independent draws.")
        ),
        tabPanel("Target and sample",
          plotOutput("target_plot", height = 440)
        ),
        tabPanel("Running diagnostics",
          plotOutput("running_plot", height = 420)
        ),
        tabPanel("Algorithm",
          wellPanel(withMathJax(
            h3("General Metropolis–Hastings"),
            tags$ol(
              tags$li("At iteration \\(t\\), propose \\(Y\\sim q(\\cdot\\mid X_t)\\)."),
              tags$li("Compute \\(\\alpha(X_t,Y)=\\min\\{1,\\pi(Y)q(X_t\\mid Y)/[\\pi(X_t)q(Y\\mid X_t)]\\}\\)."),
              tags$li("Set \\(X_{t+1}=Y\\) with probability \\(\\alpha\\); otherwise retain \\(X_{t+1}=X_t\\).")) ,
            h3("Random-walk special case"),
            p("If \\(Y=X_t+\\varepsilon\\) and the increment distribution is symmetric, then \\(q(Y\\mid X_t)=q(X_t\\mid Y)\\). The proposal terms cancel, leaving \\(\\min\\{1,\\pi(Y)/\\pi(X_t)\\}\\). They do not cancel for the asymmetric log-normal or Beta independence proposals in this app."),
            h3("Why dependent draws are still useful"),
            p("The transition kernel preserves the target distribution through detailed balance. After convergence, ergodic averages such as \\(n^{-1}\\sum_t h(X_t)\\) converge to \\(E_\\pi[h(X)]\\), even though consecutive draws are dependent."),
            h3("Effective sample size"),
            p("For an approximately stationary scalar chain, \\(ESS \\approx n/[1+2\\sum_{k\\ge1}\\rho_k]\\), where \\(\\rho_k\\) is the lag-\\(k\\) autocorrelation.")))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$target, {
    spec <- target_spec(input$target)
    updateNumericInput(session, "start_x", value = spec$default_start[1])
    updateNumericInput(session, "start_y", value = spec$default_start[2])
    updateSliderInput(session, "proposal_scale", value = spec$default_scale)
  }, ignoreInit = FALSE)

  output$target_note <- renderUI({
    p(class = "note", target_spec(input$target)$description)
  })

  output$proposal_name <- renderText(target_spec(input$target)$proposal_name)

  simulation <- eventReactive(input$run, {
    spec <- target_spec(input$target)
    validate(need(input$burn < input$iterations - 10, "Burn-in must leave at least 10 draws."))
    initial <- if (spec$dim == 1) input$start_x else c(input$start_x, input$start_y)
    validate(need(is.finite(spec$log_density(initial)), "The initial value must lie inside the support of the target."))
    result <- run_mh(spec, initial, input$proposal_scale, input$iterations)
    keep <- seq.int(input$burn + 1, input$iterations)
    result$post <- result$chain[keep, , drop = FALSE]
    result$keep <- keep
    result$spec <- spec
    result$target_name <- input$target
    result$proposal_scale <- input$proposal_scale
    result$burn <- input$burn
    result
  }, ignoreInit = FALSE)

  output$acceptance <- renderText({
    z <- simulation()
    sprintf("%.1f%%", 100 * mean(z$accepted[-1]))
  })

  output$retained <- renderText(format(nrow(simulation()$post), big.mark = ","))

  ess_values <- reactive({
    apply(simulation()$post, 2, effective_sample_size)
  })

  output$ess_metric <- renderText({
    e <- ess_values()
    paste(sprintf("%s: %.0f", names(e), e), collapse = " | ")
  })

  output$trace_plot <- renderPlot({
    z <- simulation(); d <- ncol(z$chain)
    old_par <- par(mfrow = c(d, 1), mar = c(3.5, 4, 2, 1))
    on.exit(par(old_par))
    for (j in seq_len(d)) {
      plot(z$chain[, j], type = "l", col = adjustcolor("#0f766e", .75), lwd = .8,
           xlab = if (j == d) "Iteration" else "", ylab = colnames(z$chain)[j],
           main = if (j == 1) paste("Trace plot:", z$target_name) else "")
      if (z$burn > 0) {
        rect(0, par("usr")[3], z$burn, par("usr")[4], col = adjustcolor("#f59e0b", .12), border = NA)
        abline(v = z$burn, col = "#d97706", lty = 2)
      }
    }
  })

  output$mixing_message <- renderUI({
    rate <- mean(simulation()$accepted[-1]); ratio <- min(ess_values()) / nrow(simulation()$post)
    if (rate < 0.1) {
      p(class = "warning", "Very few proposals are accepted. The proposal scale may be too large, or the target geometry may be difficult.")
    } else if (rate > 0.75 && ratio < 0.15) {
      p(class = "warning", "Acceptance is high, but ESS is low. The chain is probably taking many small, highly correlated steps.")
    } else if (ratio < 0.05) {
      p(class = "warning", "The ESS is less than 5% of the retained chain. Inspect the trace for sticking, slow exploration, or failure to cross modes.")
    } else {
      p(class = "good", "The diagnostics are reasonably healthy for this demonstration, but always inspect all plots rather than relying on acceptance alone.")
    }
  })

  observe({
    d <- target_spec(input$target)$dim
    updateSelectInput(session, "acf_coordinate", choices = paste0("X", seq_len(d)), selected = "X1")
  })

  output$acf_plot <- renderPlot({
    z <- simulation(); coordinate <- match(input$acf_coordinate, colnames(z$post))
    if (is.na(coordinate)) coordinate <- 1
    n <- nrow(z$post); lag_max <- min(100, floor(n / 4))
    acf(z$post[, coordinate], lag.max = lag_max, col = "#0f766e", lwd = 2,
        main = paste("ACF of", colnames(z$post)[coordinate]), xlab = "Lag")
    abline(h = c(-1, 1) * 1.96 / sqrt(n), col = "#b42318", lty = 2)
  })

  output$target_plot <- renderPlot({
    z <- simulation(); s <- z$spec
    if (s$dim == 1) {
      grid <- seq(s$limits[1], s$limits[2], length.out = 700)
      hist(z$post[, 1], probability = TRUE, breaks = "FD", col = "#99f6e4", border = "white",
           xlim = s$limits, xlab = "x", main = "Retained sample against the target")
      lines(grid, s$density_1d(grid), col = "#111827", lwd = 3)
      rug(z$post[, 1], col = adjustcolor("#0f766e", .2))
      legend("topright", "Target density", col = "#111827", lwd = 3, bty = "n")
    } else {
      grid <- seq(s$limits[1], s$limits[2], length.out = 120)
      surface <- outer(grid, grid, Vectorize(s$contour))
      contour(grid, grid, surface, nlevels = 8, drawlabels = FALSE, col = "#64748b",
              xlab = "X1", ylab = "X2", main = "Post-burn-in draws and target contours")
      points(z$post[, 1], z$post[, 2], pch = 19, cex = .45, col = adjustcolor("#0f766e", .28))
    }
  })

  output$running_plot <- renderPlot({
    z <- simulation(); x <- z$post[, 1]; n <- length(x)
    running_mean <- cumsum(x) / seq_len(n)
    block <- pmax(1, floor(seq_len(n) / 50))
    running_accept <- cumsum(z$accepted[z$keep]) / seq_len(n)
    old_par <- par(mfrow = c(2, 1), mar = c(3.5, 4, 2, 1))
    on.exit(par(old_par))
    plot(running_mean, type = "l", lwd = 2, col = "#0f766e", xlab = "Retained iteration",
         ylab = "Running mean", main = "Running mean of X1")
    abline(h = 0, lty = 2, col = "#64748b")
    plot(running_accept, type = "l", lwd = 2, col = "#7c3aed", ylim = c(0, 1),
         xlab = "Retained iteration", ylab = "Acceptance rate", main = "Cumulative acceptance rate")
  })
}

shinyApp(ui, server)
