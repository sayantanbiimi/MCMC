library(shiny)

effective_sample_size <- function(x) {
  n <- length(x)
  if (n < 4 || !is.finite(var(x)) || var(x) == 0) return(1)
  rho <- as.numeric(acf(x, lag.max = min(n - 1, max(20, floor(10 * log10(n)))),
                        plot = FALSE)$acf)[-1]
  if (length(rho) < 2) return(n)
  paired <- rho[seq(1, length(rho) - 1, 2)] + rho[seq(2, length(rho), 2)]
  stop_at <- which(paired <= 0)[1]
  if (!is.na(stop_at)) paired <- paired[seq_len(stop_at - 1)]
  tau_int <- if (length(paired)) 1 + 2 * sum(paired) else 1
  min(n, max(1, n / max(1, tau_int)))
}

gibbs_bivariate <- function(iterations, rho, mean_x, mean_y, sd_x, sd_y) {
  chain <- matrix(NA_real_, iterations, 2, dimnames = list(NULL, c("X", "Y")))
  chain[1, ] <- c(mean_x + 3 * sd_x, mean_y - 3 * sd_y)
  for (t in 2:iterations) {
    mx <- mean_x + rho * sd_x / sd_y * (chain[t - 1, 2] - mean_y)
    chain[t, 1] <- rnorm(1, mx, sd_x * sqrt(1 - rho^2))
    my <- mean_y + rho * sd_y / sd_x * (chain[t, 1] - mean_x)
    chain[t, 2] <- rnorm(1, my, sd_y * sqrt(1 - rho^2))
  }
  list(chain = chain, data = list(rho = rho, mean_x = mean_x, mean_y = mean_y,
                                   sd_x = sd_x, sd_y = sd_y), kind = "bivariate")
}

gibbs_hierarchical <- function(iterations, groups, observations_per_group,
                               true_mu, true_tau, observation_sd,
                               prior_family, m0, prior_scale, a0, b0) {
  theta_true <- rnorm(groups, true_mu, true_tau)
  ybar <- rnorm(groups, theta_true, observation_sd / sqrt(observations_per_group))

  theta <- ybar
  mu <- m0
  tau2 <- max(var(ybar), 0.25)
  chain <- matrix(NA_real_, iterations, groups + 2)
  colnames(chain) <- c("mu", "tau", paste0("theta", seq_len(groups)))

  for (t in seq_len(iterations)) {
    theta_variance <- 1 / (observations_per_group / observation_sd^2 + 1 / tau2)
    theta_mean <- theta_variance *
      (observations_per_group * ybar / observation_sd^2 + mu / tau2)
    theta <- rnorm(groups, theta_mean, sqrt(theta_variance))

    if (prior_family == "Independent Normal + inverse-gamma") {
      mu_variance <- 1 / (groups / tau2 + 1 / prior_scale^2)
      mu_mean <- mu_variance * (sum(theta) / tau2 + m0 / prior_scale^2)
      mu <- rnorm(1, mu_mean, sqrt(mu_variance))
      shape <- a0 + groups / 2
      rate <- b0 + 0.5 * sum((theta - mu)^2)
    } else {
      kappa0 <- prior_scale
      mu <- rnorm(1, (sum(theta) + kappa0 * m0) / (groups + kappa0),
                  sqrt(tau2 / (groups + kappa0)))
      shape <- a0 + (groups + 1) / 2
      rate <- b0 + 0.5 * sum((theta - mu)^2) + 0.5 * kappa0 * (mu - m0)^2
    }

    tau2 <- 1 / rgamma(1, shape = shape, rate = rate)
    chain[t, ] <- c(mu, sqrt(tau2), theta)
  }

  list(chain = chain,
       data = list(theta_true = theta_true, ybar = ybar, groups = groups,
                   observations_per_group = observations_per_group,
                   true_mu = true_mu, true_tau = true_tau,
                   observation_sd = observation_sd, prior_family = prior_family,
                   m0 = m0, prior_scale = prior_scale, a0 = a0, b0 = b0),
       kind = "hierarchical")
}

inverse_gamma_sd_density <- function(tau, shape, rate) {
  log_density <- log(2) + shape * log(rate) - lgamma(shape) -
    (2 * shape + 1) * log(tau) - rate / tau^2
  ifelse(tau > 0, exp(log_density), 0)
}

ui <- fluidPage(
  tags$head(tags$style(HTML("body{background:#f6f7fb;color:#172033}.hero{background:linear-gradient(120deg,#312e81,#7c3aed);color:white;border-radius:14px;padding:22px 28px;margin:18px 0}.well{background:white;border:0;box-shadow:0 2px 12px #dfe3ec}.metric{font-size:26px;font-weight:750;color:#6d28d9}.note{color:#5f6b7a}.warning{color:#b42318;font-weight:650}.good{color:#137333;font-weight:650}.btn-primary{background:#6d28d9;border-color:#6d28d9}.tab-content{padding-top:18px}"))),
  div(class = "hero", h2("Gibbs sampling demo"),
      p("Conditional simulation, hierarchical shrinkage and prior sensitivity")),
  sidebarLayout(
    sidebarPanel(width = 4,
      selectInput("example", "Choose an example",
                  c("Bivariate normal", "Conjugate normal hierarchical model")),

      conditionalPanel(
        condition = "input.example == 'Bivariate normal'",
        h5("Target parameters"),
        fluidRow(column(6, numericInput("mean_x", "Mean of X", 0, step = .25)),
                 column(6, numericInput("mean_y", "Mean of Y", 0, step = .25))),
        fluidRow(column(6, numericInput("sd_x", "SD of X", 1, min = .1, step = .1)),
                 column(6, numericInput("sd_y", "SD of Y", 1, min = .1, step = .1))),
        sliderInput("rho", "Correlation", -0.995, 0.995, 0.9, step = .005)
      ),

      conditionalPanel(
        condition = "input.example == 'Conjugate normal hierarchical model'",
        h5("Data-generating parameters"),
        fluidRow(column(6, numericInput("groups", "Number of groups", 8, min = 3, max = 20)),
                 column(6, numericInput("group_n", "Observations/group", 10, min = 2, max = 100))),
        fluidRow(column(6, numericInput("true_mu", "True grand mean", 2, step = .25)),
                 column(6, numericInput("true_tau", "True between-group SD", 1, min = .1, step = .1))),
        numericInput("observation_sd", "Known observation SD", 2, min = .1, step = .1),
        hr(), h5("Conjugate prior"),
        selectInput("prior_family", "Prior structure",
                    c("Independent Normal + inverse-gamma",
                      "Normal–inverse-gamma")),
        numericInput("m0", "Prior centre m₀", 0, step = .25),
        conditionalPanel(
          condition = "input.prior_family == 'Independent Normal + inverse-gamma'",
          numericInput("s0", "Prior SD s₀ for μ", 5, min = .05, step = .25)
        ),
        conditionalPanel(
          condition = "input.prior_family == 'Normal–inverse-gamma'",
          numericInput("kappa0", "Prior precision multiplier κ₀", 1, min = .01, step = .25)
        ),
        fluidRow(column(6, numericInput("a0", "IG shape a₀", 2, min = .1, step = .25)),
                 column(6, numericInput("b0", "IG rate b₀", 2, min = .1, step = .25)))
      ),

      hr(),
      sliderInput("iterations", "Iterations", 500, 15000, 5000, step = 500),
      sliderInput("burn", "Burn-in", 0, 4000, 500, step = 100),
      numericInput("seed", "Random seed", 123, min = 1, step = 1),
      actionButton("run", "Run Gibbs sampler", class = "btn-primary", width = "100%"),
      p(class = "note", "Changing a hyperparameter generates a new posterior only after you click Run Gibbs sampler.")
    ),
    mainPanel(width = 8,
      fluidRow(
        column(4, wellPanel(h4("Retained draws"), div(class = "metric", textOutput("retained")))),
        column(4, wellPanel(h4("Smallest ESS"), div(class = "metric", textOutput("min_ess")))),
        column(4, wellPanel(h4("ESS efficiency"), div(class = "metric", textOutput("efficiency"))))
      ),
      tabsetPanel(
        tabPanel("Trace plot",
                 selectInput("trace_parameter", "Parameter", "X", width = "220px"),
                 plotOutput("trace_plot", height = 400), uiOutput("diagnostic_note")),
        tabPanel("Autocorrelation",
                 selectInput("acf_parameter", "Parameter", "X", width = "220px"),
                 plotOutput("acf_plot", height = 400)),
        tabPanel("Target or recovery", plotOutput("recovery_plot", height = 450)),
        tabPanel("Prior versus posterior", uiOutput("prior_message"),
                 plotOutput("prior_posterior_plot", height = 430)),
        tabPanel("ESS summary", tableOutput("ess_table"), plotOutput("ess_plot", height = 300)),
        tabPanel("Model and full conditionals",
          uiOutput("model_explanation")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  result <- eventReactive(input$run, {
    validate(need(input$burn < input$iterations - 20, "Burn-in must leave at least 20 draws."))
    set.seed(input$seed)
    fit <- if (input$example == "Bivariate normal") {
      gibbs_bivariate(input$iterations, input$rho, input$mean_x, input$mean_y,
                      input$sd_x, input$sd_y)
    } else {
      prior_scale <- if (input$prior_family == "Independent Normal + inverse-gamma") input$s0 else input$kappa0
      gibbs_hierarchical(input$iterations, as.integer(input$groups), as.integer(input$group_n),
                         input$true_mu, input$true_tau, input$observation_sd,
                         input$prior_family, input$m0, prior_scale, input$a0, input$b0)
    }
    fit$burn <- input$burn
    fit
  }, ignoreInit = FALSE)

  post <- reactive({
    fit <- result()
    fit$chain[(fit$burn + 1):nrow(fit$chain), , drop = FALSE]
  })
  ess <- reactive(apply(post(), 2, effective_sample_size))

  observe({
    choices <- colnames(result()$chain)
    updateSelectInput(session, "trace_parameter", choices = choices, selected = choices[1])
    updateSelectInput(session, "acf_parameter", choices = choices, selected = choices[1])
  })

  output$retained <- renderText(format(nrow(post()), big.mark = ","))
  output$min_ess <- renderText(sprintf("%.0f", min(ess())))
  output$efficiency <- renderText(sprintf("%.1f%%", 100 * min(ess()) / nrow(post())))

  output$trace_plot <- renderPlot({
    chain <- result()$chain; j <- match(input$trace_parameter, colnames(chain)); if (is.na(j)) j <- 1
    plot(chain[, j], type = "l", lwd = .75, col = adjustcolor("#6d28d9", .8),
         xlab = "Iteration", ylab = colnames(chain)[j], main = paste("Trace of", colnames(chain)[j]))
    if (result()$burn > 0) {
      rect(0, par("usr")[3], result()$burn, par("usr")[4],
           col = adjustcolor("#f59e0b", .12), border = NA)
      abline(v = result()$burn, col = "#d97706", lty = 2)
    }
  })

  output$diagnostic_note <- renderUI({
    ratio <- min(ess()) / nrow(post())
    if (ratio < .05) p(class = "warning", "The smallest ESS is below 5% of the retained chain. Gibbs updates are accepted, but the chain is moving very slowly.")
    else if (ratio < .2) p(class = "warning", "The chain is noticeably autocorrelated. Inspect the slowest parameter and compare alternative hyperparameters or target correlations.")
    else p(class = "good", "The ESS is reasonably healthy for this run, although the trace and ACF should still be inspected.")
  })

  output$acf_plot <- renderPlot({
    x <- post(); j <- match(input$acf_parameter, colnames(x)); if (is.na(j)) j <- 1
    acf(x[, j], lag.max = min(100, floor(nrow(x) / 4)), col = "#6d28d9", lwd = 2,
        main = paste("ACF of", colnames(x)[j]))
    abline(h = c(-1, 1) * 1.96 / sqrt(nrow(x)), col = "#b42318", lty = 2)
  })

  output$recovery_plot <- renderPlot({
    fit <- result(); draws <- post(); dat <- fit$data
    if (fit$kind == "bivariate") {
      xg <- seq(dat$mean_x - 4*dat$sd_x, dat$mean_x + 4*dat$sd_x, length.out = 100)
      yg <- seq(dat$mean_y - 4*dat$sd_y, dat$mean_y + 4*dat$sd_y, length.out = 100)
      z <- outer(xg, yg, function(x, y) {
        u <- (x-dat$mean_x)/dat$sd_x; v <- (y-dat$mean_y)/dat$sd_y
        exp(-(u^2 - 2*dat$rho*u*v + v^2)/(2*(1-dat$rho^2)))
      })
      contour(xg, yg, z, drawlabels = FALSE, nlevels = 8, col = "#64748b",
              xlab = "X", ylab = "Y", main = "Gibbs draws and target contours")
      points(draws[,1], draws[,2], pch = 19, cex = .4, col = adjustcolor("#6d28d9", .22))
      points(dat$mean_x, dat$mean_y, pch = 4, cex = 1.5, lwd = 3, col = "#b42318")
    } else {
      groups <- dat$groups; theta_draws <- draws[, paste0("theta", seq_len(groups)), drop = FALSE]
      estimate <- colMeans(theta_draws)
      interval <- apply(theta_draws, 2, quantile, probs = c(.025, .975))
      limits <- range(interval, dat$theta_true, dat$ybar)
      plot(seq_len(groups), estimate, ylim = limits, pch = 19, col = "#6d28d9",
           xlab = "Group", ylab = expression(theta[j]),
           main = "Group-level recovery and hierarchical shrinkage")
      segments(seq_len(groups), interval[1,], seq_len(groups), interval[2,], col = "#6d28d9", lwd = 2)
      points(seq_len(groups), dat$ybar, pch = 1, col = "#64748b", cex = 1.2)
      points(seq_len(groups), dat$theta_true, pch = 4, col = "#b42318", lwd = 2, cex = 1.2)
      abline(h = mean(draws[,"mu"]), col = "#111827", lty = 2, lwd = 2)
      legend("topleft", c("Posterior mean", "95% interval", "Observed group mean", "True group effect"),
             col = c("#6d28d9", "#6d28d9", "#64748b", "#b42318"),
             pch = c(19, NA, 1, 4), lty = c(NA, 1, NA, NA), lwd = c(NA, 2, NA, 2), bty = "n")
    }
  })

  output$prior_message <- renderUI({
    if (result()$kind == "bivariate")
      p(class = "note", "This example samples from a specified bivariate target; it is not a Bayesian model and therefore has no prior distribution.")
    else p(class = "note", paste("Prior used:", result()$data$prior_family,
                                  ". Change its hyperparameters in the left panel and rerun."))
  })

  output$prior_posterior_plot <- renderPlot({
    fit <- result(); if (fit$kind == "bivariate") { plot.new(); text(.5,.5,"No prior in the bivariate target example."); return() }
    dat <- fit$data; draws <- post()
    old <- par(mfrow = c(1,2), mar = c(4,4,2,1)); on.exit(par(old))

    mu_grid <- seq(min(draws[,"mu"], dat$m0) - 2, max(draws[,"mu"], dat$m0) + 2, length.out = 500)
    hist(draws[,"mu"], probability = TRUE, breaks = "FD", col = "#ddd6fe", border = "white",
         xlab = expression(mu), main = "Grand mean")
    if (dat$prior_family == "Independent Normal + inverse-gamma") {
      lines(mu_grid, dnorm(mu_grid, dat$m0, dat$prior_scale), col = "#b42318", lwd = 3)
    } else {
      prior_scale_mu <- sqrt(dat$b0 / (dat$a0 * dat$prior_scale))
      lines(mu_grid, dt((mu_grid-dat$m0)/prior_scale_mu, df = 2*dat$a0)/prior_scale_mu,
            col = "#b42318", lwd = 3)
    }
    abline(v = dat$true_mu, col = "#111827", lty = 2, lwd = 2)

    tau_max <- max(quantile(draws[,"tau"], .995), dat$true_tau * 1.5)
    tau_grid <- seq(.001, tau_max, length.out = 500)
    hist(draws[,"tau"], probability = TRUE, breaks = "FD", col = "#ddd6fe", border = "white",
         xlim = c(0,tau_max), xlab = expression(tau), main = "Between-group SD")
    lines(tau_grid, inverse_gamma_sd_density(tau_grid, dat$a0, dat$b0), col = "#b42318", lwd = 3)
    abline(v = dat$true_tau, col = "#111827", lty = 2, lwd = 2)
    legend("topright", c("Prior", "Truth"), col = c("#b42318", "#111827"),
           lty = c(1,2), lwd = c(3,2), bty = "n")
  })

  output$ess_table <- renderTable({
    e <- ess(); data.frame(Parameter = names(e), ESS = round(e),
      `ESS / retained` = sprintf("%.1f%%", 100*e/nrow(post())), check.names = FALSE)
  })

  output$ess_plot <- renderPlot({
    e <- ess(); barplot(e, col = "#8b5cf6", border = NA, ylab = "ESS",
                        main = "Effective sample size by parameter")
    abline(h = nrow(post()), col = "#64748b", lty = 2)
  })

  output$model_explanation <- renderUI({
    if (result()$kind == "bivariate") {
      withMathJax(wellPanel(
        h3("Bivariate normal full conditionals"),
        p("For means \\(\\mu_X,\\mu_Y\\), standard deviations \\(\\sigma_X,\\sigma_Y\\), and correlation \\(\\rho\\),"),
        p("\\[X\\mid Y=y\\sim N\\!\\left(\\mu_X+\\rho\\frac{\\sigma_X}{\\sigma_Y}(y-\\mu_Y),\\;\\sigma_X^2(1-\\rho^2)\\right),\\]"),
        p("with the analogous expression for \\(Y\\mid X\\). Large \\(|\\rho|\\) produces slow component-wise movement.")))
    } else {
      dat <- result()$data
      withMathJax(wellPanel(
        h3("Hierarchical model"),
        p("\\[\\bar Y_j\\mid\\theta_j\\sim N(\\theta_j,\\sigma^2/n_j),\\qquad \\theta_j\\mid\\mu,\\tau^2\\sim N(\\mu,\\tau^2).\\]"),
        h3("Selectable conjugate prior"),
        if (dat$prior_family == "Independent Normal + inverse-gamma")
          p("\\[\\mu\\sim N(m_0,s_0^2),\\qquad \\tau^2\\sim IG(a_0,b_0).\\]")
        else
          p("\\[\\mu\\mid\\tau^2\\sim N(m_0,\\tau^2/\\kappa_0),\\qquad \\tau^2\\sim IG(a_0,b_0).\\]"),
        p("Every full conditional is Normal or inverse-gamma. The group effects borrow strength through the common \\(\\mu\\) and \\(\\tau\\).")))
    }
  })
}

shinyApp(ui, server)
