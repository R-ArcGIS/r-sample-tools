#
# Belarus thyroid cancer in children case study
#

tool_exec <- function(in_params, out_params)
{
### handle input gp tool paramaters
  dataset_path <- in_params[[1]]
  case_field <- in_params[[2]]
  pop_field <- in_params[[3]]
  covar1_field <- in_params[[4]]
  covar2_field <- in_params[[5]]

  df = regional_bayesian(dataset_path, case_field, pop_field, covar1_field, covar2_field)

### handle gp tool outputs
  output_dataset_name = out_params[[1]]
  arc.write(output_dataset_name, df)
  return(out_params)
}

#stat analysis with WINBUGS
regional_bayesian <- function(dataset_path, case_field, pop_field, covar1_field, covar2_field)
{
  if (!requireNamespace("sp", quietly = TRUE))
    install.packages("sp")
  if (!requireNamespace("spdep", quietly = TRUE))
    install.packages("spdep")
  if (!requireNamespace("R2WinBUGS", quietly = TRUE))
    install.packages("R2WinBUGS")
  if (!requireNamespace("rgeos", quietly = TRUE))
    install.packages("rgeos")
  require(sp)
  require(spdep)

  message("Loading dataset...")
  ds <- arc.open(dataset_path)
  df <- arc.select(ds, fields = list(case_field, pop_field, covar1_field, covar2_field))
  sp.df <- arc.data2sp(df)

  # convert shapefile data to a list of polygons and a data frame
  centroids <- rgeos::gCentroid(sp.df, byid = TRUE)
  centroids <- centroids@coords

  kk <- nrow(sp.df)

  # calculate the expected counts
  pmap <- probmap(df[[case_field]], df[[pop_field]])

  # make neighbours list
  message("Making Neighbours...")
  five.nn <- knn2nb(knearneigh(centroids, k=5), sym=TRUE)

  #plotpolys(th_polys, bbs, border="grey")
  #plot(five.nn, centroids, add=TRUE)

#####4 "GWR"

  ## data preparation
  dlist <- nbdists(five.nn, centroids)
  num <- sapply(five.nn, length)

  ll <- sum(num)
  # equal Weights:
  weights <- rep(1, ll)

  adj <- rep(0, ll)
  mm <- 1
  for (i in 1:kk)
  {
   for (j in 1:num[i])
   {
     adj[mm] <- five.nn[[i]][j]
     mm <- mm + 1
   }
  }

  input.data <- list(n=kk,
    cases = df[[case_field]],
    expected = pmap$expCount,
    weights=weights,
    a1 = df[[covar1_field]],
    a2 = df[[covar2_field]]/min(df[[covar2_field]]),
    num=num,
    adj=adj)

  ### inits
  b_init <- rep(0.0, kk)
  model.inits <- function() {list(alpha = 0, tau.v = 1, tau.u = 1, b1 = b_init, b2 = b_init)}
  # output parameters
  model.parameters <- c("theta", "mu", "v", "u", "alpha", "PP", "b1", "b2", "RR_exp", "RR_het", "RR_clust")

  model_file <- tempfile("BRCM-model", fileext='.bug')
  R2WinBUGS::write.model(BRCM.model, model_file)

  ### run BUGS
  message("Running WinBUGS model...")
  result.4 <- R2WinBUGS::bugs(input.data, model.inits, model.parameters, model_file, n.chains=1, n.iter=100000, debug=F, working.directory=tempdir())

  unlink(model_file)

#append results columns
  df$expected <- pmap$expCount
  df$rel_risk       <- result.4$mean$theta
  df$rel_risk_se    <- result.4$sd$theta
  df$extra_risk       <- result.4$mean$PP
  df$extra_risk_se <- result.4$sd$PP
  df$mean       <- result.4$mean$mu
  df$mean_se    <- result.4$sd$mu
  df$u_cor    <- result.4$mean$u
  df$u_cor_se <- result.4$sd$u
  df$u_heter  <- result.4$mean$v
  df$u_heter_se <- result.4$sd$v
  df$covariates   <- result.4$mean$RR_exp
  df$heterogeneous   <- result.4$mean$RR_het
  df$residuals <- result.4$mean$RR_clust
  df$r1       <- result.4$mean$b1
  df$r1_se    <- result.4$sd$b1
  df$r2       <- result.4$mean$b2
  df$r2_se    <- result.4$sd$b2

  return(df)
}

BRCM.model <- function()
{
  for (i in 1:n)
  {
    # Poisson likelihood for observed counts
    cases[i] ~ dpois(mu[i])
    log(mu[i]) <- log(expected[i]) + alpha + u[i] + v[i] + b1[i]*a1[i] + b2[i]*a2[i]
    # Relative Risk
    theta[i] <- exp(alpha + u[i] + v[i] + b1[i]*a1[i] + b2[i]*a2[i])
    # Posterior probability of RR[i]>1
    PP[i] <- step(theta[i] - 1 + eps)
    # Prior distribution for the uncorrelated heterogeneity
    v[i] ~ dnorm(0,tau.v)
    # Relative Risk decomposition
    RR_exp[i]<-exp(b1[i]*a1[i] + b2[i]*a2[i])
    RR_het[i]<-exp(v[i])
    RR_clust[i]<-exp(u[i])
  }

  eps <- 1.0E-6

  # CAR prior distribution for spatial correlated heterogeneity
  u[1:n] ~ car.normal(adj[], weights[], num[], tau.u)
  b1[1:n] ~ car.normal(adj[], weights[], num[], tau.b1)
  b2[1:n] ~ car.normal(adj[], weights[], num[], tau.b2)

  # Improper prior distribution for the mean relative risk in the study region
  alpha ~ dflat()
  mean <- exp(alpha)

  # Hyperprior distributions on inverse variance parameter of random effects
  tau.u ~ dgamma(0.5,0.0005)
  tau.v ~ dgamma(0.5,0.0005)
  tau.b1 ~ dgamma(0.5,0.0005)
  tau.b2 ~ dgamma(0.5,0.0005)
}
