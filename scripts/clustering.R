# case study using mclust package
tool_exec <- function(in_params, out_params)
{
  if (!requireNamespace("sp", quietly = TRUE))
    install.packages("sp")
  if (!requireNamespace("mclust", quietly = TRUE))
    install.packages("mclust")
  require(mclust)
  require(sp)

  source_dataset = in_params[[1]]
  #nclust = in_params[[2]]
  out_table = out_params[[1]]
  out_ellipses = out_params[[2]]
  out_dens = out_params[[3]]
  out_sim = out_params[[4]]
  ### read data
  arc.progress_label("Loading Dataset")
  d <- arc.open(source_dataset)
  # read only OID field
  data <- arc.select(d, names(d@fields[d@fields == "OID"]))
  data_shp <- arc.shape(data)

  data.xy <- cbind(data_shp$x, data_shp$y) # or data.frame(x=data_shp$x, y=data_shp$y)

  #########
  arc.progress_label("Bayesian Information Criterion for Model-Based Clustering...")
  patternBIC <- mclustBIC(data.xy) # G = nclust)

  patternModel <- summary(patternBIC, data.xy)
  print(patternModel)

  ### save to shapefile
  # patternModel$classification # cluster number
  # patternModel$uncertainty    # uncertainty of classification
  # patternModel$z              # conditional probabilities
  n <- patternModel$G
  cond_probs <- lapply(1:n, function(i) patternModel$z[,i])
  names(cond_probs) <- paste0("cond_prob", 1:n)

  q <- quantile(patternModel$uncertainty, probs = c(0.75, 0.95))
  q7595 <- sapply(patternModel$uncertainty, function(u){
      r <- if (u <= q[1]) 1 else 0
      if(u <= q[2] && r == 0) r <- 2
      if(u > q[2] && r == 0) r <- 3
      stopifnot(r != 0)
      r - 1
    })

  result <- data.frame(data,
              cluster = patternModel$classification,
              class_uncert = patternModel$uncertainty,
              sym_ucert=q7595,
              cond_probs)

  #print(head(result))
  arc.progress_label("Writing result dataset...")
  if (!is.null(out_table) && out_table != "NA")
    arc.write(out_table, result, coords = data_shp)

  ### automatic mapping
  arc.progress_label("Best model based on Bayesian Information Criterion...")
  bestModel <- mclustModel(data.xy, patternBIC)

  ### create ellipses
  if (!is.null(out_ellipses) && out_ellipses != "NA")
  {
    arc.progress_label("Writing ellipses dataset...")
    polygons <- create.ellipses(bestModel)
    mu <- bestModel$parameters$mean
    sp.df <- SpatialPolygonsDataFrame(polygons, data=data.frame(clast_n=1:n, mu_x=mu[1,1:n], mu_y=mu[2,1:n]), match.ID = F)
    ## keep original SpatialReference
    ## following is not always correct WKT -> P4 -> WKT' == WKT
    shape_info <- list(type="Polygon", WKT=arc.shapeinfo(data_shp)$WKT)
    #shape_info <- arc.shapeinfo(sp.df)
    #shape_info$WKT <- arc.shapeinfo(data_shp)$WKT
    arc.write(out_ellipses, sp.df, shape_info=shape_info)
  }

  ### create density grid
  if (!is.null(out_dens) && out_dens != "NA")
  {
    arc.progress_label("Calculating density...")
    tmp_file = tempfile("density", fileext = ".pdf")
    options(onefile=FALSE)
    pdf(tmp_file)

    grid.n <- 100
    grid <- surfacePlot(verbose=T, data.xy, parameters = bestModel$parameters, type = "image", what = "density", grid=grid.n, nlevels = 8, transformation = "none", ask=F)
    dev.off()
    options(device="windows")
    options(onefile=NULL)
    #file.remove(tmp_file)
    xy <-list(x = unlist(lapply(grid$x, function(x) rep(x, grid.n))),
              y = rep(grid$y, grid.n))
    arc.progress_label("Writing density dataset...")
    arc.write(out_dens, data=list(density = as.vector(grid$z)), xy, list(type='Point', WKT=arc.shapeinfo(data_shp)$WKT))
  }

  ### simulation, fitting model
  if (!is.null(out_sim) && out_sim != "NA")
  {
    sim.xy = simulation(bestModel, nrow(data))
    #print(length(sim.xy$x))
    arc.progress_label("Saving simulated dataset...")
    arc.write(out_sim, data=list(id=1:length(sim.xy$x)), sim.xy, list(type='Point', WKT=arc.shapeinfo(data_shp)$WKT))
  }
  arc.progress_label("Done")
  return(out_params)
}

simulation <- function(bestModel, n)
{
  n <- round(round(rnorm(1, mean=n, sd=sqrt(n))))
  sim <- sim(modelName = bestModel$modelName, parameters = bestModel$parameters, n, seed = round(sqrt(n+1)))
  list(x=sim[,2], y=sim[,3])
}

create.ellipses <- function(bestModel)
{
  n <- bestModel$G
  mu <- bestModel$parameters$mean
  sigma <- bestModel$parameters$variance$sigma
  cls.polygons <- lapply(1:n, function(i)
  {
    xy <- make.ellipse(mu = mu[, i], sigma = sigma[, , i])
    name<-paste0("ellipse", i)
    Polygons(list(Polygon(xy)), name)
  })
  #SpatialPolygons(cls.polygons, 1:n, proj4string=sp::CRS(p4str))
  SpatialPolygons(cls.polygons, 1:n, proj4string=sp::CRS())
}

make.ellipse <- function(mu, sigma, k = 60)
{
  p <- length(mu)
  if (p != 2) 
    stop("only two-dimensional case is available")
  if (any(unique(dim(sigma)) != p)) 
    stop("mu and sigma are incompatible")
  ev <- eigen(sigma, symmetric = TRUE)
  s <- sqrt(rev(sort(ev$values)))
  V <- t(ev$vectors[, rev(order(ev$values))])
  theta <- (0:k) * (2*pi/k)
  x <- s[1] * cos(theta) 
  y <- s[2] * sin(theta)
  xy <- cbind(x, y)
  xy <- xy %*% V
  cbind(xy[,1] + mu[1], xy[,2] + mu[2])
}

## standalone R testing
foo <- function()
{
  library(arcgisbinding)
  arc.check_product()
  tool_exec(
      list("d:\\Data\\R-staff\\samps\\data.gdb\\spruce_trees", NULL),
      list(NULL, NULL, NULL, "in_memory\\sim")
   )
}
