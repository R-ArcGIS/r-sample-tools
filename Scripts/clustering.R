# case study using mclust package
tool_exec <- function(in_params, out_params) {
  message(getOption("repos"))
  if (!requireNamespace("sp", quietly = TRUE))
    install.packages("sp")
  if (!requireNamespace("mclust", quietly = TRUE))
    install.packages("mclust")
  require(mclust)
  require(sp)

  source_dataset <- in_params[[1]]
  nclust <- in_params[[2]]
  out_table <- out_params[[1]]
  out_ellipses <- out_params[[2]]
  out_dens <- out_params[[3]]

  ### read data
  arc.progress_label("Loading Dataset")
  data <- arc.select(arc.open(source_dataset))
  data_shp <- arc.shape(data)

  #data.xy <- data.frame(x=data_shp$x - x_min, y=data_shp$y - y_min)
  # alternatively: data.frame(x=data_shp$x, y=data_shp$y)
  data.xy <- cbind(data_shp$x, data_shp$y)

  #########
  arc.progress_label("BIC for Model-Based Clustering...")
  treesBIC <- mclustBIC(data.xy, G = nclust)
  treesModel <- summary(treesBIC, data.xy)
  print(treesModel)

  ### save to shapefile
  # treesModel$classification # cluster number
  # treesModel$uncertainty    # uncertainty of classification
  # treesModel$z              # conditional probabilities
  n <- treesModel$G
  cond_probs <- lapply(1:n, function(i) treesModel$z[,i])
  names(cond_probs) <- paste0("cond_prob", 1:n)

  q <- quantile(treesModel$uncertainty, probs = c(0.75, 0.95))
  q7595 <- sapply(treesModel$uncertainty, function(u) {
      r <- if (u <= q[1]) 1 else 0
      if (u <= q[2] && r == 0) r <- 2
      if (u > q[2] && r == 0) r <- 3
      stopifnot(r != 0)
      r - 1
    })

  result <- data.frame(data,
                       cluster = treesModel$classification,
                       class_uncert = treesModel$uncertainty,
                       sym_ucert=q7595,
                       cond_probs)

  #print(head(result))
  arc.progress_label("Writing result dataset...")
  arc.write(out_table, result, coords = data_shp)

  ### automatic mapping
  arc.progress_label("Best model based on BIC...")
  bestModel <- mclustModel(data.xy, treesBIC)

  ### create ellipses
  if (!is.null(out_ellipses)) {
    arc.progress_label("Writing ellipses dataset...")
    polygons <- create.elipses(bestModel)
    mu <- bestModel$parameters$mean
    sp.df <- SpatialPolygonsDataFrame(polygons,
                                      data=data.frame(clust_n=1:n,
                                                      mu_x=mu[1, 1:n],
                                                      mu_y=mu[2, 1:n]),
                                      match.ID = F)
    ## keep original SpatialReference
    shape_info <- list(type="Polygon", WKT=arc.shapeinfo(data_shp)$WKT)
    #shape_info <- arc.shapeinfo(sp.df)
    arc.write(out_ellipses, sp.df, shape_info=shape_info)
  }

  ### create density grid
  if (!is.null(out_dens)) {
    arc.progress_label("Calculating density...")
    options(device="pdf")
    options(onefile=FALSE)
    grid.n <- 100
    grid <- surfacePlot(verbose=T, data.xy, parameters = bestModel$parameters,
                        type = "image", what = "density", grid = grid.n,
                        nlevels = 8, transformation = "none", ask=F)
    dev.off()
    options(device="windows")
    options(onefile=NULL)
    xy <- list(x = unlist(lapply(grid$x, function(x) rep(x, grid.n))),
              y = rep(grid$y, grid.n))
    arc.progress_label("Writing density dataset...")
    arc.write(out_dens,
              data=list(density = as.vector(grid$z)),
              xy,
              list(type="Point", WKT=arc.shapeinfo(data_shp)$WKT))
  }
  arc.progress_label("Done")
  return(out_params)
}

create.elipses <- function(bestModel) {
  n <- bestModel$G
  mu <- bestModel$parameters$mean
  sigma <- bestModel$parameters$variance$sigma
  pp <- lapply(1:n, function(i) {
    xy <- make.ellipse(mu = mu[, i], sigma = sigma[, , i])
    name <- paste0("ellipse", i)
    Polygons(list(Polygon(xy)), name)
  })
  #p4str <- arc.fromWktToP4(arc.shapeinfo(trees_shp)$WKT)
  #SpatialPolygons(pp, 1:n, proj4string=sp::CRS(p4str))
  SpatialPolygons(pp, 1:n, proj4string=sp::CRS())
}

make.ellipse <- function(mu, sigma, k = 60) {
  p <- length(mu)
  if (p != 2) 
    stop("only two-dimensional case is available")
  if (any(unique(dim(sigma)) != p)) 
    stop("mu and sigma are incompatible")
  ev <- eigen(sigma, symmetric = TRUE)
  s <- sqrt(rev(sort(ev$values)))
  V <- t(ev$vectors[, rev(order(ev$values))])
  theta <- (0:k) * (2 * pi / k)
  x <- s[1] * cos(theta) 
  y <- s[2] * sin(theta)
  xy <- cbind(x, y)
  xy <- xy %*% V
  cbind(xy[,1] + mu[1], xy[,2] + mu[2])
}
