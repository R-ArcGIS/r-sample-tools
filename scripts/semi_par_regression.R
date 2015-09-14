
make_knots <-function (x1, x2, num.knots)
{
  if (missing(num.knots)) 
  num.knots <- max(10, min(50, round(length(x1)/4)))
  X <- cbind(x1, x2)
  dup.inds <- (1:nrow(X))[dup.matrix(X) == T]
  if (length(dup.inds) > 0) 
      X <- X[-dup.inds, ]
  knots <- cluster::clara(X, num.knots)$medoids
  return(knots)
}

tool_exec <- function(in_params, out_params)
{
  #### Load Library for Analysis ####
  if (!requireNamespace("SemiPar", quietly = TRUE))
    install.packages("SemiPar")
  require(SemiPar)

  #### Get Input Parameters ####
  input_features <- in_params[[1]]
  input_predictions <- in_params[[2]]
  dep_variable <- in_params[[3]]
  lin_variables <- in_params[[4]]
  nonlin_variables <- in_params[[5]]
  input_knots <- in_params[[6]]
  output_features <- out_params[[1]]
  output_graph_pdf <- out_params[[2]]


  #### Import Dataset to Dataframe ####
  fc <- arc.open(input_features)
  df <- arc.select(fc, c(dep_variable, nonlin_variables, lin_variables))
  df['x'] <- arc.shape(df)$x
  df['y'] <- arc.shape(df)$y
  #### Import Knots to DataFrame ####
  if (is.null(input_knots))
  {
    message("Creating default knots")
    knots.est <- make_knots(df$x, df$y)
  }
  else
  {
    knots_df <- arc.select(arc.open(input_knots))
    knots.est <- make_knots(arc.shape(knots_df)$x, arc.shape(knots_df)$y)
  }

  #### Create Spatial Effect ####
  fxy <- "f(x,y,knots=knots.est)"

  #### Create Formula and Fit SemiPar ####
  e <- as.list(df)
  e$knots.est <- knots.est
  all_params <- paste0(dep_variable, "~",  fxy)

  #### Create Non-Linear Params ####
  if (!is.null(lin_variables))
  {
    nonlin_params <- paste(paste0("f(", nonlin_variables, ")"), collapse = "+")
    all_params <- paste0(all_params, "+", nonlin_params)
  }
  #### Create Linear Params ####
  if (!is.null(nonlin_variables ))
    all_params <- paste0(all_params, "+", paste(lin_variables , collapse = "+"))

  message(paste0("formula = ", all_params))
  tryCatch(
  {
    attach(e)
    form <- as.formula(all_params)
    fit <- spm(form, family="binomial")
  },finally = detach(e))

  summary(fit)
  #### Prediction ####
  pred_fc <- arc.open(input_predictions)
  oid_field <- pred_fc@fields
  oid_field <- names(oid_field[oid_field == 'OID'])[[1]]
  pred_df = arc.select(pred_fc, c(oid_field, nonlin_variables, lin_variables))
  
  if (!is.null(output_graph_pdf))
  {
   tryCatch(
     {
       pdf(output_graph_pdf)
       #nn <- (length(lin_variables) + length(nonlin_variables) + 2)/2
       #par(mfrow=c(ceiling(nn), 2))
       suppressWarnings(
          plot(fit, bdry=default.bdry(bdry = knots.est))
       )
     }, finally = { dev.off() })
  }

  pred_df_xy <- data.frame(
                    x = arc.shape(pred_df)$x,
                    y = arc.shape(pred_df)$y,
                    pred_df)
  pred = with(pred_df_xy, predict.spm(fit, newdata = pred_df_xy, se = TRUE))

  #### Write Output ####
  pred_df['prediction'] = 1 / (1+exp(-pred$fit))
  #pred_df['std_error'] = pred$se
  pred_df['LCL_95'] = 1 / (1+exp(-pred$fit+1.96*pred$se))
  pred_df['UCL_95'] = 1/(1+exp(-pred$fit-1.96*pred$se))

  arc.write(output_features, pred_df)

  return(out_params)
}

