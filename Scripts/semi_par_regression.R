make_knots <- function(x1, x2, num.knots) {
    if (missing(num.knots)) 
    num.knots <- max(10, min(50, round(length(x1)/4)))
    X <- cbind(x1, x2)
    dup.inds <- (1:nrow(X))[dup.matrix(X) == T]
    if (length(dup.inds) > 0) 
        X <- X[-dup.inds, ]
    knots <- cluster::clara(X, num.knots)$medoids
return(knots)
}

tool_exec <- function(in_params, out_params) {
    #### Load Library for Analysis ####
    if (!requireNamespace("SemiPar", quietly = TRUE))
      install.packages("SemiPar")
 
    library(SemiPar)

    #### Get Input Parameters ####
    input_features <- in_params[[1]]
    input_predictions <- in_params[[2]]
    dep_variable <- in_params[[3]]
    nonlin_variables <- in_params[[4]]
    lin_variables <- in_params[[5]]
    input_knots <- in_params[[6]]


    #### Import Dataset to Dataframe ####
    fc <- arc.open(input_features)
    df <- arc.select(fc, c(dep_variable, nonlin_variables, lin_variables))
    df['x'] <- arc.shape(df)$x
    df['y'] <- arc.shape(df)$y

    #### Import Knots to DataFrame ####
    if (is.null(input_knots))
    {
        knots.est <- make_knots(df$x, df$y)
    }
    else
    {
        knots_df <- arc.select(arc.open(input_knots))
        knots.est <- make_knots(arc.shape(knots_df)$x, arc.shape(knots_df)$y)
    }

    #### Create Spatial Effect ####
    fxy <- 'f(x, y, knots=knots.est)'

    #### Create Non-Linear Params ####
    nonlin_params <- nonlin_variables
    for (i in 1:length(nonlin_params)){
        nonlin_params[i] <- paste("f(", nonlin_params[i], ")", sep = "")
    }
    nonlin_params <- paste(nonlin_params, collapse = " + ")

    #### Create Linear Params ####
    lin_params <- paste(lin_variables, collapse = " + ")

    #### Create Formula and Fit SemiPar ####
    e <- as.list(df)
    e$knots.est <- knots.est
    all_params <- paste(fxy, nonlin_params, lin_params, sep = " + ")
    all_params <- paste(dep_variable, " ~ ", all_params)
    tryCatch({
      attach(e)
      form <- as.formula(all_params)
      fit <- spm(form, family="binomial")
    }, finally <- detach(e))

    #### Prediction ####
    pred_fc <- arc.open(input_predictions)
    pred_df <- arc.select(pred_fc, c(nonlin_variables, lin_variables))
    pred_df['x'] <- arc.shape(pred_df)$x
    pred_df['y'] <- arc.shape(pred_df)$y
    pred <- with(pred_df,  predict.spm(fit, newdata = pred_df, se = TRUE))

    #### Write Output ####
    pred_df['pred_1'] <- pred$fit
    pred_df['back_pred_1'] <- 1 / (1+exp(-pred$fit))
    pred_df['stderr_1'] <- pred$se
    pred_df['back_se1l'] <- 1 / (1+exp(-pred$fit+1.96*pred$se))
    pred_df['back_se1u'] <- 1/(1+exp(-pred$fit-1.96*pred$se))

    output_features <- out_params[[1]]
    arc.write(output_features, pred_df)

    return(out_params)
}
