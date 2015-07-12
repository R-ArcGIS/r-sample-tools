tool_exec <- function(in_params, out_params)
{
### handle input gp tool paramaters
  input_dataset_name <- in_params[[1]]
  primary_field <- in_params[[2]]
  fields <- in_params[[3]]

# --- import dataset to dataframe ---
  fc <- arc.open(input_dataset_name)
  df <- arc.select(fc, c(primary_field, fields))

#------- begin common R script ----------
# Make Formula
  form <- as.formula(paste(primary_field, "~", paste(fields, collapse="+")))
  cat("formula = ")
  print(form, showEnv=FALSE)

  m <- glm(form, data=df, family=binomial("logit"))

# display summary as GP Messages
  print(summary(m))

# append predicted and residual

  df["predicted"] <- fitted(m)
  r <- fitted(m) - df[primary_field]
  df["residual"] <- r

#------- end

### handle gp tool outputs
  output_dataset_name <- out_params[[1]]
# export result dataframe to FeatureClass
  arc.write(output_dataset_name, df)

# update derived parameter
  out_params[[2]] <- m$aic
  return(out_params)
}
