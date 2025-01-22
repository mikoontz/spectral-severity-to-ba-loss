# unloadNamespace("cpi")
# remotes::install_github(repo = "mikoontz/cpi@resampling-fix") # install from my dev
library(ggplot2)
library(cpi)
library(caret)
library(rPref)

# set.seed(20240724)
# 
# # read data and set up spatial folds
# # https://spatialsample.tidymodels.org/articles/spatialsample.html
# ard = sf::st_read("data/ARD_08262024.gpkg") |>
#   dplyr::filter(!is.na(pcnt_ba_mo)) |> 
#   dplyr::select(-lat, -aspectRad)

set.seed(20250121)

# read data and set up spatial folds
# https://spatialsample.tidymodels.org/articles/spatialsample.html
ard = sf::st_read("data/ARD_01212025.gpkg") |>
  dplyr::filter(!is.na(pcnt_ba_mo)) |> 
  dplyr::select(-lat, -aspectRad)

ard_for_task_by_ecoregion = ard |>
  dplyr::filter(!(ecoregion %in% c("California interior chaparral and woodlands", "Great Basin shrub steppe"))) |> 
  dplyr::group_by(ecoregion) |> 
  dplyr::group_split() |> 
  purrr::map(.f = spatialsample::spatial_clustering_cv, v = 5, .progress = TRUE)

x <- ard_for_task_by_ecoregion[[1]]

ard_with_spatial_folds = ard_for_task_by_ecoregion |>
  purrr::map(
    .f = \(x) {
      out <- x |> 
        purrr::pmap(
          .f = \(id, splits) {
            
            spatial_fold <- id
            assessment_data =
              splits |>
              rsample::assessment() |>
              sf::st_drop_geometry()
            
            return(cbind(assessment_data, spatial_fold))
          }
        ) |>
        data.table::rbindlist() |>
        dplyr::mutate(spatial_fold = factor(spatial_fold)) |>
        tibble::as_tibble()
      
      out
    }
  )

features = ard |> 
  sf::st_drop_geometry() |> 
  dplyr::select(-c(PlotID, YrFireName, Dataset, pcnt_ba_mo, UniqueID)) |> 
  colnames()

target = "pcnt_ba_mo"

full_rf_formula = glue::glue("{target} ~ {paste(features, collapse = ' + ')}")

hyperparameters_full =
  expand.grid(
    variablesPerSplit = 3:12, 
    bagFraction = c(0.4, 0.5, (1 - 1/exp(1)), 0.7, 0.8, 0.9),
    minLeafPopulation = c(1, 5, 10, 25, 50, 60, 70, 80, 90, 100, 125, 150)
  )

# hyperparameters = hyperparameters_full[sample(x = 1:nrow(hyperparameters_full), size = 10), ]
hyperparameters = hyperparameters_full

tune_validate_varselect_assess = function(variablesPerSplit, 
                                          bagFraction, 
                                          minLeafPopulation, 
                                          resampling_approach,
                                          ecoregion) {
  
  library(mlr3verse)
  
  # Set up the leaner with the (currently) 3 hyperparameters
  learner_sev_biomass <- mlr3::lrn(
    .key = "regr.ranger",
    mtry = variablesPerSplit,
    num.trees = 300,
    sample.fraction = bagFraction,
    replace = FALSE,
    min.node.size = minLeafPopulation,
    num.threads = 1,
    keep.inbag = TRUE
  )
  
  # Set up the task using the formula notation with the full set of predictors
  task_sev_biomass =
    mlr3::as_task_regr(
      x = as.formula(full_rf_formula),
      data = ard_with_spatial_folds[, c(target, features)],
      id = target
    )
  
  # Set up and instantiate the resampler using the known spatial folds as the
  # groups
  resampler_sev_biomass = rsmp("custom_cv")
  resampler_sev_biomass$instantiate(
    task_sev_biomass, 
    f = ard_with_spatial_folds$spatial_fold
  )
  
  if (resampling_approach == "resampler") {
    # Calculate conditional predictive impact
    cpi_results = cpi::cpi(
      task = task_sev_biomass,
      learner = learner_sev_biomass,
      measure = "regr.mse",
      resampling = resampler_sev_biomass,
      test = "t"
    )
  } else if (resampling_approach == "oob") {
    cpi_results = cpi::cpi(
      task = task_sev_biomass,
      learner = learner_sev_biomass,
      measure = "regr.mse",
      resampling = "oob",
      test = "t"
    )
  }
  
  # Spatially cross validated model assessment the {mlr3} way
  assessment_full = resample(
    task = task_sev_biomass, 
    learner = learner_sev_biomass, 
    resampling = resampler_sev_biomass
  )
  
  # Pull out a specific model skill metric aggregated across the spatial folds
  obs_preds_full = assessment_full$predictions(predict_sets = "test") |> 
    purrr::map(.f = \(x) tibble::tibble(
      obs = getElement(x, "truth"), 
      preds = getElement(x, "response")
    )
    ) |> 
    data.table::rbindlist()
  
  pred_full = obs_preds_full$preds
  obs_full = obs_preds_full$obs
  
  rmse_full = assessment_full$aggregate(measures = msr("regr.rmse"))
  r2_full = assessment_full$aggregate(measures = msr("regr.rsq"))
  mae_full = assessment_full$aggregate(measures = msr("regr.mae"))
  mse_full = assessment_full$aggregate(measures = msr("regr.mse"))
  
  rmse_full_overall = caret::RMSE(pred = pred_full, obs = obs_full)
  r2_full_overall = caret::R2(pred = pred_full, obs = obs_full)
  mae_full_overall = caret::MAE(pred = pred_full, obs = obs_full)
  mse_full_overall = rmse_full_overall^2
  
  # Initial pass at finding important variables
  important_variables = cpi_results |> 
    dplyr::filter(ci.lo > 0) |> 
    dplyr::pull(Variable)
  
  # Create the formula that could be used for the next iteration of the 
  # spatial cross validation (using the reduced set of only important predictors)
  important_variable_rf_formula = glue::glue(
    "{target} ~ {paste(important_variables, collapse = ' + ')}"
  )
  
  if (variablesPerSplit <= length(important_variables)) {
    # Set up the task using the formula notation with the full set of predictors
    task_sev_biomass_important_variables =
      mlr3::as_task_regr(
        x = as.formula(important_variable_rf_formula),
        data = ard_with_spatial_folds[, c(target, important_variables)],
        id = target
      )
    
    resampler_sev_biomass_important_variables = rsmp("custom_cv")
    resampler_sev_biomass_important_variables$instantiate(
      task_sev_biomass_important_variables, 
      f = ard_with_spatial_folds$spatial_fold
    )
    
    # Spatially cross validated model assessment the {mlr3} way
    assessment_important_variables = resample(
      task = task_sev_biomass_important_variables, 
      learner = learner_sev_biomass, 
      resampling = resampler_sev_biomass_important_variables
    )
    
    obs_preds_important_variables = assessment_important_variables$predictions(predict_sets = "test") |> 
      purrr::map(.f = \(x) tibble::tibble(
        obs = getElement(x, "truth"), 
        preds = getElement(x, "response")
      )
      ) |> 
      data.table::rbindlist()
    
    pred_important_variables = obs_preds_important_variables$preds
    obs_important_variables = obs_preds_important_variables$obs
    
    # Pull out a specific model skill metric aggregated across the spatial folds
    rmse_important_variables = assessment_important_variables$aggregate(
      measures = msr("regr.rmse")
    )
    
    r2_important_variables = assessment_important_variables$aggregate(
      measures = msr("regr.rsq")
    )
    
    mae_important_variables = assessment_important_variables$aggregate(
      measures = msr("regr.mae")
    )
    
    mse_important_variables = assessment_important_variables$aggregate(
      measures = msr("regr.mse")
    )
    
    rmse_important_variables_overall = caret::RMSE(pred = pred_important_variables, obs = obs_important_variables)
    r2_important_variables_overall = caret::R2(pred = pred_important_variables, obs = obs_important_variables)
    mae_important_variables_overall = caret::MAE(pred = pred_important_variables, obs = obs_important_variables)
    mse_important_variables_overall = rmse_important_variables_overall^2
    
    
  } else {
    
    rmse_important_variables = NA
    r2_important_variables = NA
    mae_important_variables = NA
    mse_important_variables = NA
    rmse_important_variables_overall = NA
    r2_important_variables_overall = NA
    mae_important_variables_overall = NA
    mse_important_variables_overall = NA
    
  }
  
  # Build the final output table that includes the unique set of hyperparameters
  # the CPI resutls, the new formula to use only important variables, the 
  # number of important variables (so we don't bother trying to re-run the 
  # reduced set of variables when the mtry/variablesPerSplit hyperparameter is greater than the
  # number of variables in the model)
  out = tibble::tibble(
    mtry = variablesPerSplit,
    sample.fraction = bagFraction,
    min.node.size = minLeafPopulation,
    rmse_full = rmse_full,
    r2_full = r2_full,
    mae_full = mae_full,
    mse_full = mse_full,
    rmse_full_overall = rmse_full_overall,
    r2_full_overall = r2_full_overall,
    mae_full_overall = mae_full_overall,
    mse_full_overall = mse_full_overall,
    n_important_variables = length(important_variables),
    important_variable_rf_formula = important_variable_rf_formula,
    rmse_important_variables = rmse_important_variables,
    r2_important_variables = r2_important_variables,
    mae_important_variables = mae_important_variables,
    mse_important_variables = mse_important_variables,
    rmse_important_variables_overall = rmse_important_variables_overall,
    r2_important_variables_overall = r2_important_variables_overall,
    mae_important_variables_overall = mae_important_variables_overall,
    mse_important_variables_overall = mse_important_variables_overall,
    cpi_results = list(cpi_results),
  ) |> 
    tidyr::unnest(cols = cpi_results)
  
  return(out)
}


variablesPerSplit = hyperparameters$variablesPerSplit[1]
bagFraction = hyperparameters$bagFraction[1]
minLeafPopulation = hyperparameters$minLeafPopulation[1]
resampling_approach = "resampler"

test_out = tune_validate_varselect_assess(
  variablesPerSplit = variablesPerSplit,
  bagFraction = bagFraction,
  minLeafPopulation = minLeafPopulation,
  resampling_approach = resampling_approach
)

# # Define the learner to be a {ranger} regression and give it the tuned hyperparameters
tictoc::tic()
future::plan(future::multisession, workers = 10)
# future::plan("sequential")
results_list = furrr::future_pmap(
  .l = hyperparameters, 
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE),
  .f = tune_validate_varselect_assess,
  resampling_approach = "resampler"
)
tictoc::toc()

results = data.table::rbindlist(results_list)

dir.create("data/processed")
data.table::fwrite(x = results, file = "data/processed/conditional-predictive-impact-results_v4.0.csv")

cpi_results_full = data.table::fread("data/processed/conditional-predictive-impact-results_v4.0.csv")

cpi_results = cpi_results_full |>
  dplyr::select(-Variable, -CPI, -SE, -test, -statistic, -estimate, -p.value, -ci.lo) |>
  unique()

# Plot Pareto frontier of cross-fold mean R2 and overall R2 of important variable reduced model
r2_pareto_front = rPref::psel(
  df = cpi_results, 
  pref = rPref::high(r2_important_variables_overall) * rPref::high(r2_important_variables)
)

ggplot2::ggplot(cpi_results, ggplot2::aes(x = r2_important_variables_overall, y = r2_important_variables)) + 
  ggplot2::geom_point() +
  ggplot2::geom_point(data = r2_pareto_front, color = "red", size = 3) +
  ggplot2::theme_bw()

r2_pareto_front |> na.omit()

top_results = r2_pareto_front |> na.omit()
top_results = structure(
  list(
    mtry = c(4L, 3L, 4L), 
    sample.fraction = c(0.7, 0.632120558828558, 0.7), 
    min.node.size = c(50L, 60L, 60L), 
    rmse_full = c(0.268083261313512, 0.26872065739475, 0.267812714273246), 
    r2_full = c(-0.303325414602873, -0.336928485033941, -0.306773573957064), 
    mae_full = c(0.205434762569152, 0.206793084624741, 0.205591477394943), 
    mse_full = c(0.0758998325166776, 0.0764617624306148, 0.0758349130378365), 
    rmse_full_overall = c(0.283676467249458, 0.284965447853962, 0.284733322894045), 
    r2_full_overall = c(0.543345881715101, 0.539332194225563, 0.539985005772275), 
    mae_full_overall = c(0.213296590576411, 0.215123118254161, 0.214940048108721), 
    mse_full_overall = c(0.0804723380711328, 0.0812053064706093, 0.0810730651662846), 
    n_important_variables = c(5L, 3L, 5L), 
    important_variable_rf_formula = c("pcnt_ba_mo ~ dndvi + dred + northness + post_swir2swir1 + zScorePrecip1", "pcnt_ba_mo ~ dswir2swir1 + northness + zScorePrecip1", "pcnt_ba_mo ~ dndvi + dswir2 + northness + post_swir2nir + zScorePrecip1"),
    rmse_important_variables = c(0.249324154440508, 0.254920503135365, 0.252569880773709), 
    r2_important_variables = c(0.118508630120102, 0.267309439297053, -0.290180212311028), 
    mae_important_variables = c(0.178943386446459, 0.181787332462959, 0.18018098713846), 
    mse_important_variables = c(0.0661507334466124, 0.0700145686587712, 0.0662784542817479), 
    rmse_important_variables_overall = c(0.274822298181919, 0.282639313677022, 0.271625889583909), 
    r2_important_variables_overall = c(0.572848080518024, 0.549886935546503, 0.582282453974439), 
    mae_important_variables_overall = c(0.194134380637982, 0.198871633350667, 0.190291467955942), 
    mse_important_variables_overall = c(0.0755272955779917, 0.0798849816358181, 0.0737806238922502)
  ), 
  row.names = c(NA, -3L), 
  class = c("data.table", "data.frame")
)









results = data.table::fread("data/processed/conditional-predictive-impact-results_v4.0.csv")

results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique() |> 
  dplyr::filter(rmse_full== min(rmse_full, na.rm = TRUE))

results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique() |> 
  dplyr::filter(r2_full == max(r2_full, na.rm = TRUE))

results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique() |> 
  dplyr::filter(mae_full == min(mae_full, na.rm = TRUE))

results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique() |> 
  dplyr::filter(mse_full == min(mse_full, na.rm = TRUE))

### Reduced
best_fit = results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique() |> 
  dplyr::filter(rmse_important_variables == min(rmse_important_variables, na.rm = TRUE))

best_fit = results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique() |> 
  dplyr::filter(r2_important_variables == max(r2_important_variables, na.rm = TRUE))

best_fit = results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique() |> 
  dplyr::filter(mae_important_variables == min(mae_important_variables, na.rm = TRUE))

best_fit = results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique() |> 
  dplyr::filter(mse_important_variables == min(mse_important_variables, na.rm = TRUE))

model_skill_results = results |> 
  dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, mae_full, mse_full, n_important_variables, important_variable_rf_formula, rmse_important_variables, r2_important_variables, mae_important_variables, mse_important_variables) |> 
  unique()

model_skill_results |> dplyr::arrange(dplyr::desc(r2_important_variables))

fm1 = ranger::ranger(
  formula = as.formula(best_fit$important_variable_rf_formula), 
  data = sf::st_drop_geometry(ard[, c(target, features)]), 
  num.trees = 1000, 
  mtry = best_fit$mtry, 
  min.node.size = best_fit$min.node.size, 
  sample.fraction = best_fit$sample.fraction
)

spatial_folds = unique(ard_with_spatial_folds$spatial_fold)
model_assessment_data = vector(mode = "list", length = length(spatial_folds))

for(i in seq_along(spatial_folds)) {
  
  print(i)
  train_data = ard_with_spatial_folds |> 
    sf::st_drop_geometry() |> 
    dplyr::filter(spatial_fold != spatial_folds[i])
  
  test_data = ard_with_spatial_folds |> 
    sf::st_drop_geometry() |> 
    dplyr::filter(spatial_fold == spatial_folds[i])
  
  
  fm = ranger::ranger(
    formula = as.formula(best_fit$important_variable_rf_formula), 
    data = train_data, 
    num.trees = 1000, 
    mtry = best_fit$mtry, 
    min.node.size = best_fit$min.node.size, 
    sample.fraction = best_fit$sample.fraction
  )
  
  out = tibble::tibble(
    obs = test_data$pcnt_ba_mo,
    preds = predict(object = fm, data = test_data)$predictions
  )
  
  model_assessment_data[[i]] = out
}

model_assessment_data = dplyr::bind_rows(model_assessment_data)
caret::R2(pred = model_assessment_data$preds, obs = model_assessment_data$obs)

ggplot(model_assessment_data, aes(x = obs, y = preds)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw()


full_model = ranger::ranger(
  formula = full_rf_formula, 
  data = sf::st_drop_geometry(ard[, c(target, features)]), 
  num.trees = 1000, 
  mtry = best_fit$mtry, 
  min.node.size = best_fit$min.node.size, 
  sample.fraction = best_fit$sample.fraction
)

per_variable_cpi_results = results |> 
  dplyr::filter(mtry == best_fit$mtry & min.node.size == best_fit$min.node.size & sample.fraction == best_fit$sample.fraction) |> 
  dplyr::arrange(dplyr::desc(ci.lo)) |> 
  dplyr::filter(ci.lo > 0)

calc_plot_data = function(var_names, fitted_model, variable_order) {
  yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata)$predictions)
  plot_data = purrr::map(
    .x = var_names,
    .f = \(J) {
      ale = ALEPlot::ALEPlot(
        X = sf::st_drop_geometry(ard[, fitted_model$forest$independent.variable.names]), 
        X.model = fitted_model, 
        J = J, 
        pred.fun = yhat
      )
      
      out = tibble::tibble(var = J, x = ale$x.values, y = ale$f.values)
      return(out)
    }
  ) |> 
    data.table::rbindlist() |> 
    dplyr::mutate(var = factor(var, levels = variable_order))
  
  return(plot_data)
}

plot_data_important_vars = calc_plot_data(
  var_names = per_variable_cpi_results$Variable,
  fitted_model = fm1,
  variable_order = per_variable_cpi_results$Variable
)


plot_data_full = calc_plot_data(
  var_names = per_variable_cpi_results$Variable,
  fitted_model = full_model,
  variable_order = per_variable_cpi_results$Variable
)

ggplot(plot_data_important_vars, aes(x = x, y = y)) +
  geom_line() +
  facet_wrap(facets = "var", scales = "free") +
  theme_bw()

ggplot(plot_data_full, aes(x = x, y = y)) +
  geom_line() +
  facet_wrap(facets = "var", scales = "free") +
  theme_bw()

summary(lm(formula = pcnt_ba_mo ~ zScorePrecip1, data = ard[, c("pcnt_ba_mo", "zScorePrecip1")]))

ggplot(ard, aes(x = zScorePrecip1, y = pcnt_ba_mo)) +
  geom_point() +
  geom_smooth()

ggplot(ard, aes(x = rdnbr, y = zScorePrecip1, col = pcnt_ba_mo)) +
  geom_point()

ard |> 
  dplyr::mutate(zscoreprecip1_fct = zScorePrecip1 > 0) |> 
  ggplot(aes(x = rdnbr, y = pcnt_ba_mo)) +
  geom_point() +
  facet_wrap(facets = "zscoreprecip1_fct") +
  geom_smooth()

ggplot(ard, aes(x = dndvi, y = pcnt_ba_mo, col = zScorePrecip1 > 0)) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous(limits = c(-100, 700)) +
  theme_bw()

ggplot(ard, aes(x = dswir2swir1, y = pcnt_ba_mo, col = cut(zScorePrecip1, breaks = c(-Inf, 0, Inf)))) +
  geom_point() +
  geom_smooth() +
  # scale_x_continuous(limits = c(-1, 0.1)) +
  theme_bw() +
  labs(col = "Precipitation category (z-score positive or negative)",
       y = "Basal area loss (%)",
       x = "pre-fire SWIR2/SWIR1 minus post-fire SWIR2/SWIR1")

ggplot(ard, aes(x = dswir2swir1, y = pcnt_ba_mo, col = cut(northness, breaks = c(-Inf, 0, Inf)))) +
  geom_point() +
  geom_smooth() +
  # scale_x_continuous(limits = c(-1, 0.1)) +
  theme_bw() +
  labs(col = "Northness (negative = south; positive = north)",
       y = "Basal area loss (%)",
       x = "pre-fire SWIR2/SWIR1 minus post-fire SWIR2/SWIR1")

ggplot(ard, aes(x = dswir2swir1, y = pcnt_ba_mo, col = cut(zScorePrecip1, breaks = c(-Inf, -1, -0.5, 0, 0.5, 1, Inf)))) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous(limits = c(-1, 0.1)) +
  theme_bw()

precip_mort_plot_data = ard |> 
  dplyr::group_by(zScorePrecip1 > 0) |> 
  dplyr::summarize(pcnt_ba_mo = mean(pcnt_ba_mo))

plot(precip_mort_plot_data[, "zScorePrecip1 > 0"])


library(ggplot2)

ggplot(plot_data, aes(x = x, y = y)) +
  geom_line() +
  facet_wrap(facets = "var", scales = "free") +
  theme_bw()

rdnbr_post_nbr_two_way_ale = ALEPlot::ALEPlot(
  X = sf::st_drop_geometry(ard[, fm1$forest$independent.variable.names]), 
  X.model = fm1, 
  J = c("rdnbr", "post_nbr"), 
  pred.fun = yhat
)

image(
  rdnbr_post_nbr_two_way_ale$x.values[[1]], 
  rdnbr_post_nbr_two_way_ale$x.values[[2]], 
  rdnbr_post_nbr_two_way_ale$f.values, 
  xlab = "RdNBR", 
  ylab = "post_nbr", 
  col = hcl.colors(n = 100, palette = "Blue-Red")
)

contour(
  rdnbr_post_nbr_two_way_ale$x.values[[1]], 
  rdnbr_post_nbr_two_way_ale$x.values[[2]], 
  rdnbr_post_nbr_two_way_ale$f.values, 
  add=TRUE, 
  drawlabels=TRUE
)


zscoreprecip1_rdnbr_two_way_ale = ALEPlot::ALEPlot(
  X = sf::st_drop_geometry(ard[, fm1$forest$independent.variable.names]), 
  X.model = fm1, 
  J = c("zScorePrecip1", "rdnbr"), 
  pred.fun = yhat
)

image(
  zscoreprecip1_rdnbr_two_way_ale$x.values[[1]], 
  zscoreprecip1_rdnbr_two_way_ale$x.values[[2]], 
  zscoreprecip1_rdnbr_two_way_ale$f.values, 
  xlab = "zScorePrecip1", 
  ylab = "RdNBR", 
  col = hcl.colors(n = 100, palette = "Blue-Red")
)

contour(
  rdnbr_post_nbr_two_way_ale$x.values[[1]], 
  rdnbr_post_nbr_two_way_ale$x.values[[2]], 
  rdnbr_post_nbr_two_way_ale$f.values, 
  add=TRUE, 
  drawlabels=TRUE
)

# results |> 
#   dplyr::select(mtry, sample.fraction, min.node.size, rmse_full, r2_full, rmse_important_variables, r2_important_variables)

task = tsk("penguins")
learner = lrn("classif.rpart")
resampling = rsmp("cv")

# Explicitly instantiate the resampling for this task for reproduciblity
set.seed(123)
resampling$instantiate(task)

rr = resample(task, learner, resampling)
print(rr)

# Retrieve performance
rr$score(msr("classif.ce"))
rr$aggregate(msr("classif.ce"))

fm1 = ranger::ranger(
  formula = as.formula(best_fit$important_variable_rf_formula), 
  data = sf::st_drop_geometry(ard[, c(target, features)]), 
  num.trees = 1000, 
  mtry = best_fit$mtry, 
  min.node.size = best_fit$min.node.size, 
  sample.fraction = best_fit$sample.fraction
)

learner_sev_biomass <- mlr3::lrn(
  .key = "regr.ranger",
  mtry = best_fit$mtry,
  num.trees = 300,
  sample.fraction = best_fit$sample.fraction,
  replace = FALSE,
  min.node.size =  best_fit$min.node.size,
  num.threads = 1,
  keep.inbag = TRUE
)

# Full model task
task_sev_biomass =
  mlr3::as_task_regr(
    x = as.formula(full_rf_formula),
    data = ard_with_spatial_folds[, c(target, features)],
    id = target
  )

# Set up the task using the formula notation with the full set of predictors
task_sev_biomass_important_variables =
  mlr3::as_task_regr(
    x = as.formula(best_fit$important_variable_rf_formula),
    data = ard_with_spatial_folds[, c(target, features)],
    id = target
  )

# Set up and instantiate the resampler using the known spatial folds as the
# groups
resampler_sev_biomass = rsmp("custom_cv")
resampler_sev_biomass$instantiate(
  task_sev_biomass, 
  f = ard_with_spatial_folds$spatial_fold
)

# Spatially cross validated model assessment the {mlr3} way
assessment_important_variables = resample(
  task = task_sev_biomass, 
  learner = learner_sev_biomass, 
  resampling = resampler_sev_biomass
)

test = assessment_important_variables$predictions() |> 
  purrr::map(.f = \(x) {
    out = tibble::tibble(obs = x$truth, pred = x$response)
    return(out)
  }) |> 
  data.table::rbindlist()

caret::R2(pred = test$pred, obs = test$obs)


|> purrr::map(.f = \(x) {
  out = tibble::tibble(x[[1]][, c("row_id", "truth", "response")])
  return(out)
})

resampler_sev_biomass2 = rsmp("custom_cv")
resampler_sev_biomass2$instantiate(
  task_sev_biomass_important_variables, 
  f = ard_with_spatial_folds$spatial_fold
)

assessment_important_variables3 = resample(
  task = task_sev_biomass_important_variables, 
  learner = learner_sev_biomass, 
  resampling = resampler_sev_biomass2
)

assessment_important_variables2 = resample(
  task = task_sev_biomass_important_variables, 
  learner = learner_sev_biomass, 
  resampling = resampler_sev_biomass
)

assessment_important_variables$score(measures = msr("regr.rsq"))
mean(assessment_important_variables$score(measures = msr("regr.rsq"))$regr.rsq)
assessment_important_variables$aggregate(measures = msr("regr.rsq"))

mean(assessment_important_variables2$score(measures = msr("regr.rsq"))$regr.rsq)
assessment_important_variables2$aggregate(measures = msr("regr.rsq"))

mean(assessment_important_variables3$score(measures = msr("regr.rsq"))$regr.rsq)
assessment_important_variables3$score(measures = msr("regr.rsq"))
assessment_important_variables3$aggregate(measures = msr("regr.rsq"))


ard_with_spatial_folds_sf = ard |> 
  dplyr::left_join(ard_with_spatial_folds[, c("UniqueID", "spatial_fold")]) |> 
  sf::st_transform(5070)

test = ard_with_spatial_folds_sf |> dplyr::filter(spatial_fold == "Fold01") |> sf::st_transform(5070)
mapview::mapview(test)

usa = USAboundaries::us_states() |> dplyr::filter(jurisdiction_type == "state" & !(state_name %in% c("Alaska", "Hawaii"))) |> sf::st_transform(sf::st_crs(5070))
plot(ard_with_spatial_folds_sf[, "spatial_fold"])
plot(usa$geometry, add = TRUE)
mapview::mapview(ard_with_spatial_folds_sf[, "spatial_fold"])
plot(ard_with_spatial_folds_sf |> dplyr::filter(spatial_fold == "Fold01") |> sf::st_geometry(), add = TRUE, col = "red")
plot(ard_with_spatial_folds_sf |> dplyr::filter(spatial_fold == "Fold05") |> sf::st_geometry(), add = TRUE, col = "blue")
plot(ard_with_spatial_folds_sf |> dplyr::filter(spatial_fold == "Fold02") |> sf::st_geometry(), add = TRUE, col = "purple")
plot(ard_with_spatial_folds_sf |> dplyr::filter(spatial_fold == "Fold07") |> sf::st_geometry(), add = TRUE, col = "darkgreen")

plot(usa$geometry)
plot(ard_with_spatial_folds_sf[, "spatial_fold"], add = TRUE)
