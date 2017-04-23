#' @examples
#' est_var(dat, 'OrthoP_P')
#' est_var(dat, 'TP_P')
est_var <- function(dat, nut, alpha=0.05) {
  col_y <- paste0(nut, '_Flux_ty')
  col_ymin <- paste0(nut, '_Flux_lo_ty')
  col_ymax <- paste0(nut, '_Flux_hi_ty')
  
  qci <- qt(1 - alpha/2, df=Inf) # use df=Inf because we don't know the df for each load estimate, and it varies
  select_(dat, .dots=c('Water_Year', y=col_y, ymin=col_ymin, ymax=col_ymax)) %>% 
    filter(., complete.cases(.)) %>%
    mutate(
      Solute = nut,
      # rloadest CIs are based on lognormal distributions, but I think we want
      # linear space obs error, so approximating the CIs that way. They're not
      # too bad; i.e., CI_dn is only a little smaller than CI_up
      CI_up = ymax - y,
      CI_dn = y - ymin,
      CI_95 = ymax - ymin,
      SD = CI_95 / (2*qci),
      Var = SD^2
    ) %>%
    select(Solute, Water_Year, Value=y, Var, SD)
}

#' @examples
#' est_sum_var(dat_A=est_var(dat, 'NO23_N'), dat_B=est_var(dat, 'TKN_N'), sum_name='TN_N')
#' est_sum_var(dat_A=est_var(dat, 'NO23_N'), dat_B=est_var(dat, 'TKN_N'), sum_name='TN_N', assume_cor=0)
est_sum_var <- function(dat_A, dat_B, sum_name, assume_cor=1, alpha=0.05) {
  vars <- 
    full_join(
      select(dat_A, Water_Year, Value, Var), select(dat_B, Water_Year, Value, Var),
      by='Water_Year',
      suffix=c('.A','.B')) %>% 
    filter(complete.cases(.))
  vars %>%
    mutate(
      Solute=sum_name,
      Value=Value.A + Value.B,
      Var=Var.A + Var.B + 2*assume_cor*sqrt(Var.A*Var.B),
      SD=sqrt(Var)) %>%
    select(Solute, Water_Year, Value, Var, SD)
}

#' @examples
#' est_ratio_var(dat_num=est_var(dat, 'OrthoP_P'), dat_denom=est_var(dat, 'TP_P'), ratio_name='OrthoP_TP')
#' est_ratio_var(dat_num=est_var(dat, 'OrthoP_P'), dat_denom=est_var(dat, 'TP_P'), ratio_name='OrthoP_TP', assume_cor=0)
est_ratio_var <- function(dat_num, dat_denom, ratio_name, assume_cor=1, alpha=0.05) {
  vars <- 
    full_join(
      select(dat_num, Water_Year, Value, Var), select(dat_denom, Water_Year, Value, Var),
      by='Water_Year',
      suffix=c('.num','.den')) %>% 
    filter(complete.cases(.))
  vars %>%
    mutate(
      Solute=ratio_name,
      Value=Value.num/Value.den,
      Var=(Value.num/Value.den)*(Var.num/(Value.num^2) + Var.den/(Value.den^2) - 2*assume_cor*sqrt(Var.num*Var.den)/(Value.num*Value.den)),
      SD=sqrt(Var)) %>%
    select(Solute, Water_Year, Value, Var, SD)
}

#' @examples
#' est_obs_err(est_var(dat, 'OrthoP_P'))
#' est_obs_err(est_ratio_var(dat_num=est_var(dat, 'OrthoP_P'), dat_denom=est_var(dat, 'TP_P'), ratio_name='OrthoP_TP'))
#' est_obs_err(
#'   est_ratio_var(
#'     dat_num=est_var(dat, 'NO23_N'), 
#'     dat_denom=est_sum_var(dat_A=est_var(dat, 'NO23_N'), dat_B=est_var(dat, 'TKN_N'), sum_name='TN_N', assume_cor=1),
#'     ratio_name='NO23_TN', assume_cor=1))
est_obs_err <- function(dat, alpha=0.05) {
  if(exists('Var_iid', dat)) {
    # dat is ratios
    dat %>%
      summarize(
        Solute = Solute[1],
        Value = mean(Value),
        Var_iid = mean(Var_iid),
        Var_cor = mean(Var_cor),
        SD_iid = sqrt(Var_iid),
        SD_cor = sqrt(Var_cor))
  } else {
    # dat is the usual
    dat %>%
      summarize(
        Solute = Solute[1],
        Value = mean(Value),
        Var = mean(Var),
        SD = sqrt(Var))
  }
}
