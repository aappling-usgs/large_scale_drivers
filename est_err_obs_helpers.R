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
    select(Water_Year, Value=y, Var, SD)
}

#' @examples
#' est_sum_var(dat_A=est_var(dat, 'NO23_N'), dat_B=est_var(dat, 'TKN_N'))
#' est_sum_var(dat_A=est_var(dat, 'NO23_N'), dat_B=est_var(dat, 'TKN_N'), assume_cor=0)
est_sum_var <- function(dat_A, dat_B, assume_cor=1, alpha=0.05) {
  vars <- 
    full_join(
      select(dat_A, Water_Year, Value, Var), select(dat_B, Water_Year, Value, Var),
      by='Water_Year',
      suffix=c('.A','.B')) %>% 
    filter(complete.cases(.))
  vars %>%
    mutate(
      Value=Value.A + Value.B,
      Var=Var.A + Var.B + 2*assume_cor*sqrt(Var.A*Var.B),
      SD=sqrt(Var)) %>%
    select( Water_Year, Value, Var, SD)
}

#' @examples
#' est_ratio_var(dat_num=est_var(dat, 'OrthoP_P'), dat_denom=est_var(dat, 'TP_P'))
#' est_ratio_var(dat_num=est_var(dat, 'OrthoP_P'), dat_denom=est_var(dat, 'TP_P'), assume_cor=0)
est_ratio_var <- function(dat_num, dat_denom, assume_cor=1, alpha=0.05) {
  vars <- 
    full_join(
      select(dat_num, Water_Year, Value, Var), select(dat_denom, Water_Year, Value, Var),
      by='Water_Year',
      suffix=c('.num','.den')) %>% 
    filter(complete.cases(.))
  vars %>%
    mutate(
      Value=Value.num/Value.den,
      Var=(Value.num/Value.den)*(Var.num/(Value.num^2) + Var.den/(Value.den^2) - 2*assume_cor*sqrt(Var.num*Var.den)/(Value.num*Value.den)),
      SD=sqrt(Var)) %>%
    select(Water_Year, Value, Var, SD)
}

#' @examples
#' est_obs_err(est_var(dat, 'OrthoP_P'))
#' est_obs_err(est_ratio_var(dat_num=est_var(dat, 'OrthoP_P'), dat_denom=est_var(dat, 'TP_P')))
#' est_obs_err(
#'   est_ratio_var(
#'     dat_num=est_var(dat, 'NO23_N'), 
#'     dat_denom=est_sum_var(dat_A=est_var(dat, 'NO23_N'), dat_B=est_var(dat, 'TKN_N'), assume_cor=1),
#'     assume_cor=1))
est_obs_err <- function(dat, alpha=0.05) {
  dat %>%
    summarize(
      Value = mean(Value),
      Var = mean(Var),
      SD = sqrt(Var))
}

#' @examples
#' MARB <- get_obs_errs(dat2)
get_obs_errs <- function(dat) {
  # get the straight solutes
  solutes <- grep('_Flux_ty', names(dat), value=TRUE) %>% gsub('_Flux_ty','',.)
  sol_vars <- lapply(setNames(nm=solutes), function(sol) mutate(est_var(dat, sol), Solute=sol))
  
  # add TN if needed
  if(!('TN_N' %in% solutes) && all(c('NO23_N','TKN_N') %in% solutes)) {
    sol_vars$TN_N0 <- mutate(est_sum_var(sol_vars$NO23_N, sol_vars$TKN_N, assume_cor=0), Solute='TN_N', TN_Cor=0)
    sol_vars$TN_N1 <- mutate(est_sum_var(sol_vars$NO23_N, sol_vars$TKN_N, assume_cor=1), Solute='TN_N', TN_Cor=1)
  }
  solutes <- names(sol_vars)
  
  # add any ratios available
  all_ratios <- data_frame(
    Ratio = c('NO23_TN','NO23_TN' ,'NO23_TN' ,'OrthoP_TP','Si_TN','Si_TN' ,'Si_TN' ,'TP_Si','TN_TP','TN_TP' ,'TN_TP'),
    RatioC= c('NO23_TN','NO23_TN0','NO23_TN1','OrthoP_TP','Si_TN','Si_TN0','Si_TN1','TP_Si','TN_TP','TN0_TP','TN1_TP'),
    Num =   c('NO23_N', 'NO23_N',  'NO23_N',  'OrthoP_P', 'SiO2', 'SiO2',  'SiO2',  'TP_P', 'TN_N', 'TN_N0', 'TN_N1'), 
    Den =   c('TN_N',   'TN_N0',   'TN_N1',   'TP_P',     'TN_N', 'TN_N0', 'TN_N1', 'SiO2', 'TP_P', 'TP_P',  'TP_P')) %>%
    mutate(TN_Cor = ifelse(grepl('TN0', RatioC), 0, ifelse(grepl('TN1', RatioC), 1, NA)))
  ratios <- all_ratios %>% filter(mapply(function(s1, s2) all(c(s1,s2) %in% solutes), Num, Den))
  for(r in seq_len(nrow(ratios))) {
    sol_vars[[paste0(ratios$RatioC[r], '_0')]] <- mutate(
      est_ratio_var(sol_vars[[ratios$Num[r]]], sol_vars[[ratios$Den[r]]], assume_cor=0), 
      Solute=ratios$Ratio[r], TN_Cor=ratios$TN_Cor[r], Ratio_Cor=0)
    sol_vars[[paste0(ratios$RatioC[r], '_1')]] <- mutate(
      est_ratio_var(sol_vars[[ratios$Num[r]]], sol_vars[[ratios$Den[r]]], assume_cor=1), 
      Solute=ratios$Ratio[r], TN_Cor=ratios$TN_Cor[r], Ratio_Cor=1)
  }

  bind_rows(sol_vars) %>%
    group_by(Solute, TN_Cor, Ratio_Cor) %>%
    do(est_obs_err(.)) %>%
    ungroup() %>%
    return()
}
