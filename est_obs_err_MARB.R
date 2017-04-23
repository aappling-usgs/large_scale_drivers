# much of this code is duplicated from get_mississippi_loads, but edited to
# match the format of the total-MARB site and file

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# create directories to store the raw and reformatted data files
rawdir <- '01_import/raw'
if(!dir.exists(rawdir)) dir.create(rawdir)
mungedir <- '01_import/munged'
if(!dir.exists(mungedir)) dir.create(mungedir)

# use rvest to find & download the annual data files from the website
library(rvest)
loads_index_url <- 'https://toxics.usgs.gov/hypoxia/mississippi/flux_ests/delivery/index.html'
site_page <- html_session(loads_index_url)
annual_href <- site_page %>% html_nodes('li') %>% html_nodes('a') %>% html_attr('href') %>%
  .[grep('[A-Z]*-Annual-2015.xlsx', basename(.))]
annual_url <- jump_to(site_page, annual_href)$url
flux_files <- localfile <- file.path(rawdir, basename(annual_url))
if(!file.exists(localfile))
  download.file(annual_url, destfile=localfile,  mode='wb')

# read and parse the file

# preliminary read
i <- 1
dat1 <- read_excel(flux_files[i], col_names=FALSE)
# chop off any empty columns
dat_cols_last <- 30 #max(which(sapply(dat1, function(col) length(which(!is.na(col))) > 0)))
dat_rows_last <- max(sapply(dat1, function(col) if(any(!is.na(col))) max(which(!is.na(col)), na.rm=TRUE) else NA), na.rm=TRUE)
dat1 <- dat1[seq_len(dat_rows_last), seq_len(dat_cols_last)]
# extract the headers and locate the start of the data
header_rows_top <- grep('Water Year1', dat1[[1]])
header_rows_bot <- header_rows_top + which(sapply(header_rows_top + 1:10, function(row) length(grep("Lower Confidence Interval", dat1[row, ])) > 0))
header_rows <- header_rows_top:header_rows_bot
headers <- dat1[header_rows,]

# identify the header row that identifies the nutrient in that column and 
# subsequent columns
nutrient_row <- unique(na.omit(sapply(headers, function(h) {
  nut_label <- grep('Metric Tons', h)
  if(length(nut_label) == 0) NA else nut_label
})))
if(length(nutrient_row) != 1) stop('expecting exactly 1 nutrient row')

# create concise column names
header <- sapply(seq_along(headers[1:30]), function(j) {
  hj <- headers[[j]]
  # the column type comes from the bottommost non-NA value in the column
  col_type <- hj[max(which(!is.na(hj)))] %>% trimws()
  # identify a short name for the column type
  colname <- switch(
    col_type,
    'Water Year1'='Water_Year',
    'Average Flow (m3/s)'='Flow_m3s',
    'LOADEST AMLE Predicted Load'='Flux_ty', # ty == metric tons per year
    'LOADEST AMLE Predicted Flux'='Flux_ty', # ty == metric tons per year
    'Composite Method Predicted Load'='FluxCmp_ty', # ty == metric tons per year
    'Lower Confidence Interval'='Flux_lo_ty',
    'Upper Confidence Interval'='Flux_hi_ty',
    {
      paste('Flow_m3s', switch(
        substr(col_type, 1, 3),
        '(A)'='A',
        '(B)'='B',
        '(C)'='C',
        '(D)'='D',
        'Tot'='E',
        stop('unrecognized coltype/site: ', col_type)
      ), sep='_site')
    })
  if(grepl('Flux_|FluxCmp_', colname)) {
    # find the long nutrient name in this or one of the preceding 3 columns
    nut_long <- c(na.omit(unname(unlist(headers[nutrient_row, (j-3):j]))))
    # identify a short name for the nutrient
    nut <- switch(
      nut_long,
      'NO2+NO3 (Metric Tons as N)'='NO23_N',
      'TKN  (Metric Tons as N)'='TKN_N',
      'NH3  (Metric Tons as N)'='NH3_N',
      'TP  (Metric Tons as P)'='TP_P',
      'OrthoP (Metric Tons as P)'='OrthoP_P',
      'SiO2 (Metric Tons as SiO2)'='SiO2',
      stop('unrecognized nutrient: ', nut_long))
    colname <- paste0(nut, '_', colname)
  }
  colname
})

# second pass at reading the file, this time with simpler headers
dat2 <- read_excel(flux_files[i], skip=header_rows_bot, col_names=FALSE) %>%
  select(seq_len(dat_cols_last)) %>%
  lapply(function(col) if(is.numeric(col)) col else as.numeric(gsub('na','',col))) %>% as_data_frame() %>%
  setNames(header) %>%
  mutate(
    Abbr = paste(strsplit(basename(flux_files[i]), '\\.|-')[[1]][1:2], collapse='-')
  ) %>%
  select(Abbr, everything())

est_obs_err <- function(dat, nut, alpha=0.05) {
  col_y <- paste0(nut, '_Flux_ty')
  col_ymin <- paste0(nut, '_Flux_lo_ty')
  col_ymax <- paste0(nut, '_Flux_hi_ty')
  
  qci <- qt(1 - alpha/2, df=Inf) # use df=Inf because we don't know the df for each load estimate, and it varies
  select_(dat, .dots=c('Water_Year', y=col_y, ymin=col_ymin, ymax=col_ymax)) %>% 
    filter(., complete.cases(.)) %>%
    mutate(
      # rloadest CIs are based on lognormal distributions, but I think we want
      # linear space obs error, so approximating the CIs that way. They're not
      # too bad; i.e., CI_dn is only a little smaller than CI_up
      CI_up = ymax - y,
      CI_dn = y - ymin,
      CI_95 = ymax - ymin,
      SD = CI_95 / (2*qci),
      Var = SD^2
    ) %>%
    summarize(
      solute = nut,
      SD = mean(SD),
      Var = mean(Var))
}
SDs <- bind_rows(lapply(c('NO23_N','TKN_N','NH3_N','TP_P','OrthoP_P','SiO2'), function(nut)
  est_obs_err(dat=dat2, nut)
)) %>%
  mutate(sqrtVar=sqrt(Var))
write.csv(SDs, file='obs_err_MARB.csv', row.names=FALSE)


# show that df=Inf isn't so bad
data.frame(df=exp(0:10), qci=sapply(exp(0:10), function(x) qt(1 - alpha/2, df=x))) %>%
  ggplot(aes(x=df, y=qci)) + geom_point() + 
  scale_x_log10(breaks=c(1,10,100,1000,10000)) + scale_y_log10(limits=c(1.5,13), breaks=c(1.5,2,3,5,10)) +
  theme_bw()


library(ggplot2)
SDs %>% 
  select(-Var) %>%
  gather(variable, value, SD, sqrtVar) %>%
  ggplot(aes(x=solute, y=value, color=variable)) + geom_point(alpha=0.5) +
  theme_bw() + theme(axis.text.x = element_text(angle=90))
ggsave(filename='obs_err_MARB.png', height=3, width=7)
