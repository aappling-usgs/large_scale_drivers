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
loads_index_url <- 'http://toxics.usgs.gov/hypoxia/mississippi/flux_ests/five_basins/index.html'
s <- html_session(loads_index_url)
site_xlsxes <- s %>% html_nodes('li') %>% html_nodes('a') %>% html_attr('href') %>%
  .[grep('\\.xlsx', basename(.))]

flux_files <- sapply(1:length(site_xlsxes), function(i) {
  annual_url <- jump_to(s, site_xlsxes[i])$url
  localfile <- file.path(rawdir, basename(annual_url))
  if(!file.exists(localfile))
    download.file(annual_url, destfile=localfile,  mode='wb')
  localfile
})

# read and reformat each data file
library(readxl)
flux_data_list <- lapply(seq_along(flux_files), function(i) {
  message(paste0('reading site ', i, ' from ', flux_files[i]))
  
  # preliminary read
  dat1 <- read_excel(flux_files[i], col_names=FALSE)
  # chop off any empty columns
  dat_cols_last <- max(which(sapply(dat1, function(col) length(which(!is.na(col))) > 0)))
  dat_rows_last <- max(sapply(dat1, function(col) if(any(!is.na(col))) max(which(!is.na(col)), na.rm=TRUE) else NA), na.rm=TRUE)
  dat1 <- dat1[seq_len(dat_rows_last), seq_len(dat_cols_last)]
  # extract the headers and locate the start of the data
  header_rows_top <- grep('Time Period \\(Water Years(1)*\\)', dat1[[1]])
  header_rows_bot <- grep('Lower Confidence Interval', dat1[[4]])
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
  header <- sapply(seq_along(headers), function(j) {
    hj <- headers[[j]]
    # the column type comes from the bottommost non-NA value in the column
    col_type <- hj[max(which(!is.na(hj)))] %>% trimws()
    # identify a short name for the column type
    colname <- switch(
      col_type,
      'Time Period (Water Years)'='Water_Year',
      'Time Period (Water Years1)'='Water_Year',
      'Average Flow (m3/s)'='Flow_m3s',
      'LOADEST AMLE Predicted Flux'='Flux_ty', # ty == metric tons per year
      'Lower Confidence Interval'='Flux_lo_ty',
      'Upper Confidence Interval'='Flux_hi_ty',
      stop('unrecognized column type: ', col_type))
    if(grepl('Flux_', colname)) {
      # find the long nutrient name in this or one of the preceding 2 columns
      nut_long <- c(na.omit(unname(unlist(headers[nutrient_row, (j-2):j]))))
      # identify a short name for the nutrient
      nut <- switch(
        nut_long,
        'NO2+NO3 (Metric Tons as N)'='NO23_N',
        'TKN  (Metric Tons as N)'='TKN_N',
        'NH3  (Metric Tons as N)'='NH3_N',
        'TN  (Metric Tons as N)'='TN_N',
        'TP  (Metric Tons as P)'='TP_P',
        'OrthoP (Metric Tons as P)'='OrthoP_P',
        'SiO2 (Metric Tons as SiO2)'='SiO2',
        stop('unrecognized nutrient: ', nut_long))
      colname <- paste0(nut, '_', colname)
    }
    colname
  })
  
  dat2 <- read_excel(flux_files[i], skip=header_rows_bot, col_names=FALSE) %>%
    select(seq_len(dat_cols_last)) %>%
    lapply(function(col) if(is.numeric(col)) col else as.numeric(gsub('na','',col))) %>% as_data_frame() %>%
    setNames(header) %>%
    filter(!is.na(Water_Year)) %>%
    mutate(
      Abbr = paste(strsplit(basename(flux_files[i]), '\\.|-')[[1]][1:2], collapse='-'),
      Flux_Src = headers %>% unlist %>% unname %>% grep('Fluxes$', ., value=TRUE) %>% {ifelse(length(.)==1,.,NA)},
      Flow_Src = headers %>% unlist %>% unname %>% grep('^Streamflow', ., value=TRUE) %>% {ifelse(length(.)==1,.,NA)},
      Flux_Site = ifelse(is.na(Flux_Src), NA, strsplit(Flux_Src, '\\(')[[1]][1] %>% trimws()),
      Flux_SiteID = ifelse(is.na(Flux_Src), NA, strsplit(Flux_Src, ' ')[[1]] %>% grep(')$', ., value=TRUE) %>% { substring(., 1, nchar(.)-1) }),
      Flow_Site = ifelse(is.na(Flow_Src), NA, strsplit(Flow_Src, 'Streamflow data from |\\(')[[1]][2] %>% trimws()),
      Flow_SiteID = ifelse(is.na(Flow_Src), NA, strsplit(Flow_Src, ' ')[[1]] %>% grep(')$', ., value=TRUE) %>% { substring(., 1, nchar(.)-1) })
    ) %>%
    select(-Flux_Src, -Flow_Src) %>%
    select(Abbr, everything())
  
  dat2
})
flux_data <- bind_rows(flux_data_list)

source('est_err_obs_helpers.R')
obs_err_4basins <- bind_rows(lapply(unique(flux_data$Abbr), function(site) {
  mutate(get_obs_errs(filter(flux_data, Abbr==site)), Site=site)
})) %>% select(Site, everything())
write.csv(obs_err_4basins, file='obs_err_4basins.csv', row.names=FALSE)

library(ggplot2)
obs_err_4basins %>% 
  mutate(
    Cor_Ratio=factor(ifelse(is.na(Ratio_Cor), '-', Ratio_Cor)),
    Cor_TN=factor(ifelse(is.na(TN_Cor), '-', TN_Cor))) %>%
  ggplot(aes(x=Site, y=Var, color=Cor_Ratio, shape=Cor_TN)) + geom_point(alpha=0.5) + facet_wrap(~ Solute, scales='free_y') +
  theme_bw() + theme(axis.text.x = element_text(angle=90))
ggsave(filename='obs_err_4basins.png', height=8, width=7)
