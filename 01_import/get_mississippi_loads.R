# Download and format load data from http://toxics.usgs.gov/hypoxia/mississippi/flux_ests/sub_basins/index.html

library(dplyr)

# create directories to store the raw and reformatted data files
rawdir <- '01_import/raw'
if(!dir.exists(rawdir)) dir.create(rawdir)
mungedir <- '01_import/munged'
if(!dir.exists(mungedir)) dir.create(mungedir)

# use rvest to find & download the annual data files from the website
library(rvest)
loads_index_url <- 'http://toxics.usgs.gov/hypoxia/mississippi/flux_ests/sub_basins/index.html'
s <- html_session(loads_index_url)
site_pages <- s %>% html_nodes('li') %>% html_nodes('a') %>% html_attr('href') %>%
  .[grep('[A-Z]{2,4}-[A-Z]{2,4}', basename(.))]

flux_files <- sapply(1:length(site_pages), function(i) {
  site_page <- jump_to(s, site_pages[i])
  site_hrefs <- site_page %>% html_nodes('a')
  annual_href <- site_hrefs[grep('Annual Nutrient Flux and Concurrent Streamflow', site_hrefs)] %>% html_attr('href')
  annual_url <- jump_to(site_page, annual_href)$url
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
  header_rows_top <- grep('Water Year1', dat1[[1]])
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
      'Water Year1'='Water_Year',
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

# spot-check plots to make sure these are the data we expect
library(ggplot2)
ggplot(flux_data, aes(x=Water_Year, y=Flow_m3s, color=Abbr)) + geom_line()
ggplot(flux_data, aes(x=Water_Year, y=NO23_N_Flux_ty, color=Abbr)) + geom_line()
ggplot(flux_data, aes(x=Flow_m3s, y=NO23_N_Flux_ty, color=Abbr)) + geom_point()

# reshape data
library(tidyr)
# extract, distill, and elaborate on the site info
library(dataRetrieval)
site_info <- flux_data %>%
  select(Abbr, Flux_Site, Flow_Site, Flux_SiteID, Flow_SiteID) %>%
  distinct() %>%
  gather(Type, Value, Flux_Site, Flow_Site, Flux_SiteID, Flow_SiteID) %>%
  mutate(Data_Type = sapply(strsplit(Type, '_'), `[`, 1),
         Label_Type = sapply(strsplit(Type, '_'), `[`, 2)) %>%
  select(-Type) %>%
  spread(Label_Type, Value) %>%
  filter(!is.na(Site) & !is.na(SiteID)) %>%
  {left_join(., readNWISsite(.$SiteID), by=c(SiteID='site_no'))}

# add N:P ratios. Jay advocates NO3+NO2 : SRP as long as the SRPs aren't too
# low. what's too low? 100 tonnes/year looks relatively low among these sites,
# though maybe they all have plenty of P (need to ask Jay, others). ranges:
sapply(unique(flux_data$Abbr), function(abv) {
  range(filter(flux_data, Abbr==abv)$OrthoP_P_Flux_ty, na.rm=TRUE)
})
# look up atomic masses to get molar N:P ratio. the tonne:gram conversion will 
# cancel out (so we don't need to do it) because it's the same for both N and P
molN_gN <- 1/14.0067
molP_gP <- 1/30.973762
# compute molar N:P ratio from NO2+NO3 and OrthoP whenever OrthoP doesn't get
# too low for that site
flux_data_stoich <- flux_data %>%
  group_by(Abbr) %>%
  mutate(OrthoP_P_too_low = min(OrthoP_P_Flux_ty, na.rm=TRUE) < 100) %>% # mostly arbitrary for now
  ungroup() %>%
  mutate(NO23OrthoP_NP_Flux_ty = 
           ifelse(OrthoP_P_too_low, NA, 
                  NO23_N_Flux_ty * molN_gN / (OrthoP_P_Flux_ty * molP_gP)))

# reshape core data cols into one col each for Abbr, year, flow, value, lowerCI,
# upperCI
flux_data_tidy <- flux_data_stoich %>%
  select(Abbr, Water_Year, Flow_m3s, ends_with('_ty')) %>%
  gather(Variable, Value, ends_with('_ty'), na.rm=TRUE) %>%
  mutate(
    Nutrient=substring(Variable, 1, unlist(regexec('_Flux_', Variable))-1),
    Bound=substring(Variable, unlist(regexec('_Flux_', Variable))+1)) %>%
  select(-Variable) %>%
  spread(Bound, Value) %>%
  select(Abbr, Water_Year, Flow_m3s, Nutrient, Flux_ty, Flux_lo_ty, Flux_hi_ty)

# split and save by nutrient (or nutrient ratio)
lapply(unique(flux_data_tidy$Nutrient), function(nutrient) {
  flux_data_1nut <- filter(flux_data_tidy, Nutrient == nutrient)
  write.csv(flux_data_1nut, file=file.path(mungedir, paste0(nutrient, '.csv')))
  saveRDS(flux_data_1nut, file=file.path(mungedir, paste0(nutrient, '.Rds')))
})
# save the site data
write.csv(site_info, file=file.path(mungedir, 'Site_Info.csv'))
saveRDS(site_info, file=file.path(mungedir, 'Site_Info.Rds'))
