library(plyr); library(dplyr); library(tidyr)
myGDwheat = read.csv(file = 'RawDataShantelCompiled/myGD.csv') %>% rename(taxa = X)
myGMwheat = read.csv(file = 'RawDataShantelCompiled/myGM.csv')
WheatPHSraw = read.csv(file = 'RawDataShantelCompiled/PHSAll_Comb_20181222.csv',
                       head = TRUE,na.string=c(""," ","NA","na","NaN"), stringsAsFactors=FALSE) %>%
  mutate(taxa = gsub("cuGSMSU", "MSU",gsub("cuGSOH", "OH",paste('cuGS',GID, sep ='')))) %>%
  filter(!is.na(RawMean)) %>% mutate(Year = as.character(Year),
                                     RawMean = as.numeric(RawMean))

WheatPHSraw %>% group_by(taxa) %>% summarise(n = n()) %>% ungroup() %>% select(n) %>% table()
# from Shantel:
# QUALITY CHECK: Red kernels were after-ripened for 7 days while white kernels were 
# after-ripened for 5 days. Omit lines that have conflicting AR time points. 
# However, keep the NAs because that is a lack of harvest and misting date recording, 
# not necessarily the wrong AR length. NOTE: the 8 days AR red kernel color days was
# one year with a reason. That years environment induced more dormancy than usual, 
# so the sames were after-ripened for one more day longer.
Red = WheatPHSraw %>% filter(KC=='R'|KC == 'r') %>% filter(AR !=5 | is.na(AR)) %>%
  filter(!is.na(Harvest))
Red %>% group_by(taxa) %>% summarise(n = n()) %>% select(n) %>% table()
White = WheatPHSraw %>% filter(KC =='W' |KC == 'w') %>% filter(AR !=7 | is.na(AR)) %>%
  filter(!is.na(Harvest))
White %>% group_by(taxa) %>% summarise(n = n()) %>% select(n) %>% table()

#since the white dataset is much larger than the other data set lets start with it.
#Game Farm weather Data
GarmFarmWeather = read.csv(file = 'TemperatureAnalysis/Weather2007_2020.txt') %>% 
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
         Year = format(Date,'%Y'),
         JulianDate = format(Date,'%j'))


