library(Hmisc);library(plyr); library(dplyr); library(tidyr); 
myGDwheat = read.csv(file = 'RawDataShantelCompiled/myGD.csv') %>% rename(taxa = X)
myGMwheat = read.csv(file = 'RawDataShantelCompiled/myGM.csv')
WheatPHSraw = read.csv(file = 'RawDataShantelCompiled/PHSAll_Comb_20181222.csv',
                       head = TRUE,na.string=c(""," ","NA","na","NaN"), stringsAsFactors=FALSE) %>%
  mutate(taxa = gsub("cuGSMSU", "MSU",gsub("cuGSOH", "OH",paste('cuGS',GID, sep ='')))) %>%
  filter(!is.na(RawMean)) %>% mutate(Year = as.character(Year),
                                     RawMean = as.numeric(RawMean))

WheatPHSraw %>% group_by(taxa) %>% summarise(n = n()) %>% ungroup() %>% select(n) %>% table()
WheatPHSraw %>% select(X) %>% table()
NoGenotypes = WheatPHSraw %>% select(taxa) %>% unique() %>% filter(taxa %nin% myGDwheat$taxa)
WheatPHSraw %>% filter(taxa %nin% NoGenotypes$taxa)%>% group_by(taxa) %>% summarise(n = n()) %>%
  ungroup() %>% select(n) %>% table() %>% 

# from Shantel:
# QUALITY CHECK: Red kernels were after-ripened for 7 days while white kernels were 
# after-ripened for 5 days. Omit lines that have conflicting AR time points. 
# However, keep the NAs because that is a lack of harvest and misting date recording, 
# not necessarily the wrong AR length. NOTE: the 8 days AR red kernel color days was
# one year with a reason. That years environment induced more dormancy than usual, 
# so the sames were after-ripened for one more day longer.
Red = WheatPHSraw %>% filter(KC=='R'|KC == 'r') %>% filter(AR !=5 | is.na(AR)) %>%
  filter(!is.na(Harvest)) %>% filter(taxa %nin% NoGenotypes$taxa)
Red %>% group_by(taxa) %>% summarise(n = n()) %>% select(n) %>% table()
White = WheatPHSraw %>% filter(KC =='W' |KC == 'w') %>% filter(AR !=7 | is.na(AR)) %>%
  filter(!is.na(Harvest)) %>% filter(taxa %nin% NoGenotypes$taxa)
White %>% group_by(taxa) %>% summarise(n = n()) %>% select(n) %>% table() 

taxaGT12obswhite = White %>% group_by(taxa) %>% summarise(n = n()) %>% arrange(-n) %>% filter(n > 12)
length(unique(taxaGT12obswhite$taxa))
View(White %>%filter(taxa %in%taxaGT12obswhite$taxa)%>% group_by(Env) %>% summarise(n = n()) )
View(White %>%filter(taxa %in%taxaGT12obswhite$taxa)%>% group_by(Env) %>% mutate(n = n()) %>% 
       select(taxa,Env, RawMean,n))



#since the white data set is much larger than the other data set lets start with it.
#Game Farm weather Data
GarmFarmWeather = read.csv(file = 'TemperatureAnalysis/Weather2007_2020.txt') %>% 
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
         Year = format(Date,'%Y'),
         JulianDate = format(Date,'%j'))

 filter(taxa %in% c('cuGS264','cuGS35', 'cuGS262'))
White %>% group_by(Env) %>% select(taxa, RawMean, Env, Harvest, Mist, Score) %>%
  mutate(envmean = mean(RawMean)) %>%
  ggplot(aes(x = envmean, y = RawMean, group = taxa)) +geom_point()+geom_smooth(method = 'lm', se = FALSE)

White %>% group_by(Env) %>% select(taxa, RawMean, Env, Harvest, Mist, Score) %>%
  mutate(envmean = mean(RawMean)) %>% ungroup() %>% group_by(taxa) %>%
  group_modify(~broom::tidy(lm(RawMean ~ envmean, data = .x)))

White%>%filter(taxa %in% c('cuGS264','cuGS35', 'cuGS262')) %>% group_by(Env)  %>% 
  select(taxa, RawMean, Env, Harvest, Mist, Score) %>%
  mutate(envmean = mean(RawMean)) %>% ungroup() %>% 
  ggplot(aes(x = envmean, y = RawMean, group = taxa)) +geom_point() +geom_line()+
    geom_smooth(method = 'lm')
    
White %>% filter(taxa %in% taxaGT12obswhite$taxa) %>% group_by(Env)  %>% 
  select(taxa, RawMean, Env, Harvest, Mist, Score) %>%
  mutate(envmean = mean(RawMean)) %>% ungroup() %>% 
  ggplot(aes(x = envmean, y = RawMean, group = taxa)) +geom_point() +geom_line()+
  geom_smooth(method = 'lm', se = FALSE)

# IT seems like there may be some that are insensitive to the environment
# Lets use the groups of PHS means with those taxaGT12obs to find the environmental means. 
# With 12 obs. There are 28 unique environments. 
WhiteFilt = White%>%filter(taxa %in% taxaGT12obswhite$taxa) %>% group_by(taxa, Env)  %>% 
  select(taxa, RawMean, Env, Harvest, Mist, Score,Year) %>% 
  summarise(RawMean = mean(RawMean), 
            Harvest = mean(Harvest), 
            Mist = mean(Mist)) %>%
  mutate(envmean = mean(RawMean)) %>% ungroup()

View(WhiteFilt %>% group_by(taxa) %>% dplyr::summarize(n = n()) %>%arrange(-n))

test = WhiteFilt %>% filter(taxa == 'cuGS901') %>% select(Env)
View(WhiteFilt %>% filter(Env %in% test$Env) %>% group_by(taxa) %>% dplyr::summarize(n = n()) %>%arrange(-n))

