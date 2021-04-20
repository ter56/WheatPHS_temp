
library(GAPIT3)
library(rrBLUP)
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
  mutate(Date1 = as.Date(Date, format = "%Y-%m-%d"),
         Year = format(Date1,'%Y'),
         JulianDate = format(Date1,'%j'),
         Date = gsub(x = as.character(Date1),pattern = '-',replacement = ''))

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
            Mist = mean(Mist)) %>% ungroup() %>% group_by(Env)%>%
  mutate(envmean = mean(RawMean)) %>% ungroup()

View(WhiteFilt %>% group_by(taxa) %>% dplyr::summarize(n = n()) %>%arrange(-n))

test = WhiteFilt %>% filter(taxa == 'cuGS901') %>% select(Env)

LineForenvMean = WhiteFilt %>% filter(Env %in% test$Env) %>% group_by(taxa) %>% dplyr::summarize(n = n()) %>%arrange(-n) %>%
  filter(n > 8)
  
WhiteFilt %>% filter(taxa %in% LineForenvMean$taxa) %>% group_by(Env) %>% mutate(envmean = mean(RawMean)) %>%
  ggplot(aes(x = envmean, y = RawMean, group = taxa)) +geom_point() +geom_line()+
  geom_smooth(method = 'lm', se = FALSE)
# So lets use these lines, 34 lines I think to figure out the correlation with environmental variables. 

WhiteFilt %>% filter(taxa %in% LineForenvMean$taxa) %>% group_by(Env) %>% 
  mutate(envmean = mean(RawMean)) %>% ungroup() %>% group_by(taxa) %>%
  group_modify(~broom::tidy(lm(RawMean ~ envmean, data = .x)))


#get an environmental file with col names of env_code per year
  
# From Jianming Yus group:
#################
Exhaustive_search <- function(env_mean_trait, 
                              env_paras, 
                              searching_daps, 
                              exp_trait_dir, 
                              FTdaps, 
                              trait, 
                              p, 
                              dap_x, 
                              dap_y, 
                              LOO, 
                              Paras, 
                              pop_cor_file) {
  # env_paras <- PTT_PTR; FTdaps <- exp_traits$FTdap; p <- 1; dap_x <- searching_daps; dap_y <- searching_daps;
  
  exs_png_file <- paste(exp_trait_dir, 'MaxR_',trait, '_', nrow(env_mean_trait), 'Envs_', LOO, 'LOO.png', sep = ''); 
  
  nParas <- length(Paras);
  if (!file.exists(pop_cor_file)) {
    dap_win <- searching_daps * searching_daps  / 2;
    
    pop_cors_matrix <- matrix(ncol = 4 + (2 * nParas), nrow = dap_win * 1);
    colnames(pop_cors_matrix) <- c("pop_code", 'Day_x', 'Day_y', 'window', paste('R_', Paras, sep = ''), paste('nR_', Paras, sep = ''));
    n <- 0;
    for (d1 in 1:(dap_y - 6)) {
      for (d2 in (d1 + 6):dap_y) {
        days <- c(d1:d2); 
        env_facts_matrix <- matrix(nrow = nrow(env_mean_trait), ncol = nParas);
        for (e_i in 1:nrow(env_mean_trait)) {
          e <- env_mean_trait$env_code[e_i];
          env_para <- subset(env_paras, env_paras$env_code == e);
          env_mean <- colMeans(env_para[days, (1:nParas) + 2]); ### DL, GDD, DTR, PTT, PTR, PTD, PTD2, PTS
          env_facts_matrix[e_i,] <- env_mean;
        }
        n <- n + 1;
        ### leave one environment out and get the median correlation
        Ymean_envPara <- cbind(env_facts_matrix, env_mean_trait$meanY);
        rs <- c();
        if (LOO == 0) {
          for (k in 1:nParas) {
            rs[k] <- round(cor(Ymean_envPara[,nParas + 1], Ymean_envPara[,k]), digits = 4)
            #           rs[k] <- round(-log(cor.test(Ymean_envPara[,8], Ymean_envPara[,k])$p.value, 10), digits = 4)
          }
        } else {
          loo_rs_matrix <- matrix(nrow = nrow(Ymean_envPara)+ 0, ncol = nParas);
          for (k in 1:nParas) { ## 8 environment parameters
            for (e_x in c(1:nrow(Ymean_envPara))) {
              t_matrix <- Ymean_envPara[-e_x,];
              loo_rs_matrix[e_x, k] <- round(cor(t_matrix[,nParas + 1], t_matrix[,k]), digits = 4)
            }
          }
          rs <- apply(loo_rs_matrix, 2, median);
        }
        pop_cors_matrix[n, ] <- c(p, d1, d2, d2 - d1, rs, 0 - rs);
      }
    }
    pop_cors_matrix <- pop_cors_matrix[1:n,]
    # write.table(pop_cors_matrix, file = pop_cor_file, sep = "\t", row.names = F, quote = F);
    
  }
  return(pop_cors_matrix)
  # pop_cors <- read.table(pop_cor_file, header = T, sep = "\t");
  # pop_cor <- subset(pop_cors, pop_cors$pop_code == p);
  # # dev.off();
  # # pdf(exs_pdf_file,width= nParas,height= 2,pointsize=6)
  # png(exs_png_file, width = nParas * 1.5, height = 2 * 1.5, pointsize = 12, unit = "in", res = 600)
  # layout(matrix(c(1:(2*nParas)), 2, nParas, byrow = T))
  # 
  # for (k in 1:(2*nParas)) {
  #   pop_cor_0 <- subset(pop_cors, pop_cors$pop_code == p); 
  #   pop_cor <- pop_cor_0[,c(1:4, k + 4)];
  #   colnames(pop_cor)[5] <- 'R';
  #   pop_cor <- pop_cor[order(pop_cor$R),];
  #   
  #   xs <- pop_cor$Day_x;  ys <-  pop_cor$Day_y;  mid_R <- median(pop_cor$R);
  #   #   pop_cor_L <- subset(pop_cor, pop_cor$R <= mid_R); cor_range_L <- range(pop_cor_L$R);
  #   #   cell_col_L <- floor((pop_cor_L$R - min(pop_cor_L$R)) / diff(cor_range_L) * col_wdw / 2) + 1;
  #   #  
  #   #   pop_cor_G <- subset(pop_cor, pop_cor$R > mid_R); cor_range_G <- range(pop_cor_G$R);
  #   #   cell_col_G <- ceiling((pop_cor_G$R - min(pop_cor_G$R)) / diff(cor_range_G) * col_wdw / 2) + 12;
  #   #   cell_col <- c(cell_col_L, cell_col_G);
  #   
  #   cell_col <- floor(pop_cor$R * 12) + 13; ### the same color scale
  #   
  #   pop_cor$cell_col <- cell_col; 
  #   
  #   #   pop_cor_m <- subset(pop_cor, pop_cor$window > 4 & pop_cor$window < 50 ); ##& pop_cor$Day_x > (window_ref$dap_s - 10) & pop_cor$Day_y < (window_ref$dap_e + 10));
  #   #   max_R <- pop_cor_m[which.max(pop_cor_m$R)[1], ];
  #   pop_cor_6 <- subset(pop_cor, pop_cor$window > 6); max_R <- pop_cor_6[which.max(pop_cor_6$R)[1], ];
  #   
  #   par(mar = c(0.5, 1.0, 1, 0.5) , mgp = c(0.05, 0.1, 0), tck = -0.01, bty = "n");
  #   plot(-50, -50, xlim = c(0, dap_x), ylim = c(0, dap_y), col = "white",  xlab = '', xaxt = "n", yaxt = "n", ylab = 'Days after planting', bty = "n", cex.lab = .4);
  #   arrows(-1, 10, -1, dap_y - 10, length = 0.05, angle = 15, lwd = .5,  col = "grey59");
  #   mtext(c(1, 50, 100, dap_y), side = 2, at = c(1,50, 100, dap_y), line = -1, cex = .6)
  #   
  #   rect(xs - 0.5, ys - 0.5, xs + 0.5, ys + 0.5, col = col_palette[pop_cor$cell_col], border = "NA")
  #   rect(max(pop_cor$Day_x) - 0.5, max(pop_cor$Day_y) - 0.5, max(pop_cor$Day_x) + 0.5, max(pop_cor$Day_y) + 0.5, border = "NA", col = "white", lwd = 0.001)
  #   #   legend("bottom", Div_Fnd_lab, bty = "n", cex = .6);
  #   
  #   arrows(10, dap_y + 4, dap_x - 10, dap_y + 4, angle = 15, length = 0.05, lwd = .5, col = "grey59")
  #   mtext("Days after planting", side = 3, at = dap_x / 2, line = -0.1, cex = .4)
  #   mtext(c(1, 50, 100, dap_y), side = 3, at = c(1, 50, 100, dap_y), line = -1.1, cex = .6)
  #   arrows(max_R$Day_x + 4,  max_R$Day_y - 4,  max_R$Day_x,  max_R$Day_y, length = 0.05, angle = 15, lwd = .5, col = "grey59")
  #   
  #   box_ys <- seq(1, 50, by = 2); box_xs <- rep(dap_x - 15, 25); 
  #   rect(box_xs - .5 * 2, box_ys - 0.5 * 2, box_xs + 0.5 * 2, box_ys + 0.5 * 2, border = "NA", col = col_palette)
  #   text(dap_x - 10 - 5, 52, 'r', cex = .5);
  #   r_lab_top <- 1; r_lab_mid <- 0; r_lab_bottom <- -1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", max_R$R), sep = '');
  #   if (k > nParas) { r_lab_top <- -1; r_lab_bottom <- 1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", 0 - max_R$R), sep = ''); mtext(side = 1, Paras[k - nParas ], line= -0.5,  cex = .5, bty = "n")}
  #   legend(max_R$Day_x - 4 , max_R$Day_y - 4 , c(paste( max_R$Day_x, ' to ', max_R$Day_y, ' DAP', sep = ''), max_r_lab),  cex = .6, bty = "n");
  #   text(dap_x - 10 + 3, 50, r_lab_top, cex = .5)
  #   text(dap_x - 10 + 3, 27, r_lab_mid, cex = .5);
  #   text(dap_x - 10 + 3, 1,  r_lab_bottom, cex = .5)
  #   if (length(FTdaps) > 1) {
  #     boxplot(FTdaps,   at = 145,  add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
  #     boxplot(FTdaps,   at = 1, horizontal = T, add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
  #     text(mean(FTdaps), 5, 'Days to anthesis', cex = .5)
  #     text(mean(FTdaps), 10, paste('Trait: ', trait, sep = ''), cex = .6)
  #   }
  # }
  # dev.off()
  # 
  
}

#  white PHS Earliest harvest date was 180 juilian days so lets use that as the searching_daps####
EnviromentalMeans = WhiteFilt %>% filter(taxa %in% LineForenvMean$taxa) %>% group_by(Env) %>% 
  mutate(envmean = mean(RawMean)) %>% ungroup() %>% select(Env, envmean) %>% unique() %>%
  rename(env_code = Env, meanY = envmean)
env_parasWhite = GarmFarmWeather %>% 
  join(.,EnviromentalMeans %>% mutate(Year = paste0('20',substr(env_code,4,5))),type = 'right',by = 'Year') %>%
  select(env_code, Date, MaxTemperature, MinTemperature, AvgTemperature, HDD, CDD, GDD)

WhiteWheatCors = Exhaustive_search(env_mean_trait = EnviromentalMeans,
                  env_paras = env_parasWhite,
                  searching_daps = 180,
                  exp_trait_dir = '.',
                  FTdaps = 'test',
                  trait = 'PHS',
                  p = 1,
                  dap_x = 180,
                  dap_y = 180,
                  LOO = 0,
                  Paras = c('MaxTemperature', 'MinTemperature',
                            'AvgTemperature', 'HDD', 'CDD', 'GDD'),
                  pop_cor_file = 'test')
  
search.results <- WhiteWheatCors %>% as.data.frame() %>%
  pivot_longer(c(-Day_x, -Day_y, -window, -pop_code), names_to = 'Parameter', values_to = 'Corr')%>%
  arrange(-Corr)%>%
  group_by(Parameter) %>%
  top_n(5, Corr)    

View(search.results)
# white HD filtered lets just take a look at heading date and its correlations here as well. #####
EnviromentalMeansHD = White %>% group_by(Env,taxa) %>% summarise(HD = mean(HD)) %>%
 ungroup()  %>% filter(taxa %in% LineForenvMean$taxa) %>% group_by(Env) %>% 
  mutate(envmean = mean(HD, na.rm = T)) %>% ungroup() %>% select(Env, envmean) %>% unique() %>%
  rename(env_code = Env, meanY = envmean) %>% filter(!is.nan(meanY))
env_parasWhiteHD = GarmFarmWeather %>% 
  join(.,EnviromentalMeansHD %>% mutate(Year = paste0('20',substr(env_code,4,5))),type = 'right',by = 'Year') %>%
  select(env_code, Date, MaxTemperature, MinTemperature, AvgTemperature, HDD, CDD, GDD)

WhiteWheatCorsHD = Exhaustive_search(env_mean_trait = EnviromentalMeansHD,
                                   env_paras = env_parasWhiteHD,
                                   searching_daps = 130,
                                   exp_trait_dir = '.',
                                   FTdaps = 'test',
                                   trait = 'PHS',
                                   p = 1,
                                   dap_x = 130,
                                   dap_y = 130,
                                   LOO = 0,
                                   Paras = c('MaxTemperature', 'MinTemperature',
                                             'AvgTemperature', 'HDD', 'CDD', 'GDD'),
                                   pop_cor_file = 'test')
HDfiltresults <- WhiteWheatCorsHD %>% as.data.frame() %>%
  pivot_longer(c(-Day_x, -Day_y, -window, -pop_code), names_to = 'Parameter', values_to = 'Corr')%>%
  arrange(-Corr)%>%
  group_by(Parameter) %>%
  top_n(5, Corr)    

View(HDfiltresults)



# White PHS with all data #######

EnviromentalMeans2 = White %>% group_by(Env,taxa) %>% summarise(RawMean = mean(RawMean)) %>%
  ungroup() %>% group_by(Env) %>%
  mutate(envmean = mean(RawMean)) %>% ungroup() %>% select(Env, envmean) %>% unique() %>%
  rename(env_code = Env, meanY = envmean)
env_parasWhite2 = GarmFarmWeather %>% 
  join(.,EnviromentalMeans2 %>% mutate(Year = paste0('20',substr(env_code,4,5))),type = 'right',by = 'Year') %>%
  select(env_code, Date, MaxTemperature, MinTemperature, AvgTemperature, HDD, CDD, GDD)

WhiteWheatCorsall = Exhaustive_search(env_mean_trait = EnviromentalMeans2,
                                   env_paras = env_parasWhite2,
                                   searching_daps = 180,
                                   exp_trait_dir = '.',
                                   FTdaps = 'test',
                                   trait = 'PHS',
                                   p = 1,
                                   dap_x = 180,
                                   dap_y = 180,
                                   LOO = 0,
                                   Paras = c('MaxTemperature', 'MinTemperature',
                                             'AvgTemperature', 'HDD', 'CDD', 'GDD'),
                                   pop_cor_file = 'test')

search.results2 <- WhiteWheatCorsall %>% as.data.frame() %>%
  pivot_longer(c(-Day_x, -Day_y, -window, -pop_code), names_to = 'Parameter', values_to = 'Corr')%>%
  arrange(-Corr)%>%
  group_by(Parameter) %>%
  top_n(5, Corr)    

View(search.results2)

# white HD correlations with all data! ######
EnviromentalMeansHD2 = White %>% group_by(Env,taxa) %>% summarise(HD = mean(HD)) %>%
  ungroup() %>% group_by(Env) %>% 
  mutate(envmean = mean(HD, na.rm = T)) %>% ungroup() %>% select(Env, envmean) %>% unique() %>%
  rename(env_code = Env, meanY = envmean) %>% filter(!is.nan(meanY))
env_parasWhiteHD2 = GarmFarmWeather %>% 
  join(.,EnviromentalMeansHD2 %>% mutate(Year = paste0('20',substr(env_code,4,5))),type = 'right',by = 'Year') %>%
  select(env_code, Date, MaxTemperature, MinTemperature, AvgTemperature, HDD, CDD, GDD)
WhiteWheatCorsallHD = Exhaustive_search(env_mean_trait = EnviromentalMeansHD2,
                                      env_paras = env_parasWhiteHD2,
                                      searching_daps = 130,
                                      exp_trait_dir = '.',
                                      FTdaps = 'test',
                                      trait = 'PHS',
                                      p = 1,
                                      dap_x = 130,
                                      dap_y = 130,
                                      LOO = 0,
                                      Paras = c('MaxTemperature', 'MinTemperature',
                                                'AvgTemperature', 'HDD', 'CDD', 'GDD'),
                                      pop_cor_file = 'test')

HDresults <- WhiteWheatCorsallHD %>% as.data.frame() %>%
  pivot_longer(c(-Day_x, -Day_y, -window, -pop_code), names_to = 'Parameter', values_to = 'Corr')%>%
  arrange(-Corr)%>%
  group_by(Parameter) %>%
  top_n(5, Corr)    
View(HDresults)

# Red+white HD correlation with all data since HD is measured the same way.  #######
head(WheatPHSraw)
AllHDEnvMeans = WheatPHSraw %>% group_by(Env, taxa) %>% summarize(HD = mean(HD)) %>% select(taxa, Env, HD) %>%
  ungroup() %>% group_by(Env) %>% summarize(envmean = mean(HD, na.rm = T)) %>% 
    rename(env_code = Env, meanY = envmean) %>% filter(!is.nan(meanY))

WheatPHSraw %>% group_by(Env, taxa) %>% 
  summarize(HD = mean(HD), PHS = mean(RawMean),
            Yield = mean(as.numeric(PlotYield), na.rm =T),
            TW = mean(as.numeric(TW), na.rm = T)) %>% ungroup() %>% pivot_longer(cols = c(HD,PHS,Yield,TW)) %>%
  ggplot(aes(x = value, fill = Env)) +geom_histogram() +facet_wrap(vars(name), scales = 'free')

env_parasWhite = GarmFarmWeather %>% 
  join(.,AllHDEnvMeans %>% mutate(Year = paste0('20',substr(env_code,4,5))),type = 'right',by = 'Year') %>%
  select(env_code, Date, MaxTemperature, MinTemperature, AvgTemperature, HDD, CDD, GDD)

AllHD_corrs = Exhaustive_search(env_mean_trait = AllHDEnvMeans,
                                   env_paras = env_parasWhite,
                                   searching_daps = 120,
                                   exp_trait_dir = '.',
                                   FTdaps = 'test',
                                   trait = 'PHS',
                                   p = 1,
                                   dap_x = 120,
                                   dap_y = 120,
                                   LOO = 0,
                                   Paras = c('MaxTemperature', 'MinTemperature',
                                             'AvgTemperature', 'HDD', 'CDD', 'GDD'),
                                   pop_cor_file = 'test')


allHDresults <- AllHD_corrs %>% as.data.frame() %>%
  pivot_longer(c(-Day_x, -Day_y, -window, -pop_code), names_to = 'Parameter', values_to = 'Corr')%>%
  arrange(-Corr)%>%
  group_by(Parameter) %>%
  top_n(5, Corr)    
View(allHDresults)






# Well for the white wheats day 55-61 window seemed to give the highest correlation with the ####
# Avg_temperature, This may be entirely spurious, or if the avg temp is high in this time 
# winter kill is more likely, which will delay the heading and cause hotter temps during grain fill
# in the summer. On second thought lets use this: avg GDD from 150 to 155 
# so for each env lets get this value. 
EnvPHSweather = GarmFarmWeather %>% 
  join(.,White %>% select(Env)%>%unique() %>% mutate(Year = paste0('20',substr(Env,4,5))),type = 'right',by = 'Year') %>%
  select(Env, Date, JulianDate, GDD, MaxTemperature, MinTemperature, AvgTemperature) %>% 
  group_by(Env) %>% filter(JulianDate %in% c(150:155)) %>% 
  summarize(MeanGDD = mean(GDD),
            MeanMaxTemp = mean(MaxTemperature),
            MeanMinTemp = mean(MinTemperature),
            MeanAvgTemp = mean(AvgTemperature),
            CumGDD = sum(GDD)) %>% mutate(MeanGDDc = MeanGDD - mean(MeanGDD))
EnvPHSweather %>% summary()
White %>% group_by(Env,taxa) %>% summarize(PHS = mean(RawMean, na.rm = T)) %>%
  ungroup() %>% group_by(taxa) %>% summarize(n = n()) %>% select(n)%>%table()

WithSuffceintObs = White %>% group_by(Env,taxa) %>% summarize(PHS = mean(RawMean, na.rm = T)) %>%
  ungroup() %>% group_by(taxa) %>% summarize(n = n()) %>% filter(n >6) 

White %>% filter(taxa %in% WithSuffceintObs$taxa) %>% mutate(PHS = RawMean) %>%
  join(., EnvPHSweather) %>%
  ggplot(aes(x = MeanGDDc, y= PHS, group = taxa)) +geom_point()+
  geom_smooth(method = 'lm',se = F)

White %>% filter(taxa %in% WithSuffceintObs$taxa) %>% mutate(PHS = RawMean)%>%
  join(., EnvPHSweather) %>% group_by(taxa) %>%
  group_modify(~broom::tidy(lm(PHS ~ MeanGDDc, data = .x))) %>% 
  ggplot(aes(x = estimate)) +geom_histogram()+ facet_wrap(vars(term), scales = 'free')

WhitePHSbyEnv = White %>% filter(taxa %in% WithSuffceintObs$taxa) %>% mutate(PHS = RawMean)%>%
  join(., EnvPHSweather) %>%  group_by(taxa) %>%
  group_modify(~broom::tidy(lm(PHS ~ MeanGDDc, data = .x)))

# wheatRealtions = A.mat(X =myGDwheat[,-1]-1,impute.method = 'EM' )
# PCAvalues = eigen(wheatRealtions)
plot(c(1:10),PCAvalues$values[1:10]/sum(PCAvalues$values))
# 2 PCs look good to me
myGDwheat[1:5,1:5]
myGMwheatGAPIT = myGMwheat %>% dplyr::select(rs, chrom, pos) %>%
  rename(SNP = rs, Chromosome = chrom, Position = pos)

Intercept.gwas.mlm = WhitePHSbyEnv %>% filter(term == '(Intercept)') %>%
  dplyr::select(taxa, estimate) %>% data.frame() %>%
  GAPIT(Y =.,GD = myGDwheat, GM = myGMwheatGAPIT, model = 'MLM', PCA.total = 2,
        Geno.View.output=F, Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)  
  
slope.gwas.mlm = WhitePHSbyEnv %>% filter(term == 'MeanGDDc') %>%
  dplyr::select(taxa, estimate) %>% data.frame() %>%
  GAPIT(Y =.,GD = myGDwheat, GM = myGMwheatGAPIT, model = 'MLM', PCA.total = 2,
        Geno.View.output=F, Major.allele.zero = F, file.output=F,SNP.MAF = 0.05)  

Intercept.gwas.mlm$GWAS %>% filter( maf>0.05) %>%arrange(P.value) %>% slice_head(n = 5)
slope.gwas.mlm$GWAS %>% filter( maf>0.05) %>% arrange(P.value) %>% slice_head(n = 5)

  
  
  
  