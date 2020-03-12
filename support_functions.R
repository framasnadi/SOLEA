
import<-function(i){
  read_csv(i) %>%dplyr::select(-X1)%>%dplyr::mutate_if(is.character, as.factor)%>%dplyr::mutate(Catch_weight=factor(Catch_weight, levels = c("Low", "Medium_Low", "Medium_High", "High")))%>%dplyr::mutate(status=as.factor(ifelse(Vitality_class=="Dead", 0, 1))) # Import Dataset  
 }

panel.smooth2<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                         cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = 1, ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor ,...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  #r <- abs(cor(x, y))
  r <- (cor(x, y))
  #txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- format(c(r, 0.123456789), digits = 1)[1]
  txt <- paste(prefix, txt, sep = "")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * r)
  text(0.5, 0.5, txt, cex = 3*abs(r))
}
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
}

Collinearity<-function(i){
  Z<-i %>%dplyr::select(-Vitality_class,-Survivability_days, -status) ### Remove not relevant info
 colnames(Z) <- str_replace(names(Z), "_", " ")
 tiff(paste0(plotdir,"Coll_plot.tif"),width = 85, height = 85, units = "mm", res = 1200, pointsize = 6) # graph collinearity
 colplot<-pairs(Z, lower.panel = panel.smooth2, upper.panel = panel.cor, diag.panel = panel.hist)
 dev.off()
 print("Check coplot in Figures folder")
 }

Boruta_screen<-function(i){
  db <- i %>% dplyr::select(-Survivability_days, -Vitality_class) # remove not relevant information for this phase
 boruta.db <- Boruta(status ~., data = db, doTrace = 2)
 decision<-enframe(boruta.db$finalDecision)
 return(decision)
 print(boruta.db$finalDecision)
}

RF_screen<-function(i){
  db <- i %>% dplyr::select(-Survivability_days, -Vitality_class) # remove not relevant information for this phase
  
  # parameters selection
  rf_ranges <- list(ntree = c(500, 1000, 1500, 2000), mtry = 2:5) 
  #rf_tune <- tune(randomForest, as.factor(status) ~ ., data = db, ranges = rf_ranges)
  #tree_best<-as.numeric(rf_tune$best.parameters[1])
  #try_best<-as.numeric(rf_tune$best.parameters[2])
  
  # run the model
  rf <- randomForest(status ~ ., data = db, importance = TRUE, ntree = tree_best, mtry = try_best)
  round(randomForest::importance(rf), 2)
  print(rf)
  # Dotchart of variable importance as measured by a Random Forest
  tiff(paste0(plotdir,"RF.tif"),width = 85, height = 85, units = "mm", res = 1200, pointsize = 5)
  varImpPlot(rf,type=2, main = "Variable importance in Random Forest")
  dev.off()
  return(rf)
} 

PD_plot<-function(i){ # Partial dependence plot gives a graphical depiction of the marginal effect of a variable on the classprobability (classification) or response
  tiff(paste0(plotdir,"PartialDep_RF.tif"),width = 85, height = 85, units = "mm", res = 1200, pointsize = 5)
  op<-par(mfrow=c(2, 3))
  for (j in seq_along(rownames(i))) {
    partialPlot(rf, as.data.frame(db), rownames(i)[j], xlab=rownames(i)[j],ylab ="Centered Log Odds of Dead",main=paste("Partial Dependence on", rownames(i)[j]))
  }
  par(op)
  dev.off()
  return(rf)
}

CT_create<-function(i){
  db <- i %>% dplyr::select(-Survivability_days, -Vitality_class) 
  par(mfrow=c(1, 1))
  tree.rpart <- rpart(status ~ . , data = db ,control = rpart.control(minsplit = as.integer(nrow(db)/as.integer(n_scenarios)), cp = 0.001)  )
  tiff(paste0(plotdir,"CT.tif"),width = 170, height = 85, units = "mm", res = 1200, pointsize = 5)
  rpart.plot(tree.rpart, type = 2, extra = 1)
  dev.off()
  return(tree.rpart)
}

SI_calculation<- function(i){
  name<-as.character(unique(i$scen_set))
  plot_scenario<-i%>% dplyr::count(Vitality_class) %>%dplyr::mutate(Percentage= n/(sum(n))) %>%dplyr::mutate(Scenar=rep(name,nrow(.)))
  return(plot_scenario)
}

Cox_model<-function(i){
  i<-i  %>% na.omit(.) %>%dplyr::filter(Vitality_class != "Dead")
  name<-as.character(unique(i$scen_set))
  i<-i %>%dplyr::mutate(Vitality_class=factor(Vitality_class))
  i<-i%>%dplyr::mutate(Survivability_hours= Survivability_days*24) %>%dplyr::mutate(cens=ifelse(Survivability_hours< censor, 1,0))
  Vitality_class<-relevel(i$Vitality_class, ref="A") ### setup reference level
  modsel <- coxph(Surv(i$Survivability_hours, i$cens) ~ Vitality_class )
  f<-summary(modsel)#;print(modsel)
  ci<-f$conf.int
  ci<-as.tibble(ci) %>% dplyr::mutate(Indicator=str_remove(rownames(ci), "Vitality_class"))%>%dplyr::rename(ul='upper .95', ll='lower .95')%>%dplyr::mutate(ul=1/ul, ll=1/ll)%>%dplyr::select(Indicator, ll,ul)
  p_value<-f$coefficients ## significatività
  p_value<-as.tibble(p_value) %>%dplyr::rename("sign" = `Pr(>|z|)`, "se"='se(coef)') %>% dplyr::mutate(Indicator=str_remove(rownames(p_value), "Vitality_class"), sign=ifelse(sign >= 0.05, "ns", "s"))%>%dplyr::inner_join(., ci, by="Indicator") %>%dplyr::select(Indicator,sign, coef,se, ul, ll) 
  plot_cox<-ggforest(modsel, data = i, main = "Hazard ratio",   cpositions = c(0.02, 0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
  ggsave(paste0(name, ".tiff"), plot_cox,path = plotdir )
  return(p_value)
}

KM_model<-function(i){
  i<-i  %>% na.omit(.) 
  name<-as.character(unique(i$scen_set))
  i$Vitality_class<-as.factor(i$Vitality_class)
  i<-i%>%dplyr::mutate(Survivability_hours= Survivability_days*24) %>%dplyr::mutate(cens=ifelse(Survivability_hours< censor, 1,0))
  Model<-survfit(Surv(i$Survivability_hours, i$cens) ~ Vitality_class, data = i);print(Model)
  # Pairwise comparisons using Peto & Peto test (Peto & Peto, 1972)
  Post_Hoc_KM<-pairwise_survdiff(Surv(Survivability_hours, cens) ~ Vitality_class ,data =i, rho = 1);print(Post_Hoc_KM)
  confAB<-as.data.frame(Post_Hoc_KM$p.value)[1,1]
  confBC<-as.data.frame(Post_Hoc_KM$p.value)[2,2]
  confAC<-as.data.frame(Post_Hoc_KM$p.value)[2,1]
  significance<-tibble(Indicator=c("A", "B", "C"),sign=c("s",NA,NA))
  if(confAB< 0.05){
    significance[2,2]<-"s"
  }else{
    significance[2,2]<-"ns"
  }
  if(confBC< 0.05){
    significance[3,2]<-"s"
  } else {
    significance[3,2]<-"ns"
  }
  z<-summary(Model)
  values<-tibble(Indicator= str_remove(z$strata, "Vitality_class="), time=z$time, coef=z$surv, se=z$std.err, ul=z$upper, ll=z$lower)%>%arrange(Indicator, desc(time))%>%distinct(Indicator,.keep_all=T)%>%dplyr::inner_join(., significance, by="Indicator")%>%dplyr::select(Indicator,sign, coef,se, ul, ll)
    return(values)
}

SR<-function(i,j){
  ##### SI
  name<-as.character(unique(i$Scenar))
  SI<-i %>% dplyr::mutate(Indicator= ifelse(Vitality_class== "Dead","0","1"))%>% dplyr::group_by(Indicator)%>%dplyr::mutate(Percentage= sum(Percentage))%>%dplyr::distinct(Indicator, Percentage,.keep_all=T) %>%dplyr::select(-n, -Vitality_class) %>%dplyr::select(Indicator, Percentage, Scenar)%>%dplyr::filter(Indicator ==1)
  print(SI)
  
  ##### SD
  # weights
  SD_weights<-i %>% dplyr::filter(Vitality_class != "Dead") %>%dplyr::mutate(Percentage=(n/sum(n)))%>%dplyr::rename("Indicator" = "Vitality_class") %>%dplyr::select(-n);print(SD_weights)
  # sp
  
  
  coef_b<-as.numeric()
  coef_c<-as.numeric()
  se_b<-as.numeric()
  se_c<-as.numeric()
  
  if(surv_data=="Absolute"){
    
    j<-j%>%dplyr::mutate(se=ifelse(is.nan(se)==T,0,se))
    
    coef_a<-j[j$Indicator=="A",]$coef
    se_a<-j[j$Indicator=="A",]$se
    
    if(j[j$Indicator=="B",]$sign=="s"){
      coef_b<-j[j$Indicator=="B",]$coef
      se_b<-j[j$Indicator=="B",]$se
    }else{
      coef_b<-coef_a
      se_b<-se_a
    }
    if(j[j$Indicator=="C",]$sign=="s"){
      coef_c<-j[j$Indicator=="C",]$coef
      se_c<-j[j$Indicator=="C",]$se
    }else if (j[j$Indicator=="C",]$sign=="ns"){
      coef_c<-coef_b
      se_c<-se_b
    } else {
      coef_c<-coef_a
      se_c<-se_a
    }
    
    qq<-rep(0,10000)
    for(k in 1:10000){
      qq[k]<-SI$Percentage *    ((as.numeric(SD_weights[1,2])* (rnorm(1,as.numeric(coef_a),as.numeric(se_a))))+
                                   (as.numeric(SD_weights[2,2])*(rnorm(1,as.numeric(coef_b),as.numeric(se_b))))+
                                   (as.numeric(SD_weights[3,2])*(rnorm(1,as.numeric(coef_c),as.numeric(se_c)))))
    }
    
  }else{
  
  if(j[j$Indicator=="B",]$sign=="s"){
    coef_b<-j[j$Indicator=="B",]$coef
    se_b<-j[j$Indicator=="B",]$se
  }else{
    coef_b<-1
    se_b<-0
  }
  if(j[j$Indicator=="C",]$sign=="s"){
    coef_c<-j[j$Indicator=="C",]$coef
    se_c<-j[j$Indicator=="C",]$se
  }else{
    coef_c<-1
    se_c<-0
  }
  
  qq<-rep(0,10000)
  for(k in 1:10000){
    qq[k]<-SI$Percentage *    ((as.numeric(SD_weights[1,2])* 1)+
                              (as.numeric(SD_weights[2,2])*(1/(exp(rnorm(1,as.numeric(coef_b),as.numeric(se_b))))))+
                              (as.numeric(SD_weights[3,2])*(1/(exp(rnorm(1,as.numeric(coef_c),as.numeric(se_c)))))))
  }
  
  }
  
  
  SR_mean<-as.numeric(mean(qq))
  SR_CI<-as.numeric(quantile(qq, c(0.05, 0.95)))
  SR<-tibble(SR = SR_mean, upper_ci = SR_CI[2], low_ci=SR_CI[1], Scenario = name)
  print(SR)
  return(SR)
}

SR_propagazione<-function(i,j){
  ##### SI
  name<-as.character(unique(i$Scenar))
  SI<-i %>% dplyr::mutate(Indicator= ifelse(Vitality_class== "Dead","0","1"))%>% dplyr::group_by(Indicator)%>%dplyr::mutate(Percentage= sum(Percentage))%>%dplyr::distinct(Indicator, Percentage,.keep_all=T) %>%dplyr::select(-n, -Vitality_class) %>%dplyr::select(Indicator, Percentage, Scenar)%>%dplyr::filter(Indicator ==1)
  print(SI)
  
  ##### SD
  # weights
  SD_weights<-i %>% dplyr::filter(Vitality_class != "Dead") %>%dplyr::mutate(Percentage=(n/sum(n)))%>%dplyr::rename("Indicator" = "Vitality_class") %>%dplyr::select(-n);print(SD_weights)
  # sp
  
  
  coef_b<-as.numeric()
  coef_c<-as.numeric()
  se_b<-as.numeric()
  se_c<-as.numeric()
  
  if(surv_data=="Absolute"){
    
    j<-j%>%dplyr::mutate(se=ifelse(is.nan(se)==T,0,se))
    
    coef_a<-j[j$Indicator=="A",]$coef
    se_a<-j[j$Indicator=="A",]$se
    ul_a<-j[j$Indicator=="A",]$ul
    ll_a<-j[j$Indicator=="A",]$ll
    
    if(j[j$Indicator=="B",]$sign=="s"){
      coef_b<-j[j$Indicator=="B",]$coef
      se_b<-j[j$Indicator=="B",]$se
      ul_b<-j[j$Indicator=="B",]$ul
      ll_b<-j[j$Indicator=="B",]$ll
    }else{
      coef_b<-coef_a
      se_b<-se_a
      ul_b<-ul_a
      ll_b<-ll_a
    }
    if(j[j$Indicator=="C",]$sign=="s"){
      coef_c<-j[j$Indicator=="C",]$coef
      se_c<-j[j$Indicator=="C",]$se
      ul_c<-j[j$Indicator=="C",]$ul
      ll_c<-j[j$Indicator=="C",]$ll
    }else if (j[j$Indicator=="C",]$sign=="ns"){
      coef_c<-coef_b
      se_c<-se_b
      ul_c<-ul_b
      ll_c<-ll_b
    } else {
      coef_c<-coef_a
      se_c<-se_a
      ul_c<-ul_a
      ll_c<-ll_a
    }
    
    ######################################
    # calcolo propagazione degli errori  #
    ######################################
    
    var <- c("si", "wa", "sa", "wb", "sb", "wc", "sc") # variabili
    dd <- vector(mode = "list", length = length(var)) # creazione lista
    vv <- rep(0, length(var)) # per i valori eval()
    for (i in 1:length(var)) {
      dd[[i]] <- D(expression(si*(sa*wa + sb*wb + sc*wc)), var[i]) # calcolo derivate della funzione expression() 
      si <- SI$Percentage; wa <- as.numeric(SD_weights[1,2]); sa <- coef_a; wb <-as.numeric(SD_weights[2,2]); sb <- coef_b; wc <- as.numeric(SD_weights[3,2]); sc <- coef_c; # valori medi variabili
      vv[i] <- eval(dd[[i]]) # valutazione numerica delle derivate
      err <- c(0, 0, (ul_a - ll_a)/2, 0, (ul_b - ll_b)/2, 0, (ul_c - ll_c)/2) # errori variabili, nello stesso ordine delle variabili
      med <- eval(expression(si*(sa*wa + sb*wb + sc*wc))) # valore medio della funzione
      xx <- sqrt(t(vv^2)%*%err^2) # percentuale errore medio della funzione
    }
    
    
  }else{
    
    if(j[j$Indicator=="B",]$sign=="s"){
      coef_b<-j[j$Indicator=="B",]$coef
      se_b<-j[j$Indicator=="B",]$se
      ul_b<-j[j$Indicator=="B",]$ul
      ll_b<-j[j$Indicator=="B",]$ll
      
    }else{
      coef_b<-1
      se_b<-0
      ul_b<-1
      ll_b<-1
    }
    if(j[j$Indicator=="C",]$sign=="s"){
      coef_c<-j[j$Indicator=="C",]$coef
      se_c<-j[j$Indicator=="C",]$se
      ul_c<-j[j$Indicator=="C",]$ul
      ll_c<-j[j$Indicator=="C",]$ll
    }else{
      coef_c<-1
      se_c<-0
      ul_c<-1
      ll_c<-1
    }
    
    var <- c("si", "wa", "sa", "wb", "sb", "wc", "sc") # variabili
    dd <- vector(mode = "list", length = length(var)) # creazione lista
    vv <- rep(0, length(var)) # per i valori eval()
    for (i in 1:length(var)) {
      dd[[i]] <- D(expression(si*(sa*wa + sb*wb + sc*wc)), var[i]) # calcolo derivate della funzione expression() 
      si <- SI$Percentage; wa <- (as.numeric(SD_weights[1,2])); sa <- 1; wb <-(as.numeric(SD_weights[2,2])); sb <- exp(-coef_b); wc <- as.numeric(SD_weights[3,2]); sc <- exp(-coef_c); # valori medi variabili
      vv[i] <- eval(dd[[i]]) # valutazione numerica delle derivate
      err <- c(0, 0, 0, 0, (ul_b - ll_b)/2, 0, (ul_c - ll_c)/2) # errori variabili, nello stesso ordine delle variabili
      med <- eval(expression(si*(sa*wa + sb*wb + sc*wc))) # valore medio della funzione
      xx <- sqrt(t(vv^2)%*%err^2) # percentuale errore medio della funzione
    }
    
    
  }
  
  
  SR_mean<-as.numeric(med)
  SR<-tibble(SR = SR_mean, upper_ci = as.numeric(med)+xx[[1]], low_ci=as.numeric(med)-xx[[1]], Scenario = name)
  print(SR)
  return(SR)
}


