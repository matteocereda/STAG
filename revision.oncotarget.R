library(xlsx)
library(scales)
library(plyr)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# fin = "Input/Figure2_CAL51_test.xlsx"
proliferation_norm_24h_log2 = function(fin
                                       , pout
                                       , lvs=c("NTC","STAG1","STAG2","STAG1_2")
                                       , clrs=c('grey','red','blue','green')
                                       , tgene=c("STAG1","STAG2","STAG1_2")
                                       , tcont="NTC"
                                       , tpoints=c(72,96)){
  require(xlsx)
  x = read.xlsx(fin, 1)

  x = melt(x, id.vars = c("Time_point","date"))
  
  # is_na=which(is.na(x$value))
  
  tmp = strsplit(as.character(x$variable),"\\.")
  
  x$gene = sapply( tmp, "[[", 1)
  x$gene = factor( x$gene, levels=lvs)
  
  tmp[ which(sapply(tmp,length)==1)] = lapply(tmp[ which(sapply(tmp,length)==1)], function(x) c(0,0))
  
  x$rep = factor(sapply(tmp, "[[", 2), levels=c(0,1,2) )
  
  # normalized each replicate by the 24 value
  
  x = ddply(x, .(date, gene, rep), mutate, value_norm_24h=value/value[Time_point==24])
  
  # get the log2
  
  x$log2_value_norm_24h = log2(x$value_norm_24h)
  
  # measure mean and se per experiemnt
  
  sy =  summarySE(x, measurevar = "log2_value_norm_24h", groupvars = c("gene","Time_point","date"), na.rm = T)
  
  # measure mean and se across experiemnt

  s =  summarySE(sy, measurevar = "log2_value_norm_24h", groupvars = c("gene","Time_point"), na.rm = T)
  

  write.xlsx(x, file=paste0(pout,"_table.xlsx"),'Proliferation', row.names = F)
  write.xlsx(sy[,c('gene','date','Time_point',"N","log2_value_norm_24h",'se')], file=paste0(pout,"_table.xlsx"),'Summary', append = T, row.names = F)
  write.xlsx(s[,c('gene','Time_point',"N","log2_value_norm_24h",'se')],         file=paste0(pout,"_table.xlsx"),'Summary_exp', append = T, row.names = F)
  
  # test ALL experiment
  
  table=cbind()
  for(z in tpoints){
    f = subset(sy, Time_point==z)
    p = rep(z, length(tgene)*length(tcont)+1); names(p) = c('Time_point',paste0(tgene,"_",tcont))
    
    for(i in tgene) for(j in tcont) p[paste0(i,"_",j)] =  with(f, t.test(  log2_value_norm_24h[gene==i]
                                                                           , log2_value_norm_24h[gene==j],
                                                                           alternative="less" )$p.value) #, alt="l"
    table=cbind(table,p)
  }
  
  # # f = subset(x, Time_point==72)
  # f = subset(sy, Time_point==72)
  # p = rep(72, length(tgene)+1); names(p) = c('Time_point',tgene)
  # 
  # for(i in tgene) p[i] =  with(f, t.test(  log2_value_norm_24h[gene==i]
  #                                          , log2_value_norm_24h[gene=="NTC"], 
  #                                          alternative="less" )$p.value) #, alt="l"
  # table=p
  # # f = subset(x, Time_point==96)
  # f = subset(sy, Time_point==96)
  # 
  # p = rep(96, length(tgene)+1); names(p) = c('Time_point',tgene)
  # 
  # for(i in tgene) p[i] =  with(f, t.test(  log2_value_norm_24h[which(gene==i)]
  #                                          , log2_value_norm_24h[which(gene=="NTC")], 
  #                                          alternative="less" )$p.value)
  # 
  # table=cbind(table,p)
  m = table[2:nrow(table),1:ncol(table)];
  if( ncol(table)==1) m = as.matrix(m)
  
  colnames(m) = table[1,]
  write.csv(m,file=paste0(pout,".T_test.csv"))  
  
  
  pdf(file=paste0(pout,".pdf"), h=6, w=6)
  
  pl =
    ggplot(s,aes(x=Time_point, y=log2_value_norm_24h,color=gene))+
    geom_line()+
    geom_errorbar(aes(ymin=log2_value_norm_24h-se, ymax=log2_value_norm_24h+se),width=.5)+
    geom_point()+
    theme_classic()+
    scale_color_manual(values=clrs)+
    scale_x_continuous(breaks =c(24,48,72,96))+
    ggtitle("AVERAGE across experiments")+
    annotation_custom(tableGrob(scientific_format(2)(m), theme=ttheme_minimal(base_size=8)), xmin=24, xmax=48, ymin=0.5, ymax=1)
  print(pl)    
  
  for(gg in unique(sy$date)){
    
    table=cbind()
    for(z in tpoints){
      f = subset(x, date==gg & Time_point==z)
      p = rep(z, length(tgene)*length(tcont)+1); names(p) = c('Time_point',paste0(tgene,"_",tcont))
      
      for(i in tgene) for(j in tcont) p[paste0(i,"_",j)] =  with(f, t.test(  log2_value_norm_24h[gene==i]
                                                                             , log2_value_norm_24h[gene==j],
                                                                             alternative="less" )$p.value) #, alt="l"
      table=cbind(table,p)
    }
    # f = subset(x, date==gg & Time_point==96)
    # p = rep(96, length(tgene)+1); names(p) = c('Time_point',tgene)
    # 
    # p = rep(96, length(tgene)*length(tcont)+1); names(p) = c('Time_point',paste0(tgene,"_",tcont))
    # 
    # for(i in tgene) for(j in tcont) p[paste0(i,"_",j)] =  with(f, t.test(  log2_value_norm_24h[gene==i]
    #                                                                        , log2_value_norm_24h[gene==j],
    #                                                                        alternative="less" )$p.value) #, alt="l"
    # 
    # 
    # table=cbind(table,p)
    m = table[2:nrow(table),1:ncol(table)];
    if( ncol(table)==1) m = as.matrix(m)
    colnames(m) = table[1,]
    
    pl =
      ggplot(subset(sy,date==gg),aes(x=Time_point, y=log2_value_norm_24h,color=gene))+
      geom_line()+
      geom_errorbar(aes(ymin=log2_value_norm_24h-se, ymax=log2_value_norm_24h+se),width=.5)+
      geom_point()+
      theme_classic()+
      scale_color_manual(values=clrs)+
      scale_x_continuous(breaks =c(24,48,72,96))+
      ggtitle(gg)+
      annotation_custom(tableGrob(scientific_format(2)(m), theme=ttheme_minimal(base_size=8)), xmin=24, xmax=48, ymin=0.5, ymax=1)
    print(pl)    
  }
  # grid.newpage()
  # grid.draw(tableGrob(scientific_format(2)(m), theme=ttheme_minimal(base_size=10)))
  dev.off()

  print("done")
}

# Fig 2B
proliferation_norm_24h_log2("Input/Figure2_CAL51_test.xlsx","Paper/Revision/Original_data/Figure_2/2B_CAL51.25_11_16")
# proliferation_norm_24h_log2("Input/Figure2B_CAL51_10_05_15.csv","Paper/Revision/Original_data/Figure_2/2B_CAL51")
# Fig 2E

proliferation_norm_24h_log2("Input/Figure2_MCF7.xlsx","Paper/Revision/Original_data/Figure_2/2E_MCF7.25_11_16")
# proliferation_norm_24h_log2("Input/Figure2E_MCF7_26_10_15.csv","Paper/Revision/Original_data/Figure_2/2E_MCF7")

# Fig 3H
proliferation_norm_24h_log2("Input/Proliferation_CAL51_L161fs_STAG2.xlsx"
                            ,"Paper/Revision/Original_data/Figure_3/3H_CAL51_L161fs_STAG2.25_11_16"
                            ,c("NTC","STAG1"), tgene=c("STAG1"))
# proliferation_norm_24h_log2("Input/Figure3H_CAL51_L161fs_STAG2_26_01_16.csv"
#                             ,"Paper/Revision/Original_data/Figure_3/3H_CAL51_L161fs_STAG2"
#                             ,c("NTC","STAG1"), tgene=c("STAG1"))

# Fig 4E
proliferation_norm_24h_log2("Input/Figure_4E_STAG1aL_3_exp_1104.xlsx"
                            ,"Paper/Revision/Original_data/Figure_4/4E_CAL51_cdLenti_STAG1.25_11_16"
                            ,c("NTC","STAG2"),c("grey",'blue'), tgene = 'STAG2')
# proliferation_norm_24h_log2("Input/Figure4E_CAL51_cdLenti_STAG1_26_10_15.csv"
#                             ,"Paper/Revision/Original_data/Figure_4/4E_CAL51_cdLenti_STAG1"
#                             ,c("NTC","STAG2"),c("grey",'blue'), tgene = 'STAG2')
proliferation_norm_24h_log2("Input/Figure_4E_CAL51_cdLenti_STAG1.25_11_16.xlsx"
                            ,"Paper/Revision/Original_data/Figure_4/4E_CAL51_cdLenti_STAG1.02_12_16"
                            ,c("Cas9_NTC","STAG1L_NTC","STAG1L_STAG2")
                            ,c("grey","black",'blue')
                            ,tgene = 'STAG1L_STAG2', tcont=c("Cas9_NTC","STAG1L_NTC"), tpoints = 72)

# EXPRESSION ==================

expression_norm_NTC = function(fin
                                       , pout
                                       , lvs=c("NTC","STAG1","STAG2","STAG1_2")
                                       , clrs=c('red','blue')
                                       , tgene=c("STAG1","STAG2")
                                       , tcont="B2M"
                                       , ref="NTC"
                                       , lmts=c(0,2)
                                       , bb = c(0,1,2)
                                       ){
  require(xlsx)
  x = read.xlsx(fin, 1)
  y = x[,tgene]-x[,tcont]
  # if( length(tgene)==1) y = as.matrix(y)
  y = cbind(x[,1],y)
  
  if(length(tgene)==1){
    y = as.data.frame(y)
    colnames(y) =c('ref', tgene)
    y[,2] =  as.numeric(y[,2])
  }else{
    colnames(y)[1]='ref'
    
  }
  
  
  x = melt(y, id.vars = 'ref'); 
  colnames(x) = c('condition','gene','CT')
  x$condition = factor(x$condition,lvs)
  x = ddply(x, .(condition, gene), mutate, rep=1:length(CT))  

  x$rep =factor(x$rep, levels=c(1,2,3) )
  # normalized each replicate by the value in NTC
  
  x = ddply(x, .( gene, rep), mutate, CT_norm_to_NTC=CT[condition==ref]-CT)
  
  x = ddply(x, .( gene, rep), mutate, ratio=2^CT_norm_to_NTC)
  
   # measure mean and se per replicate
  
  sy =  summarySE(x, measurevar = "ratio", groupvars = c("gene","condition"), na.rm = T)
  
  write.xlsx(x, file=paste0(pout,"_table.xlsx"),'Expression', row.names = F)
  write.xlsx(sy[,c('gene','condition',"N","ratio",'se')], file=paste0(pout,"_table.xlsx"),'Summary', append = T, row.names = F)
  # write.xlsx(s[,c('gene','Time_point',"N","log2_value_norm_24h",'se')],         file=paste0(pout,"_table.xlsx"),'Summary_exp', append = T, row.names = F)
  
  # test ALL experiment
  
  # table=cbind()
  # f = subset(sy, Time_point==z)
  # p = rep(z, length(tgene)*length(tcont)+1); names(p) = c('Time_point',paste0(tgene,"_",tcont))
  #   
  # for(i in tgene) for(j in tcont) p[paste0(i,"_",j)] =  with(f, t.test(  log2_value_norm_24h[gene==i]
  #                                                                          , log2_value_norm_24h[gene==j],
  #                                                                          alternative="less" )$p.value) #, alt="l"
  # table=cbind(table,p)
  # 
  #  m = table[2:nrow(table),1:ncol(table)];
  # if( ncol(table)==1) m = as.matrix(m)
  # 
  # colnames(m) = table[1,]
  # write.csv(m,file=paste0(pout,".T_test.csv"))  
  # 
  
  pdf(file=paste0(pout,".pdf"), h=6, w=6)
  
  pl =
    ggplot(sy,aes(x=condition, y=ratio, fill=gene))+
    geom_bar(stat = 'identity', position=position_dodge())+
    geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),width=.5,position=position_dodge(.9))+
    theme_classic()+
    scale_fill_manual(values=clrs)+
    scale_y_continuous(limits = lmts,breaks =bb)+
    ylab('Normalized ratio')+xlab('')
        # annotation_custom(tableGrob(scientific_format(2)(m), theme=ttheme_minimal(base_size=8)), xmin=24, xmax=48, ymin=0.5, ymax=1)
  print(pl)    
  
  dev.off()
  
  print("done")
}
#fig 2A
expression_norm_NTC("Input/Figure2A_CAL51_KD_expression.xlsx","Paper/Revision/Original_data/Figure_2/2A_CAL51_KD_expression.02_12_16")

#fig 2D
expression_norm_NTC("Input/Figure2D_MCF7_KD_expression.xlsx","Paper/Revision/Original_data/Figure_2/2D_MCF7_KD_expression.02_12_16"
                    , lmts = c(0,3), bb=c(0,1,2,3))

#fig 2H
expression_norm_NTC("Input/Figure2H_SKES1_KD_expression_exp1.xlsx","Paper/Revision/Original_data/Figure_2/2H_SKES1_KD_expression.exp1.02_12_16"
                    , lvs=c("NTC","STAG1","STAG2")
                    , clrs=c('red')
                    , tgene=c("STAG1")
                    , tcont="B2M"
                    , ref="NTC"
                    , lmts=c(0,2)
                    , bb = c(0,1,2)

                    )
expression_norm_NTC("Input/Figure2H_SKES1_KD_expression_paper.xlsx","Paper/Revision/Original_data/Figure_2/2H_SKES1_KD_expression.paper.02_12_16"
                    , lvs=c("NTC","STAG1","STAG2")
                    , clrs=c('red')
                    , tgene=c("STAG1")
                    , tcont="B2M"
                    , ref="NTC"
                    , lmts=c(0,10)
                    , bb = c(0,1,2)
)

#fig4c
expression_norm_NTC("Input/Figure4C_STAG1_KD_expression.xlsx","Paper/Revision/Original_data/Figure_4/4C_STAG1_KD_expression.02_12_16"
                    , lvs=c("NTC","STAG1aL")
                    , clrs=c('red')
                    , tgene=c("STAG1")
                    , tcont="B2M"
                    , ref="NTC"
                    , lmts=c(0,2)
                    , bb = c(0,1,2)
)


