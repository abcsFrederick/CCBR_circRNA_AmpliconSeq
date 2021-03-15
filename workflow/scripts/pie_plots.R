
rm(list=ls())
library("tidyverse")
library("RColorBrewer")
n <- 101
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

read_data_table <-function(fname,circRNA){
  d=read.csv(fname,sep="\t",header=TRUE)
  samplelist=colnames(d)[2:length(colnames(d))]
  d %>% separate(col="Name",into=c("circRNA","ASorSS","start","end","fiveDA","threeDA","DA","NT"),sep="##") -> d
  d$fiveDA=NULL
  d$threeDA=NULL
  x=list()
  x[["table"]]=d
  x[["samplelist"]]=samplelist
  x[["circRNAlist"]]=unique(d$circRNA)
  return(x)
}
filter_table_by_cname <-function(d,circRNA){
  k=d$circRNA==circRNA
  d=d[k,]
  return(d)
}
filter_table_by_sname <-function(d,samplename){
  d=data.frame(DA=d$DA,counts=d[[samplename]])
  d=d[order(d$counts,decreasing = TRUE),]
  d=d[!(d$counts==0),]
  return(d)
}
add_perc_sign <- function(s){
  result=""
  if (as.character(s)!=""){
    result=paste0(as.character(s),"%")
  }
  return(result)
}

create_perc_table<-function(df){
  table_percent <- df %>% mutate(DA=DA,perc = round((counts/ sum(counts)) * 100, 1)) %>%
    mutate(DA=DA,labels=perc,y_text=cumsum(perc)-perc/2)
  nr=nrow(table_percent)
  if (nr>4){
    table_percent[5:nrow(table_percent),][["labels"]]=""
  } 
  # min_labels=min(5,nrow(table_percent))
  # if (sum(!table_percent$labels<20)<min_labels){
  #   nr=nrow(table_percent)
  #   if ((min_labels+1)<nr){
  #     table_percent[(min_labels+1):nrow(table_percent),][["labels"]]=""
  #   }
  # }else{
  #   if (sum(table_percent$labels<20)>0){
  #     table_percent[table_percent$labels<20,][["labels"]]=""
  #   }
  # }
  table_percent[table_percent$labels!="",]$labels=paste(table_percent[table_percent$labels!="",]$DA,table_percent[table_percent$labels!="",]$labels)
  table_percent$labels=lapply(table_percent$labels,add_perc_sign)
  table_percent[table_percent$labels!="",]$labels=paste0(table_percent[table_percent$labels!="",]$labels,"(",table_percent[table_percent$labels!="",]$counts,")")
  return(table_percent)
}


args = commandArgs(trailingOnly=TRUE)
args[1] = "mergedquant.filtered.tsv"
d=read_data_table(args[1])
samplelist=d$samplelist
circRNAlist=d$circRNAlist
d=d$table

for(i in 1:length(circRNAlist)) {        
  cname=circRNAlist[[i]][1]
  for (j in 1:length(samplelist)){
    sname=samplelist[[j]][1]
    df=filter_table_by_sname(filter_table_by_cname(d,cname),sname)
    if(nrow(df)>0){
      df=create_perc_table(df)
      fname=paste0(cname,"-",sname,".png")
      # print(fname)
      png(fname)
      # print(ggplot(df, aes(x = "", y = perc, fill = DA)) +
      #   geom_bar(width = 1,stat = "identity") +
      #   geom_label_repel(aes(label = labels, y = y_text)) +
      #   scale_fill_manual(values=sample(col_vector,100,replace=TRUE))+
      #   coord_polar(theta = "y", start = 0) +
      #   theme_light()+
      #   theme(legend.position="none")+
      #   labs(x = "", y = "", title = paste(cname,sname,sep="-")))
      pie(df$perc,labels = df$labels,main= paste0(cname,"-",sname))
      dev.off()
    }
  }
}

