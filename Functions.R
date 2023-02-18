## NMF functions ####
library(ggplot2)
mytheme <- theme_classic() + 
  theme(axis.text = element_text(colour="black"),
        panel.background = element_blank(),
        panel.grid = element_blank())
read_NMF <- function(NMF_INPUT,selected_rank){
  NMF <- readRDS(NMF_INPUT)
  nmf_fit <- NMF$fit[[selected_rank-1]]
  
  G <- basis(nmf_fit) %>% as.matrix()
  G_c <- coefficients(nmf_fit) %>% as.matrix() 
  
  rownames(G_c)=paste0("Metagene_",c(1:selected_rank))
  colnames(G)=paste0("Metagene_",c(1:selected_rank))
  
  return(list(nmf_fit=nmf_fit,basis_matrix=G,coef_matrix=G_c))
}

selectBestFitNMF <- function(NMF_INPUT,byw=c("cophenetic")){
  library(NMF)
  library(patchwork)
  cop_data <- getMetricsNMF(NMF_INPUT)
  cop_data %>% filter(variable==byw) %>% 
    pull(3) %>% .[2:length(.)] %>% which.max()+2 ##except for 2
}

getMetricsNMF <-  function(NMF_INPUT){
  library(NMF)
  library(patchwork)
  NMF <- readRDS(NMF_INPUT) #All NMF runs
  cop_data <- as.data.frame(plot(NMF)[1]$data) 
  return(cop_data)
}


plotRanksNMF <- function(NMF_INPUT,cohort){
  library(NMF)
  library(patchwork)
  
  best_rank_cop <- selectBestFitNMF(NMF_INPUT, "cophenetic") 
  best_rank_sillheoutte <- selectBestFitNMF(NMF_INPUT, "silhouette.consensus")
  
  message(paste0("The best rank according to cophenetic correlation is: ", best_rank_cop,
                 "\nThe best rank according to silhouette is: ",best_rank_sillheoutte))  
  cop_data <- getMetricsNMF(NMF_INPUT)
  cop_cor <- ggplot(cop_data %>% filter(variable=="cophenetic"),aes(rank,value))+
    geom_point(col="hotpink3",size=3)+geom_line()+
    ylab("Cophenetic correlation")+
    scale_x_continuous(breaks = c(2:max(cop_data$rank %>% unique)),
                       labels = c(2:max(cop_data$rank %>% unique)))+
    mytheme+
    ggtitle(paste0(cohort," cophenetic correlation"))+
    geom_vline(xintercept = best_rank_cop,lty=2)
  
  sillhouette <- ggplot(cop_data %>% filter(variable==c("silhouette.consensus")),aes(rank,value))+
    geom_point(col="mediumseagreen",size=3)+geom_line()+
    ylab("Consensus silhouette")+
    scale_x_continuous(breaks = c(2:max(cop_data$rank %>% unique)),
                       labels = c(2:max(cop_data$rank %>% unique)))+
    mytheme+
    ggtitle(paste0(cohort," silhouette value across ranks"))+
    geom_vline(xintercept = best_rank_sillheoutte,lty=2)
  
  cop_cor+sillhouette
  
}

plotConsensusSelected <- function(NMF_INPUT, desired_fits){
  library(NMF)
  NMF <- readRDS(NMF_INPUT)
  listOfPlots <- list()
  for (i in desired_fits){
    listOfPlots[[i]] <- NMF::consensusmap(NMF$fit[[i-1]],tracks = c(),
                                          main = paste0("rank = ",i))
  }
  
}

getSubtypesNMF <- function(NMF_INPUT,selected_rank){
  NMF <- readRDS(NMF_INPUT)
  selected_fit <- NMF$fit[[selected_rank-1]]
  subtypes <- apply(coefficients(selected_fit),2,which.max) %>%
    as.data.frame(.) %>% tibble::rownames_to_column("sample") %>%
    rename(subtype=2) %>% mutate(subtype=as.character(subtype))
  return(subtypes)
}

getfitNMF <-  function(NMF_INPUT,selected_rank){
  library(NMF)
  library(patchwork)
  NMF <- readRDS(NMF_INPUT) #All NMF runs
  selected_fit <- NMF$fit[[selected_rank-1]]
  return(selected_fit)
}

survPlot <- function(data,Time,Event,var,palette,p_val=F,risktable=F){
  library(survival);library(survminer)
  f <- as.formula(paste0("Surv(",Time,",",Event,")~",paste(var,collapse="+")))
  ff=surv_fit(f, data=data)
  plot=ggsurvplot(ff,data=data,pval=T,risk.table = risktable,
                  palette = palette,
                  legend.labs=sort(unique(data[,var])),
                  ggtheme = theme_bw()+
                    theme(panel.grid = element_blank(),
                          panel.background = element_blank(),
                          legend.title = element_blank(),
                          axis.text = element_text(colour="black")))+
    xlab("Time")
  if(p_val){return(list=c(plot=plot,p_val=surv_pvalue(ff)))}
  else{return(plot)}
}
surv_p_val_univariate=function(data,Time,Event,var,test=c("logrank","wald")){
  library(survival)
  f <- as.formula(paste0("Surv(",Time,",",Event,")~",paste(var,collapse="+")))
  if(test=="wald"){
    p=summary(coxph(f, data=data))$waldtest[3] %>% as.numeric()
  }
  else if (test=="logrank"){
    p=summary(coxph(f, data=data))$sctest[3] %>% as.numeric()
  }
  return(p)
}


### Pathway activity score calculation ####
calc_pathway_activity_subtypes=function(expression_full,subtype_df,col=c(brewer.pal(4,"Set1"))){
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(fgsea)
  library(data.table)
  library(openxlsx)
  ##'calculate aggregate score for KNOWN HCC subtypes and choose the best fit for each sample 
  ##'Comprehensive molecular and immunological characterization of hepatocellular carcinoma
  ##'method adapted
  ##'subtype df : columns "sample", "subtype"
  genesets=gmtPathways("data/c2.cgp.v7.0.symbols.gmt")
  Hoshida=genesets[grepl("HOSHIDA",names(genesets))][5:7]
  Chiang_genesets=genesets[grepl("CHIANG",names(genesets))][c(1,3,5,7,9)]
  Lee_genesets=genesets[grepl("LEE",names(genesets))][c(32,33)]
  Boyault_gs=genesets[grepl("BOYAULT",names(genesets))][c(9,4,13)]
  Lee_hepatoblast=genesets[grepl("LEE",names(genesets))][c(42)]
  Roessler_met=genesets[grepl("ROESSLER",names(genesets))][1]
  YAMASHITA=genesets[grepl("YAMASHITA",names(genesets))][c(3,5)]
  ImmuneClass=fread("data/Sia_et_al_immune_class_signature.txt") %>% 
    as.list()
  ##calculate activitiy score
  
  Boyault_scores=sapply(1:length(Boyault_gs),function(i) apply(apply(expression_full[Boyault_gs[[i]],],1,percent_rank) %>% 
                                                                 t() %>% as.data.frame() %>% 
                                                                 na.omit(),2,function(x) (mean(x,na.rm=T)-0.5))) %>% as.data.frame()
  colnames(Boyault_scores)=names(Boyault_gs)
  
  Boyault_subtypes=Boyault_scores%>% 
    apply(.,1,which.max) %>% as.data.frame(.) %>% 
    rename(Boyault=1) %>% tibble::rownames_to_column("sample") %>% 
    mutate(Boyault=lapply(names(Boyault_gs)[.$Boyault], function(x) strsplit(x,split = "_")[[1]][5]) %>% unlist())
  
  
  Hoshida_scores=sapply(1:length(Hoshida),function(i) apply(apply(expression_full[Hoshida[[i]],],1,percent_rank) %>% 
                                                              t() %>% as.data.frame() %>% 
                                                              na.omit(),2,function(x) (mean(x,na.rm=T)-0.5))) %>% as.data.frame()
  colnames(Hoshida_scores)=names(Hoshida)
  Hoshida_subtypes=Hoshida_scores%>% 
    apply(.,1,which.max) %>% as.data.frame(.) %>% 
    rename(Hoshida=1) %>% tibble::rownames_to_column("sample") %>% mutate(Hoshida=paste0("S",Hoshida))
  
  
  
  Chiang_scores=sapply(1:length(Chiang_genesets),function(i) apply(apply(expression_full[Chiang_genesets[[i]],],1,percent_rank) %>% 
                                                                     t() %>% as.data.frame() %>% 
                                                                     na.omit(),2,function(x) (mean(x,na.rm=T)-0.5))) %>% as.data.frame()
  colnames(Chiang_scores)=names(Chiang_genesets)
  Chiang_subtypes=Chiang_scores%>% 
    apply(.,1,which.max) %>% as.data.frame(.) %>% 
    rename(Chiang=1) %>% tibble::rownames_to_column("sample") %>% 
    mutate(Chiang=lapply(names(Chiang_genesets)[.$Chiang], function(x) strsplit(x,split = "_")[[1]][5]) %>% unlist())
  
  Lee_scores=sapply(1:length(Lee_genesets),function(i) apply(apply(expression_full[Lee_genesets[[i]],],1,percent_rank) %>% 
                                                               t() %>% as.data.frame() %>% 
                                                               na.omit(),2,function(x) (mean(x,na.rm=T)-0.5))) %>% as.data.frame()
  
  colnames(Lee_scores)=names(Lee_genesets)
  yamashita_scores=sapply(1:length(YAMASHITA),function(i) apply(apply(expression_full[YAMASHITA[[i]],],1,percent_rank) %>% 
                                                                  t() %>% as.data.frame() %>% 
                                                                  na.omit(),2,function(x) (mean(x,na.rm=T)-0.5))) %>% as.data.frame()
  
  colnames(yamashita_scores)=names(YAMASHITA)
  
  
  Lee_subtypes=Lee_scores%>% 
    apply(.,1,which.max) %>% as.data.frame(.) %>% 
    rename(Lee=1) %>% tibble::rownames_to_column("sample") %>% 
    mutate(Lee=lapply(names(Lee_genesets)[.$Lee], function(x) paste0("SURVIVAL_",strsplit(x,split = "_")[[1]][5])) %>% unlist())
  
  Lee_hepatoblast_score=sapply(1:length(Lee_hepatoblast),function(i) apply(apply(expression_full[Lee_hepatoblast[[i]],],1,percent_rank) %>% 
                                                                             t() %>% as.data.frame() %>% 
                                                                             na.omit(),2,function(x) (mean(x,na.rm=T)-0.5))) %>% as.data.frame()
  
  colnames(Lee_hepatoblast_score)=names(Lee_hepatoblast)
  
  Roessler_met_score=sapply(1:length(Roessler_met),function(i) apply(apply(expression_full[Roessler_met[[i]],],1,percent_rank) %>% 
                                                                       t() %>% as.data.frame() %>% 
                                                                       na.omit(),2,function(x) (mean(x,na.rm=T)-0.5))) %>% as.data.frame()
  
  colnames(Roessler_met_score)=names(Roessler_met)
  
  Immune_score=sapply(1:length(ImmuneClass),function(i) apply(apply(expression_full[ImmuneClass[[i]],],1,percent_rank) %>% 
                                                                t() %>% as.data.frame() %>% 
                                                                na.omit(),2,function(x) (mean(x,na.rm=T)-0.5))) %>% as.data.frame()
  
  colnames(Immune_score)=names(ImmuneClass)
  
  to_merge=Lee_hepatoblast_score %>% tibble::rownames_to_column("sample")
  to_merge2=Roessler_met_score %>% tibble::rownames_to_column("sample")
  to_merge3=yamashita_scores %>% tibble::rownames_to_column("sample")
  to_merge4=Immune_score %>% tibble::rownames_to_column("sample")
  
  other_studies=inner_join(Hoshida_subtypes,Chiang_subtypes) %>% inner_join(.,Lee_subtypes) %>% inner_join(.,subtype_df) %>% 
    inner_join(.,Boyault_subtypes)%>% inner_join(.,to_merge) %>% inner_join(.,to_merge2) %>% inner_join(.,to_merge3) %>% inner_join(.,to_merge4) %>% 
    arrange(subtype)  %>% mutate(subtype=as.character(subtype))  %>% rename(Lee_HB=7,Roessler_met=8,Yamashita_EpCAM=9,Yamashita_Stem=10)
  
  
  ##compare with my assignment
  
  #associations=sapply(colnames(other_studies)[c(2,3,4,6,7,8,9,10)],function(x) calculate_p_value(other_studies,"subtype",x))
  col_fun = circlize::colorRamp2(c(min(other_studies$Lee_HB), 0, max(other_studies$Lee_HB)), c("#4d9221", "white", "#c51b7d"))
  col_fun2 = circlize::colorRamp2(c(min(other_studies$Roessler_met), 0, max(other_studies$Roessler_met)), c("#4d9221", "white", "#c51b7d"))
  col_fun3 = circlize::colorRamp2(c(min(other_studies$Yamashita_EpCAM), 0, max(other_studies$Yamashita_EpCAM)), c("#4d9221", "white", "#c51b7d"))
  col_fun4 = circlize::colorRamp2(c(min(other_studies$Yamashita_Stem), 0, max(other_studies$Yamashita_Stem)), c("#4d9221", "white", "#c51b7d"))
  col_fun5 = circlize::colorRamp2(c(min(other_studies$ImmuneClassGenes), 0, max(other_studies$ImmuneClassGenes)), c("#4d9221", "white", "#c51b7d"))
  ha=HeatmapAnnotation(Lee_2004=other_studies$Lee,
                       Lee_HB_2006=other_studies$Lee_HB,
                       Boyault_2007=other_studies$Boyault,
                       Chiang_2008=other_studies$Chiang,
                       Yamahita_Epcam_2008=other_studies$Yamashita_EpCAM,
                       Yamahita_Stem_2008=other_studies$Yamashita_Stem,
                       Hoshida_2009=other_studies$Hoshida,
                       Roessler_met_2010=other_studies$Roessler_met,
                       Sia_et_al_2017=other_studies$ImmuneClassGenes,
                       show_annotation_name = T,border=T,show_legend = T,
                       annotation_name_side = "left",
                       col = list(Lee_2004=c("SURVIVAL_DN"="hotpink4","SURVIVAL_UP"="khaki3"),
                                  Hoshida_2009=c("S1"=brewer.pal(8,"Set1")[3],"S2"=brewer.pal(8,"Set1")[4],"S3"=brewer.pal(8,"Set1")[5]),
                                  Chiang_2008=c("CTNNB1"=brewer.pal(5,"Accent")[1],"INTERFERON"=brewer.pal(5,"Accent")[2],
                                                "POLYSOMY7"=brewer.pal(5,"Accent")[3],"PROLIFERATION"=brewer.pal(5,"Accent")[4],
                                                "UNANNOTATED"=brewer.pal(5,"Accent")[5]),
                                  Boyault_2007=c("G12"="pink","G3"="#999999","G56"="dodgerblue"),
                                  Lee_HB_2006=col_fun,Sia_et_al_2017=col_fun5,
                                  Roessler_met_2010=col_fun2,Yamahita_Epcam_2008=col_fun3,Yamahita_Stem_2008=col_fun4))
  
  hm=Heatmap(other_studies %>% dplyr::select(sample,subtype) %>% tibble::column_to_rownames("sample") %>% t(),
             col = col,name="Subtype",
             show_column_dend = F,show_column_names = F,bottom_annotation = ha,heatmap_height = unit(6,"cm"))
  
  
  plt=draw(hm,annotation_legend_side="bottom")
  
  return(list(plot= plt,out=other_studies))
}


## PCA plot function ####
pca <- function(countdata,A,B=NULL,title="PCA",nameA="subtypeA",nameB="subtypeB",colpal=brewer.pal(5,"Set1"),scale=T){
  library(ggplot2)
  library(ggfortify)
  library(RColorBrewer)
  pr <- prcomp(t(countdata),center = T,scale. =scale)
  shapes=c(21,24,19,18,17,11)
  if(is.null(B)){
    anno <- data.frame(A)
    colnames(anno)[1] <- nameA
    autoplot(pr,data=anno,colour=nameA,main=title,
             size=4,label.size=6,alpha=0.7)+
      scale_color_manual(values=unname(colpal)) + 
      scale_shape_manual(values=shapes)+theme_bw()+
      theme(panel.background = element_blank(),
            panel.grid = element_blank())
  }
  else{
    anno <- data.frame(as.character(A),as.character(B))
    colnames(anno) <- c(nameA,nameB)
    autoplot(pr,data=anno,colour=nameA,shape=nameB,main=title,
             size=4,label.size=6,alpha=0.7)+
      scale_color_manual(values=unname(colpal)) +
      scale_shape_manual(values=shapes)+theme_bw()+
      theme(panel.background = element_blank(),
            panel.grid = element_blank())
  }
}

## plot Submap ####
plot_submap <- function(submapMatrix,col.title="",row.title=""){
  ComplexHeatmap::Heatmap(submapMatrix,col=colorRampPalette(c("#d73027","#4575b4"))(10),name = "p-value",
                          column_title = col.title,
                          column_title_side="bottom",
                          row_title=row.title,
                          row_title_side="right",
                          cluster_rows =F,cluster_columns = F,
                          heatmap_legend_param = list(color_bar = "continuous"),
                          show_heatmap_legend=F,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            #if(submapRes[i,j]<0.1) 
                            grid::grid.text(paste0("p=",signif(submapMatrix[i, j],2)), 
                                            x, y, gp = grid::gpar(fontsize = 12))
                          })
}


## Calculate p-value taken from cjbmisc https://github.com/brightchan/cjbmisc/blob/master/R/getPvalue.R
{
  
# Only return p values:
p_fish.chi.t <- function(df,v1,v2,alt="two.sided",p.test="chi",ws=2e6){
  df <- as.data.frame(df)
  if(p.test=="both"){
    tryCatch(fish.t(df,v1,v2,alt,ws)$p.value,
             error=function(e){
               warning(e)
               warning("Using ChiSQ test instead.")
               chisq.t(df,v1,v2)$p.value})
  }else if (p.test=="fisher"){
    fish.t(df,v1,v2,alt,ws)$p.value
  }else if (p.test=="chi"){
    chisq.t(df,v1,v2)$p.value
  }else stop('p.test must be "fisher" or "chi"')
}
#' @rdname getPvalue
#' @export
p_aov.t <- function(df,v1,v2){
  df <- as.data.frame(df)
  if (plyr::is.discrete(df[,v1,drop=T]))
    f <- as.formula(paste0(v2,"~",v1)) else
      f <- as.formula(paste0(v1,"~",v2))
    aov <- aov(f, df)
    summary(aov)[[1]][["Pr(>F)"]][1]
}
#' @rdname getPvalue
#' @export
p_ContCont <- function(df,v1,v2,method="spearman"){
  # method to choose:
  # spearman,pearson,kendall,lm
  df <- as.data.frame(df)
  
  if (method=="lm"){
    return(p_lm(df,v1,v2))
  }else{
    return(cor.test(df[,v1],df[,v2],method = method)$p.value)
  }
}
#' @rdname getPvalue
#' @export
p_ContDisc <- function(df,v1,v2,method="kruskal.test"){
  # remove NAs to find out how many groups
  df <- df[!is.na(df[[v1]]) & !is.na(df[[v2]]), ]
  
  if (plyr::is.discrete(df[[v1]])) {
    df[[v1]] <- as.factor(df[[v1]])
    if (length(unique(na.omit(df[[v1]]))) == 1) {
      warning(paste0("Only one group found for ", v1,
                     ", returning NA for ", method))
      return(NA)
    }
    f <- as.formula(paste0(v2, "~", v1))
  }
  else {
    df[[v2]] <- as.factor(df[[v2]])
    if (length(unique(na.omit(df[[v2]]))) == 1) {
      warning(paste0("Only one group found for ", v1,
                     ", returning NA for ", method))
      return(NA)
    }
    f <- as.formula(paste0(v1, "~", v2))
  }
  get(method)(f, df)$p.value
}
#' @rdname getPvalue
#' @export
p_lm <- function(df,v1,v2){
  x <- df[[v1]]
  y <- df[[v2]]
  pvalue <- coef(summary(lm(y~x)))[2,4]
}
#' @rdname getPvalue
#' @export
p_xVsAll <- function(df,x.coln,y.coln=NULL,
                     num.num.test="spearman",
                     cat.num.test="kruskal.test",
                     cat.cat.test="chi"){
  
  df <- as.data.frame(df)
  # use all columns if y.coln not provided
  if(is.null(y.coln)) y.coln <- setdiff(colnames(df),x.coln)
  
  df <- df[,c(x.coln,y.coln)]
  
  feat.cate <- colnames(df)[sapply(df,function(x)!is.numeric(x))]
  
  # when x is categorical:
  
  if(x.coln %in% feat.cate){
    vec.p <- sapply(y.coln,function(y){
      if(y %in% feat.cate){
        return(p_fish.chi.t(df,x.coln,y,p.test=cat.cat.test))
      }else{
        return(p_ContDisc(df,x.coln,y,cat.num.test))
      }
    })
    
  }else{
    vec.p <- sapply(y.coln,function(y){
      if(y %in% feat.cate){
        return(p_ContDisc(df,x.coln,y,cat.num.test))
      }else{
        return(p_ContCont(df,x.coln,y,num.num.test))
      }
    })
  }
  
  vec.p <- setNames(vec.p,y.coln)
  
  return(vec.p)
}

#' @rdname getPvalue
#' @export
p_dfAll <- function(df,x.coln,y.coln=NULL,...){
  if(is.null(y.coln)) y.coln <- setdiff(colnames(df),x.coln)
  
  lst.pval <- lapply(setNames(x.coln,x.coln),function(x)
    p_xVsAll(df,x,y.coln,...) %>% as.data.frame())
  
  
  df.pval <- lapply(names(lst.pval), function(n) {
    colnames(lst.pval[[n]]) <- n
    lst.pval[[n]]$Features <- row.names(lst.pval[[n]])
    return(lst.pval[[n]])
  }) %>% Reduce(full_join, .) %>% select(Features, everything())
  
  return(df.pval)
}

#' @rdname getPvalue
#' @export
p_uniCox <- function(survdf,
                     feat,
                     surv.time = "PFS.days",
                     surv.status = "PFS.status",
                     signif.cutoff=0.05,
                     keep.all.in.barplot=F,
                     plot.surv=T,
                     survp.xlim=c(0,1400),
                     survp.timebreak= 365){
  feat.missing <- setdiff(feat,colnames(survdf))
  if(length(feat.missing)!=0) message(paste0(feat.missing,collapse = ", "),
                                      " are not in the provided dataframe!")
  feat <- intersect(feat,colnames(survdf))
  
  # fit each feat into a Cox model and return the minimal pvalue.
  lst.fit <- lapply(feat %>% setNames(.,.),
                    function(x){
                      message(x," ",appendLF =F)
                      if(length(unique(na.omit(pull(survdf,!!x))))<2) return(NA)
                      fit <- coxph(as.formula(paste0("Surv(",surv.time,",",surv.status,")~`",x,"`")),
                                   data = survdf)
                    })
  pval <- sapply(lst.fit,function(x)if(!is.na(x))min(summary(x)$coefficients[,5]) else return(NA))
  names(pval) <- gsub("\\.pvalue","",names(pval))
  
  # Plot the pvalues
  
  if(keep.all.in.barplot){
    par(mar=c(max(nchar(feat))/3+0.5,5,1,1))
    barplot(-log(sort(pval)),
            las=2,cex.names =0.8,ylab="-log(Pvalue)",
            main="Pvalues of Univariate Cox")
    abline(h=-log(signif.cutoff),lty=2)
  }else if(sum(pval<signif.cutoff,na.rm=T)>0){
    par(mar=c(1+(max(nchar(feat))/2.8),5,1,1))
    barplot(-log(sort(pval[pval<signif.cutoff])),
            las=2,cex.names =0.8,ylab="-log(Pvalue)",
            main="Pvalues of Univariate Cox")
    abline(h=-log(signif.cutoff),lty=2)
    
  }
  
  ## plot the Kaplan-Meier plots
  if(plot.surv & (sum(pval<signif.cutoff,na.rm=T)>0) ){
    lst.plot <- lapply(names(sort(pval[pval<signif.cutoff])) %>%
                         setNames(.,.),function(x){
                           p <- cjbmisc::ggsurvpWrap(survdf %>% as.data.frame(),
                                                     x,surv.time,surv.status,
                                                     risk.table = T,xlim = survp.xlim,
                                                     break.time.by = survp.timebreak)$plot+
                             theme(legend.key.size=unit(0.2,"cm"))
                           return(p)
                         }
    )
    
    p <-ggpubr::ggarrange(plotlist = lst.plot)
    plot(p)
    
    return(list(fit=lst.fit,pval=pval,plot=p))
  }else
    
    return(list(fit=lst.fit,pval=pval))
  
}


#' @rdname getPvalue
#' @export
p_adjust_mat <- function(pvaldf, p.adjust.method = "BH") {
  
  pvec <- as.vector(pvaldf)
  nai <- is.na(pvec)
  qvec <- rep(NA, length(pvec))
  qvec[!nai] <- p.adjust(pvec[!nai], method = p.adjust.method)
  qmat <- matrix(qvec, nrow = nrow(pvaldf))
  dimnames(qmat) <- dimnames(pvaldf)
  qmat
  
}
}
## calculate p value for combinatioon of features

# Fisher's exact test with columns in dataframe
fish.t <- function(df,v1,v2,alt="two.sided",ws=2e7){#v1,v2:name of col
  f <- as.formula(paste0("~",v1,"+",v2))
  cont.table <- xtabs(f,data=df)
  fisher.test(cont.table,alternative=alt,workspace=ws)
}
#' @rdname stattests
#' @export
# Chi Square Test with columns in dataframe
chisq.t <- function(df,v1,v2){#v1,v2:name of col
  f <- as.formula(paste0("~",v1,"+",v2))
  cont.table <- xtabs(f,data=df)
  chisq.test(cont.table)
}
#' @rdname stattests
#' @export
# One versus all test of contingency table
one.v.all.t <- function(df,test=fisher.test,...){
  sapply(setNames(colnames(df),colnames(df)),
         function(s){
           tmp <- data.frame(df[,s],rowSums(df)-df[,s])
           pvalue <- test(tmp,...)$p.value
           if (pvalue<0.01) message(s," p=",pvalue,"**")
           else if (pvalue<0.05) message(s," p=",pvalue,"*")
           else message(s," p=",pvalue)
         })
}
calculate_p_value=function(df,x,y){
  #library(cjbmisc)
  if (class(df[,x]) %in% c("factor","character") && class(df[,y])%in% c("factor","character")){
    message("using fisher's exact test")
    p=p_fish.chi.t(df,x,y)
  }
  else if (class(df[,x]) %in% c("integer","numeric","matrix") && class(df[,y])%in% c("integer","numeric","matrix")){
    message("using linear model p_value")
    p=p_lm(df,x,y)[1]
  }
  else if (class(df[,x]) %in% c("integer","numeric") && class(df[,y])%in% c("factor","character")){
    message("using wilcoxon test")
    p=p_ContDisc(df,x,y,method="kruskal.test") ###for more than 2 level data
  }
  else if (class(df[,x]) %in% c("factor","character") && class(df[,y])%in% c("integer","numeric")){
    message("using wilcoxon test")
    p=p_ContDisc(df,x,y,method="kruskal.test") 
  }
  return(p)
}

## plotting 
stacked_bar <- function(df,x,y,title,col=ggthemes::tableau_color_pal("Classic 10 Medium")(10),ldir=2){
  library(tidyverse)
  library(ggthemes)
  x=quo_name(x)
  y=quo_name(y)
  plt <- df[complete.cases(df[,as.character(c(x,y))]),] %>%
    mutate_at(c(y),as.factor) %>% 
    group_by_at(vars(c(x,y))) %>%
    summarise(n=n())%>%
    mutate(percent = (n / sum(n)), cumsum = cumsum(percent), 
           label=paste0("n=", sum(n)))
  fp=fisher.test(table(df[complete.cases(df[,as.character(c(x,y))]),][,x],df[complete.cases(df[,as.character(c(x,y))]),][,y]))$p.value
  ggplot(plt,aes_string(x=x, y="percent", fill=y)) +
    scale_y_continuous(labels = scales::percent) +
    labs ( y = "Percentage") +
    geom_bar(position = 'fill',color = "black", stat="identity") +
    geom_text(aes(y=-0.01, label=label),vjust=1,check_overlap = T) +
    guides(fill = guide_legend(title = NULL, reverse = F,nrow=ldir))+
    scale_fill_manual(values = col )+
    theme_bw()+theme(panel.grid = element_blank(),
                     legend.position = "bottom",
                     legend.text = element_text(color = "black"),
                     legend.key.size = unit(3,"mm"),
                     plot.title = element_text(hjust=0.5),
                     axis.text = element_text(colour="black"),
                     axis.ticks.x =element_blank())+
    labs(title = title, subtitle =paste0("p-value = ",formatC(fp,format="e",digits = 2)) ,fill="")+ylab(NULL)+xlab(NULL)
  
}
DESeq2_DEG <- function(Count.a,Count.b,name.a,name.b,test=c("Wald")){
  library(DESeq2)
  library(EdgeR)
  count = cbind(Count.a,Count.b)
  isexpr <-  rowMeans(cpm(count))>2
  message(paste0("Removing ",table(isexpr)[1]," low expresed genes"))
  count <- count[isexpr,]
  samples=c(colnames(Count.a),colnames(Count.b))
  coldata = c(rep(name.a,ncol(Count.a)),rep(name.b,ncol(Count.b))) %>% as.data.frame() %>% dplyr::rename(condition=1)
  rownames(coldata)=samples
  # Normalize data using DeSeq2
  dds = DESeqDataSetFromMatrix(countData=round(count), colData=coldata, design=~condition) 
  dds = DESeq(dds,test = test)
  #to get A vs B, define in the result function
  res = DESeq2::results(dds,contrast = c("condition",name.a,name.b))
  return(res)
}

