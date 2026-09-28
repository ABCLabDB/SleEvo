############################################################
## 3기준(NJ 거리법 · MP 절약법 · ML 우도법) 재현성 비교
############################################################
.libPaths(c("/disk2/bijsy/Test/Evolution_pressure/0.Code/Rlib", .libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(ggrepel); library(patchwork)})

ROOT <- "/disk2/bijsy/Test/Sleep_association_Test"
MP   <- file.path(ROOT, "Result/Maximum_Parsimony_for_Cladistics")
OUT  <- file.path(MP, "3method_comparison")
TAB  <- file.path(OUT,"tables"); FIG <- file.path(OUT,"figures")
dir.create(TAB, recursive=TRUE, showWarnings=FALSE); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

SPEC <- list(
 list(lab="NREM ratio", ord=1, crit="BH",
   nj=file.path(ROOT,"Result/NREM_optimalK_ANOVA_selectedP_1M_FAST.tsv"), njp="P.anova.nominal", njq="adj.P.anova.nominal",
   mp=file.path(MP,"v2/tables/MP_v2_NREM_ratio_optimalK.tsv"), mpp="P.anova", mpq="adj.P",
   ml=file.path(MP,"ML_based/tables/ML_NREM_ratio_optimalK.tsv"), mlp="P.anova", mlq="adj.P"),
 list(lab="Sleep duration", ord=2, crit="BH",
   nj=file.path(ROOT,"Result/99.Previous/TST_key_12_optimal_Kruskal.tsv"), njp="P.anova", njq="adj.P",
   mp=file.path(MP,"v2/tables/MP_v2_Sleep_duration_optimalK.tsv"), mpp="P.anova", mpq="adj.P",
   ml=file.path(MP,"ML_based/tables/ML_Sleep_duration_optimalK.tsv"), mlp="P.anova", mlq="adj.P"),
 list(lab="Sleep timing", ord=3, crit="raw",
   nj=file.path(ROOT,"Result/Sleeptiming/Sleeptiming_Fisher_Result.tsv"), njp="P_Value", njq=NA,
   mp=file.path(MP,"v2/tables/MP_v2_Sleep_timing_Fisher_perm20k.tsv"), mpp="P.fisher.2xK.perm", mpq="adj.P",
   ml=file.path(MP,"ML_based/tables/ML_Sleep_timing_Fisher_perm20k.tsv"), mlp="P.fisher.2xK.perm", mlq="adj.P"),
 list(lab="Sleep frequency", ord=4, crit="raw",
   nj=file.path(ROOT,"Result/Sleepfrequency/Sleep_frequency_Fisher.tsv"), njp="P_Value", njq=NA,
   mp=file.path(MP,"v2/tables/MP_v2_Sleep_frequency_Fisher_perm20k.tsv"), mpp="P.fisher.2xK.perm", mpq="adj.P",
   ml=file.path(MP,"ML_based/tables/ML_Sleep_frequency_Fisher_perm20k.tsv"), mlp="P.fisher.2xK.perm", mlq="adj.P"))

grab <- function(f,pc,qc) {
  d <- as.data.frame(fread(f)); d <- d[!duplicated(d$Gene),]
  data.table(Gene=d$Gene, P=suppressWarnings(as.numeric(d[[pc]])),
             Q=if (is.na(qc)) NA_real_ else suppressWarnings(as.numeric(d[[qc]])))
}
ALL <- rbindlist(lapply(SPEC, function(s){
  a<-grab(s$nj,s$njp,s$njq); b<-grab(s$mp,s$mpp,s$mpq); c_<-grab(s$ml,s$mlp,s$mlq)
  setnames(a,c("Gene","P_NJ","Q_NJ")); setnames(b,c("Gene","P_MP","Q_MP")); setnames(c_,c("Gene","P_ML","Q_ML"))
  m <- merge(merge(a,b,by="Gene"), c_, by="Gene")
  m[, `:=`(Phenotype=s$lab, ord=s$ord, Criterion=s$crit)]
  sig <- function(P,Q) { v <- if (s$crit=="BH") Q else P; v <- as.numeric(v); !is.na(v) & v < .05 }
  m[, `:=`(Sig_NJ=sig(P_NJ,Q_NJ), Sig_MP=sig(P_MP,Q_MP), Sig_ML=sig(P_ML,Q_ML))]
  m[, N_methods_significant := as.integer(Sig_NJ)+as.integer(Sig_MP)+as.integer(Sig_ML)]
  m[, Support := c("0/3","1/3","2/3","3/3")[N_methods_significant+1]]
  m
}), fill=TRUE)
ALL$Phenotype <- factor(ALL$Phenotype, levels=sapply(SPEC[order(sapply(SPEC,`[[`,"ord"))],`[[`,"lab"))
setcolorder(ALL, c("Phenotype","Gene","Criterion","P_NJ","Q_NJ","Sig_NJ",
                   "P_MP","Q_MP","Sig_MP","P_ML","Q_ML","Sig_ML",
                   "N_methods_significant","Support"))
ALL[, ord:=NULL]
fwrite(ALL, file.path(TAB,"SupplementaryTable_three_method_comparison.tsv"), sep="\t")

SUM <- ALL[, .(N_genes=.N, Criterion=Criterion[1],
  NJ=sum(Sig_NJ), MP=sum(Sig_MP), ML=sum(Sig_ML),
  All_three=sum(N_methods_significant==3L),
  Two_of_three=sum(N_methods_significant==2L),
  One_only=sum(N_methods_significant==1L),
  None=sum(N_methods_significant==0L),
  NJ_MP=sum(Sig_NJ&Sig_MP), NJ_ML=sum(Sig_NJ&Sig_ML), MP_ML=sum(Sig_MP&Sig_ML),
  Union=sum(N_methods_significant>0L),
  Min_possible_overlap=max(0L, sum(Sig_NJ)+sum(Sig_MP)-.N),
  Spearman_NJ_MP=round(cor(-log10(pmax(P_NJ,1e-300)),-log10(pmax(P_MP,1e-300)),method="spearman",use="complete.obs"),3),
  Spearman_NJ_ML=round(cor(-log10(pmax(P_NJ,1e-300)),-log10(pmax(P_ML,1e-300)),method="spearman",use="complete.obs"),3),
  Spearman_MP_ML=round(cor(-log10(pmax(P_MP,1e-300)),-log10(pmax(P_ML,1e-300)),method="spearman",use="complete.obs"),3)
), by=Phenotype]
fwrite(SUM, file.path(TAB,"SupplementaryTable_three_method_summary.tsv"), sep="\t")
print(SUM)

############################################################
## 피겨
############################################################
PCOL <- c(`NREM ratio`="#7FB3D9", `Sleep duration`="#8FC98F",
          `Sleep timing`="#F2A65A", `Sleep frequency`="#BFBFBF")
MCOL <- c(NJ="#5B8FA8", MP="#C77B4E", ML="#6B8E5A")
LAB_GENES <- c("PER1","PER2","PER3","CRY1","CRY2","CLOCK","ARNTL","ARNTL2","NFIL3",
               "HDAC1","ATF4","ATF5","PPARA","CAVIN3","MAGEL2","ADRB1","CPT1A","CARTPT")
nlog <- function(p) -log10(pmax(as.numeric(p),1e-300))

## ── Fig 1: 방법별 -log10(p) 분포 (원 논문 스타일) ──
L <- melt(ALL[, .(Phenotype,Gene,NJ=P_NJ,MP=P_MP,ML=P_ML)],
          id.vars=c("Phenotype","Gene"), variable.name="Method", value.name="P")
L <- L[is.finite(P)]; L[, y := nlog(P)]
Lc <- melt(ALL[, .(Phenotype,Gene,NJ=Sig_NJ,MP=Sig_MP,ML=Sig_ML)],
           id.vars=c("Phenotype","Gene"), variable.name="Method", value.name="Sig")
L <- merge(L, Lc, by=c("Phenotype","Gene","Method"))
cnt <- L[, .(n=sum(Sig), N=.N), by=.(Phenotype,Method)]
lab <- L[Gene %in% LAB_GENES]
set.seed(1)
f1 <- ggplot(L, aes(Method, y)) +
  geom_jitter(aes(fill=Phenotype), shape=21, colour="white", stroke=.15,
              size=1.7, width=.26, alpha=.8) +
  geom_boxplot(width=.30, fill=NA, colour="black", linewidth=.5, outlier.shape=NA) +
  geom_hline(yintercept=-log10(.05), colour="red", linetype="dashed", linewidth=.5) +
  geom_text_repel(data=lab, aes(label=Gene), size=1.9, fontface="italic",
                  max.overlaps=14, segment.size=.15, segment.colour="grey55", seed=1) +
  geom_text(data=cnt, aes(x=Method, y=Inf, label=sprintf("%d/%d", n, N)),
            vjust=1.5, size=2.5, fontface="bold", inherit.aes=FALSE) +
  facet_wrap(~Phenotype, nrow=1, scales="free_y") +
  scale_fill_manual(values=PCOL, guide="none") +
  labs(title="Gene-level association under three tree-reconstruction criteria",
       subtitle="NJ = distance, MP = parsimony (v2), ML = likelihood (IQ-TREE);  dashed line = p 0.05",
       x=NULL, y=expression(-log[10]*"(p value)")) +
  theme_classic(base_size=9) +
  theme(plot.title=element_text(face="bold",size=11),
        plot.subtitle=element_text(size=7.4,colour="grey35"),
        strip.background=element_rect(fill="grey94",colour=NA),
        strip.text=element_text(face="bold",size=8.4),
        axis.text.x=element_text(face="bold",size=8))

## ── Fig 2: 방법 간 p값 일치도 ──
PR <- rbindlist(list(
  ALL[, .(Phenotype,Gene,x=nlog(P_NJ),y=nlog(P_MP),pair="NJ vs MP")],
  ALL[, .(Phenotype,Gene,x=nlog(P_NJ),y=nlog(P_ML),pair="NJ vs ML")],
  ALL[, .(Phenotype,Gene,x=nlog(P_MP),y=nlog(P_ML),pair="MP vs ML")]))
PR$pair <- factor(PR$pair, levels=c("NJ vs MP","NJ vs ML","MP vs ML"))
f2 <- ggplot(PR, aes(x,y)) +
  geom_abline(slope=1,intercept=0,colour="grey70",linetype="22",linewidth=.35) +
  geom_hline(yintercept=-log10(.05),colour="red",linetype="dashed",linewidth=.3) +
  geom_vline(xintercept=-log10(.05),colour="red",linetype="dashed",linewidth=.3) +
  geom_point(aes(fill=Phenotype), shape=21, colour="white", stroke=.13, size=1.5, alpha=.8) +
  facet_grid(Phenotype ~ pair, scales="free") +
  scale_fill_manual(values=PCOL, guide="none") +
  labs(title="Pairwise agreement between reconstruction criteria",
       subtitle="grey = identity line; red = p 0.05",
       x=expression(-log[10]*"(p)"), y=expression(-log[10]*"(p)")) +
  theme_classic(base_size=8) +
  theme(plot.title=element_text(face="bold",size=11),
        plot.subtitle=element_text(size=7.2,colour="grey35"),
        strip.background=element_rect(fill="grey94",colour=NA),
        strip.text=element_text(face="bold",size=7.4))

## ── Fig 3: 재현 수준 ──
SP <- ALL[, .N, by=.(Phenotype, Support)]
SP$Support <- factor(SP$Support, levels=c("3/3","2/3","1/3","0/3"))
f3 <- ggplot(SP[Support!="0/3"], aes(Phenotype, N, fill=Support)) +
  geom_col(width=.62, colour="white", linewidth=.3) +
  geom_text(aes(label=N), position=position_stack(vjust=.5), size=2.6, colour="white", fontface="bold") +
  scale_fill_manual(values=c(`3/3`="#2C5F8A", `2/3`="#7FA8C4", `1/3`="#C9D6E0"),
                    name="Methods significant") +
  labs(title="Cross-criterion support", x=NULL, y="genes") +
  theme_classic(base_size=9) +
  theme(plot.title=element_text(face="bold",size=11),
        legend.position="right", legend.title=element_text(size=7.5),
        legend.text=element_text(size=7), legend.key.size=unit(9,"pt"),
        axis.text.x=element_text(face="bold",size=8))

## ── Fig 4: 지정 유전자 ──
HL <- list(`Sleep duration`=c("ATF4","ADRB1","MAGEL2"),
           `NREM ratio`=c("ARNTL2","PPARA","CPT1A"),
           `Sleep timing`=c("PER1","CRY2","CAVIN3","ATF5"),
           `Sleep frequency`=c("ATF5","HDAC1","NFIL3","CARTPT"))
HD <- rbindlist(lapply(names(HL), function(p) ALL[Phenotype==p & Gene %in% HL[[p]]]))
HM <- melt(HD[, .(Phenotype,Gene,NJ=P_NJ,MP=P_MP,ML=P_ML)],
           id.vars=c("Phenotype","Gene"), variable.name="Method", value.name="P")
HS <- melt(HD[, .(Phenotype,Gene,NJ=Sig_NJ,MP=Sig_MP,ML=Sig_ML)],
           id.vars=c("Phenotype","Gene"), variable.name="Method", value.name="Sig")
HM <- merge(HM,HS,by=c("Phenotype","Gene","Method"))
HM[, y := nlog(P)]
HM[, GeneL := factor(paste0(Gene,"\n(",Phenotype,")"))]
ordv <- HD[order(-N_methods_significant, Gene), paste0(Gene,"\n(",Phenotype,")")]
HM$GeneL <- factor(HM$GeneL, levels=rev(unique(ordv)))
f4 <- ggplot(HM, aes(Method, GeneL, fill=y)) +
  geom_tile(colour="white", linewidth=.6) +
  geom_point(data=HM[Sig==TRUE], shape=8, size=1.3, colour="white") +
  geom_text(aes(label=ifelse(P<1e-3, sprintf("%.0e",P), sprintf("%.3f",P))),
            size=2.1, colour=ifelse(HM$y>2.2,"white","grey15")) +
  scale_fill_gradientn(colours=c("#F7F9FA","#BBD8E4","#5C93AE","#1F4E68"),
                       name=expression(-log[10]*"(p)")) +
  labs(title="Candidate genes across the three criteria",
       subtitle="white asterisk = significant;  value = p", x=NULL, y=NULL) +
  theme_minimal(base_size=8.5) +
  theme(plot.title=element_text(face="bold",size=11),
        plot.subtitle=element_text(size=7.2,colour="grey35"),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=6.6, lineheight=.85),
        axis.text.x=element_text(face="bold",size=8),
        legend.key.width=unit(7,"pt"), legend.key.height=unit(16,"pt"),
        legend.title=element_text(size=6.6), legend.text=element_text(size=6.2))

save2 <- function(p,n,w,h) for (fm in c("jpg","pdf")) {
  f <- file.path(FIG, paste0(n,".",fm))
  if (fm=="pdf") ggsave(f,p,width=w,height=h,device=cairo_pdf,bg="white")
  else ggsave(f,p,width=w,height=h,dpi=500,device="jpeg",quality=97,bg="white")
  }
save2(f1,"Fig1_three_method_distribution", 11.5, 4.2)
save2(f2,"Fig2_pairwise_agreement",         8.5, 9.5)
save2(f3,"Fig3_cross_criterion_support",    6.4, 3.4)
save2(f4,"Fig4_candidate_genes",            5.2, 5.6)
save2((f1/(f3|f4)) + plot_annotation(tag_levels="A"), "Figure_three_method_summary", 11.5, 9.0)
cat("\n피겨 저장:", FIG, "\n[DONE]\n")
