############################################################
## ML-02. IQ-TREE(ML) 트리 기반 수면 표현형 연관 분석
##   연속형 : K=2:12 -> ANOVA -> max F -> BH        (v1 과 동일)
##   범주형 : 기존 optimal-K 의 K 고정, Fisher 2xK + 순열
##            순열 귀무모형 = 원본 그대로
##              clp <- u[sample.int(K, n, replace = TRUE)]
##            NSIM 은 환경변수로 지정 (20000 또는 1000000)
##   추가   : 동률(equal-best) MP 트리 간 연관 안정성
############################################################
.libPaths(c("/disk2/bijsy/Test/Evolution_pressure/0.Code/Rlib", .libPaths()))
suppressMessages({library(ape); library(data.table); library(parallel)})
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

ROOT <- "/disk2/bijsy/Test/Sleep_association_Test"
V2   <- file.path(ROOT, "Result/Maximum_Parsimony_for_Cladistics/ML_based")
TR   <- file.path(V2, "trees")
EQ   <- file.path(V2, "trees_equalbest")
TAB  <- file.path(V2, "tables")
META <- file.path(ROOT, "METADATA_fixedVersion.txt")
N.CORES <- as.integer(Sys.getenv("N_CORES","64"))
NSIM    <- as.integer(Sys.getenv("NSIM","20000"))
SUF     <- if (NSIM >= 1e6) "_perm1M" else "_perm20k"
EQ.MAX  <- as.integer(Sys.getenv("EQ_MAX","20"))   # 안정성 검사에 쓸 동률 트리 수
SEED <- 319L; FISHER2K.MAXK <- 20L
K.MIN <- 2L; K.MAX <- 12L
cat(sprintf("NSIM = %s | 접미사 %s | 코어 %d\n", format(NSIM, big.mark=","), SUF, N.CORES))

meta0 <- as.data.frame(fread(META))
meta0$sp <- gsub("[ ']","_", meta0$Species_name_ensembl)

rd_tree <- function(f) read.tree(text=gsub("'","_", paste(readLines(f, warn=FALSE), collapse="")))
hc_of <- function(tr) hclust(as.dist(cophenetic.phylo(tr)), method="average")
tree_file <- function(tag,g) file.path(TR, tag, sprintf("%s_ML.nwk", g))
genes_of <- function(tag) sub("_ML\\.nwk$", "", list.files(file.path(TR, tag), "_ML\\.nwk$"))

############################################################
## A. 연속형
############################################################
run_cont <- function(tag, col) {
  m <- meta0[!is.na(suppressWarnings(as.numeric(meta0[[col]]))), ]
  gs <- genes_of(tag)
  one <- function(g) tryCatch({
    hc <- hc_of(rd_tree(tree_file(tag,g)))
    y <- as.numeric(m[[col]])[match(hc$labels, m$sp)]; names(y) <- hc$labels
    y <- y[is.finite(y)]
    if (length(y) < 4) stop("종<4")
    kmax <- min(K.MAX, length(y)-1L); if (kmax < K.MIN) stop("K없음")
    rbindlist(lapply(K.MIN:kmax, function(k) {
      cl <- cutree(hc, k=k)[names(y)]
      if (length(unique(cl)) < 2) return(NULL)
      fit <- aov(y ~ as.factor(cl)); s <- summary(fit)[[1]]
      data.table(Gene=g, Cluster=k, F_value=as.numeric(s$`F value`[1]),
                 P.anova=as.numeric(s$`Pr(>F)`[1]),
                 Shapiro=tryCatch(shapiro.test(residuals(fit))$p.value, error=function(e) NA_real_),
                 N_species=length(y))
    }))
  }, error=function(e) data.table(Gene=g, Cluster=NA_integer_, F_value=NA_real_,
       P.anova=NA_real_, Shapiro=NA_real_, N_species=NA_integer_))
  allK <- rbindlist(mclapply(gs, one, mc.cores=N.CORES), fill=TRUE)
  fwrite(allK, file.path(TAB, sprintf("ML_%s_allK.tsv", tag)), sep="\t")
  best <- allK[!is.na(F_value), .SD[which.max(F_value)], by=Gene]
  best[, adj.P := p.adjust(P.anova, "BH")]
  setorder(best, P.anova)
  fwrite(best, file.path(TAB, sprintf("ML_%s_optimalK.tsv", tag)), sep="\t")
  cat(sprintf("[%s] %d유전자 | P<0.05 %d | BH<0.05 %d\n", tag, nrow(best),
      sum(best$P.anova<0.05, na.rm=TRUE), sum(best$adj.P<0.05, na.rm=TRUE)))
}

############################################################
## B. 범주형
############################################################
mk_fish2k <- function() {
  cache <- new.env(hash=TRUE, parent=emptyenv())
  function(r1, r2) {
    key <- paste(r1, collapse=",")
    hit <- cache[[key]]; if (!is.null(hit)) return(hit)
    tb <- rbind(r1, r2); tb <- tb[, colSums(tb) > 0, drop=FALSE]
    p <- if (ncol(tb) < 2) NA_real_ else
      tryCatch(fisher.test(tb, workspace=2e8)$p.value, error=function(e)
        tryCatch(fisher.test(tb, simulate.p.value=TRUE, B=1e5)$p.value,
                 error=function(e2) NA_real_))
    cache[[key]] <- p; p
  }
}

run_cat <- function(tag, col, lev, label_map, bestk_file) {
  bk <- as.data.frame(fread(bestk_file))
  bk$K <- suppressWarnings(as.integer(gsub("[^0-9]","", as.character(bk$Cluster))))
  bk <- bk[!is.na(bk$K) & bk$K >= 2, ]; bk <- bk[!duplicated(bk$Gene), ]
  gs <- genes_of(tag); bk <- bk[bk$Gene %in% gs, ]
  m <- meta0[!is.na(meta0[[col]]) & trimws(as.character(meta0[[col]]))!="", ]
  m$Sleep_label <- unname(label_map[trimws(as.character(m[[col]]))])
  m <- m[!is.na(m$Sleep_label), ]

  tab_of <- function(tr, K_req) {
    hc <- hc_of(tr); K_req <- min(K_req, length(hc$labels)-1L)
    cl <- cutree(hc, k=K_req)
    df <- merge(data.frame(idx=names(cl), clusterInfo=as.integer(cl), stringsAsFactors=FALSE),
                m[, c("sp","Sleep_label")], by.x="idx", by.y="sp")
    if (!nrow(df)) return(NULL)
    ylab <- as.integer(factor(df$Sleep_label, levels=lev))
    list(df=df, i1=which(ylab==1L), i2=which(ylab==2L),
         u=unique(df$clusterInfo), KMAX=max(df$clusterInfo), n=nrow(df), K_req=K_req)
  }

  one <- function(i) tryCatch({
    g <- bk$Gene[i]
    tr <- rd_tree(tree_file(tag,g))
    X <- tab_of(tr, bk$K[i]); if (is.null(X)) stop("자료 없음")
    if (X$n < 5 || length(X$u) < 2 || !length(X$i1) || !length(X$i2)) stop("종/클러스터 부족")
    f2 <- mk_fish2k()
    o1 <- tabulate(X$df$clusterInfo[X$i1], nbins=X$KMAX)
    o2 <- tabulate(X$df$clusterInfo[X$i2], nbins=X$KMAX)
    obs <- f2(o1, o2)
    K <- length(X$u); n <- X$n
    n2k <- if (X$KMAX <= FISHER2K.MAXK && is.finite(obs)) NSIM else 0L
    n_ex <- 0L; d <- 0L
    if (n2k > 0L) {
      set.seed(SEED)
      for (b in seq_len(n2k)) {
        clp <- X$u[sample.int(K, n, replace = TRUE)]     ## 원본 귀무모형
        p2 <- f2(tabulate(clp[X$i1], nbins=X$KMAX), tabulate(clp[X$i2], nbins=X$KMAX))
        if (!is.na(p2)) { d <- d + 1L; if (p2 <= obs) n_ex <- n_ex + 1L }
      }
    }
    ## 동률 트리 간 관측 p 안정성
    fe <- ""
    eq_p <- NA_real_; eq_min <- NA_real_; eq_max <- NA_real_; eq_n <- 1L; eq_sig <- NA_real_
    if (file.exists(fe)) {
      tt <- tryCatch(read.tree(text=gsub("'","_", paste(readLines(fe,warn=FALSE), collapse="\n"))),
                     error=function(e) NULL)
      if (!is.null(tt)) {
        if (!inherits(tt,"multiPhylo")) tt <- structure(list(tt), class="multiPhylo")
        tt <- tt[seq_len(min(length(tt), EQ.MAX))]
        ps <- vapply(tt, function(t1) {
          Y <- tab_of(t1, bk$K[i]); if (is.null(Y)) return(NA_real_)
          f2(tabulate(Y$df$clusterInfo[Y$i1], nbins=Y$KMAX),
             tabulate(Y$df$clusterInfo[Y$i2], nbins=Y$KMAX)) }, numeric(1))
        ps <- ps[is.finite(ps)]
        if (length(ps)) { eq_n <- length(ps); eq_p <- median(ps)
                          eq_min <- min(ps); eq_max <- max(ps)
                          eq_sig <- mean(ps < 0.05) }
      }
    }
    data.table(Gene=g, Cluster=X$K_req, K_used=X$KMAX, N_species=n,
      N_grp1=length(X$i1), N_grp2=length(X$i2),
      P.fisher.2xK.observed=obs,
      P.fisher.2xK.perm=if (d>0) (n_ex+1)/(d+1) else NA_real_,
      N.extreme=n_ex, N.perm.valid=d, N_permutations=NSIM,
      EqualBest_n_tested=eq_n, EqualBest_obs_median=eq_p,
      EqualBest_obs_min=eq_min, EqualBest_obs_max=eq_max,
      EqualBest_frac_sig=eq_sig)
  }, error=function(e) data.table(Gene=bk$Gene[i], Cluster=NA_integer_, K_used=NA_integer_,
      N_species=NA_integer_, N_grp1=NA_integer_, N_grp2=NA_integer_,
      P.fisher.2xK.observed=NA_real_, P.fisher.2xK.perm=NA_real_,
      N.extreme=NA_integer_, N.perm.valid=NA_integer_, N_permutations=NSIM,
      EqualBest_n_tested=NA_integer_, EqualBest_obs_median=NA_real_,
      EqualBest_obs_min=NA_real_, EqualBest_obs_max=NA_real_,
      EqualBest_frac_sig=NA_real_, Error=conditionMessage(e)))

  R <- rbindlist(mclapply(seq_len(nrow(bk)), one, mc.cores=N.CORES), fill=TRUE)
  R[, adj.P := p.adjust(P.fisher.2xK.perm, "BH")]
  setorder(R, P.fisher.2xK.perm)
  fwrite(R, file.path(TAB, sprintf("ML_%s_Fisher%s.tsv", tag, SUF)), sep="\t")
  cat(sprintf("[%s] %d유전자 | perm P<0.05 %d | BH<0.05 %d | 동률트리 불안정(부호바뀜) %d\n",
      tag, nrow(R), sum(R$P.fisher.2xK.perm<0.05, na.rm=TRUE),
      sum(R$adj.P<0.05, na.rm=TRUE),
      sum(R$EqualBest_frac_sig > 0 & R$EqualBest_frac_sig < 1, na.rm=TRUE)))
}

if (Sys.getenv("SKIP_CONT","0") != "1") {
  run_cont("Sleep_duration", "Total_sleep_time_per_day")
  run_cont("NREM_ratio",     "Percentage_of_NREM_time_per_day")
}
run_cat("Sleep_timing", "_Sleep_timing_per_day", c("less_sleep","more_sleep"),
        c("Sleep at daytime"="more_sleep","Sleep at night"="less_sleep",
          "Sleep_at_daytime"="more_sleep","Sleep_at_night"="less_sleep",
          "nocturnal"="more_sleep","diurnal"="less_sleep"),
        file.path(ROOT,"Result/Sleeptiming/Sleep_Timing_optimal_K.tsv"))
run_cat("Sleep_frequency", "Number_of_sleep_times_per_day", c("More than twice","Once"),
        c("More than twice"="More than twice","Once"="Once",
          "more than twice"="More than twice","once"="Once"),
        file.path(ROOT,"Result/Sleepfrequency/Number_of_sleep_optimal_K.tsv"))
cat("[DONE]\n")
