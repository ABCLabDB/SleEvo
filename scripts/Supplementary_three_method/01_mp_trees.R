############################################################
## v2-01. 강화 MP 탐색 + 동률 트리 + strict consensus
##
##  v1 대비 변경
##    maxit 200 -> 1000, minit 20 -> 100  (k=10 유지)
##    all = TRUE            동률 최소 PS 트리 전부 회수
##    독립 탐색 N회          best score 재회수율 기록
##    acctran -> di2multi    zero-length 내부분지 접기
##    strict consensus       생성 + 해상도 기록
##    node-level bootstrap   prop.part 로 bipartition 매칭해 저장
##    N_variable / N_parsimony_informative / gap 기록
############################################################
.libPaths(c("/disk2/bijsy/Test/Evolution_pressure/0.Code/Rlib", .libPaths()))
suppressMessages({library(ape); library(phangorn); library(Biostrings)
                  library(data.table); library(parallel)})
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

ROOT <- "/disk2/bijsy/Test/Sleep_association_Test"
V2   <- file.path(ROOT, "Result/Maximum_Parsimony_for_Cladistics/v2")
DATA <- "/disk2/bijsy/Test/Muscle/"
META <- file.path(ROOT, "METADATA_fixedVersion.txt")
GENES<- file.path(ROOT, "supplementaryS2.tsv")
for (d in c("trees","trees_equalbest","consensus","tables","checkpoint"))
  dir.create(file.path(V2,d), recursive=TRUE, showWarnings=FALSE)

N.CORES <- as.integer(Sys.getenv("N_CORES","64"))
N.BOOT  <- as.integer(Sys.getenv("N_BOOT","2000"))
N.IND   <- as.integer(Sys.getenv("N_IND","5"))     # 독립 탐색 횟수
MAXIT   <- 1000L; MINIT <- 100L; KRAT <- 10L
MAX.KEEP<- 200L                                    # 저장할 동률 트리 상한
SEED    <- 319L

PH <- list(
  list(tag="Sleep_duration",  col="Total_sleep_time_per_day",        type="num"),
  list(tag="NREM_ratio",      col="Percentage_of_NREM_time_per_day", type="num"),
  list(tag="Sleep_timing",    col="_Sleep_timing_per_day",           type="chr"),
  list(tag="Sleep_frequency", col="Number_of_sleep_times_per_day",   type="chr"))

meta0 <- as.data.frame(fread(META))
meta0$sym <- gsub(" ","_",tolower(meta0$Species_symbol_name_ensembl))
norm_sp <- function(x) gsub("[ ']","_",x)

gl   <- as.data.frame(fread(GENES))
gene <- unique(gl$Gene_symbol)
gene <- gene[file.exists(file.path(DATA, paste0(gene,"_muscle.fasta")))]

sp_of <- function(s) {
  v <- meta0[[s$col]]
  keep <- if (s$type=="num") !is.na(suppressWarnings(as.numeric(v)))
          else (!is.na(v) & trimws(as.character(v))!="")
  meta0[keep, ]
}
for (s in PH) cat(sprintf("  %-16s %d종\n", s$tag, nrow(sp_of(s))))
cat(sprintf("유전자 %d | 트리 %d | 코어 %d | bootstrap %d | 독립탐색 %d | maxit %d minit %d\n\n",
            length(gene), length(gene)*length(PH), N.CORES, N.BOOT, N.IND, MAXIT, MINIT))

load_dat <- function(g, m) {
  cds <- readDNAStringSet(file.path(DATA, paste0(g,"_muscle.fasta")))
  ids <- vapply(strsplit(names(cds),":"), function(x) if(length(x)>=2) x[2] else x[1], character(1))
  names(cds) <- ids
  cds <- cds[names(cds) %in% m$sym]
  nm <- m$Species_name_ensembl[match(names(cds), m$sym)]
  cds <- cds[!is.na(nm)]; names(cds) <- norm_sp(nm[!is.na(nm)])
  if (length(cds) < 5) stop("종 < 5")
  dna <- as.DNAbin(cds)
  list(pd = as.phyDat(dna), dna = dna)
}

site_stats <- function(dna) {
  M <- toupper(as.character(as.matrix(dna))); base <- c("A","C","G","T")
  nv <- 0L; npi <- 0L
  for (j in seq_len(ncol(M))) {
    z <- M[,j]; z <- z[z %in% base]; tb <- table(z)
    if (length(tb) > 1L) nv <- nv + 1L
    if (sum(tb >= 2L) >= 2L) npi <- npi + 1L
  }
  list(len=ncol(M), nvar=nv, npi=npi, gap=mean(!(M %in% base)))
}

one <- function(task) {
  g <- task$g; tag <- task$tag; m <- task$m
  ck <- file.path(V2,"checkpoint", sprintf("%s__%s.rds", tag, g))
  if (file.exists(ck)) return(readRDS(ck))
  t0 <- Sys.time()
  res <- tryCatch({
    D <- load_dat(g, m); pd <- D$pd
    ss <- site_stats(D$dna)

    ## ---- 독립 탐색 N회: best score 재회수율 ----
    scores <- numeric(N.IND); pool <- list()
    for (i in seq_len(N.IND)) {
      set.seed(SEED + i)
      st <- if (i == 1L) NJ(dist.hamming(pd)) else rtree(length(pd), tip.label=names(pd))
      tr <- pratchet(pd, start=st, maxit=MAXIT, minit=MINIT, k=KRAT,
                     trace=0, all=TRUE, perturbation="ratchet")
      if (!inherits(tr,"multiPhylo")) tr <- structure(list(tr), class="multiPhylo")
      scores[i] <- parsimony(tr[[1]], pd)
      pool[[i]] <- tr
    }
    best <- min(scores)
    recov <- mean(scores == best)                       # best score 재회수율
    eq <- do.call(c, pool[scores == best])
    if (!inherits(eq,"multiPhylo")) eq <- structure(list(eq), class="multiPhylo")
    eq <- unique(eq)                                    # 동일 토폴로지 중복 제거
    n_eq <- length(eq)
    if (n_eq > MAX.KEEP) eq <- eq[seq_len(MAX.KEEP)]

    ## ---- strict consensus ----
    if (length(eq) > 1L) {
      cs  <- consensus(eq, p = 1)
      resn<- (Nnode(cs) - 1) / (length(cs$tip.label) - 2) * 100
      write.tree(cs, file.path(V2,"consensus", sprintf("%s__%s_strict.nwk", tag, g)))
    } else { cs <- eq[[1]]; resn <- 100 }

    ## ---- 대표 트리: acctran -> di2multi ----
    rep <- acctran(eq[[1]], pd)
    n_zero <- sum(rep$edge.length == 0)
    rep <- di2multi(rep)

    ps <- parsimony(rep, pd); ci <- CI(rep, pd); ri <- RI(rep, pd)

    ## ---- bootstrap ----
    set.seed(SEED)
    bsl <- bootstrap.phyDat(pd, FUN=function(x)
             pratchet(x, start=NJ(dist.hamming(x)), maxit=60, minit=10, k=5,
                      trace=0, all=FALSE), bs=N.BOOT, multicore=FALSE)
    sup <- prop.clades(rep, bsl, rooted=FALSE)/N.BOOT*100
    sup[is.na(sup)] <- 0
    rep$node.label <- round(sup)
    inode <- sup[-1]
    write.tree(rep, file.path(V2,"trees", sprintf("%s__%s_MP.nwk", tag, g)))
    if (length(eq) > 1L)
      write.tree(eq, file.path(V2,"trees_equalbest", sprintf("%s__%s_equalbest.nwk", tag, g)))

    ## ---- node-level bootstrap (bipartition 기준) ----
    pp <- prop.part(rep)
    lab <- attr(pp, "labels")
    keys <- vapply(pp, function(idx) {
      a <- sort(lab[idx]); b <- sort(setdiff(lab, lab[idx]))
      k1 <- paste(a, collapse="|"); k2 <- paste(b, collapse="|")
      if (k1 < k2) k1 else k2 }, character(1))
    nodes <- data.table(Phenotype=tag, Gene=g, Split_key=keys,
                        N_taxa_one_side=lengths(pp), MP_bootstrap=round(sup,1))

    list(status="OK", Phenotype=tag, Gene=g,
      N_tips=length(pd), Alignment_length=ss$len,
      N_variable_sites=ss$nvar, N_parsimony_informative=ss$npi,
      Gap_or_missing_proportion=round(ss$gap,4),
      MP_parsimony_score=ps, CI=ci, RI=ri, RC=ci*ri, HI=1-ci,
      Tree_length=sum(rep$edge.length),
      N_best_trees_recovered=n_eq, Best_score_recovery=round(recov,3),
      N_independent_searches=N.IND,
      Strict_consensus_resolution=round(resn,1),
      N_zero_length_edges_collapsed=n_zero,
      N_internal_nodes=length(inode),
      Bootstrap_mean=mean(inode), Bootstrap_median=median(inode),
      N_nodes_bs70=sum(inode>=70), Pct_nodes_bs70=mean(inode>=70)*100,
      N_nodes_bs95=sum(inode>=95), Pct_nodes_bs95=mean(inode>=95)*100,
      seconds=as.numeric(difftime(Sys.time(), t0, units="secs")),
      Error=NA_character_, nodes=nodes)
  }, error=function(e) list(status="ERROR", Phenotype=tag, Gene=g,
                            Error=conditionMessage(e)))
  saveRDS(res, ck)
  if (identical(res$status,"OK"))
    message(sprintf("[OK] %-16s %-12s tips=%2d PS=%6d eq=%3d recov=%.2f cons=%.0f%% bs70=%.0f%% (%.0fs)",
      tag, g, res$N_tips, res$MP_parsimony_score, res$N_best_trees_recovered,
      res$Best_score_recovery, res$Strict_consensus_resolution, res$Pct_nodes_bs70, res$seconds))
  else message(sprintf("[ERR] %-16s %-12s %s", tag, g, res$Error))
  res
}

tasks <- unlist(lapply(PH, function(s){ m <- sp_of(s)
  lapply(gene, function(g) list(g=g, tag=s$tag, m=m)) }), recursive=FALSE)

t0 <- Sys.time()
res <- mclapply(tasks, one, mc.cores=N.CORES, mc.preschedule=FALSE)
ok  <- Filter(function(x) identical(x$status,"OK"), res)
cat(sprintf("\n완료 %d / %d | %.1f분\n", length(ok), length(res),
            as.numeric(difftime(Sys.time(), t0, units="mins"))))

DF <- rbindlist(lapply(ok, function(x) as.data.frame(x[setdiff(names(x), c("status","nodes"))])), fill=TRUE)
fwrite(DF, file.path(V2,"tables","SupplementaryTable_MP_tree_statistics_v2.tsv"), sep="\t")
for (s in PH) fwrite(DF[Phenotype==s$tag],
  file.path(V2,"tables", sprintf("MP_tree_statistics_%s_v2.tsv", s$tag)), sep="\t")

ND <- rbindlist(lapply(ok, function(x) x$nodes), fill=TRUE)
fwrite(ND, file.path(V2,"tables","node_level_bootstrap_MP_v2.tsv"), sep="\t")

cat("\n=== 표현형별 요약 ===\n")
print(DF[, .(genes=.N, tips=median(N_tips), PS=round(median(MP_parsimony_score)),
             CI=round(median(CI),3), RI=round(median(RI),3),
             eq_trees=round(median(N_best_trees_recovered),1),
             recov=round(median(Best_score_recovery),2),
             cons_res=round(median(Strict_consensus_resolution),1),
             zero_edges=round(median(N_zero_length_edges_collapsed),1),
             bs70=round(median(Pct_nodes_bs70),1)), by=Phenotype])
cat(sprintf("\n노드별 부트스트랩: %d행\n", nrow(ND)))
cat("[DONE]\n")
