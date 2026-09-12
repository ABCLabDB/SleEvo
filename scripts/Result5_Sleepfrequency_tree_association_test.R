############################################################
## Result5_Sleepfrequency_tree_association_test.R
##
## Sleep-frequency tree / cluster association with permutation
## Fisher (2 x K). Output is slim:
##   Gene, Cluster, P_Value (= P.fisher.2xK.perm)
##
## Based on Sleepfrequency/Sleep_frequency_permutation.R
## Paths point at the SleEvo repository data layout.
############################################################

############################
## Paths
############################
DATA.DIR   <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Fasta"
META.FILE  <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Result1/species_sleep_metadata.txt"
BESTK.FILE <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Result5/Number_of_sleep_optimal_K.tsv"
OUT.FILE   <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Result5/Sleep_frequency_Fisher.tsv"
CHECK.DIR  <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Result5/checkpoint_literal"

############################
## Column names in METADATA
############################
COL.SPECIES.SYMBOL <- "Species_symbol_name_ensembl"
COL.SPECIES.NAME   <- "Species_name_ensembl"          # becomes idx
COL.PHENO          <- "Number_of_sleep_times_per_day"

## The original mapped "once" -> "Once", then took Sleep_label from
## "Once" and "More than twice" only; everything else became NA and the
## species was dropped. Spell that out here - add a line if this file
## uses another spelling.
LABEL.MAP <- c(
  "once"            = "Once",
  "Once"            = "Once",
  "More than twice" = "More than twice"
)

############################
## Column names in the optimal-K file
############################
COL.GENE          <- "Gene"
COL.BESTCLUSTER   <- "Cluster"              # renamed best.Cluster in the original
COL.PURITY        <- "purity_score"         # original.purity
COL.EST.COCHRAN   <- "Estimate.cochran"     # original.cochran
COL.EST.MERGE.COC <- "E.merge.cochran"      # original.merge.cochran
COL.EST.FISHER    <- "Estimate.fisher"      # original.Pvalue

############################
## Run
############################
N.CORES <- 24L
NSIM    <- 1000000L        # as in the original
SEED    <- 319L
RESUME  <- TRUE
ONLY.GENES <- NULL         # e.g. c("PER1") to test one gene

## Fisher on the RAW 2 x K table (no merging).
## The original only ran Fisher on the merged 2 x 2. A 2 x K exact test
## enumerates every table with the observed margins - measured at ~29 ms
## per call for K = 14, n = 45, with ~98% of drawn tables distinct, so
## memoisation does not help. A million of them is about 8 h per gene, so
## it runs on the FIRST FISHER2K.NSIM iterations of the SAME stream (no
## extra RNG is consumed, so every other column is unaffected).
## Set FISHER2K.NSIM <- 0L to skip it.
FISHER2K.NSIM <- 20000L
FISHER2K.MAXK <- 20L       # skip genes with a larger K

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
library(ape); library(data.table); library(Biostrings); library(parallel)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)
}
setDTthreads(1)
dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT.FILE), recursive = TRUE, showWarnings = FALSE)

must_have <- function(df, col, where) {
  if (!col %in% colnames(df))
    stop("column '", col, "' is not in ", where, ".\nColumns present:\n  ",
         paste(colnames(df), collapse = "\n  "),
         "\n\nSet the matching COL.* value at the top of the script.",
         call. = FALSE)
  invisible(TRUE)
}

############################################################
## Metadata
############################################################
meta <- as.data.frame(fread(META.FILE))
for (cc in c(COL.SPECIES.SYMBOL, COL.SPECIES.NAME, COL.PHENO))
  must_have(meta, cc, "the metadata file")

cat("raw", COL.PHENO, "values:\n")
print(table(meta[[COL.PHENO]], useNA = "ifany"))

meta$idx <- as.character(meta[[COL.SPECIES.NAME]])
meta$sym <- gsub(" ", "_", tolower(trimws(meta[[COL.SPECIES.SYMBOL]])))
meta$Sleep_label <- unname(LABEL.MAP[trimws(as.character(meta[[COL.PHENO]]))])

unmapped <- setdiff(unique(trimws(as.character(meta[[COL.PHENO]]))[
  is.na(meta$Sleep_label)]), c(NA, ""))
if (length(unmapped))
  cat("\nnot in LABEL.MAP, so dropped (as in the original):\n  ",
      paste(unmapped, collapse = " | "), "\n")

meta <- meta[!is.na(meta$Sleep_label) & !is.na(meta$idx), ]
meta <- meta[!duplicated(meta$idx), c("idx", "sym", "Sleep_label")]
rownames(meta) <- NULL
cat("\nSleep_label:\n"); print(table(meta$Sleep_label))
cat("species:", nrow(meta), "\n")
if (length(unique(meta$Sleep_label)) < 2) stop("only one Sleep_label present")

## table() sorts labels alphabetically: "More than twice" < "Once"
LEV <- c("More than twice", "Once")

############################################################
## optimal-K file
############################################################
bk <- as.data.frame(fread(BESTK.FILE))
cat("\noptimal-K file columns:\n  ", paste(colnames(bk), collapse = ", "), "\n")
for (cc in c(COL.GENE, COL.BESTCLUSTER))
  must_have(bk, cc, "the optimal-K file")
hasc  <- COL.EST.COCHRAN   %in% colnames(bk)
hasmc <- COL.EST.MERGE.COC %in% colnames(bk)
hasf  <- COL.EST.FISHER    %in% colnames(bk)
hasp  <- COL.PURITY        %in% colnames(bk)
for (nm in c(COL.EST.COCHRAN, COL.EST.MERGE.COC, COL.EST.FISHER, COL.PURITY))
  if (!nm %in% colnames(bk))
    cat("!! '", nm, "' absent - the p-value that compares against it will be NA\n",
        sep = "")

getcol <- function(nm, ok) if (ok) suppressWarnings(as.numeric(bk[[nm]])) else NA_real_
bestK <- data.frame(
  Gene = as.character(bk[[COL.GENE]]),
  K    = suppressWarnings(as.integer(gsub("[^0-9]", "",
                                          as.character(bk[[COL.BESTCLUSTER]])))),
  orig_cochran   = getcol(COL.EST.COCHRAN,   hasc),
  orig_mcochran  = getcol(COL.EST.MERGE.COC, hasmc),
  orig_fisher    = getcol(COL.EST.FISHER,    hasf),
  orig_purity    = getcol(COL.PURITY,        hasp),
  stringsAsFactors = FALSE)
bestK <- bestK[!is.na(bestK$K) & bestK$K >= 2, ]
bestK <- bestK[!duplicated(bestK$Gene), ]
bestK <- bestK[file.exists(file.path(DATA.DIR, paste0(bestK$Gene, "_muscle.fasta"))), ]
if (!is.null(ONLY.GENES)) bestK <- bestK[bestK$Gene %in% ONLY.GENES, ]
NG <- nrow(bestK)
if (NG == 0L) stop("no genes left")

N.CORES <- min(N.CORES, parallel::detectCores(logical = TRUE))
cat("\ngenes:", NG, "| workers:", N.CORES,
    "| nsim:", format(NSIM, big.mark = ","), "| seed:", SEED, "\n\n")


############################################################
## one gene
############################################################
analyze_gene <- function(g) {
  
  gene <- bestK$Gene[g]
  ck <- file.path(CHECK.DIR, paste0(gene, ".rds"))
  if (RESUME && file.exists(ck)) return(readRDS(ck))
  
  out <- tryCatch({
    t0 <- Sys.time()
    
    original.purity        <- bestK$orig_purity[g]
    original.cochran       <- bestK$orig_cochran[g]
    original.merge.cochran <- bestK$orig_mcochran[g]
    original.Pvalue        <- bestK$orig_fisher[g]
    
    ## ---- sequences ----
    cds <- readDNAStringSet(file.path(DATA.DIR, paste0(gene, "_muscle.fasta")))
    names(cds) <- gsub(" ", "_", tolower(trimws(
      vapply(strsplit(names(cds), ":"),
             function(x) if (length(x) >= 2L) x[2] else x[1], character(1)))))
    cds <- cds[names(cds) %in% meta$sym]
    sp  <- meta$idx[match(names(cds), meta$sym)]
    keep <- !is.na(sp); cds <- cds[keep]; sp <- sp[keep]
    dp <- duplicated(sp); if (any(dp)) { cds <- cds[!dp]; sp <- sp[!dp] }
    names(cds) <- sp
    if (length(cds) < 5) stop("Too few species.")
    
    ## ---- tree ----
    dm   <- dist.dna(as.DNAbin(cds), as.matrix = TRUE, pairwise.deletion = TRUE)
    tree <- njs(dm)
    dm   <- cophenetic.phylo(tree)
    average.dm <- hclust(as.dist(dm), method = "average")
    
    K_group  <- min(bestK$K[g], length(average.dm$labels) - 1L)
    clusters <- cutree(average.dm, K_group)
    clusters <- data.frame(idx = names(clusters),
                           clusterInfo = as.integer(clusters),
                           stringsAsFactors = FALSE)
    best.cluster.df <- merge(clusters, meta, by = "idx", all = FALSE, sort = FALSE)
    permu.df <- best.cluster.df[, c("idx", "clusterInfo", "Sleep_label")]
    
    n <- nrow(permu.df)
    K <- length(unique(permu.df$clusterInfo))   # what sample.int() draws from
    if (n < 5 || K < 2L) stop("too few species or clusters")
    
    ylab <- as.integer(factor(permu.df$Sleep_label, levels = LEV))
    i1 <- which(ylab == 1L)          # More than twice  (table row 1)
    i2 <- which(ylab == 2L)          # Once             (table row 2)
    n1 <- length(i1); n2 <- length(i2)
    if (n1 == 0L || n2 == 0L) stop("one phenotype missing")
    
    ## ---- Cochran-Armitage, DescTools closed form ----
    Rb2  <- (n1 + 2 * n2) / n
    s2_2 <- n1 * (1 - Rb2)^2 + n2 * (2 - Rb2)^2
    pp   <- n1 / n
    co_raw <- function(r1, r2) {
      tot <- r1 + r2
      np  <- sum(tot > 0)
      if (np == 2L) {                     # a 2-column table is not transposed
        j <- which(tot > 0)[1]
        p1 <- tot[j] / n; den <- p1 * (1 - p1) * s2_2
        if (den <= 0) return(NA_real_)
        return((r1[j] * (1 - Rb2) + r2[j] * (2 - Rb2)) / sqrt(den))
      }
      pres <- tot > 0
      S <- cumsum(pres) * pres
      Rbar <- sum(tot * S) / n
      s2 <- sum(tot * (S - Rbar)^2 * pres)
      den <- pp * (1 - pp) * s2
      if (den <= 0) return(NA_real_)
      sum(r1 * (S - Rbar)) / sqrt(den)
    }
    
    ## ---- the 2-column table the original fed to purity / merge / fisher
    ## if the permuted table already has 2 columns it is used as is,
    ## otherwise making.sum.cluster() merges it
    two_col <- function(r1, r2) {
      tot <- r1 + r2
      pres <- which(tot > 0)
      if (length(pres) == 2L) {
        j <- pres[1]; k <- pres[2]
        return(c(r1[j], r2[j], r1[k], r2[k]))     # (col1: r1,r2 ; col2: r1,r2)
      }
      tie <- (r1 == r2) & (tot > 0)
      ti  <- sum(r1[tie]) / 2
      ## bottom.res = Once = r2 ; top.res = More than twice = r1
      lessm <- r2 > r1        # "Once" dominates  -> Cluster.less
      morem <- r1 > r2        # "More than twice" -> Cluster.more
      c(sum(r1[lessm]) + ti, sum(r2[lessm]) + ti,
        sum(r1[morem]) + ti, sum(r2[morem]) + ti)
    }
    purity_of <- function(m) {
      s1 <- m[1] + m[2]; s2 <- m[3] + m[4]
      (if (s1 > 0) max(m[1], m[2]) / s1 else NaN) / 2 +
        (if (s2 > 0) max(m[3], m[4]) / s2 else NaN) / 2
    }
    ## CochranArmitageTest(table + 1) on the 2x2
    co_merge <- function(m) {
      a <- m + 1
      R1 <- a[1] + a[3]; R2 <- a[2] + a[4]; N <- R1 + R2
      Rbar <- (R1 + 2 * R2) / N
      s2 <- R1 * (1 - Rbar)^2 + R2 * (2 - Rbar)^2
      p1 <- (a[1] + a[2]) / N
      den <- p1 * (1 - p1) * s2
      if (den <= 0) return(NA_real_)
      (a[1] * (1 - Rbar) + a[2] * (2 - Rbar)) / sqrt(den)
    }
    ## fisher.test(ceiling(table) + 1): $estimate (odds ratio) and $p.value
    ocache <- rep(NA_real_, 1024L * 1024L)
    pcache <- rep(NA_real_, 1024L * 1024L)
    fish <- function(m) {
      a <- ceiling(m) + 1
      key <- as.integer(a[1]) * 1024L + as.integer(a[2]) + 1L
      if (key >= 1L && key <= length(ocache) && !is.na(ocache[key]))
        return(c(ocache[key], pcache[key]))
      ft <- tryCatch(fisher.test(matrix(a, 2)),
                     error = function(e) NULL)
      v <- if (is.null(ft)) c(NA_real_, NA_real_)
      else c(as.numeric(ft$estimate), ft$p.value)
      if (key >= 1L && key <= length(ocache)) {
        ocache[key] <<- v[1]; pcache[key] <<- v[2]
      }
      v
    }
    ## ADDED: Pearson chi-square on the full 2 x K table
    chi2 <- function(r1, r2) {
      tot <- r1 + r2
      inv <- ifelse(tot > 0, 1 / tot, 0)
      n * (sum(r1^2 * inv) / n1 + sum(r2^2 * inv) / n2) - n
    }
    
    ## ADDED: Fisher on the RAW 2 x K table (no merging, no +1)
    f2cache <- new.env(hash = TRUE, parent = emptyenv())
    fish2k <- function(r1, r2) {
      key <- paste(r1, collapse = ",")
      hit <- f2cache[[key]]; if (!is.null(hit)) return(hit)
      tb <- rbind(r1, r2)
      tb <- tb[, colSums(tb) > 0, drop = FALSE]
      p <- if (ncol(tb) < 2) NA_real_ else
        tryCatch(fisher.test(tb, workspace = 2e8)$p.value,
                 error = function(e)
                   tryCatch(fisher.test(tb, simulate.p.value = TRUE, B = 1e5)$p.value,
                            error = function(e2) NA_real_))
      f2cache[[key]] <- p; p
    }
    
    ## ---- observed values, needed only by the ADDED columns ----
    o1 <- tabulate(permu.df$clusterInfo[i1], nbins = K)
    o2 <- tabulate(permu.df$clusterInfo[i2], nbins = K)
    chi_obs      <- chi2(o1, o2)
    fisher_obs   <- fish(two_col(o1, o2))[2]     # merged 2 x 2 (+1), as before
    fisher2k_obs <- fish2k(o1, o2)               # raw 2 x K
    n2k <- if (K <= FISHER2K.MAXK && is.finite(fisher2k_obs))
      min(FISHER2K.NSIM, NSIM) else 0L
    
    ## ---- the loop: identical draws ----
    set.seed(SEED)
    n_pur <- 0L; n_coc <- 0L; n_mco <- 0L; n_odd <- 0L
    n_chi <- 0L; n_fis <- 0L; n_2k <- 0L; d_2k <- 0L
    a_pur <- abs(original.purity)
    a_coc <- abs(original.cochran)
    a_mco <- abs(original.merge.cochran)
    
    for (i in seq_len(NSIM)) {
      perm <- sample.int(K, n, replace = TRUE)
      r1 <- tabulate(perm[i1], nbins = K)
      r2 <- tabulate(perm[i2], nbins = K)
      
      z <- co_raw(r1, r2)
      if (!is.na(z) && !is.na(a_coc) && abs(z) >= a_coc) n_coc <- n_coc + 1L
      
      m  <- two_col(r1, r2)
      rr <- purity_of(m)
      if (!is.na(rr) && !is.na(a_pur) && abs(rr) >= a_pur) n_pur <- n_pur + 1L
      
      zm <- co_merge(m)
      if (!is.na(zm) && !is.na(a_mco) && abs(zm) >= a_mco) n_mco <- n_mco + 1L
      
      fv <- fish(m)
      if (!is.na(fv[1]) && !is.na(original.Pvalue) &&
          fv[1] >= original.Pvalue) n_odd <- n_odd + 1L
      
      ## ADDED
      if (chi2(r1, r2) >= chi_obs) n_chi <- n_chi + 1L
      if (!is.na(fv[2]) && fv[2] <= fisher_obs) n_fis <- n_fis + 1L
      
      ## raw 2 x K Fisher, on the first n2k iterations of this same stream
      if (i <= n2k) {
        p2 <- fish2k(r1, r2)
        if (!is.na(p2)) {
          d_2k <- d_2k + 1L
          if (p2 <= fisher2k_obs) n_2k <- n_2k + 1L
        }
      }
    }
    
    res <- data.frame(
      gene.idx = gene,
      gene.cluster.idx = paste0("cluster_", bestK$K[g]),
      ## original columns, original formulas
      p.purity        = (n_pur + 1) / NSIM,
      P.fisher.merge  = (n_odd + 1) / NSIM,
      P.cochran       = (n_coc + 1) / NSIM,
      P.cochran.merge = (n_mco + 1) / NSIM,
      ## added
      P.chisq            = (n_chi + 1) / NSIM,
      P.fisher.2x2.perm  = (n_fis + 1) / NSIM,             # merged, "shrunk"
      P.fisher.2xK.perm  = if (d_2k > 0) (n_2k + 1) / (d_2k + 1) else NA_real_,
      N.perm.fisher.2xK  = d_2k,
      Chisq_observed = chi_obs,
      P.fisher.2x2.observed = fisher_obs,
      P.fisher.2xK.observed = fisher2k_obs,
      K_used = K, N_species = n,
      N_more_than_twice = n1, N_once = n2,
      Estimate.cochran = original.cochran,
      E.merge.cochran  = original.merge.cochran,
      Estimate.fisher  = original.Pvalue,
      purity_score     = original.purity,
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      Error = NA_character_, stringsAsFactors = FALSE)
    
    message(sprintf(
      "[%s] K=%d n=%d | p.purity=%.6f P.cochran=%.6f P.cochran.merge=%.6f P.fisher.merge=%.6f | P.chisq=%.6f F2x2=%.6f F2xK=%s (%.0fs)",
      gene, K, n, res$p.purity, res$P.cochran, res$P.cochran.merge,
      res$P.fisher.merge, res$P.chisq, res$P.fisher.2x2.perm,
      ifelse(is.na(res$P.fisher.2xK.perm), "skipped",
             sprintf("%.6f", res$P.fisher.2xK.perm)), res$seconds))
    res
  }, error = function(e) {
    message("[ERROR] ", gene, " : ", conditionMessage(e))
    data.frame(gene.idx = gene, Error = conditionMessage(e), stringsAsFactors = FALSE)
  })
  
  saveRDS(out, ck)
  out
}


############################################################
## run - one gene per task (splitting would change the RNG)
############################################################
t0 <- Sys.time(); cpu0 <- proc.time()
results <- mclapply(seq_len(NG), analyze_gene, mc.cores = N.CORES,
                    mc.preschedule = FALSE)
cpu1 <- proc.time()
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cpu  <- (cpu1["user.child"] - cpu0["user.child"]) + (cpu1["sys.child"] - cpu0["sys.child"])

result.df <- rbindlist(results, fill = TRUE)

## Keep only the permutation Fisher 2xK p-value used downstream
## (Gene / Cluster / P_Value == P.fisher.2xK.perm)
out.slim <- result.df[
  ,
  .(
    Gene    = gene.idx,
    Cluster = gene.cluster.idx,
    P_Value = P.fisher.2xK.perm
  )
]
fwrite(out.slim, OUT.FILE, sep = "\t")

cat(sprintf("\nwall %.1f s | child CPU %.1f s | avg cores busy %.1f / %d\n",
            wall, cpu, cpu / max(wall, 1e-9), N.CORES))
ok <- result.df[is.na(Error)]
cat("genes done:", nrow(ok), "| failed:", nrow(result.df) - nrow(ok), "\n")
cat("P.fisher.2xK.perm < 0.05:", sum(ok$P.fisher.2xK.perm < 0.05, na.rm = TRUE),
    " (of", sum(!is.na(ok$P.fisher.2xK.perm)), "genes computed)\n")
cat("Output columns: Gene, Cluster, P_Value (= P.fisher.2xK.perm)\n")
cat("Output:", OUT.FILE, "\n")