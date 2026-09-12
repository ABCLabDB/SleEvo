############################################################
## Result4_Sleeptiming_tree_association_test.R
##
## Sleep-timing tree / cluster association with permutation
## Fisher (2 x K). Output is slim:
##   Gene, Cluster, P_Value (= P.fisher.2xK.perm)
##
## Based on Sleeptiming/Sleep_timing_permutation.R
## Paths point at the SleEvo repository data layout.
############################################################

############################
## Paths
############################
DATA.DIR   <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Fasta"
META.FILE  <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Result1/species_sleep_metadata.txt"
BESTK.FILE <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Result4/Sleep_Timing_optimal_K.tsv"
OUT.FILE   <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Result4/Sleeptiming_Fisher_Result.tsv"
CHECK.DIR  <- "C:/Users/user/Desktop/Cursor/SleEvo/data/Result4/checkpoint_literal"

############################
## Column names in METADATA
##   the original used "species symbol name ensembl" / "Sleep_timing"
############################
COL.SPECIES.SYMBOL <- "Species_symbol_name_ensembl"
COL.SPECIES.NAME   <- "Species_name_ensembl"     # becomes idx (column 1 in the original)
COL.TIMING         <- "Sleep_timing_per_day"

## the original recoded Cathemeral/diurnal/nocturnal then labelled
##   Sleep_at_daytime -> more_sleep , Sleep_at_night -> less_sleep
## spell out the same mapping for whatever spellings this file uses
LABEL.MAP <- c(
  "Sleep at daytime" = "more_sleep",
  "Sleep at night"   = "less_sleep",
  "Sleep_at_daytime" = "more_sleep",
  "Sleep_at_night"   = "less_sleep",
  "nocturnal"        = "more_sleep",
  "diurnal"          = "less_sleep"
)

############################
## Column names in the optimal-K file
############################
COL.GENE           <- "Gene"
COL.BESTCLUSTER    <- "Cluster"            # the original renamed column 2 to best.Cluster
COL.PURITY         <- "purity_score"       # original.purity
COL.EST.COCHRAN    <- "Estimate.cochran"   # original.cochran
COL.EST.FISHER     <- "Estimate.fisher"    # original.Pvalue

############################
## Run
############################
N.CORES <- 128L
NSIM    <- 2000000L
SEED    <- 319L
RESUME  <- TRUE
ONLY.GENES <- NULL          # e.g. c("PER1") to test one gene

## Fisher on the RAW 2 x K table (no merging).
## The original only ever ran Fisher on the merged 2 x 2, because a 2 x K
## exact test enumerates every table with the observed margins: measured
## at ~29 ms per call for K = 14, n = 45, and ~98% of drawn tables are
## distinct, so memoisation does not help. Two million of them would take
## about 8 h per gene, so it is applied to the FIRST FISHER2K.NSIM
## iterations of the SAME stream (no extra RNG is consumed, so every other
## column is unaffected). 20,000 is roughly 10 min for a K = 14 gene.
## Set FISHER2K.NSIM <- 0L to skip it entirely.
FISHER2K.NSIM <- 20000L
FISHER2K.MAXK <- 20L        # skip genes with a larger K

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
## Metadata - same steps as the original
############################################################
meta <- as.data.frame(fread(META.FILE))
for (cc in c(COL.SPECIES.SYMBOL, COL.SPECIES.NAME, COL.TIMING))
  must_have(meta, cc, "the metadata file")

cat("raw", COL.TIMING, "values:\n")
print(table(meta[[COL.TIMING]], useNA = "ifany"))

meta$idx <- as.character(meta[[COL.SPECIES.NAME]])
meta$sym <- gsub(" ", "_", tolower(trimws(meta[[COL.SPECIES.SYMBOL]])))
meta$Sleep_label <- unname(LABEL.MAP[trimws(as.character(meta[[COL.TIMING]]))])

unmapped <- setdiff(unique(trimws(as.character(meta[[COL.TIMING]]))[
  is.na(meta$Sleep_label)]), c(NA, ""))
if (length(unmapped))
  cat("\nnot in LABEL.MAP, so dropped:", paste(unmapped, collapse = " | "), "\n")

meta <- meta[!is.na(meta$Sleep_label) & !is.na(meta$idx), ]
meta <- meta[!duplicated(meta$idx), c("idx", "sym", "Sleep_label")]
rownames(meta) <- NULL
cat("\nSleep_label:\n"); print(table(meta$Sleep_label))
cat("species:", nrow(meta), "\n")
if (length(unique(meta$Sleep_label)) < 2) stop("only one Sleep_label present")

############################################################
## optimal-K file - the original read these columns by name
############################################################
bk <- as.data.frame(fread(BESTK.FILE))
cat("\noptimal-K file columns:\n  ", paste(colnames(bk), collapse = ", "), "\n")
for (cc in c(COL.GENE, COL.BESTCLUSTER, COL.EST.COCHRAN))
  must_have(bk, cc, "the optimal-K file")
has_fisher <- COL.EST.FISHER %in% colnames(bk)
has_purity <- COL.PURITY     %in% colnames(bk)
if (!has_fisher) cat("!! '", COL.EST.FISHER, "' absent - P.fisher (the original\n",
                     "   column) cannot be computed and will be NA\n", sep = "")
if (!has_purity) cat("!! '", COL.PURITY, "' absent - p.purity will be NA\n", sep = "")

bestK <- data.frame(
  Gene = as.character(bk[[COL.GENE]]),
  K    = suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(bk[[COL.BESTCLUSTER]])))),
  orig_cochran = suppressWarnings(as.numeric(bk[[COL.EST.COCHRAN]])),
  orig_fisher  = if (has_fisher) suppressWarnings(as.numeric(bk[[COL.EST.FISHER]])) else NA_real_,
  orig_purity  = if (has_purity) suppressWarnings(as.numeric(bk[[COL.PURITY]]))     else NA_real_,
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
    
    original.cochran <- bestK$orig_cochran[g]
    original.Pvalue  <- bestK$orig_fisher[g]
    original.purity  <- bestK$orig_purity[g]
    
    ## ---- sequences (original steps) ----
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
    
    ## ---- tree (original steps) ----
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
    u <- unique(permu.df$clusterInfo)      # what sample() draws from
    K <- length(u)
    if (n < 5 || K < 2L) stop("too few species or clusters")
    
    ## table() sorts rows alphabetically: less_sleep first
    ylab <- as.integer(factor(permu.df$Sleep_label,
                              levels = c("less_sleep", "more_sleep")))
    i1 <- which(ylab == 1L); i2 <- which(ylab == 2L)
    n1 <- length(i1); n2 <- length(i2)
    if (n1 == 0L || n2 == 0L) stop("one phenotype missing")
    KMAX <- max(u)                          # table() sorts columns too
    
    ## ---- statistics ----
    Rbar2 <- (n1 + 2 * n2) / n
    s2_2  <- n1 * (1 - Rbar2)^2 + n2 * (2 - Rbar2)^2
    pp    <- n1 / n
    
    co_raw <- function(less, more) {        # CochranArmitageTest closed form
      tot <- less + more
      if (KMAX == 2L) {                     # a 2x2 is not transposed
        p1 <- tot[1] / n; den <- p1 * (1 - p1) * s2_2
        if (den <= 0) return(NA_real_)
        return((less[1] * (1 - Rbar2) + more[1] * (2 - Rbar2)) / sqrt(den))
      }
      pres <- tot > 0                       # table() omits empty clusters
      S <- cumsum(pres) * pres
      Rbar <- sum(tot * S) / n
      s2 <- sum(tot * (S - Rbar)^2 * pres)
      den <- pp * (1 - pp) * s2
      if (den <= 0) return(NA_real_)
      sum(less * (S - Rbar)) / sqrt(den)
    }
    mergeit <- function(less, more) {       # making.sum.cluster
      tot <- less + more
      tie <- (less == more) & (tot > 0)
      ti  <- sum(less[tie]) / 2
      lm <- less > more; mm <- more > less
      c(sum(less[lm]) + ti, sum(more[lm]) + ti,
        sum(less[mm]) + ti, sum(more[mm]) + ti)
    }
    purity_of <- function(m) {              # sum(colmax/colsum)/ncol
      s1 <- m[1] + m[2]; s2 <- m[3] + m[4]
      (if (s1 > 0) max(m[1], m[2]) / s1 else NaN) / 2 +
        (if (s2 > 0) max(m[3], m[4]) / s2 else NaN) / 2
    }
    chi2 <- function(less, more) {          # NEW
      tot <- less + more
      inv <- ifelse(tot > 0, 1 / tot, 0)
      n * (sum(less^2 * inv) / n1 + sum(more^2 * inv) / n2) - n
    }
    fcache <- rep(NA_real_, 512L * 512L)    # NEW: memoised fisher.test
    fish <- function(m) {
      key <- as.integer(2 * m[1]) * 512L + as.integer(2 * m[2]) + 1L
      if (key >= 1L && key <= length(fcache) && !is.na(fcache[key])) return(fcache[key])
      p <- tryCatch(fisher.test(matrix(ceiling(m), 2))$p.value,
                    error = function(e) NA_real_)
      if (key >= 1L && key <= length(fcache)) fcache[key] <<- p
      p
    }
    
    ## NEW: Fisher on the RAW 2 x K table, memoised on the less_sleep row
    f2cache <- new.env(hash = TRUE, parent = emptyenv())
    fish2k <- function(less, more) {
      key <- paste(less, collapse = ",")
      hit <- f2cache[[key]]; if (!is.null(hit)) return(hit)
      tb <- rbind(less, more)
      tb <- tb[, colSums(tb) > 0, drop = FALSE]
      p <- if (ncol(tb) < 2) NA_real_ else
        tryCatch(fisher.test(tb, workspace = 2e8)$p.value,
                 error = function(e)
                   tryCatch(fisher.test(tb, simulate.p.value = TRUE, B = 1e5)$p.value,
                            error = function(e2) NA_real_))
      f2cache[[key]] <- p; p
    }
    
    ## ---- observed (only needed by the NEW statistics) ----
    less_o <- tabulate(permu.df$clusterInfo[i1], nbins = KMAX)
    more_o <- tabulate(permu.df$clusterInfo[i2], nbins = KMAX)
    chi_obs      <- chi2(less_o, more_o)
    fisher_obs   <- fish(mergeit(less_o, more_o))     # merged 2 x 2
    fisher2k_obs <- fish2k(less_o, more_o)            # raw 2 x K
    n2k <- if (KMAX <= FISHER2K.MAXK && is.finite(fisher2k_obs))
      min(FISHER2K.NSIM, NSIM) else 0L
    
    ## ---- the loop, same draws as the original ----
    set.seed(SEED)
    n_pur <- 0L; n_fis_orig <- 0L; n_inf <- 0L; n_coc <- 0L
    n_chi <- 0L; n_fis_new <- 0L; n_2k <- 0L; d_2k <- 0L
    a_coc <- abs(original.cochran)
    a_pur <- abs(original.purity)
    
    for (i in seq_len(NSIM)) {
      clp  <- u[sample.int(K, n, replace = TRUE)]
      less <- tabulate(clp[i1], nbins = KMAX)
      more <- tabulate(clp[i2], nbins = KMAX)
      
      z <- co_raw(less, more)               # Co_estimate[i] and P_odds[i]
      m <- mergeit(less, more)
      r <- purity_of(m)                     # res[i]
      
      if (!is.na(r) && abs(r) >= a_pur)               n_pur      <- n_pur + 1L
      if (!is.na(z)) {
        if (z >= original.Pvalue)                     n_fis_orig <- n_fis_orig + 1L
        if (abs(z) >= a_coc)                          n_coc      <- n_coc + 1L
      } else if (is.infinite(z))                      n_inf      <- n_inf + 1L
      
      if (chi2(less, more) >= chi_obs)                n_chi      <- n_chi + 1L
      pf <- fish(m)
      if (!is.na(pf) && pf <= fisher_obs)             n_fis_new  <- n_fis_new + 1L
      
      ## raw 2 x K Fisher, on the first n2k iterations of this same stream
      if (i <= n2k) {
        p2 <- fish2k(less, more)
        if (!is.na(p2)) {
          d_2k <- d_2k + 1L
          if (p2 <= fisher2k_obs) n_2k <- n_2k + 1L
        }
      }
    }
    
    res <- data.frame(
      gene.idx = gene,
      gene.cluster.idx = paste0("cluster_", bestK$K[g]),
      ## the original three columns, same formulas
      p.purity  = (n_pur + 1) / NSIM,
      P.fisher  = (n_fis_orig + n_inf) / NSIM,
      P.cochran = (n_coc + 1) / NSIM,
      ## added
      P.chisq            = (n_chi + 1) / NSIM,
      P.fisher.2x2.perm  = (n_fis_new + 1) / NSIM,          # merged, "shrunk"
      P.fisher.2xK.perm  = if (d_2k > 0) (n_2k + 1) / (d_2k + 1) else NA_real_,
      N.perm.fisher.2xK  = d_2k,
      Chisq_observed = chi_obs,
      P.fisher.2x2.observed = fisher_obs,
      P.fisher.2xK.observed = fisher2k_obs,
      K_used = KMAX, N_species = n, N_less = n1, N_more = n2,
      Estimate.cochran = original.cochran,
      Estimate.fisher  = original.Pvalue,
      purity_score     = original.purity,
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      Error = NA_character_, stringsAsFactors = FALSE)
    
    message(sprintf("[%s] K=%d n=%d | p.purity=%.6f P.fisher=%.6f P.cochran=%.6f | P.chisq=%.6f F2x2=%.6f F2xK=%s (%.0fs)",
                    gene, KMAX, n, res$p.purity, res$P.fisher, res$P.cochran,
                    res$P.chisq, res$P.fisher.2x2.perm,
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