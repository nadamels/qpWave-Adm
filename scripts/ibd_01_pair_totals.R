#!/usr/bin/env Rscript
# =============================================================================
# Step 1: stream the filtered IBD segment table once and aggregate it to
#         total cM shared per (pair of individuals, chromosome).
#
# Input  : IBD_FILE = *.final.with_cM.tsv.gz  (one IBD segment per line, tab separated, with header; already filtered
#          LOD >= 3, >= 2 cM, excess-long-IBD regions removed). Column names are auto-detected, or set COL_* below.
# Output : $OUTDIR/ibd_pair_chr.rds   list(ids, pair[a,b,chr,cM,n], roh[a,chr,cM,n], chr_cM, ...)
#          $OUTDIR/ibd_pair_totals.log, roh_by_individual.tsv
# The file is read in chunks through zcat, so memory is bounded by the (small) aggregated table, not the 4 GB input.
# =============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "ibd_common.R"))

IBD_FILE   <- getenv("IBD_FILE")
SAMPLES    <- getenv("SAMPLES")                        # from ibd_make_samples.py
OUTDIR     <- getenv("OUTDIR")
CHUNK      <- as.integer(getenv("CHUNK_LINES", 4000000))
MIN_CM     <- as.numeric(getenv("MIN_CM", 0))          # extra length filter on top of the file's own (e.g. 4)
ID_SUB     <- getenv("ID_SUB", "")                     # regex removed from IDs in the IBD file, e.g. "_[12]$" or "^.*:"
# optional column overrides: a header NAME or a 1-based INDEX
COL <- list(id1 = getenv("COL_ID1", ""), id2 = getenv("COL_ID2", ""), chr = getenv("COL_CHR", ""),
            len = getenv("COL_LEN_CM", ""), scm = getenv("COL_START_CM", ""), ecm = getenv("COL_END_CM", ""))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(OUTDIR, "ibd_pair_totals.log")
cat("", file = logf)
logm <- function(...) { msg(...); cat(paste0(..., "\n"), file = logf, append = TRUE) }

samples <- fread(SAMPLES, na.strings = "NA", colClasses = list(character = "IID"))
ids <- samples$IID
stopifnot(!anyDuplicated(ids))

decomp <- if (grepl("\\.(gz|bgz)$", IBD_FILE)) { if (nzchar(Sys.which("pigz"))) "pigz -dc" else "gzip -dc" } else "cat"
cmd <- paste(decomp, shQuote(IBD_FILE))
con <- pipe(cmd, "r")
hdr <- readLines(con, n = 1)
hnames <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
logm("IBD file: ", IBD_FILE); logm("header (", length(hnames), " columns): ", paste(hnames, collapse = " | "))

pick <- function(override, patterns) {
  if (nzchar(override)) {
    i <- if (grepl("^[0-9]+$", override)) as.integer(override) else match(override, hnames)
    if (is.na(i) || i < 1 || i > length(hnames)) stop("column '", override, "' not found in header")
    return(i)
  }
  for (p in patterns) { i <- grep(p, hnames, ignore.case = TRUE); if (length(i) >= 1) return(i[1]) }
  NA_integer_
}
ci <- list(
  id1 = pick(COL$id1, c("^(id|iid|ind|indiv|individual|sample)_?1$", "^(id|iid|ind|sample)_?a$")),
  id2 = pick(COL$id2, c("^(id|iid|ind|indiv|individual|sample)_?2$", "^(id|iid|ind|sample)_?b$")),
  chr = pick(COL$chr, c("^(chr|chrom|chromosome)$")),
  len = pick(COL$len, c("^(length|len)_?(in_)?cm$", "^cm_?(length|len)$", "^(ibd_)?(length|len)_?(\\(cm\\)|cm)", "^cm$", "^length_?genetic")),
  scm = pick(COL$scm, c("^(start|begin|from)_?cm$", "^cm_?(start|begin)", "^genetic_?start")),
  ecm = pick(COL$ecm, c("^(end|stop|to)_?cm$", "^cm_?(end|stop)", "^genetic_?end")))
if (anyNA(ci[c("id1", "id2", "chr")]))
  stop("Could not identify the id1/id2/chr columns from the header. Set COL_ID1, COL_ID2, COL_CHR (name or 1-based index).")
if (is.na(ci$len) && (is.na(ci$scm) || is.na(ci$ecm)))
  stop("Could not find a cM length column (or start/end cM). Set COL_LEN_CM, or COL_START_CM and COL_END_CM.")
logm("columns used: ", paste(sprintf("%s=%s", names(ci), vapply(ci, function(i) if (is.na(i)) "-" else hnames[i], "")), collapse = "  "))
use <- sort(unique(unlist(ci[!is.na(unlist(ci))])))
cc <- list(character = c(ci$id1, ci$id2, ci$chr))

acc <- NULL; roh <- NULL; chr_max <- numeric(0)
n_lines <- 0; n_kept <- 0; n_unmatched <- 0; unm <- character(0); n_chunk <- 0
while (length(lines <- readLines(con, n = CHUNK)) > 0) {
  n_chunk <- n_chunk + 1; n_lines <- n_lines + length(lines)
  d <- fread(text = lines, sep = "\t", header = FALSE, select = use, colClasses = cc, showProgress = FALSE)
  setnames(d, hnames[use])
  i1 <- d[[hnames[ci$id1]]]; i2 <- d[[hnames[ci$id2]]]
  if (nzchar(ID_SUB)) { i1 <- sub(ID_SUB, "", i1); i2 <- sub(ID_SUB, "", i2) }
  a <- match(i1, ids); b <- match(i2, ids)
  bad <- is.na(a) | is.na(b)
  if (any(bad)) { n_unmatched <- n_unmatched + sum(bad); unm <- unique(c(unm, head(unique(c(i1[is.na(a)], i2[is.na(b)])), 20))) }
  chr <- chr_to_int(d[[hnames[ci$chr]]])
  if (!is.na(ci$len)) len <- as.numeric(d[[hnames[ci$len]]]) else len <- as.numeric(d[[hnames[ci$ecm]]]) - as.numeric(d[[hnames[ci$scm]]])
  if (!is.na(ci$ecm)) { m <- tapply(as.numeric(d[[hnames[ci$ecm]]]), chr, max, na.rm = TRUE)
    chr_max[names(m)] <- pmax(chr_max[names(m)], m, na.rm = TRUE) }
  ok <- !bad & !is.na(chr) & !is.na(len) & len >= MIN_CM
  n_kept <- n_kept + sum(ok)
  s <- data.table(a = pmin(a, b)[ok], b = pmax(a, b)[ok], chr = chr[ok], cM = len[ok])
  roh_c <- s[a == b, .(cM = sum(cM), n = .N), by = .(a, chr)]
  agg   <- s[a != b, .(cM = sum(cM), n = .N), by = .(a, b, chr)]
  acc <- rbind(acc, agg); roh <- rbind(roh, roh_c)
  if (n_chunk %% 5 == 0) { acc <- acc[, .(cM = sum(cM), n = sum(n)), by = .(a, b, chr)]; roh <- roh[, .(cM = sum(cM), n = sum(n)), by = .(a, chr)] }
  logm(sprintf("chunk %d: %s lines read, %s kept, %s pair-chr rows", n_chunk, format(n_lines, big.mark = ","), format(n_kept, big.mark = ","), format(nrow(acc), big.mark = ",")))
}
acc <- acc[, .(cM = sum(cM), n = sum(n)), by = .(a, b, chr)]
roh <- if (is.null(roh)) data.table(a = integer(), chr = integer(), cM = numeric(), n = integer()) else roh[, .(cM = sum(cM), n = sum(n)), by = .(a, chr)]
setorder(acc, a, b, chr)

if (n_unmatched > 0) {
  logm(sprintf("WARNING: %s segments dropped because an ID is not in the sample table (e.g. %s)", format(n_unmatched, big.mark = ","), paste(head(unm, 8), collapse = ", ")))
  if (n_unmatched > 0.5 * n_lines) stop("More than half of the segments have unknown IDs: check the ID format in the IBD file vs ", SAMPLES, " (use ID_SUB to strip a suffix/prefix).")
}
if (length(chr_max) == 0) { chr_cM <- FALLBACK_CHR_CM; logm("no cM end column: using approximate chromosome lengths (FALLBACK_CHR_CM)")
} else { chr_cM <- chr_max[order(as.integer(names(chr_max)))]; chr_cM <- chr_cM[names(chr_cM) %in% as.character(1:22)] }
if (nzchar(Sys.getenv("GENOME_CM"))) { g <- as.numeric(Sys.getenv("GENOME_CM")); chr_cM <- chr_cM * g / sum(chr_cM) }
logm("chromosome lengths used (cM): ", paste(sprintf("%s=%.0f", names(chr_cM), chr_cM), collapse = " "), "  total=", round(sum(chr_cM)))
logm(sprintf("done: %s segment lines, %s kept (%s), %s pair-chromosome rows, %s pairs sharing IBD",
             format(n_lines, big.mark = ","), format(n_kept, big.mark = ","), if (MIN_CM > 0) paste0(">= ", MIN_CM, " cM") else "no extra filter",
             format(nrow(acc), big.mark = ","), format(uniqueN(acc[, .(a, b)]), big.mark = ",")))

saveRDS(list(ids = ids, pair = acc, roh = roh, chr_cM = chr_cM, file = IBD_FILE, min_cM = MIN_CM), file.path(OUTDIR, "ibd_pair_chr.rds"))
r <- roh[, .(roh_cM = sum(cM), roh_n = sum(n)), by = a]
r <- data.table(IID = ids[r$a], roh_cM = r$roh_cM, roh_n = r$roh_n)
fwrite(r, file.path(OUTDIR, "roh_by_individual.tsv"), sep = "\t")
logm("within-individual (ROH) segments in this file: ", nrow(r), " individuals (file-level cut-off only; use IBDseq output >= 1 cM for ROH analyses)")
close(con)
