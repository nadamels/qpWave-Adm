#!/usr/bin/env Rscript
# =============================================================================
# Step 3: IBD "chromosome painting" tables.
#
# For every individual i and donor group g, per chromosome c:
#     S[i, g, c] = total cM shared IBD between i and the members of g, EXCLUDING i itself (leave-one-out: an individual
#                  cannot be its own donor, so a member of g sees g with n_g - 1 donors)
#     N[i, g]    = number of donors in g available to i (n_g, or n_g - 1 if i is in g)
# The painting profile of i (built in step 4 for the donor set of each model) is  S[i, g, .] / N[i, g]  (mean per donor
# individual; SIZE_NORM=0 uses the plain sums), normalised to sum to 1 over the donors ("[0,1]" in the methods).
# Keeping S per chromosome is what lets step 4 do the leave-one-chromosome-out (weighted block) jackknife.
#
# Env: PAIR_RDS GROUPS_TSV(ibd_groups.tsv) OUTDIR [DONOR_MIN_N=2]
# Out: painting.rds, group_mean_ibd_cM.tsv (mean cM shared between groups), painting_profiles_by_group.tsv
# =============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "ibd_common.R"))

PAIR_RDS <- getenv("PAIR_RDS"); GROUPS <- getenv("GROUPS_TSV"); OUTDIR <- getenv("OUTDIR")
DONOR_MIN_N <- as.integer(getenv("DONOR_MIN_N", 2))        # groups with fewer members are not used as donors (they can still be modelled)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

P <- readRDS(PAIR_RDS)
G <- fread(GROUPS, colClasses = list(character = c("IID", "group", "type")))
G <- G[type != "none"]
stopifnot(!anyDuplicated(G$IID))
ids <- P$ids
grp_of <- G$group[match(ids, G$IID)]                     # NA = individual not used
groups <- sort(unique(G$group))
gsize <- table(factor(grp_of, levels = groups))
n_g <- as.numeric(gsize)
donors <- names(gsize)[gsize >= DONOR_MIN_N]
chrs <- as.integer(names(P$chr_cM)); nI <- length(ids); nG <- length(groups); nC <- length(chrs)
msg(nI, " individuals in IBD table, ", sum(!is.na(grp_of)), " assigned to ", nG, " groups; ", length(donors), " groups have >= ", DONOR_MIN_N, " members (donors)")
if (length(donors) < 3) stop("fewer than 3 donor groups")

gi <- match(grp_of, groups)                               # group index per individual
pr <- P$pair[chr %in% chrs]
pr <- pr[!is.na(gi[a]) & !is.na(gi[b])]
long <- rbind(pr[, .(i = a, g = gi[b], c = match(chr, chrs), cM)], pr[, .(i = b, g = gi[a], c = match(chr, chrs), cM)])
long <- long[, .(cM = sum(cM)), by = .(i, g, c)]
S <- array(0, dim = c(nI, nG, nC), dimnames = list(ids, groups, chrs))
S[cbind(long$i, long$g, long$c)] <- long$cM
Nn <- matrix(0, nI, nG, dimnames = list(ids, groups))
for (k in seq_len(nG)) { Nn[, k] <- gsize[k]; Nn[gi == k & !is.na(gi), k] <- gsize[k] - 1 }
rm(long, pr)

# ---- group x group mean IBD (cM per pair of individuals), for a first look ---------------------------------
Stot <- apply(S, c(1, 2), sum)
mean_ibd <- t(sapply(groups, function(a) { m <- which(grp_of == a); colSums(Stot[m, , drop = FALSE]) / (length(m) * pmax(1, n_g - ifelse(groups == a, 1, 0))) }))
dimnames(mean_ibd) <- list(groups, groups)
fwrite(data.table(group = groups, n = as.integer(gsize), mean_ibd), file.path(OUTDIR, "group_mean_ibd_cM.tsv"), sep = "\t")

# ---- all-donor painting profile per group (for inspection; models re-normalise over their own donor set) --------
prof_ind <- function(rows, don) { M <- Stot[rows, don, drop = FALSE] / Nn[rows, don, drop = FALSE]; s <- rowSums(M); M[s > 0, , drop = FALSE] / s[s > 0] }
prof <- t(sapply(groups, function(a) { m <- which(grp_of == a); don <- setdiff(donors, a); p <- prof_ind(m, don); out <- setNames(rep(NA_real_, length(donors)), donors); if (nrow(p)) out[don] <- colMeans(p); out }))
fwrite(data.table(group = groups, n = as.integer(gsize), prof), file.path(OUTDIR, "painting_profiles_by_group.tsv"), sep = "\t")

zero <- ids[!is.na(grp_of) & rowSums(Stot[, donors, drop = FALSE]) == 0]
if (length(zero)) msg("WARNING: ", length(zero), " assigned individuals share no IBD with any donor group (excluded from profiles), e.g. ", paste(head(zero, 5), collapse = ", "))
saveRDS(list(S = S, N = Nn, ids = ids, grp_of = grp_of, groups = groups, donors = donors, gsize = gsize, chr_cM = P$chr_cM, types = unique(G[, .(group, type)])),
        file.path(OUTDIR, "painting.rds"))
msg("wrote ", file.path(OUTDIR, "painting.rds"))
