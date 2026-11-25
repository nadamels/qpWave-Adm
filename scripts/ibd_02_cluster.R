#!/usr/bin/env Rscript
# =============================================================================
# Step 2: genetic clustering of the ANCIENT individuals by hierarchical (iterative) Leiden community detection
#         on the network of pairwise IBD sharing, then collapse the hierarchy into final clusters using
#         spatiotemporal annotations.
#
#  network  : nodes = ancient individuals, edge weight = fraction of the genome shared IBD (total cM / genome cM)
#  Leiden   : start at resolution r = 0.5. > 1 community found -> split the network into the communities and repeat
#             (each sub-network restarts at r = 0.5). Only 1 community -> r += 0.5, up to R_MAX (3). Leaf = nothing left
#             to split at R_MAX in any sub-community.
#  collapse : top-down, a node becomes ONE final cluster if its members are spatiotemporally coherent
#             (date span <= MAX_DATE_SPAN yr AND max pairwise distance <= MAX_KM km); otherwise descend.
#             config/ibd_cluster_overrides.tsv (node_path <TAB> cluster_name | SPLIT) overrides this per node.
#
# Env: PAIR_RDS SAMPLES OUTDIR [R_START=0.5 R_STEP=0.5 R_MAX=3 MIN_COMMUNITY=1 MAX_DATE_SPAN=1000 MAX_KM=500
#      MIN_PAIR_FRAC=0 OVERRIDES=<file> SEED=1]
# =============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "ibd_common.R"))
suppressPackageStartupMessages(library(igraph))

PAIR_RDS <- getenv("PAIR_RDS"); SAMPLES <- getenv("SAMPLES"); OUTDIR <- getenv("OUTDIR")
R_START <- as.numeric(getenv("R_START", 0.5)); R_STEP <- as.numeric(getenv("R_STEP", 0.5)); R_MAX <- as.numeric(getenv("R_MAX", 3))
MIN_COMMUNITY <- as.integer(getenv("MIN_COMMUNITY", 1))        # a split is only accepted if every community has >= this many members
MAX_DATE_SPAN <- as.numeric(getenv("MAX_DATE_SPAN", 1000))     # years
MAX_KM        <- as.numeric(getenv("MAX_KM", 500))
MIN_PAIR_FRAC <- as.numeric(getenv("MIN_PAIR_FRAC", 0))        # drop edges with weight below this fraction of the genome
OVERRIDES     <- getenv("OVERRIDES", "")
SEED          <- as.integer(getenv("SEED", 1))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

P <- readRDS(PAIR_RDS)
S <- fread(SAMPLES, na.strings = "NA", colClasses = list(character = "IID"))
stopifnot(identical(S$IID, P$ids))
genome_cM <- sum(P$chr_cM)
anc <- which(S$is_ancient)
msg("ancient individuals: ", length(anc), "   genome length used: ", round(genome_cM), " cM   leiden: ",
    if (requireNamespace("leidenAlg", quietly = TRUE)) "leidenAlg" else "igraph::cluster_leiden (leidenAlg not installed)")

# ---- network -------------------------------------------------------------------------------------------------
E <- P$pair[a %in% anc & b %in% anc, .(cM = sum(cM)), by = .(a, b)]
E[, weight := cM / genome_cM]
E <- E[weight > MIN_PAIR_FRAC]
connected <- sort(unique(c(E$a, E$b)))
isolated <- setdiff(anc, connected)
if (length(isolated)) msg("WARNING: ", length(isolated), " ancient individuals share no IBD with any other ancient individual: left unclustered")
g <- graph_from_data_frame(data.frame(from = S$IID[E$a], to = S$IID[E$b], weight = E$weight), directed = FALSE,
                           vertices = data.frame(name = S$IID[connected]))
msg("network: ", vcount(g), " nodes, ", ecount(g), " edges")

# ---- iterative Leiden ----------------------------------------------------------------------------------------
nodes <- new.env()      # path -> list(members, r)  (r = resolution at which this node was split; NA for leaves)
split_node <- function(members, path) {
  nodes[[path]] <- list(members = members, r = NA_real_)
  if (length(members) < 2 * MIN_COMMUNITY || length(members) < 2) return(invisible())
  sub <- induced_subgraph(g, members)
  for (r in seq(R_START, R_MAX, by = R_STEP)) {
    memb <- leiden_membership(sub, r, SEED)
    if (length(unique(memb)) > 1 && min(table(memb)) >= MIN_COMMUNITY) {
      nodes[[path]]$r <- r
      tab <- table(memb)
      comms <- as.integer(names(sort(tab, decreasing = TRUE)))             # children numbered by decreasing size
      for (k in seq_along(comms)) split_node(members[memb == comms[k]], paste0(path, ".", k))
      return(invisible())
    }
  }
  invisible()
}
split_node(V(g)$name, "1")
paths <- sort(ls(nodes))
is_leaf_path <- vapply(paths, function(p) is.na(nodes[[p]]$r), TRUE)
msg("hierarchy: ", length(paths), " nodes, ", sum(is_leaf_path), " leaves")

# ---- node annotation -----------------------------------------------------------------------------------------
setkey(S, IID)
haversine_max <- function(lat, lon) {
  ok <- !is.na(lat) & !is.na(lon); if (sum(ok) < 2) return(NA_real_)
  la <- lat[ok] * pi / 180; lo <- lon[ok] * pi / 180
  d <- outer(seq_along(la), seq_along(la), function(i, j)
    2 * 6371 * asin(pmin(1, sqrt(sin((la[j] - la[i]) / 2)^2 + cos(la[i]) * cos(la[j]) * sin((lo[j] - lo[i]) / 2)^2))))
  max(d)
}
top_pop <- function(x) { t <- sort(table(x), decreasing = TRUE); names(t)[1] }
rng <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
node_info <- rbindlist(lapply(paths, function(p) {
  m <- nodes[[p]]$members; s <- S[m]
  data.table(path = p, depth = lengths(regmatches(p, gregexpr("\\.", p))), n = length(m), split_r = nodes[[p]]$r,
             top_population = top_pop(s$population), n_populations = uniqueN(s$population),
             date_min_BP = rng(s$date_mid_BP, min), date_max_BP = rng(s$date_mid_BP, max),
             date_median_BP = rng(s$date_mid_BP, median), max_km = haversine_max(s$lat, s$lon))
}))
node_info[, date_span := date_max_BP - date_min_BP]
node_info[, coherent := n == 1 | ((is.na(date_span) | date_span <= MAX_DATE_SPAN) & (is.na(max_km) | max_km <= MAX_KM))]

ov <- if (nzchar(OVERRIDES) && file.exists(OVERRIDES)) read_tsv_nocomment(OVERRIDES) else
  data.table(node_path = character(), cluster_name = character())
if (nrow(ov)) msg("overrides: ", nrow(ov), " node(s) from ", OVERRIDES)
bad <- setdiff(ov$node_path, paths); if (length(bad)) stop("override node(s) not in hierarchy: ", paste(bad, collapse = ", "))

# ---- collapse (top-down) -------------------------------------------------------------------------------------
final <- list()   # path -> cluster name
descend <- function(p) {
  info <- node_info[path == p]; o <- ov[node_path == p, cluster_name]
  if (length(o) && o != "SPLIT") { final[[p]] <<- o; return(invisible()) }
  leaf <- is.na(info$split_r)
  if (leaf || (info$coherent && !length(o))) {
    final[[p]] <<- paste0(gsub("[^A-Za-z0-9]+", "_", info$top_population), "_c", gsub("\\.", "-", p)); return(invisible())
  }
  for (k in grep(paste0("^", gsub(".", "\\.", p, fixed = TRUE), "\\.[0-9]+$"), paths, value = TRUE)) descend(k)
}
descend("1")
cl <- rbindlist(lapply(names(final), function(p) data.table(IID = nodes[[p]]$members, cluster_node = p, cluster = final[[p]])))
leaf_of <- rbindlist(lapply(paths[is_leaf_path], function(p) data.table(IID = nodes[[p]]$members, leaf_path = p)))
cl <- merge(cl, leaf_of, by = "IID", all.x = TRUE)
if (length(isolated)) cl <- rbind(cl, data.table(IID = S$IID[isolated], cluster_node = NA_character_, cluster = "Unclustered", leaf_path = NA_character_))
cl <- merge(cl, S[, .(IID, population, date_mid_BP, lat, lon)], by = "IID", all.x = TRUE)
setorder(cl, cluster, IID)

cs <- cl[, .(n = .N, top_population = top_pop(population), n_populations = uniqueN(population),
             populations = paste(names(sort(table(population), decreasing = TRUE))[1:min(4, uniqueN(population))], collapse = ";"),
             median_date_BP = rng(date_mid_BP, median), date_min_BP = rng(date_mid_BP, min), date_max_BP = rng(date_mid_BP, max),
             median_lat = rng(lat, median), median_lon = rng(lon, median)), by = .(cluster, cluster_node)]
setorder(cs, -n)

fwrite(node_info, file.path(OUTDIR, "ibd_cluster_hierarchy_nodes.tsv"), sep = "\t")
fwrite(cl, file.path(OUTDIR, "ibd_clusters.tsv"), sep = "\t")
fwrite(cs, file.path(OUTDIR, "ibd_cluster_summary.tsv"), sep = "\t")

# ---- group table used by the painting / ancestry-modelling steps ------------------------------------------------
grp <- rbind(
  data.table(IID = cl$IID, group = cl$cluster, type = ifelse(cl$cluster == "Unclustered", "none", "ancient_cluster")),
  S[is_ancient == FALSE & !is.na(modern_group),
    .(IID, group = modern_group, type = ifelse(grepl("^Mod_", modern_group), "modern_target", "modern_donor"))])
fwrite(grp, file.path(OUTDIR, "ibd_groups.tsv"), sep = "\t")
msg("final: ", uniqueN(cl$cluster), " ancient clusters (largest sizes: ", paste(head(cs$n, 12), collapse = ","), "), ",
    uniqueN(grp[type %in% c("modern_target", "modern_donor"), group]), " modern groups -> ", file.path(OUTDIR, "ibd_groups.tsv"))
msg("Inspect ibd_cluster_hierarchy_nodes.tsv / ibd_cluster_summary.tsv; optionally edit config/ibd_cluster_overrides.tsv and rerun this step.")
