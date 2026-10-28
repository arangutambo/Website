# Builds the six panels of "Figure 1" on the home page from the real project outputs.
# Run from the Website folder:  Rscript tools/build_panels.R   (writes figures/panel-*.svg)
suppressPackageStartupMessages({library(data.table); library(DESeq2); library(airway); library(phyloseq)})
dev <- normalizePath("~/Developer"); gw <- file.path(dev, "gwas-1000g-simulation")
INK <- "#2A2521"; MUTED <- "#5C5448"; RUST <- "#8C4B34"; GREEN <- "#53624A"; BROWN <- "#5C4736"; AMBER <- "#B8923A"; SLATE <- "#7A8C9E"; PAPER <- "#EAE3D6"
f1 <- function(x) formatC(x, format = "f", digits = 1)
mono <- 'font-family="ui-monospace,Menlo,Consolas,monospace" font-size="8.5"'
wrap <- function(W, H, body, title, desc) sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" role="img" aria-labelledby="t d"><title id="t">%s</title><desc id="d">%s</desc><g %s fill="%s">%s</g></svg>', W, H, title, desc, mono, MUTED, body)
frame <- function(W, H, L = 34, R = 8, T = 16, B = 24, xlab = "", ylab = "", xt = NULL, yt = NULL) {
  s <- sprintf('<path d="M%d %d H%d M%d %d V%d" stroke="%s" fill="none"/>', L, H - B, W - R, L, T, H - B, INK)
  for (t in xt) s <- paste0(s, sprintf('<text x="%s" y="%d" text-anchor="middle">%s</text>', f1(t$x), H - B + 11, t$l))
  for (t in yt) s <- paste0(s, sprintf('<text x="%d" y="%s" text-anchor="end">%s</text>', L - 4, f1(t$y + 3), t$l))
  paste0(s, sprintf('<text x="%d" y="%d" text-anchor="end">%s</text><text x="%d" y="%d">%s</text>', W - R, H - 3, xlab, L, T - 5, ylab))
}
sx <- function(v, lo, hi, W, L = 34, R = 8) L + (v - lo) / (hi - lo) * (W - L - R)
sy <- function(v, lo, hi, H, T = 16, B = 24) H - B - (v - lo) / (hi - lo) * (H - B - T)
dots <- function(x, y, col, r = 1.3, op = 1) paste0(sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s" fill-opacity="%s"/>', f1(x), f1(y), r, col, op), collapse = "")
W <- 300; H <- 210

# ---- A. Manhattan pair (naive over PC-adjusted), scenario B replicate 1 ----
g <- file.path(gw, "results/tasks/B_1"); cz <- fread(file.path(g, "causal.tsv"))
WA <- 620; HA <- 222; ymax <- 70; L <- 34; R <- 8
pan <- function(file, top, bot, lab) {
  d <- readRDS(file.path(g, file)); d[, y := pmin(-log10(P), ymax)]; set.seed(3)
  k <- rbind(d[y > 8][sample(.N, min(.N, 500))], d[y <= 8][sample(.N, min(.N, 1100))])
  xl <- c(15.5, 51)
  yy <- function(v) bot - (v / ymax) * (bot - top)
  s <- sprintf('<path d="M%d %s V%s H%d" stroke="%s" fill="none"/>', L, f1(top), f1(bot), WA - R, INK)
  s <- paste0(s, sprintf('<path d="M%d %s H%d" stroke="%s" stroke-dasharray="3 3"/>', L, f1(yy(-log10(5e-8))), WA - R, RUST))
  s <- paste0(s, dots(sx(k$POS / 1e6, xl[1], xl[2], WA), yy(k$y), ifelse(k$y > -log10(5e-8), RUST, MUTED), 1.2))
  paste0(s, sprintf('<text x="%d" y="%s">%s</text>', L + 4, f1(top + 8), lab))
}
body <- paste0(pan("scan_naive.rds", 18, 98, "unadjusted"), pan("scan_pc.rds", 110, 188, "adjusted for 10 PCs"))
body <- paste0(body, paste0(sprintf('<path d="M%s 188 v5" stroke="%s" stroke-width="1.4"/>', f1(sx(cz$POS / 1e6, 15.5, 51, WA)), GREEN), collapse = ""))
for (v in c(20, 30, 40, 50)) body <- paste0(body, sprintf('<text x="%s" y="%d" text-anchor="middle">%d</text>', f1(sx(v, 15.5, 51, WA)), 201, v))
body <- paste0(body, sprintf('<text x="%d" y="%d" text-anchor="end">chr22 position (Mb) &#183; green ticks = true causal variants</text>', WA - R, HA - 4))
writeLines(wrap(WA, HA, body, "Manhattan plots, unadjusted and adjusted", "Two stacked Manhattan plots of a simulated association scan; the unadjusted scan is significant across the chromosome, the adjusted one only near true causal variants."), "figures/panel-a-manhattan.svg")

# ---- C. QQ plot, recomputed for the same replicate ----
tmp <- tempfile(); dir.create(tmp)
run <- function(extra, out) system(sprintf("cd %s && bash -c 'source scripts/env.sh && plink2 --pfile results/qc/chr22_common --pheno results/tasks/B_1/pheno.tsv --pheno-name y %s --out %s' > /dev/null 2>&1", shQuote(gw), extra, file.path(tmp, out)))
run("--glm allow-no-covars", "n"); run("--covar results/qc/pcs.eigenvec --covar-name PC1-PC10 --covar-variance-standardize --glm hide-covar", "p")
qq <- function(f) { d <- fread(f); setnames(d, sub("^#", "", names(d))); p <- sort(d[TEST == "ADD" & !is.na(P)]$P); n <- length(p); list(e = -log10((seq_len(n) - 0.5) / n), o = -log10(p), lam = median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(.5, 1)) }
qn <- qq(file.path(tmp, "n.y.glm.linear")); qp <- qq(file.path(tmp, "p.y.glm.linear"))
thin <- function(q) { i <- unique(c(round(seq(1, length(q$e), length.out = 700)), seq_len(min(150, length(q$e))))); lapply(q[c("e", "o")], `[`, i) }
tn <- thin(qn); tp <- thin(qp); omax <- 60; emax <- max(qn$e)
body <- frame(W, H, xlab = "expected -log10 p", ylab = "observed", xt = lapply(c(0, 2, 4), function(v) list(x = sx(v, 0, emax, W), l = v)), yt = lapply(c(0, 20, 40, 60), function(v) list(y = sy(v, 0, omax, H), l = v)))
body <- paste0(body, sprintf('<path d="M%s %s L%s %s" stroke="%s" stroke-dasharray="3 3"/>', f1(sx(0, 0, emax, W)), f1(sy(0, 0, omax, H)), f1(sx(emax, 0, emax, W)), f1(sy(emax, 0, omax, H)), INK))
body <- paste0(body, dots(sx(tn$e, 0, emax, W), sy(pmin(tn$o, omax), 0, omax, H), RUST, 1.3), dots(sx(tp$e, 0, emax, W), sy(pmin(tp$o, omax), 0, omax, H), GREEN, 1.3),
  sprintf('<text x="44" y="30" fill="%s">&#955; = %.1f unadjusted</text><text x="44" y="41" fill="%s">&#955; = %.2f adjusted</text>', RUST, qn$lam, GREEN, qp$lam))
writeLines(wrap(W, H, body, "QQ plot", "Quantile-quantile plot of association p-values: the unadjusted scan departs strongly from the diagonal, the adjusted scan stays near it."), "figures/panel-c-qq.svg")

# ---- B. volcano, D. heatmap (airway) ----
data(airway); dds <- DESeqDataSet(airway, ~ cell + dex); dds$dex <- relevel(dds$dex, "untrt")
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]; dds <- DESeq(dds, quiet = TRUE)
res <- as.data.table(as.data.frame(results(dds, alpha = 0.05))); res[, id := rownames(dds)]; res <- res[!is.na(padj)]; res[, y := pmin(-log10(padj), 120)]
set.seed(4); k <- rbind(res[padj >= 0.05][sample(.N, 800)], res[padj < 0.05][sample(.N, 1100)])
body <- frame(W, H, xlab = "log2 fold change", ylab = "-log10 adj. p", xt = lapply(c(-5, 0, 5, 10), function(v) list(x = sx(v, -6, 11, W), l = v)), yt = lapply(c(0, 60, 120), function(v) list(y = sy(v, 0, 120, H), l = v)))
body <- paste0(body, dots(sx(k$log2FoldChange, -6, 11, W), sy(k$y, 0, 120, H), ifelse(k$padj < 0.05, RUST, MUTED), 1.2))
sym <- rowData(dds)$symbol; names(sym) <- rownames(dds); lab <- res[sym[id] %in% c("CRISPLD2", "DUSP1", "KLF15", "PER1", "TSC22D3")]
body <- paste0(body, dots(sx(lab$log2FoldChange, -6, 11, W), sy(lab$y, 0, 120, H), INK, 2.1), paste0(sprintf('<text x="%s" y="%s" fill="%s">%s</text>', f1(sx(lab$log2FoldChange, -6, 11, W) + 3.5), f1(sy(lab$y, 0, 120, H) - 3), INK, sym[lab$id]), collapse = ""))
writeLines(wrap(W, H, body, "Volcano plot", "Differential expression after dexamethasone; five genes named in the original paper are labelled."), "figures/panel-b-volcano.svg")

vsd <- vst(dds, blind = TRUE); top <- res[order(padj)][1:40]$id; m <- assay(vsd)[top, ]; m <- m - rowMeans(m)
m <- m[hclust(dist(m))$order, order(dds$dex, dds$cell)]; lim <- 2
cw <- (W - 44) / 8; ch <- (H - 54) / 40
col <- function(v) { v <- max(-lim, min(lim, v)) / lim; if (v >= 0) { a <- c(234, 227, 214); b <- c(140, 75, 52) } else { a <- c(234, 227, 214); b <- c(83, 98, 74) }; w <- abs(v); sprintf("#%02x%02x%02x", round(a[1] + (b[1] - a[1]) * w), round(a[2] + (b[2] - a[2]) * w), round(a[3] + (b[3] - a[3]) * w)) }
cells <- ""
for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m))) cells <- paste0(cells, sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="%s"/>', f1(36 + (j - 1) * cw), f1(30 + (i - 1) * ch), f1(cw + .3), f1(ch + .3), col(m[i, j])))
body <- paste0(cells, sprintf('<path d="M36 26 H%s" stroke="%s" stroke-width="2"/><path d="M%s 26 H%s" stroke="%s" stroke-width="2"/>', f1(36 + 4 * cw - 2), GREEN, f1(36 + 4 * cw + 2), f1(36 + 8 * cw), RUST),
  sprintf('<text x="%s" y="20" text-anchor="middle">untreated</text><text x="%s" y="20" text-anchor="middle">dexamethasone</text>', f1(36 + 2 * cw), f1(36 + 6 * cw)),
  sprintf('<text x="36" y="%d">top 40 genes &#183; 4 paired cell lines</text>', H - 4))
writeLines(wrap(W, H, body, "Expression heatmap", "Top 40 differentially expressed genes across eight samples; treated and untreated samples separate."), "figures/panel-d-heatmap.svg")

# ---- E. PCoA, F. 1000G PCA ----
data(GlobalPatterns); gp <- subset_samples(GlobalPatterns, SampleType != "Mock"); gp <- prune_taxa(taxa_sums(gp) > 0, gp)
rf <- rarefy_even_depth(gp, sample.size = min(sample_sums(gp)), rngseed = 1, verbose = FALSE)
p <- cmdscale(phyloseq::distance(rf, "bray"), k = 2, eig = TRUE); host <- as.character(sample_data(rf)$SampleType) %in% c("Feces", "Skin", "Tongue")
xr <- range(p$points[, 1]) + c(-.05, .05); yr <- range(p$points[, 2]) + c(-.05, .05)
body <- frame(W, H, xlab = "PCo1", ylab = "PCo2") 
px <- sx(p$points[, 1], xr[1], xr[2], W); py <- sy(p$points[, 2], yr[1], yr[2], H)
body <- paste0(body, dots(px[host], py[host], RUST, 3.6, .95),
  paste0(sprintf('<circle cx="%s" cy="%s" r="3.2" fill="none" stroke="%s" stroke-width="1.5"/>', f1(px[!host]), f1(py[!host]), GREEN), collapse = ""),
  sprintf('<circle cx="%d" cy="25" r="3.4" fill="%s"/><text x="%d" y="28" fill="%s">host-associated</text><circle cx="%d" cy="37" r="3.2" fill="none" stroke="%s" stroke-width="1.5"/><text x="%d" y="40" fill="%s">free-living</text>', W - 108, RUST, W - 100, RUST, W - 108, GREEN, W - 100, GREEN))
writeLines(wrap(W, H, body, "Ordination of microbial communities", "Bray-Curtis PCoA of 23 communities; host-associated and free-living samples separate."), "figures/panel-e-pcoa.svg")

pc <- fread(file.path(gw, "results/qc/pcs.eigenvec")); smp <- fread(file.path(gw, "results/qc/samples.tsv")); pc <- merge(pc, smp, by.x = "#IID", by.y = "IID"); set.seed(8); pc <- pc[sample(.N, 900)]
cc <- c(AFR = RUST, AMR = AMBER, EAS = GREEN, EUR = BROWN, SAS = SLATE)
body <- frame(W, H, xlab = "PC1", ylab = "PC2")
body <- paste0(body, dots(sx(pc$PC1, min(pc$PC1) - .003, max(pc$PC1) + .003, W), sy(pc$PC2, min(pc$PC2) - .003, max(pc$PC2) + .003, H), cc[pc$super_pop], 1.4, .8),
  paste0(sprintf('<rect x="%d" y="%d" width="6" height="6" fill="%s"/><text x="%d" y="%d">%s</text>', W - 62, 22 + 10 * (seq_along(cc) - 1), cc, W - 53, 28 + 10 * (seq_along(cc) - 1), names(cc)), collapse = ""))
writeLines(wrap(W, H, body, "Population structure in 1000 Genomes", "First two principal components of 2,504 genotypes, coloured by continental group."), "figures/panel-f-pca.svg")
cat(paste(list.files("figures", pattern = "^panel-"), collapse = ", "), "\n")
