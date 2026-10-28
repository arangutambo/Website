# Builds the small SVG figures used on projects.html from the real analysis outputs.
# Run from the Website folder:  Rscript tools/build_figures.R
# Needs the analysis repos next to this one (~/Developer/gwas-1000g-simulation, etc.).
suppressPackageStartupMessages({library(data.table)})
dev <- normalizePath("~/Developer")
INK <- "#2A2521"; MUTED <- "#5C5448"; RUST <- "#8C4B34"; GREEN <- "#53624A"; PAPER <- "#EAE3D6"
W <- 600; H <- 250; L <- 36; R <- 12; T <- 14; B <- 26
sx <- function(v, lo, hi) L + (v - lo) / (hi - lo) * (W - L - R)
sy <- function(v, lo, hi) H - B - (v - lo) / (hi - lo) * (H - B - T)
f1 <- function(x) formatC(x, format = "f", digits = 1)

svg_wrap <- function(body, title, desc) sprintf(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" role="img" aria-labelledby="t d"><title id="t">%s</title><desc id="d">%s</desc><g font-family="ui-monospace,Menlo,Consolas,monospace" font-size="9" fill="%s">%s</g></svg>',
  W, H, title, desc, MUTED, body)
axes <- function(xlab, ylab, xt = NULL, yt = NULL) {
  s <- sprintf('<path d="M%d %d H%d M%d %d V%d" stroke="%s" stroke-width="1" fill="none"/>', L, H - B, W - R, L, T, H - B, INK)
  for (t in xt) s <- paste0(s, sprintf('<text x="%s" y="%d" text-anchor="middle">%s</text>', f1(t$x), H - B + 12, t$l))
  for (t in yt) s <- paste0(s, sprintf('<text x="%d" y="%s" text-anchor="end">%s</text>', L - 4, f1(t$y + 3), t$l))
  paste0(s, sprintf('<text x="%d" y="%d" text-anchor="end">%s</text>', W - R, H - 4, xlab),
         sprintf('<text x="%d" y="%d">%s</text>', L, T - 3, ylab))
}
dots <- function(x, y, col, r = 1.35) paste0(sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s"/>', f1(x), f1(y), r, col), collapse = "")

# ---- 1. GWAS: Manhattan, scenario B replicate 1, naive vs PC-adjusted (same y scale) ----
g <- file.path(dev, "gwas-1000g-simulation/results/tasks/B_1")
cz <- fread(file.path(g, "causal.tsv"))
ymax <- 75
manh <- function(file, caption) {
  d <- readRDS(file.path(g, file)); d[, y := pmin(-log10(P), ymax)]
  set.seed(3); hi <- d[y > 8]; lo <- d[y <= 8]
  keep <- rbind(hi[sample(.N, min(.N, 700))], lo[sample(.N, min(.N, 1500))])
  xl <- c(min(d$POS), max(d$POS)) / 1e6
  body <- axes("chr22 position (Mb)", "-log10 p", xt = lapply(c(20, 30, 40, 50), function(v) list(x = sx(v, xl[1], xl[2]), l = v)),
               yt = lapply(c(0, 25, 50, 75), function(v) list(y = sy(v, 0, ymax), l = v)))
  body <- paste0(body, sprintf('<path d="M%d %s H%d" stroke="%s" stroke-dasharray="3 3" stroke-width="1"/>', L, f1(sy(-log10(5e-8), 0, ymax)), W - R, RUST))
  body <- paste0(body, dots(sx(keep$POS / 1e6, xl[1], xl[2]), sy(keep$y, 0, ymax), ifelse(keep$y > -log10(5e-8), RUST, MUTED)))
  body <- paste0(body, paste0(sprintf('<path d="M%s %d v6" stroke="%s" stroke-width="1.4"/>', f1(sx(cz$POS / 1e6, xl[1], xl[2])), H - B, GREEN), collapse = ""))
  svg_wrap(body, caption, "Manhattan plot from a simulation on 1000 Genomes chromosome 22. Green ticks mark true causal variants; dashed line is p = 5e-8.")
}
writeLines(manh("scan_naive.rds", "Unadjusted GWAS scan"), "figures/gwas-naive.svg")
writeLines(manh("scan_pc.rds", "GWAS scan adjusted for 10 principal components"), "figures/gwas-pc.svg")

# ---- 2. RNA-seq volcano (airway, DESeq2, paired design) ----
suppressPackageStartupMessages({library(DESeq2); library(airway)})
data(airway); dds <- DESeqDataSet(airway, ~ cell + dex); dds$dex <- relevel(dds$dex, "untrt")
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]; dds <- DESeq(dds, quiet = TRUE); res <- as.data.table(as.data.frame(results(dds, alpha = 0.05)))
res[, sym := rowData(dds)$symbol]; res <- res[!is.na(padj)]; res[, y := pmin(-log10(padj), 120)]
set.seed(4); ns <- res[padj >= 0.05]; sg <- res[padj < 0.05]
keep <- rbind(ns[sample(.N, min(.N, 900))], sg[sample(.N, min(.N, 1300))])
body <- axes("log2 fold change", "-log10 adj. p", xt = lapply(c(-5, 0, 5, 10), function(v) list(x = sx(v, -6, 11), l = v)),
             yt = lapply(c(0, 40, 80, 120), function(v) list(y = sy(v, 0, 120), l = v)))
body <- paste0(body, dots(sx(keep$log2FoldChange, -6, 11), sy(keep$y, 0, 120), ifelse(keep$padj < 0.05, RUST, MUTED)))
lab <- res[sym %in% c("CRISPLD2", "DUSP1", "KLF15", "PER1", "TSC22D3")]
body <- paste0(body, dots(sx(lab$log2FoldChange, -6, 11), sy(lab$y, 0, 120), INK, 2.2),
  paste0(sprintf('<text x="%s" y="%s" fill="%s">%s</text>', f1(sx(lab$log2FoldChange, -6, 11) + 4), f1(sy(lab$y, 0, 120) - 3), INK, lab$sym), collapse = ""))
writeLines(svg_wrap(body, "Volcano plot of dexamethasone response in airway smooth muscle", "Each point is a gene. Rust points have adjusted p below 0.05. Five genes named in the original paper are labelled."), "figures/rnaseq-volcano.svg")

# ---- 3. Microbiome: Bray-Curtis PCoA coloured by habitat ----
suppressPackageStartupMessages({library(phyloseq)})
data(GlobalPatterns); gp <- subset_samples(GlobalPatterns, SampleType != "Mock"); gp <- prune_taxa(taxa_sums(gp) > 0, gp)
rf <- rarefy_even_depth(gp, sample.size = min(sample_sums(gp)), rngseed = 1, verbose = FALSE)
p <- cmdscale(phyloseq::distance(rf, "bray"), k = 2, eig = TRUE)
st <- as.character(sample_data(rf)$SampleType); host <- st %in% c("Feces", "Skin", "Tongue")
xr <- range(p$points[, 1]) + c(-0.05, 0.05); yr <- range(p$points[, 2]) + c(-0.05, 0.05)
ve <- round(100 * p$eig[1:2] / sum(p$eig[p$eig > 0]))
body <- axes(sprintf("PCo1 (%d%%)", ve[1]), sprintf("PCo2 (%d%%)", ve[2]))
body <- paste0(body, dots(sx(p$points[, 1], xr[1], xr[2]), sy(p$points[, 2], yr[1], yr[2]), ifelse(host, RUST, GREEN), 4.2),
  sprintf('<text x="%d" y="%d" fill="%s">host-associated</text><text x="%d" y="%d" fill="%s">free-living</text>', W - 110, T + 6, RUST, W - 110, T + 18, GREEN))
writeLines(svg_wrap(body, "Ordination of microbial communities", "Bray-Curtis principal coordinates of 23 samples; rust is host-associated, green is free-living."), "figures/microbiome-pcoa.svg")
cat("wrote", paste(list.files("figures"), collapse = ", "), "\n")

# ---- 4. Hero tree: UPGMA dendrogram of the 23 Global Patterns communities (Bray-Curtis) ----
hc <- hclust(phyloseq::distance(rf, "bray"), "average")
n <- length(hc$order); ypos <- numeric(n); ypos[hc$order] <- seq_len(n)
node_x <- numeric(n - 1); node_y <- numeric(n - 1)
pos <- function(i) if (i < 0) c(0, ypos[-i]) else c(node_x[i], node_y[i])
segs <- character(0)
tw <- 300; th <- 330; x0 <- 6; y0 <- 10; maxh <- max(hc$height)
X <- function(h) x0 + (1 - h / maxh) * (tw - x0) * 0.82  # root at left, tips at right
Y <- function(y) y0 + (y - 1) / (n - 1) * (th - 2 * y0)
for (k in seq_len(n - 1)) {
  a <- hc$merge[k, 1]; b <- hc$merge[k, 2]; h <- hc$height[k]
  pa <- pos(a); pb <- pos(b); node_y[k] <- (pa[2] + pb[2]) / 2; node_x[k] <- h
  xa <- if (a < 0) X(0) else X(node_x[a]); xb <- if (b < 0) X(0) else X(node_x[b]); xk <- X(h)
  segs <- c(segs, sprintf("M%s %sH%s", f1(xa), f1(Y(pa[2])), f1(xk)), sprintf("M%s %sH%s", f1(xb), f1(Y(pb[2])), f1(xk)),
            sprintf("M%s %sV%s", f1(xk), f1(Y(pa[2])), f1(Y(pb[2]))))
}
lab <- st[hc$order]; col <- ifelse(host[hc$order], RUST, GREEN)
tips <- paste0(sprintf('<rect x="%s" y="%s" width="5" height="5" fill="%s"/><text x="%s" y="%s" fill="%s">%s</text>',
  f1(X(0) + 2), f1(Y(seq_len(n)) - 2.5), col, f1(X(0) + 11), f1(Y(seq_len(n)) + 3), col, tolower(sub(" \\(.*\\)", "", lab))), collapse = "")
writeLines(sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" role="img" aria-label="Dendrogram of 23 microbial communities from the Global Patterns dataset: host-associated samples (rust) cluster apart from free-living samples (green)"><path d="%s" fill="none" stroke="%s" stroke-width="1.2" stroke-linecap="square"/><g font-family="ui-monospace,Menlo,Consolas,monospace" font-size="8.5">%s</g></svg>',
  tw, th, paste(segs, collapse = ""), INK, tips), "figures/hero-tree.svg")
cat("tree leaves:", n, "\n")

# ---- 5. Skyline: the PC-adjusted Manhattan scan as a thin divider under the hero ----
d <- readRDS(file.path(g, "scan_pc.rds")); d[, y := pmin(-log10(P), 55)]
set.seed(5); keep <- rbind(d[y > 3], d[y <= 3][sample(.N, min(.N, 1400))])
SW <- 1200; SH <- 64
xs <- (keep$POS - min(d$POS)) / (max(d$POS) - min(d$POS)) * (SW - 8) + 4
ys <- SH - 4 - keep$y / 55 * (SH - 10)
blk <- floor(keep$POS / 5e6) %% 2
colr <- ifelse(keep$y > -log10(5e-8), RUST, ifelse(blk == 0, "#5C4736", GREEN))
sky <- sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" preserveAspectRatio="none" role="img" aria-label="Manhattan plot of a simulated association scan on chromosome 22"><path d="M0 %s H%d" stroke="%s" stroke-width="1" stroke-dasharray="4 4" vector-effect="non-scaling-stroke"/>%s</svg>',
  SW, SH, f1(SH - 4 - (-log10(5e-8)) / 55 * (SH - 10)), SW, RUST,
  paste0(sprintf('<rect x="%s" y="%s" width="2.2" height="2.2" fill="%s"/>', f1(xs), f1(ys), colr), collapse = ""))
writeLines(sky, "figures/skyline.svg")
cat("skyline points:", nrow(keep), "\n")

# ---- 6. Small plot motifs for the project header bars (cream dots, drawn on the coloured bar) ----
CREAM <- "#EAE3D6"; TW <- 170; TH <- 36
motif <- function(x, y, fill, op, r = 1.1, hi = NULL) {
  x <- (x - min(x)) / diff(range(x)) * (TW - 12) + 6; y <- TH - 5 - (y - min(y)) / diff(range(y)) * (TH - 10)
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" aria-hidden="true">%s</svg>', TW, TH,
    paste0(sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s" fill-opacity="%s"/>', f1(x), f1(y), r, fill, op), collapse = ""))
}
pc <- fread(file.path(dev, "gwas-1000g-simulation/results/qc/pcs.eigenvec")); set.seed(6); pc <- pc[sample(.N, 650)]
writeLines(motif(pc$PC1, pc$PC2, CREAM, 0.7), "figures/tag-gwas.svg")
set.seed(7); r2 <- res[sample(.N, 700)]
writeLines(motif(r2$log2FoldChange, r2$y, CREAM, ifelse(r2$padj < 0.05, 0.85, 0.35)), "figures/tag-rnaseq.svg")
writeLines(motif(p$points[, 1], p$points[, 2], CREAM, ifelse(host, 0.95, 0.5), r = 2.4), "figures/tag-microbiome.svg")
cat("tag motifs written\n")
