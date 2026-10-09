# =============================================================================
# Complete reanalysis of the tomato root metabolome
# Metarhizium anisopliae (MA) x Meloidogyne incognita (MI)
# 2 x 2 factorial, completely randomised, 6 plants per treatment (N = 24)
#
# Analyses (all with the threshold FDR <= FDR_LIMIAR; results with
# FDR <= FDR_ESTRITO are flagged separately, so that the text can grade the
# strength of each claim):
#   1. Two-way ANOVA per metabolite (log2; type III SS; NA omitted)
#   2. limma: cell-means model + factorial contrasts (empirical Bayes)
#   3. Two-way PERMANOVA (global profile) + PERMDISP
#   4. Exact pairwise PERMANOVA, only for the 4 simple-effect comparisons
#   5. Distances between centroids with bootstrap intervals
#   6. PERMANOVA by chemical class (BH across classes, within each term)
#   7. Sensitivity: without MA plants 10 and 12 (ANOVA, limma, PERMANOVA)
#   8. Summary table of findings and figures for the manuscript
#
# Input   : MetabolomicData.csv (same folder; or change ARQUIVO)
# Outputs : folder DIR_SAIDA with the spreadsheet, figures (vector PDF and
#           600 dpi PNG) and sessionInfo.txt
# Packages: car, vegan, ggplot2, ggrepel, scales, writexl (CRAN); limma
#           (Bioconductor; with conda: conda install -c bioconda bioconductor-limma)
# =============================================================================


# ---- 0. Configuration -------------------------------------------------------
ARQUIVO                <- "MetabolomicData.csv"
DIR_SAIDA              <- "resultados_FDR010"
FDR_LIMIAR             <- 0.10     # requested "significant" threshold
FDR_ESTRITO            <- 0.05     # flagged separately in tables and figures
METODO_PRINCIPAL       <- "limma"  # "limma" or "anova": defines figures and summary
MIN_DET                <- 3        # minimum number of plants with a value in EACH treatment
EXCLUIR_NAO_BIOLOGICOS <- TRUE     # TRUE = 46 metabolites; FALSE = 49
SOMAR_PICOS_DUPLICADOS <- FALSE    # TRUE sums peaks of the same compound (see below)
PICOS_DUPLICADOS       <- list(Fructose = c("Fructofuranose", "Fructopyranose"),
                               Glucose  = c("Glucose", "Glucopyranose"))
PLANTAS_SENSIBILIDADE  <- c(10, 12)
ROTULAR_PLANTAS_PCA    <- TRUE
N_PERM                 <- 9999
N_BOOT                 <- 9999
SEMENTE                <- 20260930

GRUPOS  <- c("Controle", "MA", "MI", "MA_MI")
ROTULOS <- c(Controle = "Control", MA = "MA", MI = "MI", MA_MI = "MA + MI")

# categorical palette checked for colour-vision deficiency (all pairs); the
# aqua colour has low contrast against a white background, so each treatment
# also has its own shape
CORES              <- c("Control" = "#2a78d6", "MA" = "#1baf7a", "MI" = "#eb6834", "MA + MI" = "#4a3aa7")
FORMAS             <- c("Control" = 16, "MA" = 17, "MI" = 15, "MA + MI" = 18)
FORMAS_PREENCHIDAS <- c("Control" = 21, "MA" = 24, "MI" = 22, "MA + MI" = 23)

PROVAVEL_NAO_BIOLOGICO <- c(
  "Phthalic acid, di(2-propylpentyl) ester", "Phthalic acid, diisobutyl ester",
  "2-Ethylhexanoic acid", "Bisphenol A monomethyl ether",
  "4-octylphenyl 4-octylbenzoate", "2,2,4-Trimethyl-1,3-pentanediol diisobutyrate",
  "Diethylene glycol", "4-tert-Butylphenol", "2,7,10-Trimethyldodecane",
  "4,6-Dimethyldodecane", "Phosphoric acid, bis(trimethylsilyl)monomethyl ester",
  "Alanine, N-methyl-N-ethoxycarbonyl-, dodecyl ester")

# CHEMICAL CLASSES: fix this list by chemical criteria BEFORE running and do
# not change it after seeing the results. Class names appear in the figures
# (hence in English). Tested metabolites outside these classes are not
# included in this part of the analysis (they appear as "Other" in the heatmap).
CLASSES <- list(
  "Amino acids"         = c("Asparagine", "Glutamic acid", "Isoleucine", "Leucine",
                            "Phenylalanine", "Proline", "Serine", "Threonine", "Tyrosine",
                            "Valine", "4-Aminobutanoic acid", "3-Amino-2-piperidone"),
  "Free sugars"         = c("Fructofuranose", "Fructopyranose", "Glucopyranose", "Glucose",
                            "Sucrose", "Trehalose", "Melibiose", "Xylose", "Talofuranose"),
  "Polyols"             = c("Glucitol", "Mannitol", "Inositol", "Galactinol"),
  "Sugar phosphates"    = c("Glucose-6-Phosphate", "Mannose 6-phosphate", "Inositol-1- phosphate"),
  "Sugar acids"         = c("Erythronic acid", "Galactonic acid", "Galacturonic acid",
                            "Gluconic acid", "Glyceric acid"),
  "TCA intermediates"   = c("2-Butenedioic acid", "2-Ketoglutaric acid", "Malic acid"),
  "Amines and alkaloids" = c("Putrescine", "Ethanolamine", "Tropine")
)

pacotes_cran <- c("car", "vegan", "ggplot2", "ggrepel", "scales", "writexl")
faltando <- pacotes_cran[!vapply(pacotes_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltando) > 0) install.packages(faltando)
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("limma")
}
suppressPackageStartupMessages({
  library(car); library(vegan); library(limma)
  library(ggplot2); library(ggrepel); library(writexl)
})

DIR_FIG <- file.path(DIR_SAIDA, "figuras")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# symbols: with cairo, PDF and PNG show "≤" and "×"; without cairo, ASCII is used
USA_CAIRO <- isTRUE(capabilities("cairo"))
SIMB_LE   <- if (USA_CAIRO) "\u2264" else "<="
SIMB_X    <- if (USA_CAIRO) "\u00d7" else "x"
DISP_PDF  <- if (USA_CAIRO) grDevices::cairo_pdf else grDevices::pdf
ROTULO_TERMO <- c(MA = "MA", MI = "MI", MAxMI = paste("MA", SIMB_X, "MI"))
CORES_TERMO  <- setNames(c("#1baf7a", "#eb6834", "#4a3aa7"), ROTULO_TERMO)
TERMOS       <- c("MA", "MI", "MAxMI")
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)


# ---- 1. Data input ----------------------------------------------------------
if (!file.exists(ARQUIVO)) {
  stop("Arquivo '", ARQUIVO, "' não encontrado em ", getwd(),
       ". Ajuste o diretório de trabalho (Session > Set Working Directory).")
}
bruto <- read.table(ARQUIVO, sep = ";", header = FALSE, colClasses = "character",
                    quote = "", comment.char = "", strip.white = TRUE,
                    na.strings = c("NA", ""))
amostra <- as.integer(unlist(bruto[1, -1]))
grupo   <- trimws(unlist(bruto[2, -1]))
nomes   <- trimws(bruto[-(1:2), 1])
stopifnot(!anyDuplicated(nomes))
txt <- as.matrix(bruto[-(1:2), -1])
X <- suppressWarnings(matrix(as.numeric(gsub(",", ".", txt)), nrow = nrow(txt)))
if (sum(!is.na(txt) & is.na(X)) > 0) warning("Há valores não numéricos convertidos em NA.")
dimnames(X) <- list(nomes, amostra)
ordem <- order(amostra)
X <- X[, ordem]; amostra <- amostra[ordem]; grupo <- grupo[ordem]
grupo[grupo == "Teste"] <- "Controle"

# optional summing of peaks of the same compound (GC-MS derivatives/isomers).
# When only one of the peaks has a value in a plant, the sum uses only that
# peak, which underestimates the total; the option is therefore off by default
# and must be decided beforehand.
if (SOMAR_PICOS_DUPLICADOS) {
  for (novo in names(PICOS_DUPLICADOS)) {
    pk <- intersect(PICOS_DUPLICADOS[[novo]], rownames(X))
    if (length(pk) < 2) next
    soma <- colSums(X[pk, , drop = FALSE], na.rm = TRUE)
    soma[colSums(!is.na(X[pk, , drop = FALSE])) == 0] <- NA
    X <- rbind(X[setdiff(rownames(X), pk), , drop = FALSE],
               matrix(soma, nrow = 1, dimnames = list(novo, colnames(X))))
    CLASSES <- lapply(CLASSES, function(v) unique(ifelse(v %in% pk, novo, v)))
  }
}

meta <- data.frame(
  amostra = amostra,
  grupo   = factor(grupo, levels = GRUPOS),
  MA      = factor(ifelse(grupo %in% c("MA", "MA_MI"), "With", "Without"), levels = c("Without", "With")),
  MI      = factor(ifelse(grupo %in% c("MI", "MA_MI"), "With", "Without"), levels = c("Without", "With")))
meta$trat <- factor(unname(ROTULOS[as.character(meta$grupo)]), levels = unname(ROTULOS))
stopifnot(!any(is.na(meta$grupo)))
print(table(MA = meta$MA, MI = meta$MI))

det <- sapply(GRUPOS, function(g) rowSums(!is.na(X[, meta$grupo == g, drop = FALSE])))
testados <- rownames(det)[apply(det >= MIN_DET, 1, all)]
if (EXCLUIR_NAO_BIOLOGICOS) testados <- setdiff(testados, PROVAVEL_NAO_BIOLOGICO)
cat("\nMetabólitos na tabela:", nrow(X), "| testados:", length(testados), "\n")

nivel_fdr <- function(q) {
  n1 <- paste0("FDR ", SIMB_LE, " ", fmt(FDR_ESTRITO))
  n2 <- paste0(fmt(FDR_ESTRITO), " < FDR ", SIMB_LE, " ", fmt(FDR_LIMIAR))
  factor(ifelse(is.na(q) | q > FDR_LIMIAR, "n.s.", ifelse(q <= FDR_ESTRITO, n1, n2)),
         levels = c(n1, n2, "n.s."))
}


# ---- 2. Analysis functions ---------------------------------------------------
# 2.1 Two-way ANOVA per metabolite: log2(y) = mu + MA + MI + MA:MI + error
#     Type III SS with contr.sum; NA omitted; logFC = contrasts of cell means
#     (the same as in limma)
anova_2vias <- function(X, meta, metabolitos) {
  res <- do.call(rbind, lapply(metabolitos, function(m) {
    d <- data.frame(y = log2(X[m, ]), meta)
    d <- d[!is.na(d$y), ]
    fit <- lm(y ~ MA * MI, data = d, contrasts = list(MA = "contr.sum", MI = "contr.sum"))
    tab <- car::Anova(fit, type = 3)
    mu  <- tapply(d$y, d$grupo, mean)
    lfc <- c((mu[["MA"]] + mu[["MA_MI"]] - mu[["Controle"]] - mu[["MI"]]) / 2,
             (mu[["MI"]] + mu[["MA_MI"]] - mu[["Controle"]] - mu[["MA"]]) / 2,
             (mu[["MA_MI"]] - mu[["MI"]]) - (mu[["MA"]] - mu[["Controle"]]))
    linha  <- c("MA", "MI", "MA:MI")
    ss_res <- tab["Residuals", "Sum Sq"]
    data.frame(metabolito = m, termo = TERMOS, logFC = lfc,
               estatistica = tab[linha, "F value"], P = tab[linha, "Pr(>F)"],
               eta2p = tab[linha, "Sum Sq"] / (tab[linha, "Sum Sq"] + ss_res),
               n = nrow(d), gl_residuo = tab["Residuals", "Df"],
               Shapiro_P = shapiro.test(residuals(fit))$p.value,
               BrownForsythe_P = car::leveneTest(y ~ grupo, data = d, center = median)[1, "Pr(>F)"])
  }))
  res$FDR <- ave(res$P, res$termo, FUN = function(p) p.adjust(p, method = "BH"))
  res$nivel <- nivel_fdr(res$FDR)
  res
}

# 2.2 limma: cell means (~ 0 + grupo) and factorial contrasts; NA allowed
limma_fatorial <- function(X, meta, metabolitos) {
  Y <- log2(X[metabolitos, , drop = FALSE])
  design <- model.matrix(~ 0 + grupo, data = meta)
  colnames(design) <- levels(meta$grupo)
  cm <- makeContrasts(
    MA                = (MA + MA_MI - Controle - MI) / 2,
    MI                = (MI + MA_MI - Controle - MA) / 2,
    MAxMI             = (MA_MI - MI) - (MA - Controle),
    MI_vs_Controle    = MI - Controle,
    MA_MI_vs_MA       = MA_MI - MA,
    MA_vs_Controle    = MA - Controle,
    MA_MI_vs_MI       = MA_MI - MI,
    levels = design)
  fit2 <- eBayes(contrasts.fit(lmFit(Y, design), cm))
  res <- do.call(rbind, lapply(colnames(cm), function(k) {
    tt <- topTable(fit2, coef = k, number = Inf, sort.by = "none", adjust.method = "BH")
    data.frame(metabolito = rownames(tt), termo = k, logFC = tt$logFC,
               estatistica = tt$t, P = tt$P.Value, FDR = tt$adj.P.Val,
               n = rowSums(!is.na(Y[rownames(tt), , drop = FALSE])))
  }))
  res$nivel <- nivel_fdr(res$FDR)
  list(tabela = res, prior_gl = fit2$df.prior, prior_s2 = fit2$s2.prior)
}

# 2.3 multivariate matrix: NA -> half of the minimum (only here), log2, centring
#     and Pareto scaling per metabolite
preparar_matriz <- function(X, metabolitos) {
  M <- X[metabolitos, , drop = FALSE]
  M <- t(apply(M, 1, function(v) { v[is.na(v)] <- min(v, na.rm = TRUE) / 2; v }))
  L <- t(log2(M))
  L <- scale(L, center = TRUE, scale = FALSE)
  sweep(L, 2, sqrt(apply(L, 2, sd)), "/")
}

# 2.4 Two-way PERMANOVA, Euclidean distance, type III SS (-1/+1 coding),
#     unrestricted permutation of plants
permanova_t3 <- function(Z, meta, n_perm = N_PERM, semente = SEMENTE) {
  Z <- as.matrix(Z)
  a <- ifelse(meta$MA == "With", 1, -1)
  b <- ifelse(meta$MI == "With", 1, -1)
  D <- cbind(1, a, b, a * b)
  q_comp <- qr(D)
  q_red  <- lapply(2:4, function(j) qr(D[, -j, drop = FALSE]))
  gl_res <- nrow(Z) - 4
  estat <- function(Y) {
    sqr <- sum(qr.resid(q_comp, Y)^2)
    sq  <- vapply(q_red, function(q) sum(qr.resid(q, Y)^2), numeric(1)) - sqr
    list(sq = sq, F = sq / (sqr / gl_res), sqr = sqr)
  }
  obs <- estat(Z)
  set.seed(semente)
  cont <- numeric(3)
  for (i in seq_len(n_perm)) {
    p <- estat(Z[sample(nrow(Z)), , drop = FALSE])
    cont <- cont + (p$F >= obs$F - 1e-12)
  }
  sq_tot <- sum(scale(Z, scale = FALSE)^2)
  data.frame(termo    = c(TERMOS, "Residuo", "Total"),
             gl       = c(1, 1, 1, gl_res, nrow(Z) - 1),
             SQ       = c(obs$sq, obs$sqr, sq_tot),
             R2       = c(obs$sq, obs$sqr, sq_tot) / sq_tot,
             pseudo_F = c(obs$F, NA, NA),
             P        = c((cont + 1) / (n_perm + 1), NA, NA))
}

# 2.5 Exact pairwise PERMANOVA (all labellings), simple effects only
PARES <- list(c("Controle", "MI"), c("MA", "MA_MI"), c("Controle", "MA"), c("MI", "MA_MI"))
DESCR_PARES <- c("MI effect without MA", "MI effect with MA",
                 "MA effect without MI", "MA effect with MI")
pares_exatos <- function(Z, grupo) {
  res <- do.call(rbind, lapply(seq_along(PARES), function(k) {
    a <- PARES[[k]][1]; b <- PARES[[k]][2]
    sel <- grupo %in% c(a, b)
    Y   <- Z[sel, , drop = FALSE]
    rot <- as.character(grupo[sel]) == a
    n <- nrow(Y); na_ <- sum(rot)
    sq_tot <- sum(scale(Y, scale = FALSE)^2)
    estat <- function(m) {
      sq_w <- sum(scale(Y[m, , drop = FALSE], scale = FALSE)^2) +
              sum(scale(Y[!m, , drop = FALSE], scale = FALSE)^2)
      c(F = (sq_tot - sq_w) / (sq_w / (n - 2)), R2 = (sq_tot - sq_w) / sq_tot)
    }
    obs  <- estat(rot)
    comb <- combn(n, na_)
    Fs <- apply(comb, 2, function(idx) { m <- rep(FALSE, n); m[idx] <- TRUE; estat(m)[["F"]] })
    data.frame(comparacao = paste(ROTULOS[a], "vs", ROTULOS[b]), efeito = DESCR_PARES[k],
               n = paste0(na_, "+", n - na_), pseudo_F = obs[["F"]], R2 = obs[["R2"]],
               P_exato = mean(Fs >= obs[["F"]] - 1e-12), rotulagens = ncol(comb))
  }))
  res$FDR <- p.adjust(res$P_exato, method = "BH")
  res$nivel <- nivel_fdr(res$FDR)
  res
}

# 2.6 distances between centroids with bootstrap (resampling within treatment)
centroides <- function(Z, grupo, n_boot = N_BOOT, semente = SEMENTE) {
  idx <- split(seq_len(nrow(Z)), grupo)
  cent_de <- function(lista) t(sapply(lista, function(ii) colMeans(Z[ii, , drop = FALSE])))
  D_obs <- as.matrix(dist(cent_de(idx)))[GRUPOS, GRUPOS]
  pares <- combn(GRUPOS, 2)
  ij <- cbind(pares[1, ], pares[2, ])
  set.seed(semente)
  boot <- replicate(n_boot, {
    ib <- lapply(idx, function(ii) ii[sample.int(length(ii), length(ii), replace = TRUE)])
    Db <- as.matrix(dist(cent_de(ib)))
    c(Db[ij], Db["MA_MI", "Controle"] - Db["MA_MI", "MI"])
  })
  ic <- apply(boot, 1, quantile, probs = c(0.025, 0.975))
  k <- ncol(pares)
  list(pares = data.frame(par = paste(ROTULOS[pares[1, ]], "vs", ROTULOS[pares[2, ]]),
                          g1 = pares[1, ], g2 = pares[2, ], d = D_obs[ij],
                          IC_inf = ic[1, 1:k], IC_sup = ic[2, 1:k]),
       dif = D_obs["MA_MI", "Controle"] - D_obs["MA_MI", "MI"],
       dif_ic = ic[, k + 1])
}

# 2.7 PERMANOVA by chemical class, BH across classes within each term
permanova_classes <- function(Z, meta, classes) {
  res <- do.call(rbind, lapply(names(classes), function(cl) {
    ms <- intersect(classes[[cl]], colnames(Z))
    if (length(ms) < 2) return(NULL)
    pm <- permanova_t3(Z[, ms, drop = FALSE], meta)
    data.frame(classe = cl, n_metabolitos = length(ms), pm[1:3, c("termo", "SQ", "R2", "pseudo_F", "P")])
  }))
  res$FDR <- ave(res$P, res$termo, FUN = function(p) p.adjust(p, method = "BH"))
  res$nivel <- nivel_fdr(res$FDR)
  res
}


# ---- 3. Main analyses (24 plants) -------------------------------------------
cat("\nRodando as análises (alguns minutos por causa das permutações)...\n")
res_anova <- anova_2vias(X, meta, testados)
lim       <- limma_fatorial(X, meta, testados)
res_limma <- lim$tabela

Z <- preparar_matriz(X, testados)
res_perm <- permanova_t3(Z, meta)
set.seed(SEMENTE)
res_adonis <- adonis2(dist(Z) ~ MA * MI, data = meta, permutations = N_PERM, by = "terms")
sq_confere <- isTRUE(all.equal(unname(res_adonis[1:3, 2]), res_perm$SQ[1:3], tolerance = 1e-6))

bd <- betadisper(dist(Z), meta$grupo, type = "centroid")
F_disp <- anova(bd)[1, "F value"]
set.seed(SEMENTE)
perm_disp <- permutest(bd, permutations = N_PERM)
P_disp <- tryCatch(perm_disp$tab[1, "Pr(>F)"], error = function(e) NA)

res_pares <- pares_exatos(Z, meta$grupo)
res_cent  <- centroides(Z, meta$grupo)

colnames(Z) <- testados
res_classes <- permanova_classes(Z, meta, CLASSES)


# ---- 4. Sensitivity: without the specified plants ---------------------------
manter  <- !(meta$amostra %in% PLANTAS_SENSIBILIDADE)
X_s     <- X[, manter]
meta_s  <- meta[manter, ]
res_anova_s <- anova_2vias(X_s, meta_s, testados)
res_limma_s <- limma_fatorial(X_s, meta_s, testados)$tabela
Z_s         <- preparar_matriz(X_s, testados)
res_perm_s  <- permanova_t3(Z_s, meta_s)
res_cent_s  <- centroides(Z_s, meta_s$grupo)


# ---- 5. Summary of findings --------------------------------------------------
pegar <- function(df, rotulo) {
  d <- df[df$termo %in% TERMOS, c("metabolito", "termo", "logFC", "FDR")]
  names(d)[3:4] <- paste0(names(d)[3:4], "_", rotulo)
  d
}
achados <- Reduce(function(a, b) merge(a, b, by = c("metabolito", "termo")),
                  list(pegar(res_anova, "ANOVA"), pegar(res_limma, "limma"),
                       pegar(res_anova_s, "ANOVA_sens"), pegar(res_limma_s, "limma_sens")))
col_fdr <- grep("^FDR_", names(achados), value = TRUE)
achados$metodos_com_FDR_ate_limiar <- rowSums(achados[, col_fdr] <= FDR_LIMIAR)
achados$nivel_principal <- nivel_fdr(achados[[if (METODO_PRINCIPAL == "limma") "FDR_limma" else "FDR_ANOVA"]])
achados$consistente_nas_4 <- achados$metodos_com_FDR_ate_limiar == 4
achados <- achados[achados$metodos_com_FDR_ate_limiar > 0, ]
achados <- achados[order(achados$termo, -achados$metodos_com_FDR_ate_limiar, achados[[col_fdr[1]]]), ]
achados$provavel_nao_biologico <- achados$metabolito %in% PROVAVEL_NAO_BIOLOGICO

principal <- if (METODO_PRINCIPAL == "limma") res_limma else res_anova


# ---- 6. Console summary -----------------------------------------------------
cat("\n== Efeitos principais e interação (FDR <=", FDR_LIMIAR, ") ==\n")
for (nome in c("ANOVA", "limma")) {
  r <- if (nome == "ANOVA") res_anova else res_limma
  for (termo_i in TERMOS) {
    s <- r[r$termo == termo_i & r$FDR <= FDR_LIMIAR, ]
    s <- s[order(s$FDR), ]
    cat(sprintf("%-6s %-6s %s\n", nome, termo_i,
                if (nrow(s)) paste0(s$metabolito, " (", fmt(s$FDR, 3), ")", collapse = "; ") else "nenhum"))
  }
}
cat("\nlimma: gl do prior =", fmt(lim$prior_gl), "; s0^2 =", fmt(lim$prior_s2, 4), "\n")
cat("\n== PERMANOVA de duas vias ==\n"); print(res_perm, digits = 4, row.names = FALSE)
cat("SQ iguais às do vegan::adonis2:", sq_confere, "\n")
cat("PERMDISP: F =", fmt(F_disp, 3), "; P =", P_disp, "\n")
cat("\n== Pares (efeitos simples) ==\n"); print(res_pares[, c("comparacao", "R2", "P_exato", "FDR")], digits = 3, row.names = FALSE)
cat("\nd(MA+MI, Controle) - d(MA+MI, MI) =", fmt(res_cent$dif), "; IC 95%:", fmt(res_cent$dif_ic), "\n")
cat("\n== Classes ==\n"); print(res_classes[, c("classe", "termo", "R2", "P", "FDR")], digits = 3, row.names = FALSE)


# ---- 7. Figures -------------------------------------------------------------
tema <- theme_classic(base_size = 9) +
  theme(axis.line  = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        axis.text  = element_text(colour = "grey20"),
        legend.position = "top",
        legend.key.size = unit(3.5, "mm"),
        legend.text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"))

salvar <- function(p, nome, largura_mm, altura_mm) {
  ggsave(file.path(DIR_FIG, paste0(nome, ".pdf")), p, width = largura_mm, height = altura_mm,
         units = "mm", device = DISP_PDF)
  ggsave(file.path(DIR_FIG, paste0(nome, ".png")), p, width = largura_mm, height = altura_mm,
         units = "mm", dpi = 600)
}

# 7.1 PCA (same matrix as the PERMANOVA)
pca <- prcomp(Z, center = FALSE)
ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
d_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], trat = meta$trat, planta = meta$amostra)
p_pca <- ggplot(d_pca, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") +
  stat_ellipse(aes(colour = trat, fill = trat), geom = "polygon", type = "t", level = 0.95,
               alpha = 0.07, linewidth = 0.35, show.legend = FALSE) +
  geom_point(aes(colour = trat, shape = trat), size = 2.3) +
  scale_colour_manual(values = CORES, name = NULL) +
  scale_fill_manual(values = CORES, guide = "none") +
  scale_shape_manual(values = FORMAS, name = NULL) +
  labs(x = paste0("PC1 (", fmt(ve[1], 1), "%)"), y = paste0("PC2 (", fmt(ve[2], 1), "%)")) +
  tema + theme(panel.grid.major.y = element_blank())
if (ROTULAR_PLANTAS_PCA) {
  p_pca <- p_pca + geom_text_repel(aes(label = planta), size = 2.2, colour = "grey40",
                                   min.segment.length = 0.3, max.overlaps = Inf, seed = SEMENTE)
}
salvar(p_pca, "Fig_PCA", 110, 100)

# 7.2 Interaction plots of metabolites with FDR <= threshold (main method)
sig_princ <- unique(principal$metabolito[principal$termo %in% TERMOS & principal$FDR <= FDR_LIMIAR])
if (length(sig_princ) > 0) {
  d_ind <- do.call(rbind, lapply(sig_princ, function(m) data.frame(metabolito = m, meta, y = log2(X[m, ]))))
  d_ind <- d_ind[!is.na(d_ind$y), ]
  desloc <- c(Without = -0.12, With = 0.12)
  set.seed(SEMENTE)
  d_ind$x <- as.numeric(d_ind$MI) + desloc[as.character(d_ind$MA)] + runif(nrow(d_ind), -0.035, 0.035)
  chaves <- unique(d_ind[, c("metabolito", "MA", "MI", "trat")])
  d_med <- do.call(rbind, lapply(seq_len(nrow(chaves)), function(i) {
    k <- chaves[i, ]
    v <- d_ind$y[d_ind$metabolito == k$metabolito & d_ind$MA == k$MA & d_ind$MI == k$MI]
    data.frame(k, media = mean(v), ep = sd(v) / sqrt(length(v)), n_cel = length(v))
  }))
  d_med$x <- as.numeric(d_med$MI) + desloc[as.character(d_med$MA)]
  anot <- do.call(rbind, lapply(sig_princ, function(m) {
    r <- principal[principal$metabolito == m & principal$termo %in% TERMOS & principal$FDR <= FDR_LIMIAR, ]
    data.frame(metabolito = m,
               rotulo = paste(paste0(ROTULO_TERMO[r$termo], ": FDR = ", fmt(r$FDR, 3)), collapse = "\n"))
  }))
  ordem_sig <- sig_princ[order(sig_princ)]
  d_ind$metabolito <- factor(d_ind$metabolito, levels = ordem_sig)
  d_med$metabolito <- factor(d_med$metabolito, levels = ordem_sig)
  anot$metabolito  <- factor(anot$metabolito, levels = ordem_sig)
  n_col <- min(3, length(sig_princ))

  p_int <- ggplot() +
    geom_point(data = d_ind, aes(x = x, y = y, colour = trat), size = 1.1, alpha = 0.55, show.legend = FALSE) +
    geom_line(data = d_med, aes(x = x, y = media, group = MA, linetype = MA), colour = "grey25", linewidth = 0.45) +
    geom_errorbar(data = d_med, aes(x = x, ymin = media - ep, ymax = media + ep),
                  width = 0.05, colour = "grey25", linewidth = 0.35) +
    geom_point(data = d_med, aes(x = x, y = media, fill = trat, shape = trat),
               size = 2.6, colour = "grey15", stroke = 0.4) +
    geom_text(data = anot, aes(x = -Inf, y = Inf, label = rotulo),
              hjust = -0.08, vjust = 1.25, size = 2.4, lineheight = 0.9, colour = "grey20") +
    facet_wrap(~ metabolito, scales = "free_y", ncol = n_col) +
    scale_x_continuous(breaks = 1:2, labels = c("Without", "With"), expand = expansion(add = 0.35)) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.30))) +
    scale_colour_manual(values = CORES, guide = "none") +
    scale_fill_manual(values = CORES, name = NULL) +
    scale_shape_manual(values = FORMAS_PREENCHIDAS, name = NULL) +
    scale_linetype_manual(values = c(Without = "22", With = "solid"),
                          name = expression(italic("M. anisopliae"))) +
    labs(x = expression(italic("M. incognita")), y = expression(log[2] ~ "relative abundance")) +
    tema
  salvar(p_int, "Fig_interacao", 174, 18 + 55 * ceiling(length(sig_princ) / n_col))
}

# 7.3 Volcano plot of the main effects and the interaction (main method)
d_vol <- principal[principal$termo %in% TERMOS, ]
d_vol$termo_rotulo <- factor(ROTULO_TERMO[d_vol$termo], levels = ROTULO_TERMO)
cores_nivel <- setNames(c("#1f5fae", "#8fb8ea", "grey78"), levels(d_vol$nivel))
p_vol <- ggplot(d_vol, aes(x = logFC, y = -log10(P))) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_point(aes(colour = nivel), size = 1.7) +
  geom_text_repel(data = d_vol[d_vol$nivel != "n.s.", ], aes(label = metabolito),
                  size = 2.3, min.segment.length = 0, max.overlaps = Inf, seed = SEMENTE) +
  facet_wrap(~ termo_rotulo, nrow = 1) +
  scale_colour_manual(values = cores_nivel, drop = FALSE, name = NULL) +
  labs(x = expression(log[2] ~ "fold change"), y = expression(-log[10] ~ italic(P))) +
  tema
salvar(p_vol, "Fig_volcano", 174, 80)

# 7.4 Heatmap of the tested metabolites, grouped by class (requested by Reviewer 2)
Lh <- log2(X[testados, , drop = FALSE])
Zh <- t(scale(t(Lh)))                          # z-score per metabolite; NA kept
classe_de <- setNames(rep(names(CLASSES), lengths(CLASSES)), unlist(CLASSES, use.names = FALSE))
classe_met <- ifelse(testados %in% names(classe_de), classe_de[testados], "Other")
names(classe_met) <- testados
sig_qualquer <- tapply(principal$FDR[principal$termo %in% TERMOS] <= FDR_LIMIAR,
                       principal$metabolito[principal$termo %in% TERMOS], any)
lfc_mi <- setNames(principal$logFC[principal$termo == "MI"], principal$metabolito[principal$termo == "MI"])
ordem_met <- testados[order(factor(classe_met, levels = c(names(CLASSES), "Other")), -lfc_mi[testados])]
rotulo_met <- setNames(paste0(testados, ifelse(sig_qualquer[testados], " *", "")), testados)
ordem_plantas <- meta$amostra[order(meta$trat, meta$amostra)]
d_hm <- data.frame(metabolito = rep(rownames(Zh), times = ncol(Zh)),
                   amostra = rep(as.integer(colnames(Zh)), each = nrow(Zh)),
                   z = as.vector(Zh))
d_hm$trat   <- meta$trat[match(d_hm$amostra, meta$amostra)]
d_hm$classe <- factor(classe_met[d_hm$metabolito], levels = c(names(CLASSES), "Other"))
d_hm$rotulo <- factor(rotulo_met[d_hm$metabolito], levels = rev(rotulo_met[ordem_met]))
d_hm$planta <- factor(d_hm$amostra, levels = ordem_plantas)
p_hm <- ggplot(d_hm, aes(x = planta, y = rotulo, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(classe ~ trat, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f4f4f2", high = "#eb6834", midpoint = 0,
                       limits = c(-2.5, 2.5), oob = scales::squish, na.value = "white",
                       name = "z-score") +
  labs(x = "Plant", y = NULL) +
  tema +
  theme(legend.position = "right",
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6.5),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 7),
        panel.grid.major.y = element_blank(), panel.spacing = unit(0.8, "mm"))
salvar(p_hm, "Fig_heatmap", 174, 40 + 3.6 * length(testados))

# 7.5 Chemical classes: variation explained by term
d_cl <- res_classes
d_cl$termo_rotulo <- factor(ROTULO_TERMO[d_cl$termo], levels = ROTULO_TERMO)
d_cl$marca  <- ifelse(d_cl$FDR <= FDR_ESTRITO, "**", ifelse(d_cl$FDR <= FDR_LIMIAR, "*", ""))
d_cl$classe <- factor(d_cl$classe, levels = rev(names(CLASSES)))
pd <- position_dodge(width = 0.65)
p_cl <- ggplot(d_cl, aes(x = classe, y = 100 * R2, colour = termo_rotulo, shape = termo_rotulo)) +
  geom_linerange(aes(ymin = 0, ymax = 100 * R2), position = pd, linewidth = 0.5) +
  geom_point(position = pd, size = 2.2) +
  geom_text(aes(label = marca), position = pd, hjust = -0.5, vjust = 0.75, size = 3, show.legend = FALSE) +
  coord_flip() +
  scale_colour_manual(values = CORES_TERMO, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = expression("Variation explained," ~ italic(R)^2 ~ "(%)")) +
  tema + theme(panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey92"))
salvar(p_cl, "Fig_classes", 120, 20 + 11 * length(CLASSES))

# 7.6 Distances between centroids (requested by Reviewer 2)
d_ce <- res_cent$pares
subtitulo_ce <- paste0("d(MA + MI, Control) ", if (USA_CAIRO) "\u2212" else "-", " d(MA + MI, MI) = ",
                       fmt(res_cent$dif), " (95% CI ", fmt(res_cent$dif_ic[1]), " to ",
                       fmt(res_cent$dif_ic[2]), ")")
d_ce$destaque <- ifelse(d_ce$g1 == "Controle" & d_ce$g2 == "MA_MI", "ctrl",
                  ifelse(d_ce$g1 == "MI" & d_ce$g2 == "MA_MI", "mi", "outro"))
p_ce <- ggplot(d_ce, aes(x = d, y = reorder(par, d), colour = destaque)) +
  geom_linerange(aes(xmin = IC_inf, xmax = IC_sup), linewidth = 0.6) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(ctrl = "#2a78d6", mi = "#eb6834", outro = "grey60"), guide = "none") +
  labs(x = "Euclidean distance between treatment centroids (95% bootstrap interval)", y = NULL,
       subtitle = subtitulo_ce) +
  tema + theme(plot.subtitle = element_text(size = 8, colour = "grey25"),
               panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey92"))
salvar(p_ce, "Fig_centroides", 140, 70)


# ---- 8. Spreadsheet ---------------------------------------------------------
com_nome <- function(m, nome) data.frame(setNames(list(rownames(m)), nome), as.data.frame(m), check.names = FALSE)
deteccao <- data.frame(metabolito = rownames(det), det, check.names = FALSE)
deteccao$testado <- deteccao$metabolito %in% testados
deteccao$provavel_nao_biologico <- deteccao$metabolito %in% PROVAVEL_NAO_BIOLOGICO

leia_me <- data.frame(
  item = c("Desenho", "Metabólitos testados", "Limiar", "ANOVA", "limma", "PERMANOVA",
           "Pares", "Centroides", "Classes", "Sensibilidade", "Método principal", "Parâmetros"),
  descricao = c(
    "Fatorial 2 x 2 inteiramente casualizado; MA e MI presentes/ausentes; 6 plantas por tratamento",
    paste0(length(testados), " de ", nrow(X), " metabólitos com valor em >= ", MIN_DET,
           " de 6 plantas em todos os tratamentos", if (EXCLUIR_NAO_BIOLOGICOS) "; sem prováveis não biológicos" else ""),
    paste0("Significativo: FDR <= ", FDR_LIMIAR, "; coluna 'nivel' separa FDR <= ", FDR_ESTRITO,
           " de ", FDR_ESTRITO, " < FDR <= ", FDR_LIMIAR, "; BH dentro de cada termo/contraste"),
    "log2; lm(y ~ MA * MI) com contr.sum; car::Anova(type = 3); NA omitidos; logFC = contrastes de médias de célula",
    paste0("lmFit(~ 0 + grupo) + makeContrasts + eBayes; gl do prior = ", fmt(lim$prior_gl, 3),
           "; s0^2 = ", fmt(lim$prior_s2, 5)),
    paste0("Euclidiana; log2 + Pareto; NA -> metade do mínimo; SQ Tipo III; ", N_PERM,
           " permutações; conferida com vegan::adonis2; PERMDISP (betadisper, centroide)"),
    "PERMANOVA exata (todas as rotulagens) nas 4 comparações de efeito simples; BH entre as 4",
    paste0("Distâncias euclidianas entre centroides; ", N_BOOT, " bootstraps dentro de tratamento"),
    "PERMANOVA por classe (mesma matriz); BH entre classes dentro de cada termo; ver aba Classes_definicao",
    paste0("Sem as plantas ", paste(PLANTAS_SENSIBILIDADE, collapse = " e ")),
    paste0(METODO_PRINCIPAL, " (define figuras e a coluna nivel_principal em Achados)"),
    paste0("SOMAR_PICOS_DUPLICADOS = ", SOMAR_PICOS_DUPLICADOS, "; semente = ", SEMENTE)))


# ---- 9. VALIDATION against Python (default parameters) ----------------------
padrao <- EXCLUIR_NAO_BIOLOGICOS && !SOMAR_PICOS_DUPLICADOS && MIN_DET == 3 &&
          identical(sort(PLANTAS_SENSIBILIDADE), c(10, 12))
if (padrao) {
  ref <- list(n_testados = 46,
              putrescina_F_MI = 18.889587, tropina_F_MI = 21.251762,
              putrescina_FDR_MI = 0.010099, tropina_FDR_MI = 0.010099,
              permanova_SQ = c(56.301686, 94.043674, 42.306368, 819.159565),
              permdisp_F = 1.101067,
              pares_P_exato = c(0.064935, 0.082251, 0.556277, 0.036797),
              centroide_dif = 0.148833,
              sens_putrescina_F_MI = 14.251264)
  pega <- function(r, m, t, col) r[[col]][r$metabolito == m & r$termo == t]
  obt <- list(n_testados = length(testados),
              putrescina_F_MI = pega(res_anova, "Putrescine", "MI", "estatistica"),
              tropina_F_MI = pega(res_anova, "Tropine", "MI", "estatistica"),
              putrescina_FDR_MI = pega(res_anova, "Putrescine", "MI", "FDR"),
              tropina_FDR_MI = pega(res_anova, "Tropine", "MI", "FDR"),
              permanova_SQ = res_perm$SQ[1:4],
              permdisp_F = F_disp,
              pares_P_exato = res_pares$P_exato,
              centroide_dif = res_cent$dif,
              sens_putrescina_F_MI = pega(res_anova_s, "Putrescine", "MI", "estatistica"))
  validacao <- data.frame(
    item = names(ref),
    python = vapply(ref, function(v) paste(signif(v, 6), collapse = "; "), character(1)),
    R = vapply(obt, function(v) paste(signif(v, 6), collapse = "; "), character(1)),
    confere = mapply(function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-4)), ref, obt))
  validacao <- rbind(validacao, data.frame(
    item = "limma: tropina e putrescina com FDR <= 0,05 para MI (aproximado)", python = "TRUE",
    R = as.character(all(pega(res_limma, "Tropine", "MI", "FDR") <= 0.05,
                         pega(res_limma, "Putrescine", "MI", "FDR") <= 0.05)),
    confere = all(pega(res_limma, "Tropine", "MI", "FDR") <= 0.05,
                  pega(res_limma, "Putrescine", "MI", "FDR") <= 0.05)))
  cat("\n== VALIDAÇÃO (R x Python) ==\n"); print(validacao, row.names = FALSE)
  if (!all(validacao$confere)) warning("Há diferenças em relação ao Python: ver aba VALIDACAO.")
} else {
  validacao <- data.frame(item = "Validação não aplicável: parâmetros diferentes do padrão")
}

sem_fator <- function(d) { d[] <- lapply(d, function(x) if (is.factor(x)) as.character(x) else x); d }
write_xlsx(lapply(list(
  LEIA_ME              = leia_me,
  Achados              = achados,
  ANOVA                = res_anova[order(res_anova$termo, res_anova$P), ],
  limma                = res_limma[order(res_limma$termo, res_limma$P), ],
  PERMANOVA            = res_perm,
  PERMANOVA_vegan      = com_nome(res_adonis, "termo"),
  PERMDISP             = data.frame(estatistica = c("F", "P (permutação)"), valor = c(F_disp, P_disp)),
  Pares_efeito_simples = res_pares,
  Centroides           = rbind(res_cent$pares[, c("par", "d", "IC_inf", "IC_sup")],
                               data.frame(par = "d(MA+MI,Controle) - d(MA+MI,MI)", d = res_cent$dif,
                                          IC_inf = res_cent$dif_ic[1], IC_sup = res_cent$dif_ic[2])),
  Classes              = res_classes,
  Classes_definicao    = data.frame(classe = rep(names(CLASSES), lengths(CLASSES)),
                                    metabolito = unlist(CLASSES, use.names = FALSE),
                                    testado = unlist(CLASSES, use.names = FALSE) %in% testados),
  Sens_ANOVA           = res_anova_s[order(res_anova_s$termo, res_anova_s$P), ],
  Sens_limma           = res_limma_s[order(res_limma_s$termo, res_limma_s$P), ],
  Sens_PERMANOVA       = res_perm_s,
  Sens_centroides      = rbind(res_cent_s$pares[, c("par", "d", "IC_inf", "IC_sup")],
                               data.frame(par = "d(MA+MI,Controle) - d(MA+MI,MI)", d = res_cent_s$dif,
                                          IC_inf = res_cent_s$dif_ic[1], IC_sup = res_cent_s$dif_ic[2])),
  Deteccao             = deteccao,
  VALIDACAO            = validacao
), sem_fator), path = file.path(DIR_SAIDA, "resultados_completos.xlsx"))

writeLines(capture.output(sessionInfo()), file.path(DIR_SAIDA, "sessionInfo.txt"))
cat("\nConcluído. Resultados em:", normalizePath(DIR_SAIDA), "\n")
