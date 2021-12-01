### open with R studio ###

## ---- setup -----------------
if (!"tidyverse" %in% rownames(installed.packages())) {install.packages("tidyverse")}
if (!"Seurat" %in% rownames(installed.packages())) {install.packages("Seurat")}
if (!"rstudioapi" %in% rownames(installed.packages())) {install.packages("rstudioapi")}
if (!"hdf5r" %in% rownames(installed.packages())) {install.packages("hdf5r")}
if (!"scales" %in% rownames(installed.packages())) {install.packages("scales")}

library(tidyverse)
library(Seurat) # Seurat V4
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

# functions for axis text used later
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(base::pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

scientific_10 <- function(x) {
  if (all(substr(as.character(x[which(!is.na(x))]),1,1) == "1")) {
    x[which(x>=100)] <- gsub("1e\\+", "10^", scales::scientific_format()(x[which(x>=100)]))
  } else if (any(substr(as.character(x[which(!is.na(x))]),1,1) != "1")) {
    x[which(x>=100)] <- gsub("e\\+", "%*%10^", scales::scientific_format()(x[which(x>=100)]))
  }
  parse(text = x)
}


## ---- folder_structure --------------
print(list.files(wd, recursive = T)[which(grepl("^data|r_script.R", list.files(wd, recursive = T)))])

## ---- rename_gene_barcode_matrix_files_for_Kubli_et_al_(GSE130287) --------------
# renaming required to make it valid for Seurat::Read10X (Seurat4)
paths <- list.files(file.path(wd, "data", "GSE130287"), recursive = T, full.names = T)
file.rename(paths, file.path(dirname(paths), gsub("genes", "features", sapply(sapply(strsplit(basename(paths), "_"), rev, simplify = F), "[", 1))))

## ---- prepare_Seurat_object_Kubli --------------
# mm for Fcmr minus/minus, pp for Fcmr plus/plus
raw_mat <- Seurat::Read10X(data.dir = stats::setNames(list.dirs(file.path(wd, "data", "GSE130287"))[-1], nm = c("mm", "pp")), strip.suffix = T)
SO_kubli <- Seurat::CreateSeuratObject(counts = raw_mat, min.cells = 3)
SO_kubli <- Seurat::AddMetaData(SO_kubli, Seurat::PercentageFeatureSet(SO_kubli, pattern = "^mt-"), "pct.mt")
SO_kubli <- subset(SO_kubli, cells = SeuratObject::WhichCells(SO_kubli, expression = pct.mt < 30))
SO_kubli <- subset(SO_kubli, cells = SeuratObject::WhichCells(SO_kubli, expression = nFeature_RNA < 7000))
SO_kubli <- Seurat::NormalizeData(SO_kubli)
table(SO_kubli@meta.data$orig.ident)

## ---- prepare_Seurat_object_Riedel --------------
raw_mat <- do.call(cbind, lapply(list.files(file.path(wd, "data", "GSE140133"), pattern = "GSM4147195|GSM4147196", full.names = T, recursive = T), function(x) {
  y <- Seurat::Read10X_h5(x)
  colnames(y) <- paste0(strsplit(basename(x), "_")[[1]][1], "_", colnames(y))
  return(y)
}))

SO_riedel <- Seurat::CreateSeuratObject(counts = raw_mat, min.cells = 3)
SO_riedel <- Seurat::AddMetaData(SO_riedel, Seurat::PercentageFeatureSet(SO_riedel, pattern = "^mt-"), "pct.mt")
filt_pct.mt <- 
  SO_riedel@meta.data %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::filter(pct.mt < 10) %>%
  dplyr::filter(pct.mt > 2) %>%
  group_by(orig.ident) %>%
  dplyr::filter(pct.mt < quantile(pct.mt, 0.99)) %>%
  dplyr::filter(pct.mt > quantile(pct.mt, 0.01)) %>%
  ungroup() %>%
  dplyr::pull(ID)
filt_nFeat <-
  SO_riedel@meta.data %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::filter(nFeature_RNA < 4000) %>%
  dplyr::filter(nFeature_RNA > 1000) %>%
  group_by(orig.ident) %>%
  dplyr::filter(nFeature_RNA < quantile(nFeature_RNA, 0.975)) %>%
  ungroup() %>%
  dplyr::pull(ID)
SO_riedel <- subset(SO_riedel, cells = intersect(filt_pct.mt, filt_nFeat))
SO_riedel <- Seurat::NormalizeData(SO_riedel)


## ---- create_figures --------------
features <- c("Fcmr", "Gapdh", "Cd74", "Cd68", "Itgam", "Cd86")

counts_kubli <- as.data.frame(t(as.matrix(Seurat::GetAssayData(SO_kubli, slot = "counts", assay = "RNA")[features,])))
counts_kubli$type <- SO_kubli@meta.data$orig.ident
counts_kubli <- tidyr::pivot_longer(counts_kubli, cols = dplyr::all_of(features), names_to = "feature", values_to = "count")
counts_kubli$feature <- factor(counts_kubli$feature, levels = features)

counts_riedel <- as.data.frame(t(as.matrix(Seurat::GetAssayData(SO_riedel, slot = "counts", assay = "RNA")[features,])))
counts_riedel <- tidyr::pivot_longer(counts_riedel, cols = dplyr::all_of(features), names_to = "feature", values_to = "count")
counts_riedel$feature <- factor(counts_riedel$feature, levels = features)

# Fig.1a
for (i in unique(counts_kubli$type)) {
  p <- ggplot2::ggplot(counts_kubli[which(counts_kubli$type == i),], ggplot2::aes(x = count)) +
    ggplot2::geom_bar(color = "black", fill = "grey") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = element_rect(fill = "white", color = NA), strip.text = element_text(face = "bold.italic"), text = element_text(face = "bold", size = 20, color = "black")) +
    ggplot2::scale_y_log10(label = scientific_10) +
    ggplot2::scale_x_continuous(breaks = integer_breaks(n = 4)) +
    ggplot2::ylab("Number of cells") +
    ggplot2::xlab("Raw values of unique molecular identifier (UMI) per cell") +
    ggplot2::facet_wrap(ggplot2::vars(feature), scales = "free", nrow = 1)
  print(p)
  ggplot2::ggsave(p, filename = paste0("Fig.1a_", i, ".eps"), device = "eps", path = file.path(wd, "figures"), width = 16, height = 4)
}


# Fig.1b
ggplot2::ggplot(counts_riedel, ggplot2::aes(x = count)) +
  ggplot2::geom_bar(color = "black", fill = "grey") +
  ggplot2::theme_bw() +
  ggplot2::theme(strip.background = element_rect(fill = "white", color = NA), strip.text = element_text(face = "bold.italic"), text = element_text(face = "bold", size = 20, color = "black")) +
  ggplot2::scale_y_log10(label = scientific_10) +
  ggplot2::scale_x_continuous(breaks = integer_breaks(n = 4)) +
  ggplot2::ylab("Number of cells") +
  ggplot2::xlab("Raw values of unique molecular identifier (UMI) per cell") +
  ggplot2::facet_wrap(ggplot2::vars(feature), scales = "free", nrow = 1)
ggplot2::ggsave(filename = "Fig.1b.eps", device = "eps", path = file.path(wd, "figures"), width = 16, height = 4)


# Fig.1c
umi_kubli <- as.data.frame(t(as.matrix(Seurat::GetAssayData(SO_kubli, slot = "data", assay = "RNA")["Fcmr",,drop=F])))
umi_kubli$type <- as.character(SO_kubli@meta.data$orig.ident)
umi_kubli[which(umi_kubli$type == "mm"), "type"] <- "Fcmr KO TMP"
umi_kubli[which(umi_kubli$type == "pp"), "type"] <- "Fcmr WT TMP"
umi_riedel <- as.data.frame(t(as.matrix(Seurat::GetAssayData(SO_riedel, slot = "data", assay = "RNA")["Fcmr",,drop=F])))
umi_riedel$type <- "WT spl. IgG B"
umi <- rbind(umi_kubli, umi_riedel)
umi <-
  umi %>%
  dplyr::mutate(type_italic = case_when(type == "Fcmr WT TMP" ~ "bold(bolditalic('Fcmr')~'WT TMP')",
                                        type == "Fcmr KO TMP" ~ "bold(bolditalic('Fcmr')~'KO TMP')",
                                        type == "WT spl. IgG B" ~ "bold(WT~spl.~IgG~B)"))
umi$type_italic <- factor(umi$type_italic, levels = unique(umi$type_italic)[c(2,1,3)])


ggplot2::ggplot(umi, ggplot2::aes(x = Fcmr)) +
  ggplot2::geom_freqpoly(bins = 200) +
  ggplot2::theme_bw() +
  ggplot2::theme(strip.background = element_rect(fill = "white", color = NA), strip.text = element_text(face = "bold"), text = element_text(face = "bold", size = 20, color = "black")) +
  ggplot2::scale_y_log10(label = scientific_10) +
  ggplot2::scale_x_continuous(breaks = integer_breaks(n = 4)) +
  ggplot2::ylab("Number of cells") +
  ggplot2::xlab(expression(paste(bold("Normalized "), bolditalic("Fcmr"), bold(" transcripts (UMI/cell in "), bold(Log[e]), bold(")")))) +
  ggplot2::facet_wrap(ggplot2::vars(type_italic), scales = "free", nrow = 1, labeller = label_parsed)
ggplot2::ggsave(filename = "Fig.1c.eps", device = "eps", path = file.path(wd, "figures"), width = 10, height = 4)





