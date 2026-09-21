#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2026 Broad Institute
#
#  Author: Giulio Genovese <giulio.genovese@gmail.com>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
###

options(error = function() {traceback(3); q("no", 1)})

xipher_version <- "2026-09-21"

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table)) # used to read sparse matrices in MatrixMarket format
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(VGAM)) # used for beta-binomial computations to identify homozygous variants
suppressPackageStartupMessages(library(RSpectra)) # used for computing the initial state with svds

###########################################################################
## FUNCTIONS TO LOAD INPUT DATA                                          ##
###########################################################################

# by using data.table::fread(), this approach is much faster than using Matrix::readMM() which is super slow
# see http://github.com/satijalab/seurat/issues/9711 for a complete discussion
read_mm <- function(file) {
  df <- data.table::fread(file = file, header = TRUE, skip = 1, colClasses = c("integer", "integer", "integer"), quote = "", data.table = FALSE)
  ad <- Matrix::sparseMatrix(i = df[,1], j = df[,2], x = df[,3], dims = as.numeric(c(names(df)[1], names(df)[2])), repr = "T")
}

# load an scDblFinder (or scPred) output table and returns a list of cell doublets
load_scdbl <- function(file) {
  df_dbl <- read.table(file, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
  write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Read doublets file", file), stderr())
  if (all(c("cell", "class") %in% names(df_dbl))) return(df_dbl$cell[df_dbl$class == "doublet"])
  else if (all(c("CELL_BARCODE", "doublet") %in% names(df_dbl))) return(df_dbl$CELL_BARCODE[df_dbl$doublet == "TRUE"])
  else stop(paste("File", file, "is neither an scDblFinder nor an scPred table\n"))
}

# load a donor cell map table
load_donor_cell_map <- function(file, donors) {
  cell2donor <- read.table(file, sep = "\t", header = TRUE)
  if (!any(cell2donor$bestSample %in% unlist(strsplit(donors, ",")))) stop(paste("File", file, "does not contain", donors, "\n"))
  write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Read donor cell map file", file), stderr())
  return(cell2donor$cell[cell2donor$bestSample %in% unlist(strsplit(donors, ","))])
}

# read VarTrix files .SNP.vcf.gz .barcodes.tsv, .alt.mtx.gz, .ref.mtx.gz
# read cellSNP-lite files .base.vcf.gz, .samples.tsv, .tag.AD.mtx, and .tag.DP.mtx (as .tag.OTH.mtx is not required)
# see http://cellsnp-lite.readthedocs.io/en/latest/main/manual.html#output for more information
load_mtx <- function(prefix, cellsnp = FALSE, donors = NULL, doublets = NULL, name = NULL) {
  ret <- list()
  cmd <- paste("bcftools query -f \"%CHROM\\_%POS\\_%REF\\_%ALT\n\"", paste0(prefix, ifelse(cellsnp, ".base.vcf.gz", ".SNP.vcf.gz")))
  ret$vars <- data.table::fread(cmd = cmd, sep = "\t", header = FALSE, stringsAsFactors = TRUE, data.table = FALSE)$V1
  ret$cells <- readLines(paste0(prefix,  ifelse(cellsnp, ".samples.tsv", ".barcodes.tsv")))
  ret$alt <- read_mm(paste0(prefix, ifelse(cellsnp, ".tag.AD.mtx.gz", ".alt.mtx.gz")))
  ret$ref <- if (cellsnp) as(Matrix::drop0(read_mm(paste0(prefix, ".tag.DP.mtx.gz")) - ret$alt), "TsparseMatrix") else read_mm(paste0(prefix, ".ref.mtx.gz"))
  if (!is.null(donors)) {
    idx <- match(donors, ret$cells)
    ret$cells <- ret$cells[na.omit(idx)]
    ret$alt <- ret$alt[,na.omit(idx)]
    ret$ref <- ret$ref[,na.omit(idx)]
  }
  ret$doublet <- if (is.null(doublets)) rep(FALSE, length(ret$cells)) else ret$cells %in% doublets
  if (!is.null(name)) ret$cells <- paste(name, ret$cells, sep="_")
  write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Read allelic counts data from", prefix), stderr())
  return(ret)
}

# merge multiple vartrix/cellSNP-lite datasets assuming there is no intersection about the barcodes
merge_libraries <- function(lst) {
  if (length(lst) == 1) return(lst[[1]])
  ret <- list()
  ret$vars <- Reduce(union, lapply(lst, `[[`, "vars"))
  ret$cells <- unlist(lapply(lst, `[[`, "cells"), use.names = FALSE)
  idxs <- lapply(lst, function(x) match(x$vars, ret$vars))
  offsets <- cumsum(c(0L, head(unname(sapply(lst, function(x) length(x$cells))), -1)))
  ret$ref <- Matrix::sparseMatrix(
    i = unlist(Map(function(x, idx) idx[x$ref@i + 1L], lst, idxs), use.names = FALSE),
    j = unlist(Map(function(x, offset) x$ref@j + 1L + offset, lst, offsets), use.names = FALSE),
    x = unlist(lapply(lst, function(x) x$ref@x), use.names = FALSE),
    dims = c(length(ret$vars), length(ret$cells)),
    repr = "T"
  )
  ret$alt <- Matrix::sparseMatrix(
    i = unlist(Map(function(x, idx) idx[x$alt@i + 1L], lst, idxs), use.names = FALSE),
    j = unlist(Map(function(x, offset) x$alt@j + 1L + offset, lst, offsets), use.names = FALSE),
    x = unlist(lapply(lst, function(x) x$alt@x), use.names = FALSE),
    dims = c(length(ret$vars), length(ret$cells)),
    repr = "T"
  )
  ret$doublet <- unlist(lapply(lst, `[[`, "doublet"), use.names = FALSE)
  write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Merged data from", length(lst), "samples"), stderr())
  return(ret)
}

###########################################################################
## FUNCTION TO RUN THE CO-CLUSTERING ALGORITHM                           ##
###########################################################################

# co-cluster variants and cell states by iteratively multiplyng by D and t(D)
# this algorithm tries to maximize the function sign(p) %*% D %*% sign(h). To do so
# to avoid local minima it initializes the first state by computing the SVD of D
# to then converge to a suitable partition it uses deterministic annealing
# deterministic annealing could be replaced with simulated annealing:
# prob <- plogis(2 * betas[i] * as.vector(D %*% h) / rs)
# p <- ifelse(runif(length(prob)) < prob, 1, -1)
# prob <- plogis(2 * betas[i] * as.vector(p %*% D) / cs)
# h <- ifelse(runif(length(prob)) < prob, 1, -1)
# but through testing it seems that deterministic annealing works well enough
co_cluster <- function(D, rs, cs, max_iter) {
  Dn <- Matrix::Diagonal( x = ifelse(rs == 0, 0, 1 / sqrt(rs)) ) %*% D %*%
        Matrix::Diagonal( x = ifelse(cs == 0, 0, 1 / sqrt(cs)) )
  sv <- RSpectra::svds(Dn, k = 2)
  active_ratio <- (1 + c(as.numeric(sign(sv$u[,1]) %*% D %*% sign(sv$v[,1])), as.numeric(sign(sv$u[,2]) %*% D %*% sign(sv$v[,2]))) / sum(rs)) / 2
  p <- sign(as.vector(sv$u[,which.max(active_ratio)]))
  h <- sign(as.vector(sv$v[,which.max(active_ratio)]))
  write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Initial state using SVD achieved a", format(max(active_ratio), digits=4), "active ratio"), stderr())
  start <- proc.time()
  iter <- 0
  betas <- c(2, 4, 8, 16, 32, 64)
  for (i in seq_len(length(betas))) {
    repeat {
      # update phase of each variant
      old_p <- p
      p <- tanh(betas[i] * as.vector(D %*% h) / rs)
      p[rs == 0] <- 0
      if (all(abs(p - old_p) < 1e-6)) break

      # infer active haplotype in every cell
      old_h <- h
      h <- tanh(betas[i] * as.vector(p %*% D) / cs)
      h[cs == 0] <- 0
      if (all(abs(h - old_h) < 1e-6)) break
      iter <- iter + 1

      if (iter * length(betas) >= i * max_iter) break
    }
  }
  p <- sign(p)
  h <- sign(h)
  end <- proc.time()
  elapsed_time <- as.numeric((end - start)["elapsed"])
  active_ratio <- (1 + as.numeric(p %*% D %*% h) / sum(rs)) / 2
  write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), ifelse(iter == max_iter, "Did not converge" , "Converged"), "after", iter, "iterations in", format(elapsed_time, digits=4), "seconds and achieved a", format(active_ratio, digits=4), "active ratio"), stderr())
  return(list(phase = p, cell_state = h))
}

###########################################################################
## OPTION PARSER CODE                                                    ##
###########################################################################

parser <- OptionParser(paste0("xipher.R [options] --mtx <prefix>|--table <file.tsv>\nAbout: infer X-inactivation from single cell data. (version ", xipher_version, " http://github.com/freeseek/xipher)\n[Burger, S., et al. Human X-inactivation mosaicism reveals the cell-autonomous effects\nof X-linked mutations. (2023) http://nrs.harvard.edu/URN-3:HUL.INSTREPOS:37377886]"))
parser <- add_option(parser, c("--mtx"), type = "character", default = ".", help = "prefix for input files generated by VarTrix or cellsnp-lite", metavar = "<prefix>")
parser <- add_option(parser, c("--cellsnp"), action='store_true', default = FALSE, help = "whether data is in cellSNP-lite format rather than VarTrix format")
parser <- add_option(parser, c("--dbl"), type = "character", default = ".", help = "input scDblFinder/scPred file with doublet predictions", metavar = "<file.tsv>")
parser <- add_option(parser, c("--celldonor"), type = "character", default = ".", help = "input file with donor cell map predictions", metavar = "<file.tsv>")
parser <- add_option(parser, c("--table"), type = "character", default = ".", help = "input table with name, prefix, dbl, and celldonor headers to input multiple samples", metavar = "<file.tsv>")
parser <- add_option(parser, c("--mtx-data-path"), type = "character", default = ".", help = "input data path for VarTrix or cellsnp-lite data in table", metavar = "<path>")
parser <- add_option(parser, c("--dbl-data-path"), type = "character", default = ".", help = "input data path for scDblFinder/scPred files in table", metavar = "<path>")
parser <- add_option(parser, c("--celldonor-data-path"), type = "character", default = ".", help = "input data path for donor cell map files in table", metavar = "<path>")
parser <- add_option(parser, c("--donors"), type = "character", help = "comma separate list of single donor IDs to use if input files are from villages", metavar = "<donor_id>")
parser <- add_option(parser, c("--het-frac"), type = "double", default = 0.2, help = "minimum expected REF/ALT and XA/XI fraction at heterozygous sites/doublet cells [0.2]", metavar = "<float>")
parser <- add_option(parser, c("--hom-frac"), type = "double", default = 0.01, help = "maximum expected REF/ALT and XA/XI fraction at homozygous sites/singlet cells [0.01]", metavar = "<float>")
parser <- add_option(parser, c("--lod"), type = "double", default = 2, help = "LOD score for heterozygous variants and doublets selection [2]", metavar = "<float>")
parser <- add_option(parser, c("--max-iter"), type = "double", default = 1000, help = "maximum number of iterations for the EM algorithm [1000]", metavar = "<integer>")
parser <- add_option(parser, c("--max-xi-proportion"), type = "double", default = 0.2, help = "maximum xi proportion allowed [0.2]", metavar = "<float>")
parser <- add_option(parser, c("--rho"), type = "double", default = 0.001, help = "rho parameter for beta binomial distribution [0.001]", metavar = "<float>")
parser <- add_option(parser, c("--all"), action = "store_true", default = FALSE, help = "whether the output cell UMI counts at all sites")
parser <- add_option(parser, c("--out"), type = "character", default="", help = "output cells TSV file [stdout]", metavar = "<file.tsv>")
parser <- add_option(parser, c("--var"), type = "character", help = "output variants TSV file", metavar = "<file.tsv>")
parser <- add_option(parser, c("--pdf"), type = "character", help = "output PDF file", metavar = "<file.pdf>")
args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)

if ((args$mtx == ".") + (args$table == ".") != 1) {print_help(parser); stop("one and only one of options --mtx and --table is required")}
if (args$celldonor == "." && args$table == "." && !is.null(args$donors)) {print_help(parser); stop("option --donors requires option --celldonor")}
if (args$celldonor != "." && is.null(args$donors)) {print_help(parser); stop("option --celldonor requires option --donors")}
if (args$table == "." && args$mtx_data_path != ".") {print_help(parser); stop("option --mtx-data-path must be used together option --table")}
if (args$table == "." && args$dbl_data_path != ".") {print_help(parser); stop("option --dbl-data-path must be used together option --table")}
write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Started XIPHER version", xipher_version, "http://github.com/freeseek/xipher"), stderr())

if (!is.null(args$pdf)) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggrastr)) # used to plot the output of the co-clustering algorithm
}

###########################################################################
## LOAD DATA IN CELLSNP-LITE FORMAT AND CELLDONOR AND DOUBLETS DATA      ##
###########################################################################

if (args$table != ".") {
  df_tbl <- read.table(args$table, sep = "\t", header = TRUE)
  if (!("name" %in% names(df_tbl))) stop(paste("File", args$table, "does not contain a name column"))
  if (!("prefix" %in% names(df_tbl))) stop(paste("File", args$table, "does not contain a prefix column"))
  lst <- list()
  for (i in seq_len(length(df_tbl$prefix))) {
    donors <- if ("celldonor" %in% names(df_tbl) && df_tbl$celldonor[i] != "." && !is.na(args$donors)) load_donor_cell_map(paste(args$celldonor_data_path, df_tbl$celldonor[i], sep="/"), args$donors) else NULL
    doublets <- if ("dbl" %in% names(df_tbl) && df_tbl$dbl[i] != ".") load_scdbl(paste(args$dbl_data_path, df_tbl$dbl[i], sep="/")) else NULL
    lst[[i]] <- load_mtx(paste(args$mtx_data_path, df_tbl$prefix[i], sep="/"), args$cellsnp, donors, doublets, df_tbl$name[i])
  }
  ret <- merge_libraries(lst)
  rm(lst)
} else {
  ret <- load_mtx(args$mtx,
                  args$cellsnp,
                  if (args$celldonor != ".") load_donor_cell_map(args$celldonor, args$donors) else NULL,
                  if (args$dbl != ".") load_scdbl(args$dbl) else NULL)
}
ret$ref <- as(ret$ref, "CsparseMatrix")
ret$alt <- as(ret$alt, "CsparseMatrix")

###########################################################################
## REMOVE LIKELY HOMOZYGOUS VARIANTS                                     ##
###########################################################################

df_var <- data.frame(var = ret$vars)
df_cell <- data.frame(cell = ret$cells)
df_cell$doublet <- ret$doublet

df_var$ref_count <- Matrix::rowSums(ret$ref)
df_var$alt_count <- Matrix::rowSums(ret$alt)
df_var$umi_count <- df_var$ref_count + df_var$alt_count
df_var$lod <- VGAM::dbetabinom(pmin(df_var$ref_count, df_var$alt_count), df_var$ref_count + df_var$alt_count, args$het_frac, args$rho, log = TRUE) / log(10) -
              VGAM::dbetabinom(pmin(df_var$ref_count, df_var$alt_count), df_var$ref_count + df_var$alt_count, args$hom_frac, args$rho, log = TRUE) / log(10)

if (!is.null(args$pdf)) {
  p1 <- ggplot(df_var, aes(x=ref_count, y=alt_count, color=lod>args$lod)) +
    geom_function(fun = function(x) args$hom_frac/(1-args$hom_frac) * x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) args$het_frac/(1-args$het_frac) * x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) (1-args$het_frac)/args$het_frac * x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) (1-args$hom_frac)/args$hom_frac * x, color = "gray", alpha = 1/2) +
    geom_jitter(size = 1/2, alpha = 1/2) +
    scale_x_sqrt("UMIs on reference") +
    scale_y_sqrt("UMIs on alternate") +
    coord_fixed(ratio = 1, xlim = range(df_var$ref_count), ylim = range(df_var$alt_count)) +
    scale_color_manual(guide = "none", values = c("TRUE" = "black", "FALSE" = "gray")) +
    ggtitle(paste(sum(df_var$lod <= args$lod), "homozygous variants -", sum(df_var$lod > args$lod), "heterozygous variants")) +
    theme_bw(base_size = 16)
}

keep <- df_var$lod > args$lod
df_var <- df_var[keep,]
ret$ref <- ret$ref[keep,]
ret$alt <- ret$alt[keep,]
delta_mat <- ret$ref - ret$alt
write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Removed", sum(!keep), "variants while keeping", sum(keep), "variants"), stderr())

###########################################################################
## RUN CO-CLUSTERING                                                     ##
###########################################################################

df_var$use <- (df_var$ref_count > 0 & df_var$alt_count > 0) & Matrix::rowSums(ret$ref > 0) + Matrix::rowSums(ret$alt > 0) > 1
df_cell$use <- Matrix::colSums(ret$ref > 0) + Matrix::colSums(ret$alt > 0) > 1 & !df_cell$doublet

D <- delta_mat[df_var$use, df_cell$use]
rs <- Matrix::rowSums(ret$ref[df_var$use,df_cell$use]) + Matrix::rowSums(ret$alt[df_var$use,df_cell$use])
cs <- Matrix::colSums(ret$ref[df_var$use,df_cell$use]) + Matrix::colSums(ret$alt[df_var$use,df_cell$use])
res <- co_cluster(D, rs, cs, args$max_iter)

df_var$phase <- rep(NA_integer_, length(df_var$var))
df_var$phase[df_var$use] <- res$phase
df_var$phase[!df_var$use] <- sign(as.vector(delta_mat[!df_var$use,df_cell$use] %*% res$cell_state) / (Matrix::rowSums(ret$ref[!df_var$use,df_cell$use])) + Matrix::rowSums(ret$alt[!df_var$use,df_cell$use]))
df_var$phase[is.na(df_var$phase)] <- 0

df_cell$cell_state <- rep(NA_integer_, length(df_cell$cell))
df_cell$cell_state[df_cell$use] <- res$cell_state
df_cell$cell_state[!df_cell$use] <- sign(as.vector(res$phase %*% delta_mat[df_var$use,!df_cell$use]) / (Matrix::colSums(ret$ref[df_var$use,!df_cell$use])) + Matrix::colSums(ret$alt[df_var$use,!df_cell$use]))
df_cell$cell_state[is.na(df_cell$cell_state)] <- 0

###########################################################################
## GENERATE OUTPUT TABLE FOR VARIANTS                                    ##
###########################################################################

df_var$p1_count <- ( df_var$umi_count + as.vector(delta_mat[,df_cell$use] %*% df_cell$cell_state[df_cell$use])) / 2
df_var$p2_count <- df_var$umi_count - df_var$p1_count
df_var$p1_fraction <- df_var$p1_count / df_var$umi_count
df_var$p1_fraction[is.na(df_var$p1_fraction)] <- .5
df_var$xa_count <- ( df_var$umi_count + as.vector(df_var$phase * delta_mat[,df_cell$use] %*% df_cell$cell_state[df_cell$use]) ) / 2
df_var$xi_count <- df_var$umi_count - df_var$xa_count
df_var$xa_fraction <- df_var$xa_count / df_var$umi_count
df_var$xa_fraction[is.na(df_var$xa_fraction)] <- .5
if (!is.null(args$var)) {
  cols <- c("var", "ref_count", "alt_count", "p1_count", "p2_count")
  write.table(lapply(df_var[,cols], function(x) if (is.numeric(x)) sub("\\.?0+$", "", sprintf("%.4f", x)) else x), args$var, quote = FALSE, sep = "\t", row.names = FALSE)
  write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Generated output table", args$var), stderr())
}

###########################################################################
## GENERATE OUTPUT TABLE FOR CELLS                                       ##
###########################################################################

use <- pmin(df_var$p1_fraction, 1 - df_var$p1_fraction) < args$max_xi_proportion
df_cell$x1_count <- ( Matrix::colSums(ret$ref[use,]) + Matrix::colSums(ret$alt[use,]) + as.vector(df_var$phase[use] %*% delta_mat[use,]) ) / 2
df_cell$x2_count <- Matrix::colSums(ret$ref[use,]) + Matrix::colSums(ret$alt[use,]) - df_cell$x1_count
df_cell$x1_fraction <- df_cell$x1_count / ( df_cell$x1_count + df_cell$x2_count )
df_cell$x1_fraction[is.na(df_cell$x1_fraction)] <- .5
df_cell$xci_lod <- VGAM::dbetabinom(round(df_cell$x1_count), round(df_cell$x1_count + df_cell$x2_count), args$hom_frac, args$rho, log = TRUE) / log(10) -
                   VGAM::dbetabinom(round(df_cell$x2_count), round(df_cell$x1_count + df_cell$x2_count), args$hom_frac, args$rho, log = TRUE) / log(10)
df_cell$dblt_lod <- VGAM::dbetabinom(round(pmin(df_cell$x1_count, df_cell$x2_count)), round(df_cell$x1_count + df_cell$x2_count), args$het_frac, args$rho, log = TRUE) / log(10) -
                    VGAM::dbetabinom(round(pmin(df_cell$x1_count, df_cell$x2_count)), round(df_cell$x1_count + df_cell$x2_count), args$hom_frac, args$rho, log = TRUE) / log(10)
df_cell$active.x <- ifelse(df_cell$x1_fraction > .5, "X1", "X2")
df_cell$active.x[abs(df_cell$xci_lod) < args$lod] <- "NA"
df_cell$active.x[df_cell$dblt_lod > args$lod] <- "XD"
cols <- c("cell", "x1_count", "x2_count", "active.x")
if (args$all) {
  df_cell$x1_all_count <- ( Matrix::colSums(ret$ref) + Matrix::colSums(ret$alt) + as.vector(df_var$phase %*% delta_mat) ) / 2
  df_cell$x2_all_count <- Matrix::colSums(ret$ref) + Matrix::colSums(ret$alt) - df_cell$x1_all_count
  cols <- c(cols, c("x1_all_count", "x2_all_count"))
}
write.table(lapply(df_cell[,cols], function(x) if (is.numeric(x)) sub("\\.?0+$", "", sprintf("%.4f", x)) else x), args$out, quote = FALSE, sep = "\t", row.names = FALSE)
write(paste0(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), " Doublet rate estimated at ", format(100*sum(df_cell$dblt_lod>args$lod)/length(df_cell$cell), digits = 3), "%"), stderr())

active_ratio <- (1 + as.numeric(df_var$phase[df_var$use] %*% delta_mat[df_var$use,df_cell$use] %*% df_cell$cell_state[df_cell$use]) / sum(rs)) / 2
dt <- as.data.table(summary(delta_mat[order(df_var$p1_fraction), order(df_cell$x1_fraction)]))
write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Generated output table", args$out), stderr())

###########################################################################
## GENERATE QUALITY CONTROL PLOTS                                        ##
###########################################################################

if (!is.null(args$pdf)) {
  # matrix two-dimensional embedding
  sv <- RSpectra::svds(Matrix::Diagonal( x = ifelse(rs == 0, 0, 1 / sqrt(rs)) ) %*% D %*%
                       Matrix::Diagonal( x = ifelse(cs == 0, 0, 1 / sqrt(cs)) ), k = 2)
  svd_active_ratio <- (1 + c(as.numeric(sign(sv$u[,1]) %*% D %*% sign(sv$v[,1])), as.numeric(sign(sv$u[,2]) %*% D %*% sign(sv$v[,2]))) / sum(rs)) / 2
  idx <- df_cell$use
  df_cell$v1[idx] <- sv$v[,1]
  df_cell$v2[idx] <- sv$v[,2]
  p2 <- ggplot(df_cell[idx,], aes(x=v1, y=v2, color=active.x)) +
    geom_vline(xintercept = 0, color = "gray", alpha = 1/2) +
    geom_point(size = 1/2, alpha = 1/2) +
    scale_x_continuous("") +
    scale_y_continuous("") +
    scale_color_manual(guide = "none", values = c("X1" = "black", "X2" = "blue", "XD" = "red", "NA" = "gray")) +
    ggtitle(paste("SVD active ratio:", format(max(svd_active_ratio), digits = 4), "-", sum(df_cell$active.x[idx] == "X1") , "X1 -", sum(df_cell$active.x[idx] == "X2") , "X2 -", sum(df_cell$active.x[idx] == "XD"), "doublets -", sum(df_cell$active.x[idx] == "NA"), "NA")) +
    theme_bw(base_size = 16)
  
  # matrix co-clustering
  p3 <- ggplot(dt, aes(j, i, color = sign(x))) +
    rasterise(geom_point(size = 0.1, alpha = 0.2), dpi = 300) +
    scale_x_continuous(NULL, breaks = NULL, expand = FALSE) +
    scale_y_continuous(NULL, breaks = NULL, expand = FALSE) +
    scale_color_gradient2(guide = "none", low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggtitle(paste("Active ratio:", format(active_ratio, digits = 4), "-", length(df_var$var), "variants -", length(df_cell$cell), "cell barcodes")) +
    theme_bw(base_size = 16)

  # plot variants
  p4 <- ggplot(df_var, aes(x=xa_count, y=xi_count, color=pmin(p1_fraction, 1 - p1_fraction) < args$max_xi_proportion)) +
    geom_function(fun = function(x) args$max_xi_proportion * x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) x, color = "gray", alpha = 1/2) +
    geom_jitter(size = 1/2, alpha = 1/2) +
    scale_x_sqrt("UMIs on active X") +
    scale_y_sqrt("UMIs on inactive X") +
    coord_fixed(ratio = 1, xlim = range(df_var$xa_count), ylim = range(df_var$xi_count)) +
    scale_color_manual(guide = "none", values = c("TRUE" = "black", "FALSE" = "gray")) +
    ggtitle(paste(sum(pmin(df_var$p1_fraction, 1 - df_var$p1_fraction) < args$max_xi_proportion), "good variants -", sum(pmin(df_var$p1_fraction, 1 - df_var$p1_fraction) >= args$max_xi_proportion), "bad variants")) +
    theme_bw(base_size = 16)

  # plot cells
  idx <- !df_cell$doublet
  p5 <- ggplot(df_cell[idx,], aes(x=x1_count, y=x2_count, color=active.x)) +
    geom_function(fun = function(x) args$hom_frac/(1-args$hom_frac) * x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) args$het_frac/(1-args$het_frac) * x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) (1-args$het_frac)/args$het_frac * x, color = "gray", alpha = 1/2) +
    geom_function(fun = function(x) (1-args$hom_frac)/args$hom_frac * x, color = "gray", alpha = 1/2) +
    geom_jitter(size = 1/2, alpha = 1/2) +
    scale_x_continuous("UMIs on X1") +
    scale_y_continuous("UMIs on X2") +
    coord_fixed(ratio = 1, xlim = range(df_cell$x1_count), ylim = range(df_cell$x2_count)) +
    scale_color_manual(guide = "none", values = c("X1" = "black", "X2" = "blue", "XD" = "red", "NA" = "gray")) +
    ggtitle(paste(sum(df_cell$active.x[idx] == "X1") , "X1 -", sum(df_cell$active.x[idx] == "X2") , "X2 -", sum(df_cell$active.x[idx] == "XD"), "doublets -", sum(df_cell$active.x[idx] == "NA"), "NA")) +
    theme_bw(base_size = 16)

  # plot prior doublet cells
  if (any(df_cell$doublet)) {
    idx <- df_cell$doublet
    p6 <- ggplot(df_cell[idx,], aes(x=x1_count, y=x2_count, color=active.x)) +
      geom_function(fun = function(x) args$hom_frac/(1-args$hom_frac) * x, color = "gray", alpha = 1/2) +
      geom_function(fun = function(x) args$het_frac/(1-args$het_frac) * x, color = "gray", alpha = 1/2) +
      geom_function(fun = function(x) x, color = "gray", alpha = 1/2) +
      geom_function(fun = function(x) (1-args$het_frac)/args$het_frac * x, color = "gray", alpha = 1/2) +
      geom_function(fun = function(x) (1-args$hom_frac)/args$hom_frac * x, color = "gray", alpha = 1/2) +
      geom_jitter(size = 1/2, alpha = 1/2) +
      scale_x_continuous("UMIs on X1") +
      scale_y_continuous("UMIs on X2") +
      coord_fixed(ratio = 1, xlim = range(df_cell$x1_count), ylim = range(df_cell$x2_count)) +
      scale_color_manual(guide = "none", values = c("X1" = "black", "X2" = "blue", "XD" = "red", "NA" = "gray")) +
      ggtitle(paste("prior doublets:", sum(df_cell$active.x[idx] == "X1") , "X1 -", sum(df_cell$active.x[idx] == "X2") , "X2 -", sum(df_cell$active.x[idx] == "XD"), "doublets -", sum(df_cell$active.x[idx] == "NA"), "NA")) +
      theme_bw(base_size = 16)
  }

  # plot all counts for cells
  if (args$all) {
    p7 <- ggplot(df_cell, aes(x=x1_all_count, y=x2_all_count, color=active.x)) +
      geom_function(fun = function(x) x, color = "gray", alpha = 1/2) +
      geom_jitter(size = 1/2, alpha = 1/2) +
      scale_x_continuous("All UMIs on X1") +
      scale_y_continuous("All UMIs on X2") +
      coord_fixed(ratio = 1, xlim = range(df_cell$x1_all_count), ylim = range(df_cell$x2_all_count)) +
      scale_color_manual(guide = "none", values = c("X1" = "black", "X2" = "blue", "XD" = "red", "NA" = "gray")) +
      ggtitle(paste(sum(df_cell$active.x == "X1") , "X1 -", sum(df_cell$active.x == "X2") , "X2 -", sum(df_cell$active.x == "XD"), "doublets -", sum(df_cell$active.x == "NA"), "NA")) +
      theme_bw(base_size = 16)
  }

  pdf(args$pdf, width = 10, height = 7)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  if (any(df_cell$doublet)) print(p6)
  if (args$all) print(p7)
  invisible(dev.off())
  write(paste(format(.POSIXct(Sys.time(), tz="GMT"), "[%H:%M:%S]"), "Generated output pdf", args$pdf), stderr())
}
