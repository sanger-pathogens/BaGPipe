#!/usr/bin/env Rscript

# manhattan_kmers.R
# Render a Manhattan plot for bacterial GWAS k-mers from either Phandango .plot or annotated_kmers.tsv

# Set library path for user-writable directory (in case container has read-only system paths)
tmplib <- "/tmp/R_lib"
dir.create(tmplib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(tmplib, .libPaths()))

# Ensure required packages are available (install if missing inside the container)
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = tmplib, repos = "https://cran.r-project.org", quiet = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

ensure_pkg("ggplot2")
ensure_pkg("ggrepel")
suppressPackageStartupMessages(library(grid))

# ---- Simple argument parsing (no external deps) ----
args <- commandArgs(trailingOnly = TRUE)
usage <- "
Usage:
  Rscript manhattan_kmers.R --plot <mapped_kmers.plot> [--threshold-file <kmer_pattern_count.txt>] [--out-prefix <name>] [--label-top N]
or
  Rscript manhattan_kmers.R --annot <annotated_kmers.tsv> [--threshold-file <kmer_pattern_count.txt>] [--out-prefix <name>] [--label-top N]

Arguments:
  --plot            Phandango .plot file from pyseer phandango_mapper
  --annot           annotated_kmers.tsv from annotate_hits_pyseer
  --threshold-file  File that contains a line 'Threshold:\\t<value>'
  --out-prefix      Prefix for outputs (default: 'manhattan')
  --label-top       Integer: label top N points above threshold (default: 0)
"
getOpt <- function(key, default = NULL) {
  i <- which(args == key)
  if (length(i) == 0) return(default)
  if (i == length(args)) return(TRUE) # flags
  return(args[i + 1])
}
plot_file       <- getOpt("--plot", NULL)
annot_file      <- getOpt("--annot", NULL)
thresh_file     <- getOpt("--threshold-file", NULL)
gff_file        <- getOpt("--gff", NULL)  # Optional GFF for gene annotation
out_prefix      <- getOpt("--out-prefix", "manhattan")
label_top       <- as.integer(getOpt("--label-top", "0"))

if (is.null(plot_file) && is.null(annot_file)) {
  cat(usage, "\n"); stop("Provide either --plot or --annot\n")
}
if (!is.null(plot_file) && !file.exists(plot_file)) stop("Plot file not found: ", plot_file)
if (!is.null(annot_file) && !file.exists(annot_file)) stop("Annotation file not found: ", annot_file)
if (!is.null(thresh_file) && !file.exists(thresh_file)) stop("Threshold file not found: ", thresh_file)

# ---- Helpers ----
read_genome_length <- function(path) {
  # Read genome length from GFF ##sequence-region header
  # Returns total genome length in bp
  if (is.null(path) || !file.exists(path)) return(NULL)
  
  header_lines <- readLines(path, n = 100)
  seq_regions <- grep("^##sequence-region", header_lines, value = TRUE)
  if (length(seq_regions) == 0) return(NULL)
  
  # Parse ##sequence-region lines (format: ##sequence-region seqname start end)
  total_length <- 0
  for (line in seq_regions) {
    parts <- strsplit(line, "\\s+")[[1]]
    if (length(parts) >= 4) {
      end_pos <- as.numeric(parts[4])
      if (!is.na(end_pos) && end_pos > total_length) {
        total_length <- end_pos
      }
    }
  }
  
  if (total_length > 0) {
    message("Genome length from GFF: ", format(total_length, big.mark = ","), " bp (", 
            round(total_length / 1e6, 2), " Mb)")
    return(total_length)
  }
  return(NULL)
}

read_gff_genes <- function(path) {
  # Read GFF and extract genes with coordinates
  # Returns data.frame with columns: chr, start, end, gene_name
  if (is.null(path) || !file.exists(path)) return(NULL)
  
  message("Reading GFF file: ", path)
  gff <- read.delim(path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  # GFF columns: seqname, source, feature, start, end, score, strand, frame, attributes
  if (ncol(gff) < 9) return(NULL)
  
  # Extract CDS and gene features
  gff_genes <- gff[gff[[3]] %in% c("CDS", "gene"), ]
  if (nrow(gff_genes) == 0) return(NULL)
  
  # Parse gene name from attributes (e.g., 'gene="pbpX"' or 'gene=pbpX')
  genes_list <- lapply(strsplit(gff_genes[[9]], ";"), function(x) {
    # Find gene= pattern (with or without quotes)
    gene_match <- grep('gene=', x, value = TRUE)
    if (length(gene_match) > 0) {
      # Try with quotes first, then without
      gene_name <- sub('.*gene="([^"]+)".*', '\\1', gene_match[1])
      if (gene_name == gene_match[1]) {
        # No quotes, extract after gene=
        gene_name <- sub('.*gene=([^;]+).*', '\\1', gene_match[1])
      }
      return(gene_name)
    }
    return(NA_character_)
  })
  
  gene_df <- data.frame(
    chr = gff_genes[[1]],
    start = as.numeric(gff_genes[[4]]),
    end = as.numeric(gff_genes[[5]]),
    gene_name = unlist(genes_list),
    stringsAsFactors = FALSE
  )
  gene_df <- gene_df[!is.na(gene_df$gene_name), ]
  message("Loaded ", nrow(gene_df), " genes from GFF")
  return(gene_df)
}

find_overlapping_gene <- function(chr, pos, gff_genes) {
  # Find the gene that overlaps with the given position
  if (is.null(gff_genes)) return(NA_character_)
  
  matches <- gff_genes[
    gff_genes$chr == chr & 
    gff_genes$start <= pos & 
    gff_genes$end >= pos,
  ]
  
  if (nrow(matches) == 0) return(NA_character_)
  # Return first match (or combine if multiple)
  return(paste(unique(matches$gene_name), collapse = ";"))
}

read_threshold <- function(path) {
  if (is.null(path)) return(NA_real_)
  x <- readLines(path, warn = FALSE)
  hit <- grep("^\\s*Threshold:", x)
  if (length(hit) == 0) return(NA_real_)
  # Expect 'Threshold:\t<value>'
  val <- sub("^\\s*Threshold:\\s*", "", x[hit[1]])
  as.numeric(val)
}

make_cumulative_positions <- function(df, chr_col = "CHR", bp_col = "BP") {
  # Order contigs by appearance, compute per-contig offsets, midpoints for x-axis
  chr_levels <- unique(df[[chr_col]])
  df[[chr_col]] <- factor(df[[chr_col]], levels = chr_levels)
  contig_max <- tapply(df[[bp_col]], df[[chr_col]], max, na.rm = TRUE)
  offsets <- c(0, cumsum(head(contig_max, -1)))
  names(offsets) <- chr_levels
  df$BP_CUM <- df[[bp_col]] + offsets[as.character(df[[chr_col]])]
  axis_df <- data.frame(
    CHR = chr_levels,
    midpoint = offsets + contig_max / 2
  )
  list(df = df, axis_df = axis_df)
}

# ---- Ingest data ----
dat <- NULL
name_col <- NULL

if (!is.null(plot_file)) {
  # Phandango .plot: header often '#CHR SNP BP minLOG10(P) log10(p) r^2'
  # We'll be flexible to tabs/spaces and column names.
  raw <- read.delim(plot_file, header = TRUE, check.names = FALSE)
  # Keep track of original column names before normalization
  raw_colnames <- colnames(raw)
  
  # Normalise names (strip leading '#', spaces)
  cn <- gsub("^#","", colnames(raw))
  cn <- gsub("\\s+"," ", cn)            # collapse multiple spaces
  cn <- gsub("\\(.*\\)","", cn)         # drop bracket notes
  colnames(raw) <- cn

  # Identify key columns
  CHR <- if ("CHR" %in% cn) "CHR" else colnames(raw)[1]
  BP  <- if ("BP"  %in% cn) "BP"  else colnames(raw)[3]
  
  # Check if BP column contains only dots (phandango_mapper quirk)
  # If so, the actual positions are in column 3
  bp_test <- as.character(raw[[BP]][1:min(10, nrow(raw))])
  if (all(bp_test == "." | is.na(bp_test))) {
    message("BP column contains only dots; using column 3 for positions")
    BP <- colnames(raw)[3]
  }
  
  # y column: prefer minLOG10 or a -log10 column; otherwise compute -log10(p) if 'p' exists
  ycol <- NA_character_
  y <- NULL
  
  # Special case: if BP column is minLOG10, the y values are in the next column
  if (grepl("minLOG10", BP, ignore.case = TRUE)) {
    message("Position column is labeled as minLOG10; using next column for y-axis")
    # Find the actual column index in the original raw dataframe
    bp_idx_in_raw <- which(raw_colnames == raw_colnames[which(cn == BP)])
    if (bp_idx_in_raw < length(raw_colnames)) {
      # Get the actual column name from original raw dataframe
      ycol_raw <- raw_colnames[bp_idx_in_raw + 1]
      message("Using column '", ycol_raw, "' for y-axis")
      y <- as.numeric(raw[[cn[bp_idx_in_raw + 1]]])  # Use normalized name to access
    }
  }
  
  # If y wasn't set by special case, try standard column detection
  if (is.null(y)) {
    if ("minLOG10 P" %in% cn) ycol <- "minLOG10 P"
    if ("minLOG10"   %in% cn) ycol <- "minLOG10"
    if ("log10 p"    %in% cn) ycol <- "log10 p"
    if ("log10(p)"   %in% cn) ycol <- "log10(p)"

    if (!is.na(ycol)) {
      y <- as.numeric(raw[[ycol]])
      # If it's log10(p) (positive towards zero), convert to -log10(p)
      # Heuristic: if max(y) < 0, negate; if column name contains "log10 p" (no minus), negate.
      if (grepl("log10", ycol, ignore.case = TRUE) && !grepl("min", ycol, ignore.case = TRUE)) {
        y <- -y
      }
    } else {
      # Fall back: look for a 'p' column and transform
      pcol <- grep("pvalue|p-value|pval|p", cn, ignore.case = TRUE, value = TRUE)
      if (length(pcol) == 0) stop("Could not infer a p-value/log10 column in .plot")
      y <- -log10(as.numeric(raw[[pcol[1]]]))
    }
  }

  # Parse BP column - may contain ranges like "1991755..1991812" or single positions
  bp_raw <- as.character(raw[[BP]])
  bp_values <- sapply(bp_raw, function(x) {
    if (grepl("..", x, fixed = TRUE)) {
      # Range format: start..end - use midpoint
      parts <- strsplit(x, "..", fixed = TRUE)[[1]]
      start <- as.numeric(parts[1])
      end <- as.numeric(parts[2])
      return(round((start + end) / 2))
    } else {
      # Single position
      return(as.numeric(x))
    }
  }, USE.NAMES = FALSE)
  
  dat <- data.frame(
    CHR = as.character(raw[[CHR]]),
    BP  = as.numeric(bp_values),
    LOGP = as.numeric(y)
  )
  # Avoid using non-gene identifiers for labels
  name_col <- NULL

} else {
  # annotated_kmers.tsv: columns from annotate_hits_pyseer (often no header)
  # Format: variant af filter-pvalue lrt-pvalue beta beta-std-err variant_h2 notes ... mapping
  # We'll read without header and assign expected column names
  ann <- read.delim(annot_file, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # annotated_kmers.tsv typically has no header; assign standard pyseer column names
  # Columns: 1=variant, 2=af, 3=filter-pvalue, 4=lrt-pvalue, 5=beta, 6=beta-std-err, 7=variant_h2, ..., last=mapping
  if (ncol(ann) < 8) stop("annotated_kmers.tsv has fewer than expected columns")
  
  colnames(ann)[1:8] <- c("variant", "af", "filter-pvalue", "lrt-pvalue", "beta", "beta-std-err", "variant_h2", "notes")
  # Last column is the mapping (CHR:BP-BP format)
  colnames(ann)[ncol(ann)] <- "mapping"
  
  # Try to discover p-value column (prefer lrt-pvalue over filter-pvalue)
  pcol <- c("lrt-pvalue", "filter-pvalue")
  pcol <- pcol[pcol %in% colnames(ann)]
  if (length(pcol) == 0) stop("No p-value column found in annotated_kmers.tsv (expected 'lrt-pvalue' or 'filter-pvalue' at columns 3-4)")
  
  # Mapping column is the last column (CHR:start-end format)
  map_col <- "mapping"
  first_hit <- sub(",.*$", "", ann[[map_col]])
  chr <- sub(":.*$", "", first_hit)
  start <- as.numeric(sub(".*:(\\d+)-\\d+.*", "\\1", first_hit))
  end   <- as.numeric(sub(".*:\\d+-(\\d+).*", "\\1", first_hit))
  bp_mid <- suppressWarnings(round((start + end) / 2))
  
  dat <- data.frame(
    CHR = chr,
    BP  = bp_mid,
    LOGP = -log10(as.numeric(ann[[pcol[1]]]))
  )
  # Try to extract gene names from GFF if available, otherwise from mapping column
  gff_genes <- read_gff_genes(gff_file)
  
  if (!is.null(gff_genes)) {
    # Look up genes from GFF using genomic positions
    genes <- mapply(function(c, p) {
      find_overlapping_gene(c, p, gff_genes)
    }, chr, start, SIMPLIFY = TRUE)
    dat$gene <- genes
    message("Annotated ", sum(!is.na(genes)), " k-mers from GFF file")
    # Restrict to contigs present in the GFF to avoid sample-specific IDs
    valid_contigs <- unique(gff_genes$chr)
    before_n <- nrow(dat)
    dat <- dat[dat$CHR %in% valid_contigs, ]
    after_n <- nrow(dat)
    if (after_n < before_n) {
      message("Filtered out ", before_n - after_n, " points on non-GFF contigs")
    }
  } else {
    # Fall back to extracting from mapping column
    genes <- sapply(strsplit(ann$mapping, ";"), function(x) {
      if (length(x) > 1) {
        gene <- gsub('"', '', x[2])
        return(ifelse(nchar(gene) > 0, gene, NA_character_))
      }
      return(NA_character_)
    })
    dat$gene <- genes
    message("Using gene names from annotation mapping column (GFF not available)")
  }
  name_col <- "gene"
}

# Remove NAs and non-finite
initial_rows <- nrow(dat)
dat <- dat[is.finite(dat$BP) & is.finite(dat$LOGP) & !is.na(dat$CHR), ]
final_rows <- nrow(dat)
message("Filtered data: ", initial_rows, " -> ", final_rows, " rows (removed ", initial_rows - final_rows, " with NA/non-finite values)")

if (final_rows == 0) {
  stop("No valid data points after filtering. Check that BP column contains numeric values and LOGP is finite.")
}

dat <- dat[order(dat$CHR, dat$BP), ]

# Compute cumulative positions for x‑axis
cp <- make_cumulative_positions(dat, chr_col = "CHR", bp_col = "BP")
dat <- cp$df
axis_df <- cp$axis_df
message("Created plot with ", final_rows, " data points across ", nrow(axis_df), " contigs")

# Read threshold and compute -log10 line
thresh <- read_threshold(thresh_file)
thresh_line <- if (!is.na(thresh) && is.finite(thresh) && thresh > 0) -log10(thresh) else NA_real_

# Colour group alternating by contig index (blue and grey for visual distinction across contigs)
dat$grp <- as.integer(factor(dat$CHR)) %% 2

# Optionally select top labels above threshold
to_label <- NULL
if (!is.na(thresh_line)) {
  above <- dat[dat$LOGP >= thresh_line, ]
  if (label_top > 0 && nrow(above) > 0) {
    ord <- order(above$LOGP, decreasing = TRUE)
    pick <- above[ord, ][seq_len(min(label_top, nrow(above))), ]
    to_label <- pick
    # Use gene names if available; otherwise use CHR:BP
    # Prefer gene names; fallback to CHR:BP
    genes <- if ("gene" %in% colnames(above)) as.character(above$gene)[ord][seq_len(nrow(to_label))] else rep(NA_character_, nrow(to_label))
    to_label$label <- ifelse(is.na(genes) | genes == "NA",
                              paste0(to_label$CHR, ":", to_label$BP),
                              genes)
  }

  # If no mapped points are above the significance threshold, skip output plots
  # Only check for gene overlap if we're using annotated file (has gene column)
  no_sig_points <- (nrow(above) == 0)
  no_sig_genes  <- !is.null(annot_file) && ("gene" %in% colnames(dat)) && (sum(!is.na(above$gene) & above$gene != "NA") == 0)
  if (no_sig_points || no_sig_genes) {
    msg <- if (no_sig_points) {
      "No mapped k-mers above threshold; skipping plot outputs."
    } else {
      "No significant k-mers overlap annotated genes; skipping plot outputs."
    }
    message(msg)
    skip_file <- paste0(out_prefix, ".skip.txt")
    writeLines(c(
      paste0("Skipped Manhattan plotting for ", out_prefix),
      paste0("Reason: ", msg, " (threshold=", sprintf("%.2g", thresh), ")")
    ), con = skip_file)
    message("Wrote: ", skip_file)
    # Exit early to avoid generating empty/irrelevant plots
    quit(save = "no", status = 0)
  }
}

# ---- Plot ----
# Get genome length from GFF to set x-axis limits
genome_length <- read_genome_length(gff_file)

# Create full genome Manhattan plot
p <- ggplot(dat, aes(x = BP_CUM, y = LOGP, colour = factor(grp))) +
  geom_point(alpha = 0.7, size = 0.8) +
  scale_colour_manual(values = c("#3B4A5A", "#1E88E5"), guide = "none") +  # grey and blue for alternating contigs
  scale_x_continuous(
    name = "Genomic position (bp)",
    labels = scales::unit_format(unit = "Mb", scale = 1e-6),
    limits = if (!is.null(genome_length)) c(0, genome_length) else NULL
  ) +
  ylab(expression(-log[10](p))) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(colour = "lightgrey", size = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

if (!is.na(thresh_line)) {
  # Use 0 as x position if we have genome length, otherwise use min data position
  x_pos <- if (!is.null(genome_length)) 0 else min(dat$BP_CUM, na.rm = TRUE)
  p <- p + geom_hline(yintercept = thresh_line, linetype = "dashed", colour = "firebrick", size = 0.6) +
    annotate("text", x = x_pos, y = thresh_line,
             label = sprintf("Threshold = %.2g (–log10 = %.2f)", thresh, thresh_line),
             hjust = 0, vjust = -0.7, size = 3, colour = "firebrick")
}

if (!is.null(to_label) && nrow(to_label) > 0) {
    p <- p + ggrepel::geom_label_repel(
        data = to_label,
        aes(x = BP_CUM, y = LOGP, label = label),
        size = 2.6,    
        box.padding   = grid::unit(0.35, "lines"),    
        point.padding = grid::unit(0.20, "lines"),    
        label.padding = grid::unit(0.15, "lines"),    
        min.segment.length = grid::unit(0, "mm"),
        label.size = 0.15,    
        fill  = "white",    
        colour = "black"  
    )
}


# ---- Output ----
pdf_file <- paste0(out_prefix, ".manhattan.pdf")
png_file <- paste0(out_prefix, ".manhattan.png")

message("Saving plots...")
tryCatch({
  ggsave(pdf_file, p, width = 10, height = 4, units = "in", device = "pdf")
  message("Wrote: ", pdf_file)
}, error = function(e) {
  message("ERROR saving PDF: ", e$message)
})

tryCatch({
  ggsave(png_file, p, width = 10, height = 4, units = "in", dpi = 300, device = "png")
  message("Wrote: ", png_file)
}, error = function(e) {
  message("ERROR saving PNG: ", e$message)
  # Try alternative: use ragg if available
  tryCatch({
    ggsave(png_file, p, width = 10, height = 4, units = "in", dpi = 300, device = "ragg_png")
    message("Wrote: ", png_file, " (using ragg)")
  }, error = function(e2) {
    message("ERROR: Could not save PNG with ragg either: ", e2$message)
  })
})

