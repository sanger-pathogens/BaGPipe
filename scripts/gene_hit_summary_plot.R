#!/usr/bin/env Rscript

require(ggplot2)
require(ggrepel)
library(ggrepel)

gene_hits = read.table("gene_hits.tsv", stringsAsFactors=FALSE, header=TRUE, sep="\t")

pdf("Gene_Hit_Plot_small.pdf", width = 12, height = 9) 

ggplot(gene_hits, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
   geom_point(alpha=0.5) +
   geom_text_repel(aes(size=60), show.legend = FALSE, colour='black') +
   scale_size("Number of unitigs/k-mers", range=c(1,10)) +
   scale_colour_gradient('Average MAF') +
   theme_bw(base_size=14) +
   ggtitle("Your tested phenotype") +
   xlab("Average effect size") +
   ylab("Maximum -log10(p-value)")

dev.off()

pdf("Gene_Hit_Plot_medium.pdf", width = 36, height = 27)

ggplot(gene_hits, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
   geom_point(alpha=0.5) +
   geom_text_repel(aes(size=60), show.legend = FALSE, colour='black') +
   scale_size("Number of unitigs/k-mers", range=c(1,10)) +
   scale_colour_gradient('Average MAF') +
   theme_bw(base_size=18) +
   ggtitle("Your tested phenotype") +
   xlab("Average effect size") +
   ylab("Maximum -log10(p-value)")

dev.off()

pdf("Gene_Hit_Plot_large.pdf", width = 60, height = 45)

ggplot(gene_hits, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
   geom_point(alpha=0.5) +
   geom_text_repel(aes(size=60), show.legend = FALSE, colour='black') +
   scale_size("Number of unitigs/k-mers", range=c(1,10)) +
   scale_colour_gradient('Average MAF') +
   theme_bw(base_size=22) +
   ggtitle("Your tested phenotype") +
   xlab("Average effect size") +
   ylab("Maximum -log10(p-value)")

dev.off()
