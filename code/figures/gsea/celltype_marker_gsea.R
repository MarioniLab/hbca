#! make gsea of marker genes.
# C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/code/figures/gsea/celltype_marker_gsea.R
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(clusterProfiler))
suppressMessages(library(biomaRt))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(DOSE))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(readr))



##Functions

test_gsea <- function(DEGs, top100_genes, prefix) {
  organism <- "org.Hs.eg.db"
  
  # dge_pwd <- paste0(prefix, '/dge_testing/', test_var)
  # DEGs <- read.csv(paste0(dge_pwd, '/dge-block_', block_var, '-all.csv'))
  DEGs.ord <- DEGs[order(-DEGs[,"logFC"]), ]
  ranks <- DEGs.ord$logFC
  names(ranks) <- rD[DEGs.ord$X, 'gene_ids'] #want ensembl id for obtaining entrez id
  top100_genes_LFC <- DEGs[DEGs$X %in% top100_genes, 'logFC']
  names(top100_genes_LFC) <-  rD[top100_genes, 'gene_ids']
  #head(ranks)
  
  #get entrez ids from gene symbols
  ensembl <- useDataset("hsapiens_gene_ensembl",mart = useMart("ensembl"))
  gene_names <- getBM(attributes=c("ensembl_gene_id",'entrezgene_id', "hgnc_symbol"), mart = ensembl)
  ranks_entrez <- ranks
  names(ranks_entrez) <- mapvalues(names(ranks_entrez), from = gene_names$ensembl_gene_id, to = gene_names$entrezgene_id)
  top100_genes_LFC_entrez <- top100_genes_LFC
  names(top100_genes_LFC_entrez) <- mapvalues(names(top100_genes_LFC), from = gene_names$ensembl_gene_id, to = gene_names$entrezgene_id)
  
  DEGs.ord <- DEGs.ord[names(ranks) %in% names(na.omit(ranks)),]
  gene_list <- na.omit(ranks)
  gene_list_entrez <- na.omit(ranks_entrez)
  
  #GO
  enrich_GO <- enrichGO(gene = names(top100_genes_LFC),
                        universe = names(gene_list),
                        keyType = 'ENSEMBL',
                        OrgDb = org.Hs.eg.db, # organism,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05,
                        readable = TRUE,
                        minGSSize = 5)
  
  gse_GO <- gseGO(geneList=gene_list, 
                  ont ="ALL", 
                  keyType = "ENSEMBL",
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = FALSE, 
                  OrgDb = org.Hs.eg.db,#organism, 
                  pAdjustMethod = "none")
  
  #KEGG
  options(clusterProfiler.download.method = "wget")  #this fixed error with KEGG download
  enrich_KEGG <- enrichKEGG(gene = names(top100_genes_LFC_entrez),
                            universe = names(gene_list_entrez),
                            organism = 'hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.05,
                            minGSSize = 5)
  gse_KEGG <- gseKEGG(geneList=gene_list_entrez,
                      organism = 'hsa',
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = FALSE)
  
  #DO (disease ontology)
  enrich_DO <- enrichDO(gene = names(top100_genes_LFC_entrez),
                        universe = names(gene_list_entrez),
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05,
                        minGSSize = 5, 
                        readable = TRUE)
  gse_DO <- gseDO(geneList=gene_list_entrez,
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = FALSE)
  
  #DGN (DisGeNET - another disease network)
  enrich_DGN <- enrichDGN(gene = names(top100_genes_LFC_entrez),
                          universe = names(gene_list_entrez),
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05,
                          minGSSize = 5, 
                          readable = TRUE)
  gse_DGN <- gseDGN(geneList=gene_list_entrez,
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = FALSE)
  
  #save plots 
  save_pw <- paste0(prefix, '/clusterProfiler_plotting/')
  dir.create(save_pw, showWarnings = FALSE, recursive = TRUE)
  
  #GO
  if (!is.null(enrich_GO)) {
    if (dim(enrich_GO)[1] > 0) {
      pdf(paste0(save_pw, '/dotplot_GO_enrichment.pdf'), width = 8, height = 20)
      print(dotplot(enrich_GO, showCategory=30) + ggtitle("Dotplot for GO enrichment"))
      dev.off()
    }
  }
  if (dim(gse_GO)[1] > 0) {
    pdf(paste0(save_pw, '/dotplot_GO_gsea.pdf'), width = 8, height = 20)
    print(dotplot(gse_GO, showCategory=30) + ggtitle("Dotplot for GO gsea"))
    dev.off()
  }
  
  #KEGG
  if (!is.null(enrich_KEGG)) {
    if (dim(enrich_KEGG)[1] > 0) {
      pdf(paste0(save_pw, '/dotplot_KEGG_enrichment.pdf'))
      print(dotplot(enrich_KEGG, showCategory=30) + ggtitle("Dotplot for KEGG enrichment"))
      dev.off()
    }
  }
  if (dim(gse_KEGG)[1] > 0) {
    pdf(paste0(save_pw, '/dotplot_KEGG_gsea.pdf'), width = 8, height = 15)
    print(dotplot(gse_KEGG, showCategory = 30, title = "Dotplot for KEGG gsea" , split=".sign") + facet_grid(.~.sign))
    dev.off()
  }
  
  #DO
  if (!is.null(enrich_DO)) {
    if (dim(enrich_DO)[1] > 0) {
      pdf(paste0(save_pw, '/dotplot_DO_enrichment.pdf'))
      print(dotplot(enrich_DO, showCategory=30) + ggtitle("Dotplot for DO enrichment"))
      dev.off()
    }
  }
  if (dim(gse_DO)[1] > 0) {
    pdf(paste0(save_pw, '/dotplot_DO_gsea.pdf'), width = 8, height = 15)
    print(dotplot(gse_DO, showCategory=30) + ggtitle("Dotplot for DO gsea"))
    dev.off()
  }
  #DGN
  if (!is.null(enrich_DGN)) {
    if (dim(enrich_DGN)[1] > 0) {
      pdf(paste0(save_pw, '/dotplot_DGN_enrichment.pdf'))
      print(dotplot(enrich_DGN, showCategory=30) + ggtitle("Dotplot for DGN enrichment"))
      dev.off()
    }
  }
  if (dim(gse_DGN)[1] > 0){
    pdf(paste0(save_pw, '/dotplot_DGN_gsea.pdf'), width = 8, height = 15)
    print(dotplot(gse_DGN, showCategory=30) + ggtitle("Dotplot for DGN gsea"))
    dev.off()
  }
  
  # save data 
  save_pw <- paste0(prefix, '/clusterProfiler_save_S4/')
  dir.create(save_pw, showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(enrich_KEGG)) {
    enrich_GO_df <- data.frame(enrich_KEGG@result)
  }
  if (!is.null(enrich_GO)) {
    enrich_KEGG_df <- data.frame(enrich_KEGG@result)
  }
  if (!is.null(enrich_DO)) {
    enrich_DO_df <- data.frame(enrich_KEGG@result)
  }
  if (!is.null(enrich_DGN)) {
    enrich_DGN_df <- data.frame(enrich_KEGG@result)
  }
  
  df <- fortify(gse_KEGG, showCategory = 1000, split=".sign")
  
  gse_GO_df <- data.frame(gse_GO@result)
  gse_KEGG_df <- data.frame(gse_KEGG@result)
  gse_DO_df <- data.frame(gse_DO@result)
  gse_DGN_df <- data.frame(gse_DGN@result)
  
  readr::write_tsv(enrich_GO_df, paste0(save_pw, '/GO_enrich.tsv'))
  readr::write_tsv(enrich_KEGG_df, paste0(save_pw, '/KEGG_enrich.tsv'))
  readr::write_tsv(enrich_DO_df, paste0(save_pw, '/DO_enrich.tsv'))
  readr::write_tsv(enrich_DGN_df, paste0(save_pw, '/DGN_enrich.tsv'))
  
  readr::write_tsv(gse_GO_df, paste0(save_pw, '/GO_enrich.tsv'))
  readr::write_tsv(gse_KEGG_df, paste0(save_pw, '/KEGG_enrich.tsv'))
  readr::write_tsv(gse_DO_df, paste0(save_pw, '/DO_enrich.tsv'))
  readr::write_tsv(gse_DGN_df, paste0(save_pw, '/DGN_enrich.tsv'))
  
  #save(enrich_GO, ascii=TRUE, file=paste0(save_pw, '/GO_enrich'))
}




### Analysis

## IFNG+ T cells
input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/imm_2/output'
sce <- readRDS(paste0(input_path, '/input/input_sce.rds'))
rownames(sce) <- rowData(sce)$X
rD <- rowData(sce)
prefix_imm <- '/nfs/research/marioni/areed/projects/hbca/figures/src/gsea/output/imm'

DEGs_imm <- read.csv('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/nosubset/level2_nosub/dge-IFNG+ T-all.csv') 
DEGS_imm_up <- read.csv('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/nosubset/level2_nosub/dge-IFNG+ T-up.csv')
top100_genes_imm <- DEGS_imm_up[1:min(100, sum(DEGS_imm_up$FDR < 0.05)),1]

test_gsea(DEGs = DEGs_imm, 
          top100_genes = top100_genes_imm,
          prefix = prefix_imm)


## LP4
prefix_lp4 <- '/nfs/research/marioni/areed/projects/hbca/figures/src/gsea/output/lp4'

DEGs_lp4 <- read.csv('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/subset/level2/dge-LP4-all.csv') 
DEGS_lp4_up <- read.csv('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/subset/level2/dge-LP4-up.csv')
top100_genes_lp4 <- DEGS_lp4_up[1:min(100, sum(DEGS_lp4_up$FDR < 0.05)),1]

test_gsea(DEGs = DEGs_lp4, 
          top100_genes = top100_genes_lp4,
          prefix = prefix_lp4)


## HS2
prefix_hs2 <- '/nfs/research/marioni/areed/projects/hbca/figures/src/gsea/output/hs2'

DEGs_hs2 <- read.csv('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/subset/level2/dge-HS2-all.csv') 
DEGS_hs2_up <- read.csv('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/subset/level2/dge-HS2-up.csv')
top100_genes_hs2 <- DEGS_hs2_up[1:min(100, sum(DEGS_hs2_up$FDR < 0.05)),1]

test_gsea(DEGs = DEGs_hs2, 
          top100_genes = top100_genes_hs2,
          prefix = prefix_hs2)



## BSL2

prefix_bsl2 <- '/nfs/research/marioni/areed/projects/hbca/figures/src/gsea/output/bsl2'

DEGs_bsl2 <- read.csv('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/subset/level2/dge-BSL2-all.csv') 
DEGS_bsl2_up <- read.csv('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/subset/level2/dge-BSL2-up.csv')
top100_genes_bsl2 <- DEGS_bsl2_up[1:min(100, sum(DEGS_bsl2_up$FDR < 0.05)),1]

test_gsea(DEGs = DEGs_bsl2, 
          top100_genes = top100_genes_bsl2,
          prefix = prefix_bsl2)

















