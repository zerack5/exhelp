library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(plyr)
library(dplyr)
library(parallel)
library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ReactomePA)


#Marker: CIDEA, VSTM2A, VSTM2B, ITGB1, 
genelist_converter = function(msigdb_object){
  uniq_pathway = unique(msigdb_object$gs_name)
  genelist = list()
  for(pathway in uniq_pathway){
    genelist[[pathway]] = as.vector(t(msigdb_object$gene_symbol[msigdb_object$gs_name == pathway]))
  }
  return(genelist)
  
}

entrez_to_symbol = function(genes, ref){
  temp = strsplit(as.character(genes), "/")
  temp = lapply(temp, function(x) ref$SYMBOL[match(x, ref$ENTREZID)])
  temp = lapply(temp, function(x) paste(x, collapse = "/"))
  temp = as.vector(t(as.data.frame(temp)))
  return(temp)
}




KEGG_Enrichment = function(Diff_result, savename){
  
  up_regulated_genes = Diff_result$gene[Diff_result$avg_logFC > 0]
  down_regulated_genes = Diff_result$gene[Diff_result$avg_logFC < 0]
  
  up_entrez = mapIds(org.Hs.eg.db, up_regulated_genes, "ENTREZID", "SYMBOL")
  down_entrez = mapIds(org.Hs.eg.db, down_regulated_genes, "ENTREZID", "SYMBOL")
  total_entrez = mapIds(org.Hs.eg.db, (Diff_result$gene), "ENTREZID", "SYMBOL")
  gene_transform = data.frame(SYMBOL = Diff_result$gene, ENTREZID = total_entrez)
  
  up_KEGG = enrichKEGG(up_entrez, organism = "hsa")
  down_KEGG = enrichKEGG(down_entrez, organism = "hsa")
  
  up_KEGG = as.data.frame(up_KEGG@result)
  up_KEGG$geneID = entrez_to_symbol(up_KEGG$geneID, ref = gene_transform)
  
  down_KEGG = as.data.frame(down_KEGG@result)
  down_KEGG$geneID = entrez_to_symbol(down_KEGG$geneID, ref = gene_transform)
  
  sheets = list("KEGG_up" = up_KEGG, "KEGG_down" = down_KEGG)
  write.xlsx(sheets, paste(savename, ".xlsx", sep = ""))
}


GO_Enrichment = function(Diff_result, savename){
  
  up_regulated_genes = Diff_result$gene[Diff_result$avg_logFC > 0]
  down_regulated_genes = Diff_result$gene[Diff_result$avg_logFC < 0]
  
  gene_list = exp(Diff_result$avg_logFC)
  names(gene_list) = as.character(Diff_result$gene)
  gene_list = sort(gene_list, decreasing = T)
  
  BP_up = enrichGO(up_regulated_genes, "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL")
  BP_down = enrichGO(down_regulated_genes, "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL")
  
  print("Finished BP!")
  
  CC_up = enrichGO(up_regulated_genes, "org.Hs.eg.db", ont = "CC", keyType = "SYMBOL")
  CC_down = enrichGO(down_regulated_genes, "org.Hs.eg.db", ont = "CC", keyType = "SYMBOL")
  
  print("Finished CC!")
  
  MF_up = enrichGO(up_regulated_genes, "org.Hs.eg.db", ont = "MF", keyType = "SYMBOL")
  MF_down = enrichGO(down_regulated_genes, "org.Hs.eg.db", ont = "MF", keyType = "SYMBOL")
  
  print("Finished MF!")
  
  
  
  sheets =list("BP_up" = BP_up, "BP_down" = BP_down,
               "CC_up" = CC_up, "CC_down" = CC_down,  
               "MF_up" = MF_up, "MF_down" = MF_down)
  
  write.xlsx(sheets, paste(savename, ".xlsx", sep = ""))
}



MsigDB_Enrichment = function(Diff_result, savename){
  up_regulated_genes = Diff_result$gene[Diff_result$avg_logFC > 0]
  down_regulated_genes = Diff_result$gene[Diff_result$avg_logFC < 0]
  
  
  gmtfile_H = "~/Downloads/PAmission/h.all.v6.2.symbols.gmt"
  gmtfile_C2 = "~/Downloads/PAmission/c2.all.v6.2.symbols.gmt"
  
  H = read.gmt(gmtfile = gmtfile_H)
  C2 = read.gmt(gmtfile = gmtfile_C2)
  
  H_up = enricher(up_regulated_genes, TERM2GENE = H)
  H_down = enricher(down_regulated_genes, TERM2GENE = H)
  print("Finished H!")
  
  C2_up = enricher(up_regulated_genes, TERM2GENE = C2)
  C2_down = enricher(down_regulated_genes, TERM2GENE = C2)
  print("Finished C2!")
  
  sheets =list("H_up" = H_up, "H_down" = H_down,
               "C2_up" = C2_up, "C2_down" = C2_down)
  
  write.xlsx(sheets, paste(savename, ".xlsx", sep = ""))
  
}


Wiki_Pathway_Enrichment = function(Diff_result, savename){
  up_regulated_genes = Diff_result$gene[Diff_result$avg_logFC > 0]
  down_regulated_genes = Diff_result$gene[Diff_result$avg_logFC < 0]
  
  up_entrez = mapIds(org.Hs.eg.db, up_regulated_genes, "ENTREZID", "SYMBOL")
  down_entrez = mapIds(org.Hs.eg.db, down_regulated_genes, "ENTREZID", "SYMBOL")
  total_entrez = mapIds(org.Hs.eg.db, Diff_result$gene, "ENTREZID", "SYMBOL")
  gene_transform = data.frame(SYMBOL = Diff_result$gene, ENTREZID = total_entrez)
  
  wp2gene = read.gmt("~/Downloads/PAmission/wikipathways-20190610-gmt-Homo_sapiens.gmt")
  wp2gene = wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
  wpid2gene = wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name = wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
  
  Wiki_up = enricher(up_entrez, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)@result
  Wiki_down = enricher(down_entrez, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)@result
  
  Wiki_up$geneID = entrez_to_symbol(Wiki_up$geneID, ref = gene_transform)
  Wiki_down$geneID = entrez_to_symbol(Wiki_down$geneID, ref = gene_transform)
  sheets =list("Wiki_up" = Wiki_up, "Wiki_down" = Wiki_down)
  
  write.xlsx(sheets, paste(savename, ".xlsx", sep = ""))
  
}
#DEG_Heatmap = function(seurat_object, Diff_result, n_genes, clusters){

#}

Reactome_Enrichment = function(Diff_result, savename){
  up_regulated_genes = Diff_result$gene[Diff_result$avg_logFC > 0]
  down_regulated_genes = Diff_result$gene[Diff_result$avg_logFC < 0]
  
  up_entrez = mapIds(org.Hs.eg.db, up_regulated_genes, "ENTREZID", "SYMBOL")
  down_entrez = mapIds(org.Hs.eg.db, down_regulated_genes, "ENTREZID", "SYMBOL")
  total_entrez = mapIds(org.Hs.eg.db, Diff_result$gene, "ENTREZID", "SYMBOL")
  gene_transform = data.frame(SYMBOL = Diff_result$gene, ENTREZID = total_entrez)
  
  up_Reactome = enrichPathway(up_entrez)
  down_Reactome = enrichPathway(down_entrez)
  
  up_Reactome = as.data.frame(up_Reactome@result)
  up_Reactome$geneID = entrez_to_symbol(up_Reactome$geneID, ref = gene_transform)
  
  down_Reactome = as.data.frame(down_Reactome@result)
  down_Reactome$geneID = entrez_to_symbol(down_Reactome$geneID, ref = gene_transform)
  
  sheets = list("Reactome_up" = up_Reactome, "Reactome_down" = down_Reactome)
  write.xlsx(sheets, paste(savename, ".xlsx", sep = ""))
  
}

