## WGCNA_KEGGenrichment.R ##
## To test whether genes in a module enrichment in KEGG pathways ##
## Zhang Fei ##
## 2017-11-28 ##

# Rscript WGCNA_KEGGenrichment.R <dir.module> <bg.size> <dir.outsidefile>
argv <- commandArgs(T)

workingdir <- "."
setwd(workingdir)


#fl.module should write the whole path, such as "~/winhome/inMP/time-series/GCN/WGCNA/WGCNA_pipeline/result/modules/merged_modes.dissTOM.JS.pearson.deepsplit=2.MEDissThres=0.10.txt"

# path.module <- "H://inMP/time-series/DESeq2_timecourse/pathway_enrichment/different_eachtime_HZLY1.txt"
# outpath <- "H://inMP/time-series/DESeq2_timecourse/pathway_enrichment/KEGG_HZLY1_tbt.txt"
#path.module <- argv[1]

path.module <- "H://inMP/time-series/DESeq2_timecourse/pathway_enrichment/different_eachtime_JI853.txt"
outpath <- "H://inMP/time-series/DESeq2_timecourse/pathway_enrichment/KEGG_JI853_tbt.txt"

#### read kegg pathways ####

my.loadpathway <- function(){
    list.pathway <- list()
    dir.pathway <- dir("H://inMP/tools/kegg_pathway_geneget/result/")
    names.pathway <- NULL
    for(i in 1: length(dir.pathway)){
         info.pathway <- read.table(paste("H://inMP/tools/kegg_pathway_geneget/result/", dir.pathway[i], sep = ""), quote = "", header = T, sep = "\t")
         list.pathway[[i]] <- info.pathway
         names <- strsplit(dir.pathway[i], ".txt", perl = T)[[1]][1]
         names <- strsplit(names, "_", perl = T)[[1]][2]
         names.pathway <- append(names.pathway, names)
    }
    names(list.pathway) <- names.pathway
    return(list.pathway)
}


#### read modules ####
my.loadmodule <- function(path.module){
    table.module <- read.table(path.module, header = F, sep = "\t")
    colnames(table.module) <- c("Gene", "Module")
    return(table.module)
}

#### Calculate pathway gene numbers

# my.bgsummary <- function(list.pathway){
#     genelist <- NULL
#     for(i in 1: length(list.pathway)){
#         genes <- as.character(list.pathway[[i]][, 2])
#         genelist <- append(genelist, genes)
#         
#     }
#     #remove gene "-", can't trans entrz id to GRAME ID
#     gene.number <- length(genelist[!duplicated(genelist)]) - 1
#     return(gene.number)
# }

#### fisher test and FDR test ####

list.pathway <- my.loadpathway()
table.module <- my.loadmodule(path.module)
modulename <- levels(table.module$Module)
#num.bg <- as.integer(argv[2])
num.bg <- 35581

num <- length(modulename)*length(list.pathway)

result.p <- matrix(rep(NA, num), nrow = length(list.pathway))
rownames(result.p) <- names(list.pathway)
colnames(result.p) <- modulename

result.detail <- c(1, 1, 1, 1, 1, 1, 1)


for(i in 1: length(modulename)){
    module.gene <- table.module[table.module[ , 2] == modulename[i], ]
    modulesize <- dim(module.gene)[1]
    for(j in 1: length(list.pathway)){
        gene.intersect <- module.gene[module.gene[ , 1] %in% list.pathway[[j]][ , 2], ][ , 1]
        pathwaysize <- dim(list.pathway[[j]])[1]
        a <- length(gene.intersect)
        b <- modulesize - a
        c <- pathwaysize - a
        d <- num.bg - modulesize - c
        mat <- matrix(c(a, b, c, d), nrow = 2, dimnames = list(c('In pathway', 'Not in pathway'), c('Module', 'Not in module')))
        p.value <- fisher.test(mat, alternative = "greater")$p.value
        result.p[j, i] <- p.value
        result.detail <- data.frame(result.detail, c(modulename[i], names(list.pathway)[j], a, modulesize, pathwaysize, num.bg, p.value))
        # gene.intersect <- module.gene[module.gene[ , 1] %in% list.pathway[[j]][ , 2], ][ , 1]
        # size.intersect <- length(gene.intersect)
        # module.in <- size.intersect
        # module.not <- size.module - size.intersect
        # bg.in <- dim(list.pathway[[j]])[1]
        # bg.out <- num.bg - size.module - bg.in + size.intersect
        # mat <- matrix(c(module.in, module.not, bg.in, bg.out), nrow = 2, dimnames = list(c('In pahtway', 'Not in pathway'), c('Module', 'Background')))
        # p.value <- fisher.test(mat, alternative = "greater")$p.value
        # result.p[j, i] <-  p.value
        # result.detail <- data.frame(result.detail, c(modulename[i], names(list.pathway)[j], size.intersect, size.module, bg.in, num.bg, p.value))
    }
}
result.detail <- t(result.detail)
result.detail <- result.detail[-1, ]
rownames(result.detail) <- c(1: num)
fdr <- p.adjust(result.detail[ , 7], method = 'fdr')
result.detail <- data.frame(module = as.character(result.detail[ , 1]), pathway = as.character(result.detail[ , 2]), n = as.integer(result.detail[ , 3]), N = as.integer(result.detail[ , 4]), x = as.integer(result.detail[ , 5]), X = as.integer(result.detail[ , 6]), p.value = as.numeric(result.detail[ , 7]), fdr = fdr)
colnames(result.detail) <- c('module', 'pathway', 'n', 'N', 'x', 'X', 'p.value', 'fdr')
result.detail <- result.detail[order(result.detail[ , 8], decreasing =  F), ]
write.table(result.detail, outpath, quote = F, row.names = F, col.names = T, sep = "\t" )
