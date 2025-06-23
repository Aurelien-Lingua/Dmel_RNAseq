
##### Loading needed values for functions to work independently
rootdir <- getwd()
Robjectsdir <- paste0(rootdir, "/Robjects")
if (!exists("t2g")){
  t2g <- readRDS(file = paste0(Robjectsdir, "/t2g_dmel.rds"))
}
if (!exists("logCPM")){
  logCPM <- readRDS(file = paste0(Robjectsdir, "/logCPM.rds"))
  avlogCPM <- readRDS(file = paste0(Robjectsdir, "/avlogCPM.rds"))
}
if (!exists("contrasts_list")){
  contrasts_list <- readRDS(file = paste0(Robjectsdir, "/contrasts.rds"))
}
if (!exists("glmQL_res")){
  glmQL_res <- readRDS(file = paste0(Robjectsdir, "/glmQL_res_txi.rds"))
  glmQL_decided <- readRDS(file = paste0(Robjectsdir, "/glmQL_decided_txi.rds"))
  Treat_res <- readRDS(file = paste0(Robjectsdir, "/Treat_res_txi.rds"))
  Treat_decided <- readRDS(file = paste0(Robjectsdir, "/Treat_decided_txi.rds"))
  QL_res <- glmQL_res
  QL_decided <- glmQL_decided
}
if (!exists("gsego_res")){
  gsego_res <- readRDS(file = paste0(Robjectsdir, "/GSEA_results_GO_glmQL-txi.rds"))
  gsekegg_res <- readRDS(file = paste0(Robjectsdir, "/GSEA_results_K_glmQL-txi.rds"))
}
if (!exists("ORA_results_GO")){
  ORA_results_GO <- readRDS(file = paste0(Robjectsdir, "/ORA_results_GO_glmQL-txi.rds"))
  ORA_results_K <- readRDS(file = paste0(Robjectsdir, "/ORA_results_K_glmQL-txi.rds"))
}
##### Functions


.hunt <- function(x, contr, .t2g = t2g){
  found <- t2g[t2g$target_id %in% "empty",]
  # looping until the entry is found in t2g
  while(dim(found)[1]==0){
    # separating IDs if a list was given
    if(length(x)==1){
      if(grepl("/",x)|grepl(" ",x)){
        try(x <- strsplit(x, "/"))
        try(x <- strsplit(x, " "))
        x <- unlist(x)
      }
    }
    if (grepl("FBtr",x[[1]],fixed = TRUE)){
      found <- t2g[t2g$target_id %in% x,]
    } else if (grepl("FBgn",x[[1]],fixed = TRUE)){
      found <- t2g[t2g$ens_gene %in% x,]
    } else if (x[[1]] %in% t2g$entrez_gene){
      found <- t2g[t2g$entrez_gene %in% x,]
    } else if (x[[1]] %in% t2g$ext_gene){
      found <- t2g[t2g$ext_gene %in% x,]
    } else {
      cat("Please enter either a gene ID in the form of transcript id (FBtrxx), flybase gene ID (Fbgnxxx), gene symbol (e.g. pug), or entrezID (numbers):")
      x<- readline()
    }
  }
  # getting contrast values for each gene
  temp_table <- QL_decided[[contr]]$table[QL_decided[[contr]]$table$target_id %in% found$target_id,]
  # reordering the table to follow found's order, so values are correctly assigned
  temp_table<-temp_table[order(factor(temp_table, levels = found$target_id)),]
  found <- temp_table[!duplicated(temp_table$ens_gene),]
  return(found)
}

# function saves the Flybase gene IDs in found object for batch search on flyatlas
FBgns_save <- function(){
  label <- paste0(length(found$ens_gene),"_",gsub(":",".",Sys.time()))
  
  write.csv(found$ens_gene,
            quote = FALSE,
            file = paste0("C:/Users/Aurelien/Expansion/proj28/plots/tximport plots/batch_searches/",
                          label,".csv"))
  
}

# Call function to assign the contrast of interest.
wcontrast <- function(x = 404){
  x<-as.integer(x)
  num_c <- length(QL_res)
  message <- paste("Please enter a number between 1 to",num_c,"to choose your contrast. Below are the contrast options:\n")
  for (j in 0:3){  # Iterate over the 4 lines
    for (i in seq(1 + j, num_c, by = 4)){  # Print contrasts separated by 4
      message <- paste0(message, i, " - ", names(QL_res)[i], "\t")
    }
    message <- paste0(message, "\n")
  }
  while ( x > num_c| x < 1| is.na(x)){
    cat(message)
    x<- as.integer(readline())
  }
  assigned <<- names(QL_res)[x]
}
# function to use in other functions to return contrasts without assigning contrast globally
.wcontrast <- function(x = 404){
  x<-as.integer(x)
  num_c <- length(QL_res)
  message <- paste("Please enter a number between 1 to",num_c,"to choose your contrast. Below are the contrast options:\n")
  for (j in 0:3){  # Iterate over the 4 lines
    for (i in seq(1 + j, num_c, by = 4)){  # Print contrasts separated by 4
      message <- paste0(message, i, " - ", names(QL_res)[i], "\t")
    }
    message <- paste0(message, "\n")
  }
  while ( x > num_c | x < 1 | is.na(x)){
    cat(message)
    x<- as.integer(readline())
  }
  assigned <-  names(QL_res)[x]
  return(assigned)
}
# function to get all transcript IDs associated with entered gene ID (ens,ext,entrez)
.ID_to_trans <- function(x="404"){
  looking <- t2g[t2g$target_id %in% "empty",]
  # looping until the entry is found in t2g
  while(dim(looking)[1]==0){
    # separating IDs if a list was given
    if (grepl("/",x)|grepl(" ",x)){
      try(x <- strsplit(x, "/"))
      try(x <- strsplit(x, " "))
      x <- unlist(x)
    }
    if (grepl("FBgn",x[[1]],fixed = TRUE)){
      looking <- t2g[t2g$ens_gene %in% x,]
    } else if (x[[1]] %in% t2g$entrez_gene){
      looking <- t2g[t2g$entrez_gene %in% x,]
    } else if (x[[1]] %in% t2g$ext_gene){
      looking <- t2g[t2g$ext_gene %in% x,]
    } else {
      cat("Please enter gene ID in the form of flybase gene ID (Fbgnxxx), gene symbol (e.g. pug), or entrezID (numbers):")
      x<- readline()
    }
  }
  return(looking)
}
# function to retrieve the average logCPM of a list of input transcripts.
gavlogCPM <- function(){
  # retrieving list of transcript IDs with associated gene symbol
  IDs <- .ID_to_trans()
  buffer <- readline("buffer. Just press enter if required.")
  # getting contrast to extract sample groups from
  contrast <- .wcontrast()
  # getting sample groups from contrasts matrix
  sample_groups <- contrasts_list[,contrast]
  sample_groups <- names(sample_groups[sample_groups!=0])
  # getting avlogCPM values for the selected transcripts and sample groups
  avlogCPM_vals <- avlogCPM[IDs$ens_gene,sample_groups]
  # renaming rows to include the gene symbol in them
  rownames(avlogCPM_vals) <- paste0(IDs$ens_gene," (",IDs$ens_gene,")")
  # removing NA rows which come from transcript IDs pulled out of t2g that did not pass the filtering steps. These transcripts don't exist in y so they don't exist in logCPM
  avlogCPM_vals <- avlogCPM_vals[!is.na(avlogCPM_vals[,1]),]
  View(avlogCPM_vals)
  avlogCPM_vals <<- avlogCPM_vals
}

# Function finds given gene entries in t2g, saves the list globally and opens it to view. you can enter one ID or multiple IDs separated by a "/"
hunt <- function(x = "geneID"){
  looking <- t2g[t2g$target_id %in% "empty",]
  # looping until the entry is found in t2g
  while(dim(looking)[1]==0){
    # separating IDs if a list was given
    if(length(x)==1){
      if(grepl("/",x)|grepl(" ",x)){
        try(x <- strsplit(x, "/"))
        try(x <- strsplit(x, " "))
        x <- unlist(x)
      }
    }
    if (grepl("FBtr",x[[1]],fixed = TRUE)){
      looking <- t2g[t2g$target_id %in% x,]
      found <<- t2g[t2g$target_id %in% x,]
    } else if (grepl("FBgn",x[[1]],fixed = TRUE)){
      looking <- t2g[t2g$ens_gene %in% x,]
      found <<- t2g[t2g$ens_gene %in% x,]
    } else if (x[[1]] %in% t2g$entrez_gene){
      looking <- t2g[t2g$entrez_gene %in% x,]
      found <<- t2g[t2g$entrez_gene %in% x,]
    } else if (x[[1]] %in% t2g$ext_gene){
      looking <- t2g[t2g$ext_gene %in% x,]
      found <<- t2g[t2g$ext_gene %in% x,]
    } else {
      cat("Please enter either a gene ID in the form of transcript id (FBtrxx), flybase gene ID (Fbgnxxx), gene symbol (e.g. pug), or entrezID (numbers):")
      x<- readline()
    }
  }
  # I want to add to the found object a 5th column with information on its expression in the currently assigned contrast, i.e. logFC and FDR for the contrast
  # running which contrast function if not already done
  if (!exists("assigned")){
    .wcontrast()
  }
  # getting contrast values for each gene
  temp_table <- QL_decided[[assigned]]$table[QL_decided[[assigned]]$table$target_id %in% found$target_id,]
  
  
  # reordering the table to follow found's order, so values are correctly assigned
  temp_table<<-temp_table[order(factor(temp_table, levels = found$target_id)),]
  View(temp_table)
  found <<- temp_table[!duplicated(temp_table$ens_gene),]
  View(found)
}

.jhunt <- function(x="404", contrast = NULL){
  # retrieving list of transcript IDs with associated gene symbol
  if(x=="404"){
    IDs <- .ID_to_trans()
  }else{IDs <- .ID_to_trans(x)}
  # buffer <- readline("buffer. Just press enter if required.")
  # getting contrast to extract sample groups from
  contrast <- .wcontrast(contrast)
  # getting sample groups from contrasts matrix
  sample_groups <- contrasts_list[,contrast]
  sample_groups <- names(sample_groups[sample_groups!=0])
  # getting avlogCPM values for the selected transcripts and sample groups
  # if statement to check if I'm working with transcript-level or gene-level information
  avlogCPM_vals <- avlogCPM[IDs$ens_gene,sample_groups]
  # removing NA rows which come from transcript IDs pulled out of t2g that did not pass the filtering steps. These transcripts don't exist in y so they don't exist in logCPM
  avlogCPM_vals <- avlogCPM_vals[!is.na(avlogCPM_vals[,1]),]
  
  DGEvals <-.hunt(IDs$ens_gene, contrast = contrast)
  DGEvals[,6:8] <- list(NULL)
  newtable <- merge(DGEvals, avlogCPM_vals, by=0)
  newtable[,1] <- NULL
  return(newtable)
}

jhunt <- function(x="404"){
  # retrieving list of transcript IDs with associated gene symbol
  if(x=="404"){
    IDs <- .ID_to_trans()
  }else{IDs <- .ID_to_trans(x)}
  buffer <- readline("buffer. Just press enter if required.")
  # getting contrast to extract sample groups from
  contrast <- .wcontrast()
  # getting sample groups from contrasts matrix
  sample_groups <- contrasts_list[,contrast]
  sample_groups <- names(sample_groups[sample_groups!=0])
  # getting avlogCPM values for the selected transcripts and sample groups
  avlogCPM_vals <- avlogCPM[IDs$ens_gene,sample_groups]
  # removing NA rows which come from transcript IDs pulled out of t2g that did not pass the filtering steps. These transcripts don't exist in y so they don't exist in logCPM
  avlogCPM_vals <- avlogCPM_vals[!is.na(avlogCPM_vals[,1]),]
  colnames(avlogCPM_vals) <- paste0("logCPM in ",colnames(avlogCPM_vals))
  DGEvals <-.hunt(x=IDs$ens_gene, contr = contrast)
  DGEvals[,6:8] <- list(NULL)
  newtable <- merge(DGEvals, avlogCPM_vals, by=0)
  newtable[,1] <- NULL
  View(newtable)
}
