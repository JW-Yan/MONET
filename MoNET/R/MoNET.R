

#######################################################
#Funciton1. MultiOmic Network Construction
#######################################################
#' @title Multi-omic network construction
#'
#' @description Remove nodes whose gene interaction score less than the threshold from the muiti-omic network we manually integrated.
#'
#' @param threshold A constant, gives a threshhold for removing nodes whose gene interaction score less than it, with the default value 0.
#'
#' @return a dataframe of network with two columns
#' @export BuildNet
#'
#' @examples
#' MulOmicNet<-BuildNet(threshold=0.5)
#'
BuildNet <-function(threshold=0){
  if(threshold<0|threshold>1){
    stop("Invalid 'threshold' argument. It's interval is [0,1]")
    geterrmessage()
  }
  #AllSNPGene<-readRDS(file = "AllSNPGene.rds")
  #AllGeneGene<-readRDS(file = "AllGeneGene.rds")
  AllGeneGene_temp <- AllGeneGene[which(AllGeneGene$Score>=threshold),c(1,2)]
  #AllGeneReac<-readRDS(file = "AllGeneReac.rds")
  #AllReacMeta<-readRDS(file = "AllReacMeta.rds")
  #AllMetaReac<-readRDS(file = "AllMetaReac.rds")
  WholeNet <- do.call("rbind", list(AllSNPGene, AllGeneGene_temp, AllGeneReac, AllReacMeta, AllMetaReac))
  #rm(list = c('AllSNPGene','AllGG','AllGeneGene','AllGeneReac','AllReacMeta','AllMetaReac'))
  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("Package \"igraph\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  MulOmicNet <-igraph::graph_from_data_frame(d=WholeNet, directed=T)
  MulOmicNet <- igraph::simplify(MulOmicNet, remove.multiple = T, remove.loops = T,edge.attr.comb=c(weight="max", type="ignore"))
  print(sprintf("The Multi-Omic Network with gene interaction score lareger than %1.3f has been constructed.",threshold))
  rm(list = c('AllGeneGene_temp','WholeNet'))
  return(MulOmicNet)
}

##############################################################
#Function2. List the Query Format of gene/metabolite
##############################################################
#' @title Query the types of gene/metabolite ID
#'
#' @description Provide users to look up available input ID for genes or metabolites.
#'
#' @param what Character constant, 'gene'and 'metabolite' represent the querying types of gene ID and metabolite ID respectively.
#'
#' @details Gene ID includes Gene_Entrez_ID, Gene_Name, Ensembl_Gene_ID, Enzyme_IDs and RefSeq_IDs.
#' Gene_Entrez_ID: stable identifiers for genes at the NCBI (e.g.,51166).
#' Gene_Name: the official gene symbol that has been approved by the HGNC (e.g.,AADAT).
#' Ensembl_Gene_ID:curated ensembl gene ID (e.g.,ENSG00000109576).
#' Enzyme_IDs:Enzyme entries have Enzyme Commission (EC) numbers associated with them that indicate the hierarchical
#' functional classes to which they belong. This field can contain multiple values as a comma delimited list (e.g.,2.6.1.39, 2.6.1.7).
#' RefSeq_IDs:The Reference Sequence (RefSeq) identifier for that entry, provided by the NCBI. This field may contain multiple
#' values as a comma delimited list (e.g.,NM_016228).
#' Metabolite ID includes Metabolites, Names, Formulas, HMDB_ID, KEGG_ID, PubChem_ID and CHEBI_ID.
#' Metabolites:descriptive name of metabolite (e.g.,h2o[c]).
#' Names:full names of metabolite (e.g.,Water).
#' Formulas:fuomula of metabolite (e.g.,H2O).
#' HMDB_ID:identifier of metaboite from human metabolome database (e.g.,HMDB02111).
#' KEGG_ID:identifier metabolite from KEGG database (e.g.,C00001).
#' PubChem_ID:identifier of metabolite from PubChem database (e.g.,962).
#' CHEBI_ID:identifier of metabolite from Chemical Entities of Biological Interest(ChEBI)database (e.g.,15377).
#'
#' @return character vector
#' @export ListFilters
#'
#' @examples
#' ListFilters("gene")
#'
ListFilters <- function(what){
  if (what[1]=="gene"){print(sprintf("Gene_Entrez_ID, Gene_Name, Ensembl_Gene_ID, Enzyme_IDs, RefSeq_IDs"))}
  else if (what[1]=="metabolite"){print(sprintf("Metabolites, Names, Formulas, HMDB_ID, KEGG_ID, PubChem_ID, CHEBI_ID"))}
  else {
    stop("Invalid 'what' argument. It should be 'gene' or 'metabolite'.")
    geterrmessage()
    }
}

#######################################################
#Funciton3. Shortest Path Identification-SNPquery
#######################################################
#' @title Query shortest paths from SNPs
#'
#' @description Taking SNPs as input, one can get paths in our multi-omic network from SNPs to genes or metabolites.
#'
#' @param MulOmicNet A dataframe of network with two columns.
#' @param snplist A dataframe of SNP IDs with one column.
#' @param node Character constant, gives the shortest paths from SNP to which type of nodes should be calculated.
#' If 'gene',the default,then the shortest paths from SNPs to genes will be calculated.
#' If 'metabolite',then the shortest paths from SNPs to metabolites will be calculated.
#' @param minpath the lower limit of the shortest path length,the default value is 1.
#' @param maxpath the upper limit of the shortest path length,the default value is 3.
#'
#' @return A list is returned, each list element contains a shortest path from one SNP to one gene or from one SNP to one metabolite.
#' Also a text file is generated with the list written to.
#' @export SNPquery
#'
#' @examples
#' MulOmicNet<-BuildNet(threshold=0.5)
#' snplist<- as.data.frame(c("rs1000313","rs1001106"))
#' PathList<-SNPquery(MulOmicNet=MulOmicNet,snplist=snplist,node="gene",minpath=3,maxpath=4)
#'
SNPquery <- function(MulOmicNet,snplist,node=c("gene","metabolite"),minpath=1,maxpath=3){
  colnames(snplist)<-"SNP"
  #SNPAnno <- readRDS(file = "SNPAnno.rds")
  SNPNode <- as.character(unique(SNPAnno[SNPAnno$SNP %in% snplist$SNP,1]))
  if (length(SNPNode)==0){
    print(sprintf("No snps are found in MulOmicNet."))
    return()
    }
  else{
    if (node[1]=="gene"){
      #GeneAnno <- readRDS(file="GeneAnno.rds")
      GeneNode <- as.character(GeneAnno[GeneAnno$Gene_Entrez_ID %in% igraph::V(MulOmicNet)$name,1])
      print(sprintf("Identification of paths from %d snps to %d genes are about to start.",length(SNPNode),length(GeneNode)))
      PathList <-list()
      j<-0
      for (i in 1:length(SNPNode)){
        distance_m<-igraph::distances(MulOmicNet,v=SNPNode[i],to=GeneNode,mode="out")
        nodeto<-GeneNode[which(distance_m<=maxpath & distance_m>=minpath)]
        if(length(nodeto)==0){next}
        j<-j+1
        path_list <- igraph::all_shortest_paths(MulOmicNet,from=SNPNode[i],to=nodeto,mode ="out")
        PathList[[j]]<-path_list$res
        print(i)
      }
    }
    else if (node[1]=="metabolite"){
      #MetaAnno <- readRDS(file="MetaAnno.rds")
      MetaNode <- as.character(MetaAnno$Metabolite_ID)
      print(sprintf("Identification of paths from %d snps to 4901 metabolites are about to start.",length(SNPNode)))
      PathList <-list()
      j<-0
      for (i in 1:length(SNPNode)){
        distance_m<-igraph::distances(MulOmicNet,v=SNPNode[i],to=MetaNode,mode="out")
        nodeto<-MetaNode[which(distance_m<=maxpath & distance_m>=minpath)]
        if(length(nodeto)==0){next}
        j<-j+1
        path_list <- igraph::all_shortest_paths(MulOmicNet,from=SNPNode[i],to=nodeto,mode ="out")
        PathList[[j]]<-path_list$res
        print(i)
      }
    }
    else{
      stop("Invalid 'node' argument. It should be 'gene' or 'metabolite'.")
      geterrmessage()
      }
    if (length(PathList)==0){
      print(sprintf("End of searching. No paths are found."))
      return()
    }
    else{
      print(sprintf("End of searching. Writing paths with length between %d and %d into text file.",minpath,maxpath))
      for(i in (1:length(PathList))){
        temp1<-lapply(PathList[[i]], igraph::as_ids)
        temp2<-lapply(temp1,paste,collapse=",")
        lapply(temp2, write, "Paths.txt", append=TRUE)
      }
    }
  }
  return(PathList)
}

#######################################################
#Funciton4. Shortest Path Identification-Genequery
#######################################################
#' @title Query shortest paths from/to genes
#'
#' @description Taking genes as input, one can get paths in our multi-omic network from SNPs to genes or from genes to metabolites.
#'
#' @param MulOmicNet A dataframe of network with two columns
#' @param genelist A dataframe of genes with one column.
#' @param filter Character constant,the type of input gene ID.Gene ID includes Gene_Entrez_ID, Gene_Name, Ensembl_Gene_ID, Enzyme_IDs and RefSeq_IDs.
#' For the detail,please see the function ListFilters.
#' @param node Character constant, gives the shortest paths between genes to which type of nodes should be calculated.
#' If 'metabolite',the default,then the shortest paths from genes to metabolites will be calculated.
#' If 'snp',then the shortest paths from SNPs to genes will be calculated.
#' @param minpath the lower limit of the shortest path length,the default value is 1.
#' @param maxpath the upper limit of the shortest path length,the default value is 3.
#'
#' @return A list is returned, each list element contains a shortest path from one SNP to one gene or from one gene to one metabolite.
#' Also a text file is generated with the list written to.
#' @export Genequery
#'
#' @examples
#' MulOmicNet<-BuildNet(threshold=0.5)
#' genelist<- as.data.frame(c("A1BG","A1BG-AS1","A2M","A2ML1"))
#' PathList<-Genequery(MulOmicNet=MulOmicNet,genelist=genelist,filter="Gene_Name",
#'                     node="metabolite",minpath=1,maxpath=3)
#'
Genequery <- function(MulOmicNet,genelist,filter="",node=c("metabolite","snp"),minpath=1,maxpath=3){
  GeneAnno_temp <- GeneAnno[GeneAnno$Gene_Entrez_ID %in% igraph::V(MulOmicNet)$name,]
  if (filter[1]=="RefSeq_IDs"){
    colnames(genelist)<-"RefSeq_IDs"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$RefSeq_IDs %in% genelist$RefSeq_IDs,1])
    }
  else if (filter[1]=="Gene_Name"){
    colnames(genelist)<-"Gene_Name"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Gene_Name %in% genelist$Gene_Name,1])
    }
  else if (filter[1]=="Ensembl_Gene_ID"){
    colnames(genelist)<-"Ensembl_Gene_ID"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Ensembl_Gene_ID %in% genelist$Ensembl_Gene_ID,1])
    }
  else if (filter[1]=="Enzyme_IDs"){
    colnames(genelist)<-"Enzyme_IDs"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Enzyme_IDs %in% genelist$Enzyme_IDs,1])
    }
  else if(filter[1]=="Gene_Entrez_ID"){
    colnames(genelist)<-"Gene_Entrez_ID"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Gene_Entrez_ID %in% genelist$Gene_Entrez_ID,1])
    }
  else{
    stop("Invalid 'filter' argument. Please refer to function ListFilters for valid input.")
    geterrmessage()
  }
  if (length(GeneNode)==0){
    print(sprintf("No genes are found in MulOmicNet."))
    return()
  }
  else{
    if (node[1]=="snp"){
      #SNPAnno <- readRDS(file = "SNPAnno.rds")
      SNPNode <- as.character(unique(SNPAnno$SNP))
      print(sprintf("Identification of paths from 1373322 snps to %d genes are about to start.",length(GeneNode)))
      PathList <-list()
      j<-0
      for (i in 1:length(SNPNode)){
        distance_m<-igraph::distances(MulOmicNet,v=SNPNode[i],to=GeneNode,mode="out")
        nodeto<-GeneNode[which(distance_m<=maxpath & distance_m>=minpath)]
        if(length(nodeto)==0){next}
        j<-j+1
        path_list <- igraph::all_shortest_paths(MulOmicNet,from=SNPNode[i],to=nodeto,mode ="out")
        PathList[[j]]<-path_list$res
        print(i)
      }
    }
    else if (node[1]=="metabolite"){
      #MetaAnno <- readRDS(file="MetaAnno.rds")
      MetaNode <- as.character(MetaAnno$Metabolite_ID)
      print(sprintf("Identification of paths from %d genes to 4901 metabolites are about to start.",length(GeneNode)))
      PathList <-list()
      j<-0
      for (i in 1:length(GeneNode)){
        distance_m<-igraph::distances(MulOmicNet,v=GeneNode[i],to=MetaNode,mode="out")
        nodeto<-MetaNode[which(distance_m<=maxpath & distance_m>=minpath)]
        if(length(nodeto)==0){next}
        j<-j+1
        path_list <- igraph::all_shortest_paths(MulOmicNet,from=GeneNode[i],to=nodeto,mode ="out")
        PathList[[j]]<-path_list$res
        print(i)
      }
    }
    else{
      stop("Invalid 'node' argument. It should be 'snp' or 'metabolite'.")
      geterrmessage()
    }
    if (length(PathList)==0){
      print(sprintf("End of searching. No paths are found."))
      return()
    }
    else{
      print(sprintf("End of searching. Writing paths with length between %d and %d into text file.",minpath,maxpath))
      for(i in (1:length(PathList))){
        temp1<-lapply(PathList[[i]], igraph::as_ids)
        temp2<-lapply(temp1,paste,collapse=",")
        lapply(temp2, write, "Paths.txt", append=TRUE)
      }
    }
  }
  return(PathList)
}

#######################################################
#Funciton5. Shortest Path Identification-Metaquery
#######################################################
#' @title Query shortest paths to metabolites
#'
#' @description Taking metabolites as input, one can get paths in our multi-omic network from SNPs or from genes to metabolites.
#'
#' @param MulOmicNet A dataframe of network with two columns.
#' @param metalist A dataframe of metabolites with one column.
#' @param filter Character constant, the type of input metabolite ID.Metabolite ID includes Metabolites, Names, Formulas, HMDB_ID, KEGG_ID, PubChem_ID and CHEBI_ID.
#' For the detail,please see the function ListFilters.
#' @param node character constant, gives the shortest paths from which type of nodes to metabolites should be calculated.
#' If 'gene',the default,then the shortest paths from genes to metabolites will be calculated.
#' If 'snp',then the shortest paths from SNPs to metabolites will be calculated.
#' @param minpath the lower limit of the shortest path length,the default value is 1.
#' @param maxpath the upper limit of the shortest path length,the default value is 3.
#'
#' @return A list is returned, each list element contains a shortest path from one SNP to one metabolite or from one gene to one metabolite.
#' Also a text file is generated with the list written to.
#' @export Metaquery
#'
#' @examples
#' MulOmicNet<-BuildNet(threshold=0.5)
#' metalist<- as.data.frame(c("10-Formyltetrahydrofolate","Water"))
#' PathList<-Metaquery(MulOmicNet=MulOmicNet,metalist=metalist,filter="Names",
#'                     node="gene",minpath=1,maxpath=2)
#'
Metaquery <- function(MulOmicNet,metalist,filter="",node=c("gene","snp"),minpath=1,maxpath=3){
  if (filter[1]=="Names"){
    colnames(metalist)<-"Names"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Names %in% metalist$Names,1])
  }
  else if (filter[1]=="Formulas"){
    colnames(metalist)<-"Formulas"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Formulas %in% metalist$Formulas,1])
  }
  else if (filter[1]=="HMDB_ID"){
    colnames(metalist)<-"HMDB_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$HMDB_ID %in% metalist$HMDB_ID,1])
  }
  else if (filter[1]=="PubChem_ID"){
    colnames(metalist)<-"PubChem_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$PubChem_ID %in% metalist$PubChem_ID,1])
  }
  else if (filter[1]=="KEGG_ID"){
    colnames(metalist)<-"KEGG_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$KEGG_ID %in% metalist$KEGG_ID,1])
  }
  else if (filter[1]=="CHEBI_ID"){
    colnames(metalist)<-"CHEBI_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$CHEBI_ID %in% metalist$CHEBI_ID,1])
  }
  else if (filter[1]=="Metabolites"){
    colnames(metalist)<-"Metabolites"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Metabolites %in% metalist$Metabolites,1])
  }
  else{
    stop("Invalid 'filter' argument. Please refer to function ListFilters for valid input.")
    geterrmessage()
  }
  if (length(MetaNode)==0){
    print(sprintf("No metabolites are found in MulOmicNet."))
    return()
  }
  else{
    if (node[1]=="gene"){
      #GeneAnno <- readRDS(file="GeneAnno.rds")
      GeneNode <- as.character(GeneAnno[GeneAnno$Gene_Entrez_ID %in% igraph::V(MulOmicNet)$name,1])
      print(sprintf("Identification of paths from %d genes to %d metabolites are about to start.",length(GeneNode),length(MetaNode)))
      PathList <-list()
      j<-0
      for (i in 1:length(GeneNode)){
        distance_m<-igraph::distances(MulOmicNet,v=GeneNode[i],to=MetaNode,mode="out")
        nodeto<-MetaNode[which(distance_m<=maxpath & distance_m>=minpath)]
        if(length(nodeto)==0){next}
        j<-j+1
        path_list <- igraph::all_shortest_paths(MulOmicNet,from=GeneNode[i],to=nodeto,mode ="out")
        PathList[[j]]<-path_list$res
        print(i)
      }
    }
    else if (node[1]=="snp") {
      #SNPAnno <- readRDS(file="SNPAnno.rds")
      SNPNode <- as.character(unique(SNPAnno$SNP))
      print(sprintf("Identification of paths from 1373322 snps to %d metabolites are about to start.",length(MetaNode)))
      PathList <-list()
      j<-0
      for (i in 1:length(SNPNode)){
        distance_m<-igraph::distances(MulOmicNet,v=SNPNode[i],to=MetaNode,mode="out")
        nodeto<-MetaNode[which(distance_m<=maxpath & distance_m>=minpath)]
        if(length(nodeto)==0){next}
        j<-j+1
        path_list <- igraph::all_shortest_paths(MulOmicNet,from=SNPNode[i],to=nodeto,mode ="out")
        PathList[[j]]<-path_list$res
        print(i)
      }
    }
    else {
      stop("Invalid 'node' argument. It should be 'snp' or 'gene'.")
      geterrmessage()
      }
    if (length(PathList)==0){
      print(sprintf("End of searching. No paths are found"))
      return()
    }
    else{
      print(sprintf("End of searching. Writing paths with length between %d and %d into text file.",minpath,maxpath))
      for(i in (1:length(PathList))){
        temp1<-lapply(PathList[[i]], igraph::as_ids)
        temp2<-lapply(temp1,paste,collapse=",")
        lapply(temp2, write, "Paths.txt", append=TRUE)
      }
    }
  }
  return(PathList)
}

##########################################################
#Function6. Shortest Path Identification-SNP_Meta_query
##########################################################
#' @title Query shortest paths from SNPs to Metabolites
#'
#' @description Taking SNPs and metabolites as input, one can get paths in our multi-omic network from SNPs to metabolites.
#'
#' @param MulOmicNet A dataframe of network with two columns.
#' @param snplist A dataframe of SNPs with one column.
#' @param metalist a dataframe of metabolites with one column.
#' @param filter Character constant, the type of input metabolite ID.Metabolite ID includes Metabolites, Names, Formulas, HMDB_ID, KEGG_ID, PubChem_ID and CHEBI_ID.
#' For the detail,please see the function ListFilters.
#' @param minpath the lower limit of the shortest path length,the default value is 1.
#' @param maxpath the upper limit of the shortest path length,the default value is 3.
#'
#' @return A list is returned, each list element contains a shortest path from one SNP to one metabolite.
#' Also a text file is generated with the list written to.
#' @export SNP_Meta_query
#'
#' @examples
#' MulOmicNet<-BuildNet(threshold=0.5)
#' snplist<- as.data.frame(c("rs1000313","rs1001104"))
#' metalist<- as.data.frame(c("10-Formyltetrahydrofolate","Water"))
#' PathList<-SNP_Meta_query(MulOmicNet=MulOmicNet,snplist=snplist,metalist=metalist,
#'                          filter="Names",minpath=6,maxpath=7)
#'
SNP_Meta_query <- function(MulOmicNet,snplist,metalist,filter="",minpath=1,maxpath=3){
  if (filter[1]=="Names"){
    colnames(metalist)<-"Names"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Names %in% metalist$Names,1])
  }
  else if (filter[1]=="Formulas"){
    colnames(metalist)<-"Formulas"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Formulas %in% metalist$Formulas,1])
  }
  else if (filter[1]=="HMDB_ID"){
    colnames(metalist)<-"HMDB_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$HMDB_ID %in% metalist$HMDB_ID,1])
  }
  else if (filter[1]=="PubChem_ID"){
    colnames(metalist)<-"PubChem_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$PubChem_ID %in% metalist$PubChem_ID,1])
  }
  else if (filter[1]=="KEGG_ID"){
    colnames(metalist)<-"KEGG_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$KEGG_ID %in% metalist$KEGG_ID,1])
  }
  else if (filter[1]=="CHEBI_ID"){
    colnames(metalist)<-"CHEBI_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$CHEBI_ID %in% metalist$CHEBI_ID,1])
  }
  else if (filter[1]=="Metabolites"){
    colnames(metalist)<-"Metabolites"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Metabolites %in% metalist$Metabolites,1])
  }
  else{
    stop("Invalid 'filter' argument. Please refer to function ListFilters for valid input.")
    geterrmessage()
  }
  colnames(snplist)<-"SNP"
  #SNPAnno <- readRDS(file = "SNPAnno.rds")
  SNPNode <- as.character(unique(SNPAnno[SNPAnno$SNP %in% snplist$SNP,1]))
  if (length(MetaNode)==0|length(SNPNode)==0){
    print(sprintf("No snps/metabolites are found in MulOmicNet."))
    return()
  }
  else{
    print(sprintf("Identification of paths from %d snps to %d metabolites are about to start.",length(SNPNode),length(MetaNode)))
    PathList <-list()
    j<-0
    for (i in 1:length(SNPNode)){
      distance_m<-igraph::distances(MulOmicNet,v=SNPNode[i],to=MetaNode,mode="out")
      nodeto<-MetaNode[which(distance_m<=maxpath & distance_m>=minpath)]
      if(length(nodeto)==0){next}
      j<-j+1
      path_list <- igraph::all_shortest_paths(MulOmicNet,from=SNPNode[i],to=nodeto,mode ="out")
      PathList[[j]]<-path_list$res
      print(i)
    }
    if (length(PathList)==0){
      print(sprintf("End of searching. No paths are found."))
      return()
    }
    else{
      print(sprintf("End of searching. Writing paths with length between %d and %d into text file.",minpath,maxpath))
      for(i in (1:length(PathList))){
        temp1<-lapply(PathList[[i]], igraph::as_ids)
        temp2<-lapply(temp1,paste,collapse=",")
        lapply(temp2, write, "Paths.txt", append=TRUE)
      }
    }
  }
  return(PathList)
}

##########################################################
#Function7. Shortest Path Identification-SNP_Gene_query
##########################################################
#' @title Query shortest paths from SNPs to genes
#'
#' @description Taking SNPs and genes as input, one can get paths in our multi-omic network from SNPs to genes.
#'
#' @param MulOmicNet A dataframe of network with two columns.
#' @param snplist A dataframe of SNPs with one column.
#' @param genelist A dataframe of genes with one column.
#' @param filter Character constant,the type of input gene ID.Gene ID includes Gene_Entrez_ID, Gene_Name, Ensembl_Gene_ID, Enzyme_IDs and RefSeq_IDs.
#' For the detail,please see the function ListFilters.
#' @param minpath the lower limit of the shortest path length,the default value is 1.
#' @param maxpath the upper limit of the shortest path length,the default value is 3.
#'
#' @return A list is returned, each list element contains a shortest path from one SNP to one gene.
#' Also a text file is generated with the list written to.
#' @export SNP_Gene_query
#'
#' @examples
#' MulOmicNet<-BuildNet(threshold=0.5)
#' snplist<- as.data.frame(c("rs1000313","rs1001104"))
#' genelist<- as.data.frame(c("A1BG","A1BG-AS1","A2M","A2ML1"))
#' PathList<-SNP_Gene_query(MulOmicNet=MulOmicNet,snplist=snplist,genelist=genelist,
#'                          filter="Gene_Name",minpath=4,maxpath=5)
#'
SNP_Gene_query <- function(MulOmicNet,snplist,genelist,filter="",minpath=1,maxpath=3){
  GeneAnno_temp <- GeneAnno[GeneAnno$Gene_Entrez_ID %in% igraph::V(MulOmicNet)$name,]
  if (filter[1]=="RefSeq_IDs"){
    colnames(genelist)<-"RefSeq_IDs"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$RefSeq_IDs %in% genelist$RefSeq_IDs,1])
  }
  else if (filter[1]=="Gene_Name"){
    colnames(genelist)<-"Gene_Name"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Gene_Name %in% genelist$Gene_Name,1])
  }
  else if (filter[1]=="Ensembl_Gene_ID"){
    colnames(genelist)<-"Ensembl_Gene_ID"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Ensembl_Gene_ID %in% genelist$Ensembl_Gene_ID,1])
  }
  else if (filter[1]=="Enzyme_IDs"){
    colnames(genelist)<-"Enzyme_IDs"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Enzyme_IDs %in% genelist$Enzyme_IDs,1])
  }
  else if (filter[1]=="Gene_Entrez_ID"){
    colnames(genelist)<-"Gene_Entrez_ID"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Gene_Entrez_ID %in% genelist$Gene_Entrez_ID,1])
  }
  else{
    stop("Invalid 'filter' argument. Please refer to function ListFilters for valid input.")
    geterrmessage()
  }
  colnames(snplist)<-"SNP"
  #SNPAnno <- readRDS(file = "SNPAnno.rds")
  SNPNode <- as.character(unique(SNPAnno[SNPAnno$SNP %in% snplist$SNP,1]))
  if (length(GeneNode)==0|length(SNPNode)==0){
    print(sprintf("No snps/genes are found in MulOmicNet."))
    return()
  }
  else {
    print(sprintf("Identification of paths from %d snps to %d genes are about to start.",length(SNPNode),length(GeneNode)))
    PathList <-list()
    j<-0
    for (i in 1:length(SNPNode)){
      distance_m<-igraph::distances(MulOmicNet,v=SNPNode[i],to=GeneNode,mode="out")
      nodeto<-GeneNode[which(distance_m<=maxpath & distance_m>=minpath)]
      if(length(nodeto)==0){next}
      j<-j+1
      path_list <- igraph::all_shortest_paths(MulOmicNet,from=SNPNode[i],to=nodeto,mode ="out")
      PathList[[j]]<-path_list$res
      print(i)
    }
    if (length(PathList)==0){
      print(sprintf("End of searching. No paths are found."))
      return()
    }
    else{
      print(sprintf("End of searching. Writing paths with length between %d and %d into text file.",minpath,maxpath))
      for(i in (1:length(PathList))){
        temp1<-lapply(PathList[[i]], igraph::as_ids)
        temp2<-lapply(temp1,paste,collapse=",")
        lapply(temp2, write, "Paths.txt", append=TRUE)
      }
    }
  }
  return(PathList)
}

##########################################################
#Function8. Shortest Path Identification-Gene_Meta_query
##########################################################
#' @title Query shortest paths from genes to metabolites
#'
#' @description Taking genes and metabolites as input, one can get paths in our multi-omic network from genes to metabolites.
#'
#' @param MulOmicNet A dataframe of network with two columns.
#' @param genelist A dataframe of genes with one column.
#' @param metalist A dataframe of metabolites with one column.
#' @param genefilter Character constant,the type of input gene ID.Gene ID includes Gene_Entrez_ID, Gene_Name, Ensembl_Gene_ID, Enzyme_IDs and RefSeq_IDs.
#' For the detail,please see the function ListFilters.
#' @param metafilter Character constant, the type of input metabolite ID.Metabolite ID includes Metabolites, Names, Formulas, HMDB_ID, KEGG_ID, PubChem_ID and CHEBI_ID.
#' For the detail,please see the function ListFilters.
#' @param minpath The lower limit of the shortest path length,the default value is 1.
#' @param maxpath The upper limit of the shortest path length,the default value is 3.
#'
#' @return A list is returned, each list element contains a shortest path from one gene to one metabolite.
#' Also a text file is generated with the list written to.
#' @export Gene_Meta_query
#'
#' @examples
#' MulOmicNet<-BuildNet(threshold=0.5)
#' genelist<- as.data.frame(c("A1BG","AOC3"))
#' metalist<- as.data.frame(c("bamppald[c]","Water"))
#' PathList <- Gene_Meta_query(MulOmicNet=MulOmicNet,genelist=genelist,metalist=metalist,
#'                             genefilter="Gene_Name",metafilter="Names",minpath=1,maxpath=3)
#'
Gene_Meta_query <- function(MulOmicNet,genelist,metalist,genefilter="",metafilter="",minpath=1,maxpath=3){
  if (metafilter[1]=="Names"){
    colnames(metalist)<-"Names"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Names %in% metalist$Names,1])
  }
  else if (metafilter[1]=="Formulas"){
    colnames(metalist)<-"Formulas"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Formulas %in% metalist$Formulas,1])
  }
  else if (metafilter[1]=="HMDB_ID"){
    colnames(metalist)<-"HMDB_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$HMDB_ID %in% metalist$HMDB_ID,1])
  }
  else if (metafilter[1]=="PubChem_ID"){
    colnames(metalist)<-"PubChem_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$PubChem_ID %in% metalist$PubChem_ID,1])
  }
  else if (metafilter[1]=="KEGG_ID"){
    colnames(metalist)<-"KEGG_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$KEGG_ID %in% metalist$KEGG_ID,1])
  }
  else if (metafilter[1]=="CHEBI_ID"){
    colnames(metalist)<-"CHEBI_ID"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$CHEBI_ID %in% metalist$CHEBI_ID,1])
  }
  else if (metafilter[1]=="Metabolites"){
    colnames(metalist)<-"Metabolites"
    #MetaAnno <- readRDS(file="MetaAnno.rds")
    MetaNode <- as.character(MetaAnno[MetaAnno$Metabolites %in% metalist$Metabolites,1])
  }
  else{
    stop("Invalid 'metafilter' argument. Please refer to function ListFilters for valid input.")
    geterrmessage()
  }
  GeneAnno_temp <- GeneAnno[GeneAnno$Gene_Entrez_ID %in% igraph::V(MulOmicNet)$name,]
  if (genefilter[1]=="RefSeq_IDs"){
    colnames(genelist)<-"RefSeq_IDs"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$RefSeq_IDs %in% genelist$RefSeq_IDs,1])
  }
  else if (genefilter[1]=="Gene_Name"){
    colnames(genelist)<-"Gene_Name"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Gene_Name %in% genelist$Gene_Name,1])
  }
  else if (genefilter[1]=="Ensembl_Gene_ID"){
    colnames(genelist)<-"Ensembl_Gene_ID"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Ensembl_Gene_ID %in% genelist$Ensembl_Gene_ID,1])
  }
  else if (genefilter[1]=="Enzyme_IDs"){
    colnames(genelist)<-"Enzyme_IDs"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Enzyme_IDs %in% genelist$Enzyme_IDs,1])
  }
  else if (genefilter[1]=="Gene_Entrez_ID"){
    colnames(genelist)<-"Gene_Entrez_ID"
    #GeneAnno <- readRDS(file="GeneAnno.rds")
    GeneNode <- as.character(GeneAnno_temp[GeneAnno_temp$Gene_Entrez_ID %in% genelist$Gene_Entrez_ID,1])
  }
  else{
    stop("Invalid 'genefilter' argument. Please refer to function ListFilters for valid input.")
    geterrmessage()
  }
  if (length(GeneNode)==0|length(MetaNode)==0){
    print(sprintf("No genes/metabolites are found in MulOmicNet."))
    return()
  }
  else{
    print(sprintf("Identification of paths from %d genes to %d metabolites are about to start.",length(GeneNode),length(MetaNode)))
    PathList <-list()
    j<-0
    for (i in 1:length(GeneNode)){
      distance_m<-igraph::distances(MulOmicNet,v=GeneNode[i],to=MetaNode,mode="out")
      nodeto<-MetaNode[which(distance_m<=maxpath & distance_m>=minpath)]
      if(length(nodeto)==0){next}
      j<-j+1
      path_list <- igraph::all_shortest_paths(MulOmicNet,from=GeneNode[i],to=nodeto,mode ="out")
      PathList[[j]]<-path_list$res
      print(i)
    }
    if (length(PathList)==0){
      print(sprintf("End of searching. No paths are found."))
      return()
    }
    else{
      print(sprintf("End of searching. Writing paths with length between %d and %d into text file.",minpath,maxpath))
      for(i in (1:length(PathList))){
        temp1<-lapply(PathList[[i]], igraph::as_ids)
        temp2<-lapply(temp1,paste,collapse=",")
        lapply(temp2, write, "Paths.txt", append=TRUE)
      }
    }
  }
  return(PathList)
}

###########################################################
#Function9. Writing Subnetwork and Annotation File
###########################################################
#' @title Print output
#'
#' @description Print results to text files, plot subnetwork and perform gene enrichiment analysis.
#'
#' @param MulOmicNet A dataframe of network with two columns.
#' @param PathList A list of paths that generated by query function.
#' @param Anno Logical scalar,whether to annotate the nodes in the identified paths. the default is FALSE.
#' @param enrich Logical scalar,whether to perform the gene enrichment analysis, the default is FALSE.
#' @param simplot Logical scalar, whether to plot the subnetwork, the default is FALSE.
#'
#' @details A subnetwork composed of nodes involved in the identified paths is generated. The subnetwork
#' is written into a text file with each line representing one interaction.When argument Anno is set to true,
#' one can get annotation files of nodes involved in the identified paths.When argument enrich is set to true,
#' all genes involved in the identified paths are enriched in the 2019 KEGG human pathway database provided by Enrichr.
#' When argument simplot is set to true, an image of the subnetwork in PDF format is generated. The subnetwork is plotted
#' by placing vertices on the plane using the force-directed layout algorithm.
#'
#' @return text files and PDF file
#' @export SubNetRFW
#'
#' @examples
#' MulOmicNet<-BuildNet(threshold=0.5)
#' snplist<- as.data.frame(c("rs1000313","rs1001104"))
#' metalist<- as.data.frame(c("10-Formyltetrahydrofolate","Water"))
#' PathList<-SNP_Meta_query(MulOmicNet=MulOmicNet,snplist=snplist,metalist=metalist,
#'                          filter="Names",minpath=6,maxpath=7)
#' SubNetRFW(MulOmicNet=MulOmicNet,PathList=PathList,Anno=TRUE,enrich=TRUE,simplot=TRUE)
#'
SubNetRFW <- function(MulOmicNet,PathList,Anno=FALSE,enrich=FALSE,simplot=FALSE){
  if (length(PathList)==0){
    print(sprintf("PathList is empty."))
    return()
  }
  subnode_index<-unique(unlist(PathList))
  print(sprintf("Writing subnetwork file."))
  subg <- igraph::induced_subgraph(MulOmicNet,subnode_index)
  subgedges<-igraph::as_data_frame(subg,what="edges")
  utils::write.table(subgedges,file="SubNet.txt",sep="\t",quote=FALSE,row.names=FALSE)
  subnode_name<-igraph::V(MulOmicNet)$name[subnode_index]
  if(Anno){
    if(!requireNamespace("utils", quietly = TRUE)){
      stop("Package \"utils\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    print(sprintf("Writing annotation file."))
    #SNPAnno <-readRDS(file = "SNPAnno.rds")
    SubSNPAnno <-SNPAnno[SNPAnno$SNP %in% subnode_name,]
    utils::write.table(SubSNPAnno,file="SubSNPAnno.txt",sep="\t",quote=FALSE,row.names=FALSE)
    #GeneAnno <-readRDS(file = "GeneAnno.rds")
    SubGeneAnno <-GeneAnno[GeneAnno$Gene_Entrez_ID %in% subnode_name,]
    utils::write.table(SubGeneAnno,file="SubGeneAnno.txt",sep="\t",quote=FALSE,row.names=FALSE)
    #MetaAnno <-readRDS(file = "MetaAnno.rds")
    SubMetaAnno <-MetaAnno[MetaAnno$Metabolite_ID %in% subnode_name,]
    utils::write.table(SubMetaAnno,file="SubMetaAnno.txt",sep="\t",quote=FALSE,row.names=FALSE)
    #ReacAnno <- readRDS(file="ReacAnno.rds")
    SubReacAnno <-ReacAnno[ReacAnno$Reaction_ID %in% subnode_name,]
    utils::write.table(SubReacAnno,file="SubReacAnno.txt",sep="\t",quote=FALSE,row.names=FALSE)
  }
  if(enrich){
    if(!requireNamespace("enrichR", quietly = TRUE)){
      stop("Package \"enrichR\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    print(sprintf("Writing enrichment file."))
    #GeneAnno <-readRDS(file = "GeneAnno.rds")
    genes<-as.character(GeneAnno[GeneAnno$Gene_Entrez_ID %in% subnode_name,2])
    EnrichGene<-enrichR::enrichr(genes, databases = "KEGG_2019_Human")
    EnrichResult<-EnrichGene$KEGG_2019_Human[,c(1,2,4,7,9)]
    utils::write.table(EnrichResult,file="Enrichment.txt",sep="\t",quote=FALSE,row.names=FALSE)
  }
  if(simplot){
    if(!requireNamespace("grDevices", quietly = TRUE)){
      stop("Package \"grDevices\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if(!requireNamespace("graphics", quietly = TRUE)){
      stop("Package \"graphics\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    print(sprintf("Plotting subnetwork."))
    SubNode <-data.frame(Node=subnode_name)
    #SNPAnno <-readRDS(file = "SNPAnno.rds")
    SubNode[SubNode$Node %in% SNPAnno$SNP,2]<-1
    SubNode[SubNode$Node %in% SNPAnno$SNP,3]<-"SNP"
    #GeneAnno <-readRDS(file = "GeneAnno.rds")
    SubNode[SubNode$Node %in% GeneAnno$Gene_Entrez_ID,2]<-2
    SubNode[SubNode$Node %in% GeneAnno$Gene_Entrez_ID,3]<-"Gene"
    #ReacAnno <- readRDS(file="ReacAnno.rds")
    SubNode[SubNode$Node %in% ReacAnno$Reaction_ID,2]<-3
    SubNode[SubNode$Node %in% ReacAnno$Reaction_ID,3]<-"Reaction"
    #MetaAnno <-readRDS(file = "MetaAnno.rds")
    SubNode[SubNode$Node %in% MetaAnno$Metabolite_ID,2]<-4
    SubNode[SubNode$Node %in% MetaAnno$Metabolite_ID,3]<-"Metabolite"
    names(SubNode)<-c("subnode","node_attribute","node_type")
    SubNode <- SubNode[order(SubNode$node_attribute),]
    subnetwork <- igraph::graph_from_data_frame(d=subgedges, vertices=SubNode, directed=T)
    subnetwork <- igraph::simplify(subnetwork, remove.multiple = T, remove.loops = T,edge.attr.comb=c(weight="max", type="ignore"))
    l<-igraph::layout_with_fr(subnetwork)
    colrs <- c("LightGreen", "gold", "HotPink1","DeepSkyBlue")
    igraph::V(subnetwork)$color <- colrs[igraph::V(subnetwork)$node_attribute]
    grDevices::pdf(file="SubNet.pdf")
    graphics::plot(subnetwork,
         vertex.label=NA,
         vertex.color=igraph::V(subnetwork)$color,
         vertex.size=1.5,
         vertex.frame.color=NA,
         edge.width=0.01,
         egge.color="gray91",
         edge.arrow.size=0.1,
         edge.arrow.width=0.8,
         layout=l)
    graphics::legend(x=-1.5,y=-0.5,c("SNP","Gene", "Reaction","Metabolite"), pch=21,
           col=NA, pt.bg=colrs, pt.cex=1, cex=.6, bty="n", ncol=1)
    grDevices::dev.off()
  }
}

