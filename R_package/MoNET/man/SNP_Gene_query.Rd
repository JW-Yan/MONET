% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MoNET.R
\name{SNP_Gene_query}
\alias{SNP_Gene_query}
\title{Query shortest paths from SNPs to genes}
\usage{
SNP_Gene_query(
  MulOmicNet,
  snplist,
  genelist,
  filter = "",
  minpath = 1,
  maxpath = 3
)
}
\arguments{
\item{MulOmicNet}{A dataframe of network with two columns.}

\item{snplist}{A dataframe of SNPs with one column.}

\item{genelist}{A dataframe of genes with one column.}

\item{filter}{Character constant,the type of input gene ID.Gene ID includes Gene_Entrez_ID, Gene_Name, Ensembl_Gene_ID, Enzyme_IDs and RefSeq_IDs.
For the detail,please see the function ListFilters.}

\item{minpath}{the lower limit of the shortest path length,the default value is 1.}

\item{maxpath}{the upper limit of the shortest path length,the default value is 3.}
}
\value{
A list is returned, each list element contains a shortest path from one SNP to one gene.
Also a text file is generated with the list written to.
}
\description{
Taking SNPs and genes as input, one can get paths in our multi-omic network from SNPs to genes.
}
\examples{
MulOmicNet<-BuildNet(threshold=0.5)
snplist<- as.data.frame(c("rs1000313","rs1001104"))
genelist<- as.data.frame(c("A1BG","A1BG-AS1","A2M","A2ML1"))
PathList<-SNP_Gene_query(MulOmicNet=MulOmicNet,snplist=snplist,genelist=genelist,
                         filter="Gene_Name",minpath=4,maxpath=5)

}
