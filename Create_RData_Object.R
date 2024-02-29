#BiocManager::install("EnrichmentBrowser")
library(EnrichmentBrowser)
library(tidyverse)
library(limma)
library(readr)
library(dplyr)

## -------------- Get GO Terms -------------- ##
# get GO terms
# returns a list of GO terms with corresponding genes
go.gs<-getGenesets(org="hsa",db="go",onto="BP",mode="GO.db",cache=FALSE)

# list of all GO term IDs 
go.gene.ids<-unique(unlist(go.gs,use.names = FALSE))
# file to map GO gene identifiers to Uniprot ID using https://www.uniprot.org/id-mapping
write(go.gene.ids,"/Volumes/Stanton_Archive/Taub/DataObjectFiles/go_gene_ids.txt")

## -------------- Get KEGG Pathways -------------- ##
# getGeneKEGGLinks and getKEGGPathwayNames from limma package
# getGeneKEGGLinks: returns a data frame with with GeneID and PathwayID cols
kegg.gs<-getGeneKEGGLinks(species.KEGG = "hsa")
# getKEGGPathwayNames: returns a data frame with PathwayID and Description (pathway name) cols
kegg.pathnames<-getKEGGPathwayNames(species="hsa") 

# update data frame columns to have correct syntax for pathway IDs and gene IDs
kegg.gs <-kegg.gs %>% 
  mutate(PathwayID2=str_extract(PathwayID,"hsa\\d+")) %>% 
  mutate(GeneID2=paste("hsa",GeneID,sep=":")) %>% 
  dplyr::select(c("GeneID2","PathwayID2")) %>% 
  rename(PathwayID=PathwayID2,GeneID=GeneID2)

# merge pathway names, pathway IDs, and gene IDs 
kegg.pathnames.gs<-merge(kegg.gs,kegg.pathnames,by.x="PathwayID",by.y="PathwayID")
# clean Description col
kegg.pathnames.gs$Description<-sapply(kegg.pathnames.gs$Description,function(x){str_split_1(x," - Homo")[1]})

# drop encoded pathway identifier
kegg.pathnames.gs<-kegg.pathnames.gs %>% 
  dplyr::select(c("GeneID","Description")) %>% 
  rename(PathwayID=Description)

kegg.gene.ids<-unique(kegg.pathnames.gs$GeneID)
# file to map KEGG gene identifiers to Uniprot ID using https://www.uniprot.org/id-mapping
write(kegg.gene.ids,"/Volumes/Stanton_Archive/Taub/DataObjectFiles/kegg_gene_ids.txt")

## -------------- Get Uniprot Mappings -------------- ##
# GO gene IDs may need to be fed into Uniprot mapper as two separate files (length constraints)
# The downloaded results of the Uniprot mapper should be specified in the file paths below
go.filename="/Volumes/Stanton_Archive/Taub/DataObjectFiles/go_uniprot_04112023.tsv"
go.uniprot.df<-read_tsv(go.filename)

kegg.filename="/Volumes/Stanton_Archive/Taub/DataObjectFiles/kegg_uniprot_04112023.tsv"
kegg.uniprot.df<-read_tsv(kegg.filename)

## -------------- Create GO Data Structure -------------- ##
# function used to return cleaned GO term name
# utilized when creating GO term data structure
getGOTerm<-function(entry){
  entry.list<-strsplit(entry,"_")
  go.term<-entry.list[[1]][2:length(entry.list[[1]])]
  return(paste(go.term,collapse=" "))
}

# convert list to data frame in order to take advantage of merge
go.gs.df<-data.frame()
for(i in 1:length(go.gs)){
  go.terms<-rep(getGOTerm(names(go.gs[i])),length(go.gs[[i]]))
  entrez.ids<-go.gs[[i]]
  go.gs.df<-rbind(go.gs.df,data.frame("go.term"=go.terms,"From"=entrez.ids))
}

# merge GO terms, GO gene IDs and Uniprot gene IDs
go.gs.df<-go.gs.df %>% 
  mutate(From=as.numeric(From)) %>% 
  left_join(go.uniprot.df,multiple = "all")

# join GO data frame with Uniprot IDs
go.gs.df<-go.gs.df %>% 
  left_join(go.uniprot.df,multiple = "all")

# remove NAs
go.df<-go.gs.df %>% 
  drop_na(Entry) %>% 
  dplyr::select(-c("Entry Name"))

## -------------- Create KEGG Data Structure -------------- ##
# merge KEGG pathway IDs and gene IDs with Uniprot gene IDs
kegg.df<-merge(kegg.pathnames.gs,kegg.uniprot.df,by.x="GeneID",by.y="From") %>% 
  rename(From=GeneID, pathway=PathwayID) %>% 
  dplyr::select(-c("Entry Name"))

## -------------- Create Path List Data Structure -------------- ##
# this data structure is used by the app to link KEGG pathways to pathway diagrams from https://www.kegg.jp/kegg/kegg2.html
path.list<-data.frame(ID=kegg.pathnames$PathwayID,Name=kegg.pathnames$Description)
# ID and name cleaning
path.list$ID<-substr(path.list$ID,4,8)
path.list$Name<-sapply(kegg.pathnames$Description,function(x){str_split_1(x," - Homo")[1]},USE.NAMES = FALSE)

save(kegg.df,go.df,path.list,file="/Volumes/Stanton_Archive/Taub/AppData2.RData")
