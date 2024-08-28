# WikiPathway-Networking-
#My first repository on GitHub
#1 retrieving associated gene from variants.

#1.1 Set up environment
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
install.packages("tidyverse")
library(BiocManager)
library(biomaRt)
library(tidyverse)

#1.2 Connect to the BioMart database and dataset within the database. Here we chose human gene to investigate as an example.
mart <- useMart("ENSEMBL_MART_SNP") #Choose "ENSEMBL_MART_SNP" database.
data <- listDatasets(mart)
data_human <- useDataset("mmusculus_snp", mart = mart) # Choose human dataset.
listAttributes(data_human)

#1.3 Load the data of interest 
Raw_data <- read.delim("variant_list.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(Raw_data)[1] <- "variant" # We converted the txt file into a table with a column called "variant" 
                                   #and rows with the names of the variants.
listFilters(data_human)

#1.4 Retrieve the gene names from the dataset by filtering the overlapping variants from raw data and remove duplicates.
SNP_gene <- getBM(attributes = c('refsnp_id','associated_gene'),
                 filters = "snp_filter",
                  values = Raw_data$variant,
                 mart = data_human)

#1.6 Expand the columns with multiple associated genes
SNP_gene_expanded <- SNP_gene %>%
  separate_rows(associated_gene, sep = ",")

#1.7 Remove rows and columns with any NA or empty values.
SNP_gene_expanded_cleaned <- SNP_gene_expanded[!(is.na(SNP_gene_expanded$refsnp_id) | SNP_gene_expanded$refsnp_id == "" |
                                is.na(SNP_gene_expanded$associated_gene) | SNP_gene_expanded$associated_gene == ""), ]
SNP_gene_expanded_cleaned <- unique(SNP_gene_expanded_cleaned)# Remove the duplicates.
write.table(SNP_gene_expanded_cleaned$associated_gene, "Associated Gene_names", col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
print(SNP_gene_expanded_cleaned)

#2 Create a SNP-Gene network in Cytoscape. 

#2.1 Set the environment to connet RStudio to Cytoscape.
BiocManager::install("RCy3", force = TRUE)
library(RCy3)
cytoscapePing() 
if("cytargetlinker" %in% commandsHelp("")) print("Success: the CyTargetLinker app is installed") 
else print("Warning: CyTargetLinker app is not installed. Please install the CyTargetLinker app before proceeding.")
cytoscapeVersionInfo() # Checking the connection.

#2.2 Create the SNP-Gene network.
edges <- data.frame(
  source = SNP_gene_expanded_cleaned$refsnp_id, target = SNP_gene_expanded_cleaned$associated_gene, 
                    interaction = "interacts_with",
                    stringsAsFactors = FALSE)# Define the edges.
nodes <- data.frame(
  id = c(SNP_gene_expanded_cleaned$refsnp_id, SNP_gene_expanded_cleaned$associated_gene),
  stringsAsFactors = FALSE
)# Define the nodes.
createNetworkFromDataFrames(edges = edges, nodes = nodes, title = "SNP-Gene Network_mouse", collection = "SNP-Gene-pathway Analysis_mouse")

#3 Add known pathways from WikiPathways.

#3.1 Set the environment.
BiocManager::install("rWikiPathways",force = TRUE)
library(rWikiPathways)
installApp('WikiPathways')
installApp('CyTargetLinker')
installApp('stringApp')

#3.2 Extend network with pathway information
wp <- file.path(getwd(), "wikipathways_hsa_20240410.xgmml")
CTLextend.cmd = paste('cytargetlinker extend idAttribute="shared name" linkSetFiles="', wp,
                      '" network=current direction=SOURCES', sep="")
commandsRun(CTLextend.cmd)
layoutNetwork()
fitContent()
filter1.cmd = "network select edgeList=all"
filter2.cmd = "network select extendEdges=true"
filter3.cmd = "network create nodeList=selected edgeList=selected networkName=selection source=current"

commandsRun(filter1.cmd)
commandsRun(filter2.cmd)
commandsRun(filter3.cmd)

#set visual style manually

