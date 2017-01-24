# TODO: Add comment
# 
# Author: Sudeep Sahadevan
# Perform GO enrichment using topGO
# Detailed description for most of the functions are in the accompanying .Rmd and .html files
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

###############################################################################
#
# for the given input file (Excel .xls) generate a list of DEGS for each comparison
# input_param: 
#           input_path : string for file path
#                sheet :  sheet number with DEG tables in .xls/.xlsx file
# sample input .xls format:
#   BFP  	# comparision type
#       id	miRNA_subepi.	miRNA_epi	Gap	REV_epi	REV_subepi
#   30	AT1G01900	0.34005745	0.446529969	1	0.333646081	0.722220337
#   37	AT1G02230	0.105150325	0.144800866	1	1.145252733	0.190918928
#   44	AT1G02450	0.397239147	0.623130411	1	0.753432416	0.33917174
#   BY  	# comparision type				
#       id	miRNA_subepi.	miRNA_epi	Gap	REV_epi	REV_subepi
#   26	AT1G01690	0.768192857	0.876352222	1	1.819999363	1.326299012
#   34	AT1G02140	0.665887973	0.783444425	1	1.517932165	0.879594451
#   42	AT1G02370	0.565892052	0.64731211	1	1.161606912	0.552932433
#
###############################################################################
get_gene_lists <- function(input_path,sheet=1){
  library(gdata)
  # return something
  gene_list <- list()
  excel_degs <- read.xls(xls=input_path,sheet=sheet,header=FALSE,stringsAsFactors=FALSE)
  list_name <- ""
  gene_ids <- character()
  first <- TRUE
  for (r in 1:nrow(excel_degs)){
    if(excel_degs[r,2]==""){
      if(!first){
        gene_list[[list_name]] <- gene_ids
        gene_ids <- character()
      }
      list_name <- as.character(excel_degs[r,1])
      first <- FALSE
    }else{
      gene_ids <- c(gene_ids, as.character(excel_degs[r,2]))
    }
  }
  return (gene_list)
}

###############################################################################
#
# convert Gene Ontology table into a list compatabile with topGO
#	table format : as required by GOstats R package
# input_params
#                  GO_table : table format ("GO_genes","GO_ids") 
#                           : rest of the table columns are not considered even if present
#                 minThresh : minimum number of Gene annotations for a GO term
#                           : terms with less than minThresh genes will be pruned off
# 
###############################################################################
get_topGO_list <- function(GO_table,minThresh=5){
  if(dim(GO_table)[1]==0){stop("Error! empty GO_table\n")}
  go_ids <- unique(GO_table[,2])
  goLen <- c(1:length(go_ids))
  go_list <- sapply(goLen,function(goLen){
    goSub <- GO_table[GO_table[,2] %in% go_ids[goLen],]
    geneIds <- unique(goSub[,1])
    tempList <- geneIds
    tempList
  },simplify=FALSE)
  names(go_list) <- go_ids
  selector <- sapply(goLen,function(goLen) length(go_list[[goLen]])>=minThresh)
  tempNames <- names(go_list)
  go_list <- go_list[selector]
  names(go_list) <- tempNames[selector]
  return(go_list)
}

###############################################################################
#
# get Arabidopsis GO ontology table and convert it into topGO compatible list format
# input_params 
#             ontology : GO ontolgy to use (BP or MF or CC)
#
###############################################################################
get_ath_topGO_list <- function(ontology="BP"){
  library(org.At.tair.db)
  ontology <- toupper(ontology)
  ontology <- match.arg(ontology,c("BP","MF","CC"),FALSE)
  ath_annotations <- toTable(org.At.tairGO)
  ath_annotations <- ath_annotations[ ath_annotations$Ontology==ontology, ]
  ath_topGO_list <- get_topGO_list(GO_table=ath_annotations,minThresh=1)
  return(ath_topGO_list)
}

###############################################################################
#
# do topGO enrichment test for the DEG list from get_gene_list using GO list from
# function get_ath_topGO_list
# input_params
#             gene_list : DEG gene list, using the same list format as output from 
#                       : function get_gene_lists
#           topGO_list  : GO annotation list, see topGO vignette and manual for details,
#                       : output from function get_ath_topGO_list
#              ontology : MUST be one of (BP,MF,CC) (default: BP) and 
#                       : MUST be the same as the ontology type for topGO_list
#       topGO_algorithm : topGO algorithm to use, see topGO manual (default: weight01)
#       topGO_statistic : topGO statistic test to use, see topGO manual (default: fisher)
#           pval_cutoff : pvalue cut off for enriched GO terms (default: 0.1) 
#             node_size : before tests, prune GO hierarchy for terms with less than node_size genes annotated (default: 3)
#                 ncpus : whether use morethan one core at a time to run multiple tests, works only on a multicore cpu and 
#                       : works only on a multicore cpu (Duh!) and require libraries doMC and foreach to be installed           
#
###############################################################################
do_topGO_enrichment <- function(gene_lists,topGO_list,ontology="BP",topGO_algorithm="weight01",topGO_statistic="fisher",
                                pval_cutoff=0.1,node_size=3,enrichment=3,ncpus=1){
  if(ncpus>1){
    library(doMC)
    library(foreach)
    registerDoMC(cores=ncpus)
  }
  ontology <- toupper(ontology)
  ontology <- match.arg(ontology,c("BP","MF","CC"),FALSE)
  universeGenes <- as.character(unique(unlist(topGO_list)))
  tot_list <- length(gene_lists)
  message("\n\tfound ",tot_list," DEG list input gene_list\n")
  if(ncpus==1){
    library(topGO)
#     return something
    out_list <- list()
    for( g in 1:length(gene_lists)){
      list_name = names(gene_lists)[g]
      genes <- as.character(gene_lists[[g]])
      geneSet <- factor(as.integer(universeGenes %in% genes))
      names(geneSet) <- universeGenes
      topGO_session <- new("topGOdata",description="topGO enrichment analysis for DEGs",nodeSize=node_size,annot=annFUN.GO2genes,
                           ontology=ontology,allGenes=geneSet,GO2genes=topGO_list)
      topGO_test <- runTest(topGO_session,algorithm=topGO_algorithm,statistic=topGO_statistic)
      topNode_count <- sum(score(topGO_test)<=pval_cutoff)
      resultTable <- GenTable(topGO_session,testScore=topGO_test,topNodes=topNode_count)
      resultTable <- resultTable[resultTable$Significant>=enrichment,]
	  if(nrow(resultTable)>0){
		  out_list[[list_name]] <- list(enrichment=resultTable,testClass=list_name)  
	  }
      message("\t ...........finished ",g,"/",tot_list,"........\n")
	  flush.console()
    }
	message("\n")
    return(out_list)
  }else{
    out_list <- foreach(g=1:length(gene_lists))%dopar%{
      library(topGO)
      list_name = names(gene_lists)[g]
      genes <- as.character(gene_lists[[g]])
      geneSet <- factor(as.integer(universeGenes %in% genes))
      names(geneSet) <- universeGenes
      topGO_session <- new("topGOdata",description="topGO enrichment analysis for DEGs",nodeSize=node_size,annot=annFUN.GO2genes,
                           ontology=ontology,allGenes=geneSet,GO2genes=topGO_list)
      topGO_test <- runTest(topGO_session,algorithm=topGO_algorithm,statistic=topGO_statistic)
      topNode_count <- sum(score(topGO_test)<=pval_cutoff)
      resultTable <- GenTable(topGO_session,testScore=topGO_test,topNodes=topNode_count)
      resultTable <- resultTable[resultTable$Significant>=enrichment,]
      if(dim(resultTable)[1]>0){
        temp_list <- list(enrichment=resultTable,testClass=list_name)
        temp_list
      }
    }
    list_len <- c(1:length(out_list))
    names(out_list) <- sapply(list_len,function(list_len){out_list[[list_len]]$testClass})
    return(out_list)
  }
  message("\n\tFin.\n")
}

#' do topGO enrichemnt for a given table list
#' the elements of the list must have an attribute "degs", which should contain the list of DEGS
#' @param deg list : DEG table list 
#' @param topGO_list  : GO annotation list, see topGO vignette and manual for details,output from function get_ath_topGO_list
#' @param universeGenes : list of genes to be used as universal list of genes
#' @param ontology : MUST be one of (BP,MF,CC) (default: BP) and MUST be the same as the ontology type for topGO_list
#' @param topGO_algorithm : topGO algorithm to use, see topGO manual (default: weight01)
#' @param topGO_statistic : topGO statistic test to use, see topGO manual (default: fisher)
#' @param pval_cutoff : pvalue cut off for enriched GO terms (default: 0.1) 
#' @param node_size : before tests, prune GO hierarchy for terms with less than node_size genes annotated (default: 3)
#' @param enrichment : minimum cut-off threshold for the number of genes annotated to a  GO term to be selected as enriched 
#' @param ncpus : whether use morethan one core at a time to run multiple tests, works only on a multicore cpu and works only on a multicore cpu (Duh!) and require libraries doMC and foreach to be installed
#' @return a list
do_topGO_enrichment_degs <- function(deg_lists,topGO_list,universeGenes=NULL,ontology="BP",topGO_algorithm="weight01",topGO_statistic="fisher",
		pval_cutoff=0.1,node_size=3,enrichment=3,ncpus=1){
	
	if(ncpus>1){
		library(doMC)
		library(foreach)
		registerDoMC(cores=ncpus)
	}
	ontology <- toupper(ontology)
	ontology <- match.arg(ontology,c("BP","MF","CC"),FALSE)
	if(is.null(universeGenes)){
		universeGenes <- as.character(unique(unlist(topGO_list)))
	}
	tot_list <- length(deg_lists)
	message("\n\tfound ",tot_list," DEG list input gene_list\n")
	if(ncpus==1){
		library(topGO)
#     return something
		out_list <- list()
		for( g in 1:length(deg_lists)){
			list_name = names(deg_lists)[g]
			genes <- as.character(rownames(deg_lists[[g]]$degs))
			geneSet <- factor(as.integer(universeGenes %in% genes))
			names(geneSet) <- universeGenes
			topGO_session <- new("topGOdata",description="topGO enrichment analysis for DEGs",nodeSize=node_size,annot=annFUN.GO2genes,
					ontology=ontology,allGenes=geneSet,GO2genes=topGO_list)
			topGO_test <- runTest(topGO_session,algorithm=topGO_algorithm,statistic=topGO_statistic)
			topNode_count <- sum(score(topGO_test)<=pval_cutoff)
			resultTable <- GenTable(topGO_session,testScore=topGO_test,topNodes=topNode_count)
			resultTable <- resultTable[resultTable$Significant>=enrichment,]
			out_list[[list_name]] <- list(enrichment=resultTable,testClass=list_name)
			message("\t ...........finished ",g,"/",tot_list,"........\n")
		}
		return(out_list)
	}else{
		out_list <- foreach(g=1:length(deg_lists))%dopar%{
			library(topGO)
			list_name = names(deg_lists)[g]
			genes <- as.character(rownames(deg_lists[[g]]$degs))
			geneSet <- factor(as.integer(universeGenes %in% genes))
			names(geneSet) <- universeGenes
			topGO_session <- new("topGOdata",description="topGO enrichment analysis for DEGs",nodeSize=node_size,annot=annFUN.GO2genes,
					ontology=ontology,allGenes=geneSet,GO2genes=topGO_list)
			topGO_test <- runTest(topGO_session,algorithm=topGO_algorithm,statistic=topGO_statistic)
			topNode_count <- sum(score(topGO_test)<=pval_cutoff)
			resultTable <- GenTable(topGO_session,testScore=topGO_test,topNodes=topNode_count)
			resultTable <- resultTable[resultTable$Significant>=enrichment,]
			if(dim(resultTable)[1]>0){
				temp_list <- list(enrichment=resultTable,testClass=list_name)
				temp_list
			}
		}
		list_len <- c(1:length(out_list))
		names(out_list) <- sapply(list_len,function(list_len){out_list[[list_len]]$testClass})
		return(out_list)
	}
	message("\n\tFin.\n")
}

###############################################################################
#
# generate DAG igraph object of a given GO ontology
#	see here for the original source http://igraph.wikidot.com/r-recipes
# input_params:
#              ontology : MUST be BP or MF or CC for one of the ontologies 
#                       : and ALL for all the ontologies combined            
#
###############################################################################
get_GO_graph <- function(ontology="BP"){
  library(igraph)
  library(GO.db) # assuming that the latest version is installed
  ontologyFrame <- data.frame()
  if(ontology=="BP"||ontology=="bp"){
    ontologyFrame <- toTable(GOBPPARENTS)
  }else if(ontology=="CC"||ontology=="cc"){
    ontologyFrame <- toTable(GOCCPARENTS)
  }else if(ontology=="MF"||ontology=="mf"){
    ontologyFrame <- toTable(GOMFPARENTS)
  }else if(ontology=="All"||ontology=="ALL"||ontology=="all"){
    BP <- toTable(GOBPPARENTS)
    CC <- toTable(GOCCPARENTS)
    MF <- toTable(GOMFPARENTS)
    ontologyFrame <-  rbind(BP,CC,MF)
  }
  go_ontology <- graph.data.frame(ontologyFrame)
  terms <- toTable(GOTERM)[,2:3]
  terms <- terms[ !duplicated(terms[,1]), ]
  rownames(terms) <- terms[,1]
  terms <- terms[V(go_ontology)$name,]
  V(go_ontology)$goName <- as.character(gsub(pattern="\'",replacement="",x=terms[,2]))
  return(go_ontology)
}

#' given a directed Gene Ontology graph and a GO id, find all the child nodes of the given id 
#' @param go_graph: GO directed graph, igraph object
#' @param Vids: GO identifier
#' @param rel_type: GO relationship type, MUST BE one of "is_a","regulates","negatively_regulates","part_of","positively_regulates"
get_child_nodes <- function(go_graph,Vids,rel_type=c("is_a","regulates","negatively_regulates","part_of","positively_regulates"),rel_attrib="RelationshipType"){
	library(igraph)
#	function for returning child nodes of a single graph
	single_node_fun <- function(go_graph,vids,rel_type,rel_attrib){
		vid <- which( V(go_graph)$name == vids )
		child_terms <- vid
		all_child <- child_terms
		while(length(child_terms)>0){
			child_terms <- neighbors(go_graph,v=child_terms,mode="in")
			all_child <- union(all_child,child_terms)
		}
		ind_subg <- induced.subgraph(go_graph,vids=all_child)
		subg_el <- as.data.frame(get.edgelist(ind_subg,names=TRUE))
		subg_el <- data.frame(subg_el,rel=get.edge.attribute(ind_subg,name=rel_attrib))
		subg_el <- subg_el[ subg_el$rel %in% rel_type, ]
		return(union(as.character(subg_el$V1),as.character(subg_el$V2)))
	}
#	attributes, sanity checking....
	rel_type <- match.arg(rel_type,c("is_a","regulates","negatively_regulates","part_of","positively_regulates"),several.ok=TRUE)
	Vids <- intersect(toupper(Vids),V(go_graph)$name)
	if(length(Vids)==0){ stop("ERROR! the given Vids are not found in the input graph\n") }
#	do some work
	if(length(Vids)==1){
		return(single_node_fun(go_graph,Vids,rel_type,rel_attrib))
	}else{
		ret_list <- list()
		for(i in 1:length(Vids)){
			ret_list[[Vids[i]]] <- single_node_fun(go_graph,Vids[i],rel_type,rel_attrib)
		}
		return(ret_list)
	}
}

###############################################################################
#
# given a go_graph (igraph object), a GO term to gene list and a GO id return 
# all the gene annotations for the given GO term and its immediate child nodes
# input_params:
#               go_graph : GO as directed graph (igraph object, output from function get_GO_graph)
#           go_gene_list : GO term to gene annotation list (output from functions get_topGO_list or get_ath_topGO_list)
#                  go_id : GO term id (GO:\d+)
#           neighborhood : Value to return neighbors of a vertex. See igraph graph.neighborhoof function 
#
###############################################################################
get_GO_genes <- function(go_graph,go_gene_list,go_id,neighborhood=1){
  library(igraph)
  if(!is.igraph(go_graph)){stop("\tError! go_graph must be an igraph object\n")}
  if(!is.directed(go_graph)){stop("\tError! go_graph must be a directed DAG\n")}
  child_nodes <- V(graph.neighborhood(graph=go_graph,order=neighborhood,nodes=as.character(go_id),mode="in")[[1]])$name
  if(length(child_nodes)<=1){message("\tWarning! cannot find child nodes for ",go_id," returning only its annotations\n")}
  common_nodes <- intersect(child_nodes,names(go_gene_list))
  if(length(common_nodes)==0){
	  message("\tError! there are no gene annotations in go_gene_annotations for child nodes of ",go_id," removing annotations\n")
	  return(NULL)
  }else{
	  clen <- c(1:length(common_nodes))
	  all_gene_annotations <- lapply(clen,function(clen) go_gene_list[[common_nodes[clen]]])
	  all_gene_annotations <- unique(unlist(all_gene_annotations))
	  return(all_gene_annotations)  
  }
  
}

###############################################################################
#
# add genes to enriched GO terms from topGO test
# due to the nature of enrichment test done in topGO, the default test results 
# does not contain (DE) genes for an enriched GO term.
# input_params
#             go_graph : GO as an igraph object (output from function get_GO_graph)
#  go_gene_annotations : GO term to gene annotation list (output from function get_topGO_list or get_ath_topGO_list)
#     topGO_enrichment : topGO enrichment list (output from function do_topGO_enrichment)
#            gene_list : Gene list used in topGO enrichment (output from function get_gene_lists)
#
###############################################################################
get_genes_topGO_enrichment <- function(go_graph,go_gene_list,topGO_enrichment,gene_list){
  if(is.null(names(topGO_enrichment))){stop("Error! die list topGO_enrichment must have cluster names as names(list)\n")}
  if(is.null(names(gene_list))){stop("Error! die list gene_list must have cluster names as names(list)\n")}
  common_classes <- intersect(names(topGO_enrichment),names(gene_list))
  if(length(common_classes)==0){stop("Error! topGO_enrichment and gene_list no common cluster ids\n")}
  topGO_enrichment <- topGO_enrichment[common_classes]
  for(i in 1:length(topGO_enrichment)){
    list_name <- names(topGO_enrichment[i])
    class_geneIds <- as.character(gene_list[[list_name]])
    enrichmentFrame <- topGO_enrichment[[i]]$enrichment
    enrichedGenes_list <- list()
    for(g in 1:nrow(enrichmentFrame)){
      go_id <- enrichmentFrame[g,1]
	  message(go_id,"\n")
      sig_nr <- enrichmentFrame[g,4]
      if(!class(sig_nr)=="integer"){
        sig_nr <- as.integer(as.character(sig_nr))
      }
      enrichedGenes_nr <- 0
      neighborhood <- 1
#			until number of enriched genes equals "Significant" column in enrichmentFrame
     # message(names(topGO_enrichment[i])," : ",enrichedGenes_nr," : ",sig_nr,"\n")
      while(enrichedGenes_nr<sig_nr){
#				should implement a min. number of trials counter, but let's see fors now
        enrichedGenes <- get_GO_genes(go_graph = go_graph,go_gene_list = go_gene_list,
                                           go_id=go_id,neighborhood = neighborhood)
        enrichedGenes <- intersect(enrichedGenes,class_geneIds)
        enrichedGenes_nr <- length(enrichedGenes)
        enrichedGenes_list[[g]] <- as.character(enrichedGenes)
        neighborhood <- neighborhood +1
      }
    }
    eLen <- c(1:length(enrichedGenes_list))
    enrichedGenes_vec <- sapply(eLen,function(eLen) paste(enrichedGenes_list[[eLen]],collapse=", "))
    go_term_names <- get.vertex.attribute(graph=go_graph,name="goName",index=as.character(enrichmentFrame[,1]))
    enrichmentFrame[,2] <- go_term_names
    enrichmentFrame <- data.frame(enrichmentFrame,enrichedGenes=as.character(enrichedGenes_vec))
    topGO_enrichment[[i]]$enrichment <- enrichmentFrame
  }
  return(topGO_enrichment)
}

#' wrapper function to do topGO enrichment on list, and find gene annotation,
#' calls functions do_topGO_enrichment  and get_genes_topGO_enrichment
#' for details of input parameters, see the respective funcitons
topGO_enrichment_wrapper <- function(gene_lists,topGO_list,ontology="BP",topGO_algorithm="weight01",topGO_statistic="fisher",
		pval_cutoff=0.1,node_size=3,enrichment=3,ncpus=1,go_graph){
	topGO_enrichment <- do_topGO_enrichment(gene_lists,topGO_list,ontology=ontology,topGO_algorithm=topGO_algorithm,topGO_statistic=topGO_statistic,
			pval_cutoff=pval_cutoff,node_size=node_size,enrichment=enrichment,ncpus=ncpus)
	topGO_annotation <- get_genes_topGO_enrichment(go_graph,topGO_list,topGO_enrichment,gene_lists)
	return(topGO_annotation)
}

#' wrapper function to do topGO enrichment on list, and find gene annotation,
#' calls functions do_topGO_enrichment  and get_genes_topGO_enrichment
#' for details of input parameters, see the respective funcitons
topGO_enrichment_wrapper_degs <- function(gene_lists,topGO_list,ontology="BP",topGO_algorithm="weight01",topGO_statistic="fisher",
		pval_cutoff=0.1,node_size=3,enrichment=3,ncpus=1,go_graph){
	topGO_enrichment_degs <- do_topGO_enrichment_degs(gene_lists,topGO_list,ontology=ontology,topGO_algorithm=topGO_algorithm,topGO_statistic=topGO_statistic,
			pval_cutoff=pval_cutoff,node_size=node_size,enrichment=enrichment,ncpus=ncpus)
	topGO_annotation <- get_genes_topGO_degs(go_graph,topGO_list,topGO_enrichment_degs,gene_lists)
	return(topGO_annotation)
}

###############################################################################
#
# add genes to enriched GO terms from topGO test
# due to the nature of enrichment test done in topGO, the default test results 
# does not contain (DE) genes for an enriched GO term.
# input_params
#             go_graph : GO as an igraph object (output from function get_GO_graph)
#  go_gene_annotations : GO term to gene annotation list (output from function get_topGO_list or get_ath_topGO_list)
#     topGO_enrichment : topGO enrichment list (output from function do_topGO_enrichment)
#            gene_list : Gene list used in topGO enrichment (output from function get_gene_lists)
#
###############################################################################
get_genes_topGO_degs <- function(go_graph,go_gene_list,topGO_enrichment,deg_list){
	if(is.null(names(topGO_enrichment))){stop("Error! die list topGO_enrichment must have cluster names as names(list)\n")}
	if(is.null(names(deg_list))){stop("Error! die list gene_list must have cluster names as names(list)\n")}
	common_classes <- intersect(names(topGO_enrichment),names(deg_list))
	if(length(common_classes)==0){stop("Error! topGO_enrichment and gene_list no common cluster ids\n")}
	topGO_enrichment <- topGO_enrichment[common_classes]
	for(i in 1:length(topGO_enrichment)){
		list_name <- names(topGO_enrichment[i])
		class_geneIds <- as.character(rownames(deg_list[[list_name]]$degs))
		enrichmentFrame <- topGO_enrichment[[i]]$enrichment
		enrichedGenes_list <- list()
		remove_rows <- vector()
		rownames(enrichmentFrame) <- NULL
		list_id <- 1
		for(g in 1:nrow(enrichmentFrame)){
			go_id <- enrichmentFrame[g,1]
			sig_nr <- enrichmentFrame[g,4]
			if(!class(sig_nr)=="integer"){
				sig_nr <- as.integer(as.character(sig_nr))
			}
			enrichedGenes_nr <- 0
			neighborhood <- 1
#			until number of enriched genes equals "Significant" column in enrichmentFrame
			while(enrichedGenes_nr<sig_nr){
#				should implement a min. number of trials counter, but let's see fors now
				enrichedGenes <- get_GO_genes(go_graph = go_graph,go_gene_list = go_gene_list,
						go_id=go_id,neighborhood = neighborhood)
				if(is.null(enrichedGenes)){
					remove_rows <- c(remove_rows,g)
					break
				}else{
					enrichedGenes <- intersect(enrichedGenes,class_geneIds)
					enrichedGenes_nr <- length(enrichedGenes)
					enrichedGenes_list[[g]] <- as.character(enrichedGenes)
					neighborhood <- neighborhood +1
				}
			}
		}
		row_it <- c(1:nrow(enrichmentFrame))
		good_ids <- sapply(row_it,function(row_it) !any(row_it %in% remove_rows))
		enrichmentFrame <- enrichmentFrame[good_ids,]
		rownames(enrichmentFrame) <- NULL
#		now prune enrichedGenes_list
		eLen <- c(1:length(enrichedGenes_list))
		gene_selector <- sapply(eLen, function(eLen) length(enrichedGenes_list[[eLen]])>0)
		enrichedGenes_list <- enrichedGenes_list[gene_selector]
		eLen <- c(1:length(enrichedGenes_list))
		enrichedGenes_vec <- sapply(eLen,function(eLen) paste(enrichedGenes_list[[eLen]],collapse=", "))
		go_term_names <- get.vertex.attribute(graph=go_graph,name="goName",index=as.character(enrichmentFrame[,1]))
		enrichmentFrame[,2] <- go_term_names
		enrichmentFrame <- data.frame(enrichmentFrame,enrichedGenes=as.character(enrichedGenes_vec))
		topGO_enrichment[[i]]$enrichment = enrichmentFrame
	}
	return(topGO_enrichment)
}

###############################################################################
#
# function to generate Xgmml objects for enriched GO terms bapsed on input
# GO enriched terms, corresponding p-values and a GO graph igraph object
# input_params:
#              go_graph : GO graph (igraph object)
#              ontology : MUST be one of (BP,CC,MF)
#       enrichedGOTable : GO enrichment table from topGO
#               goIdCol : GO term id column (default : 1)
#               pvalCol : P-value column (default : 6)
#            beginColor : Color for the most enriched terms (default : green)
#              endColor : Color for the least enriched terms (default : yellow)
#            vertexName : go_graph vertex attribute name for GO term names (default for go_graph : goName)
#              getXgmml : Boolean, return Xgmml or igraph object (default : TRUE)
#             graphName : A string with specific name for the graph (default : enriched terms)
#
###############################################################################
generate_GO_graph <-function(go_graph,ontology="BP",enrichedGOTable,goIdCol=1,
		pvalCol=6,beginColor="green",endColor="yellow",vertexName="goName",getXgmml=TRUE,graphName="enriched terms"){
  library(igraph)
  library(gplots)
	root <- "NULL"
  ontology <- toupper(ontology)
  if(ontology=="BP"||ontology=="bp"){
		root <- "GO:0008150"
	}else if(ontology=="CC"||ontology=="cc"){
		root <- "GO:0005575"
	}else if(ontology=="MF"||ontology=="mf"){
		root <- "GO:0003674"
	}else{stop("Error! unknown ontology ", ontology," not defined\nShould be BP Or CC or MF !\n")}
#	sort enrichedGOTable on pvalues
	enrichedGOTable <- enrichedGOTable[ order(as.numeric(enrichedGOTable[,pvalCol]),decreasing=FALSE),]
#	now it can happen that the go_graph and enrichedGOTable are from different ontology versions, get common to remove inconsistencies and errors
	rownames(enrichedGOTable) <- enrichedGOTable[,goIdCol]
	goIds <- as.character(enrichedGOTable[,goIdCol])
	goIds <- intersect(goIds,V(go_graph)$name)
	if(length(goIds)==0){stop("Error! go_graph and erichedGOTable has no GO terms in common\n")}
	enrichedGOTable <- enrichedGOTable[goIds,]
	goLen <- c(1:length(goIds))
	allAncestors <- unique(unlist (sapply(goLen,function(goLen) 
								get.shortest.paths(graph=go_graph,from=goIds[goLen],
										to=root,mode="out"))))
	subGraph <- induced.subgraph(graph=go_graph,vids=allAncestors)
#	now add colors
	maxColor <- nrow(enrichedGOTable) + 5
	if(nrow(enrichedGOTable)>maxColor){
		maxColor <- nrow(enrichedGOTable) + 10
	}
	palette <- colorRampPalette(colors=c(beginColor, endColor))(maxColor)
	enrichedcolors <- rev(colorpanel(n=nrow(enrichedGOTable),low="#CCFF33",high="#006600"))
	names(enrichedcolors) <- enrichedGOTable[,1]
	graphCol <- rep("#FFFFFF",vcount(subGraph))
	names(graphCol) <- V(subGraph)$name
	graphCol[goIds] <- enrichedcolors 
	graphCol <- graphCol[V(subGraph)$name]
	V(subGraph)$color <- graphCol
	if(getXgmml){
		graph_xgmml <- generate_Xgmml_enrichedGO(go_graph = subGraph,name = graphName)
		return(graph_xgmml)
	}else{
		return(subGraph)
	}
}

#' Given a GO enrichment table, output from function: do_topGO_enrichment OR do_topGO_enrichment_degs
#' generate a directed graph (igraph object) 
generate_GO_gene_graph <- function(enrichment,GO.ID.col=1,GO.name.col=2,enriched.Genes.col=7){
	
}

###############################################################################
#
# Generate an Xgmml object for the given GO graph 
# input_params:
#             go_graph : GO graph (igraph object)
#           nameAttrib : go_graph vertex attribute with GO term names (default: goName)
#          colorAttrib : go_graph vertex attribute with vertex color
#   edgeRelationAttrib : go_graph edge attribute with relationship type between a parent and child node
#                 name : graph name
#
###############################################################################
generate_Xgmml_enrichedGO<-function(go_graph,nameAttrib="goName",colorAttrib="color",
                                    edgeRelationAttrib="RelationshipType",name=NULL){
  library(igraph)
  library(XML)
  #	always throw errors first
  if(is.null(go_graph)){stop("Error! go_graph must be given\n")}
  if(!is.igraph(go_graph)){stop("Error! go_graph not an igraph object\n")}
  #	first print some info back
  message("\n\n")
  message("node count: ",vcount(go_graph),"\n")
  message("edge count: ",ecount(go_graph),"\n\n")
  #	now assuming that the input go_graph as these attributes, or throw warnings
  #	vertices have names ?
  nodeAttrib <- data.frame(nodeId=V(go_graph)$name)
  nodeNames <- get.vertex.attribute(graph=go_graph,name=nameAttrib)
  if(!is.null(nodeNames)){
    nodeAttrib <- data.frame(nodeAttrib,names=nodeNames)
  }else{
    nodeAttrib <- data.frame(nodeAttrib,names=nodeAttrib[,1])
    message("Warning... no vertex names were found under the attribute name ",nameAttrib,"vertex Ids will be used as names\n")
  }
  #	color vertices
  nodeColors <- get.vertex.attribute(graph=go_graph,name=colorAttrib)
  if(!is.null(nodeColors)){
    nodeAttrib <- data.frame(nodeAttrib,color=nodeColors)
  }else{
    nodeAttrib <- data.frame(nodeAttrib,color=rep("#FFFFFF",vcount(go_graph)))
    message("Warning... no vertex colors were found under the attribute name ",colorAttrib," all vertices will be white\n")
  }
  #	now complain for name
  graphName="gene ontology graph"
  if(is.null(name)){
    message("Warning... no graph name was given, default gene ontology graph will be used as name\n")
  }else{
    graphName <- name
  }
  #	some XGMML thingies
  #	graph node
  graphNode <- newXMLNode("graph", attrs = c("directed" = "0","label"=graphName),namespace=c("xmlns"="http://www.cs.rpi.edu/XGMML"),
                          namespaceDefinitions = c("cy" = "http://www.cytoscape.org","dc" = "http://purl.org/dc/elements/1.1/",
                                                   "rdf" = "http://www.w3.org/1999/02/22-rdf-syntax-ns#","xlink" = "http://www.w3.org/1999/xlink"))
  #		graph attributes
  att1 <- newXMLNode("att",attrs=c(name="documentVersion", value="1.1"),parent=graphNode)
  att2 <- newXMLNode("att",attrs=c(name="networkMetadata"),parent=att1)		
  backAtt <- newXMLNode("att",attrs=c(name="backgroundColor",type="string",value="#FFFFFF"),parent=graphNode)
  desAtt <- newXMLNode("att",attrs=c(description=as.character(graphName),type="string"))
  #	add all the nodes
  for(go in 1:nrow(nodeAttrib)){
    #		go node
    goNode <- newXMLNode("node",attrs=c(id=as.character(nodeAttrib[go,1]), 
                                        label=as.character(nodeAttrib[go,1])),parent=graphNode)
    #		attributes
    geneNameAtt <- newXMLNode("att",attrs=c(name="goName",type="string",value=as.character(nodeAttrib[go,2])),
                              parent=goNode)
    #		graphics
    nodeGraphics <- newXMLNode("graphics",attrs=c(fill=as.character(nodeAttrib[go,3]),width="1",outline="#666666","cy:nodeTransparency"="1.0",
                                                  "cy:nodeLabelFont"="SansSerif.bold-0-12","cy:nodeLabel"=as.character(nodeAttrib[go,2]),type="ELLIPSE",h="40.0",w="40.0","cy:borderLineType"="solid"),
                               namespaceDefinitions = c("cy" = ""),parent=goNode)
  }
  #	make edge data.frames
  edge_el <- as.data.frame(get.edgelist(go_graph))
  edgeRelation <- get.edge.attribute(graph=go_graph,name=edgeRelationAttrib)
  if(is.null(edgeRelation)){message("Warning no edge relationship was found under", edgeRelationAttrib," graph will be undirected!\n")}
  else{
    edge_el <- data.frame(edge_el,relation=edgeRelation)
  }
  #	add edges and edge attributes
  for(el in 1:nrow(edge_el)){
    #		edge label and node
    labelString <- paste(as.character(edge_el[el,1]),"(pp)",as.character(edge_el[el,2]),sep="")
    edgeInt <- newXMLNode("edge",attrs=c(source=as.character(edge_el[el,1]),target=as.character(edge_el[el,2]), 
                                         label=labelString),parent=graphNode)
    if(is.null(edgeRelation)){
      #			edge graphics for undirected Gene Ontology graph
      cyGraphics <- newXMLNode("graphics", attrs = c("cy:edgeLineType" = "SOLID"
                                                     ,"width"=1,fill="#999999","cy:curved"="STRAIGHT_LINES",
                                                     "cy:sourceArrow"="0","cy:targetArrow"="0","cy:sourceArrowColor"="#000000",
                                                     "cy:targetArrowColor"="#000000","cy:edgeLabelFont"="Default-0-10"),
                               namespaceDefinitions = c("cy" = ""),parent=edgeInt)
    }else{
      #			give edge attribute
      eType <- newXMLNode("att",attrs=c(name="relationType",type="string", value=as.character(edge_el[el,3])),parent=edgeInt)
      relationType <- tolower(as.character(edge_el[el,3]))
      edgeColor <- "#000000" # default is_a
      targetArrowColor <- "#000000"
      if(relationType=="part_of"){
        edgeColor <- "#0000FF"
        targetArrowColor <- "#0000FF"
      }else if(relationType=="has_part"){
        edgeColor <- "#660066"
        targetArrowColor <- "#660066"
      }else if(relationType=="regulates"){
        edgeColor <- "#FF6600"
        targetArrowColor <- "#FF6600"
      }else if(relationType=="positively_regulates"){
        edgeColor <- "#006600"
        targetArrowColor <- "#006600"
      }else if(relationType=="negatively_regulates"){
        edgeColor <- "#FF0000"
        targetArrowColor <- "#FF0000"
      }else if(relationType=="occurs_in"){
        edgeColor <- "#006666"
        targetArrowColor <- "#006666"
      }
      cyGraphics <- newXMLNode("graphics", attrs = c("cy:edgeLineType" = "SOLID"
                                                     ,"width"=1,fill=as.character(edgeColor),"cy:curved"="STRAIGHT_LINES",
                                                     "cy:sourceArrow"="0","cy:targetArrow"="6","cy:sourceArrowColor"="#000000",
                                                     "cy:targetArrowColor"=as.character(targetArrowColor),"cy:edgeLabelFont"="Default-0-10"),
                               namespaceDefinitions = c("cy" = ""),parent=edgeInt)
    }
  }
  return(graphNode)
}


###############################################################################
#
# for a given list of topGO enrichment analysis tables, generate GO graph
# list wrapper function for the function generate_GO_graph()
# input_params:
#             enrichmentList : topGO enrichment list (output from function do_topGO_enrichment())
#                   go_graph : GO graph (igraph object, output from function get_go_graph())
#                    goIdCol : GO term id column (default : 1)
#                    pvalCol : P-value colmn (default : 1)
#                 beginColor : Color for the most enriched terms (default : green)
#                   endColor : Color for the least enriched terms (default : yellow)
#                 vertexName : go_graph vertex attribute name for GO term names (default for go_graph : go_Name)
#                   getXgmml : Boolean, return Xgmml or igraph object (default : TRUE)
#
###############################################################################
generate_GO_graph_enrichmentList <- function(enrichmentList,go_graph,ontology="BP",goIdCol=1,
                                             pvalCol=6,beginColor="green",endColor="yellow",vertexName="goName",getXgmml=TRUE){
  if(is.null(enrichmentList)){stop("Error! enrichmentList cannot be empty\n")}
  outList <- list()
  for(ET in 1:length(enrichmentList)){
    erichmentTable <- enrichmentList[[ET]]$enrichment
    #		graph name will be generated based on cluster id
    graphName <- paste("DEG class ",names(enrichmentList[ET])," GO enrichment",sep=" ")
    graph_visual<-generate_GO_graph(go_graph = go_graph,ontology = ontology,goIdCol = goIdCol,pvalCol = pvalCol,beginColor = beginColor,
                                    endColor = endColor,vertexName = vertexName,getXgmml = getXgmml,graphName = graphName,enrichedGOTable = erichmentTable)
    outList[[names(enrichmentList[ET])]] <- graph_visual
  }
  return(outList)
}

###############################################################################
#
# save the enrichment tables generated by the functions do_topGO_enrichment() OR get_geenes_topGO_enrichment()
#
###############################################################################
save_enrichmentList <- function(enrichmentList,fileNameHead="deg_list_",fileNameTail=".csv"){
  if(is.null(enrichmentList)){stop("Error! enrichmentList cannot be empty\n")}
  for(el in 1:length(enrichmentList)){
    file_name = gsub(pattern="\\s{1,}",replacement="\\_",x=names(enrichmentList[el]))
    saveName <- paste(fileNameHead,file_name,fileNameTail,sep="")
    saveTable <- enrichmentList[[el]]$enrichment
    write.table(x=saveTable,file=saveName,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
    message(" saved ",file_name," enrichment list: ",saveName,"\n")
  }
}

###############################################################################
#
# save the enrichment tables generated by the functions do_topGO_enrichment() OR get_geenes_topGO_enrichment()
# into a .xls file. Package WriteXLS MUST BE installed
#
###############################################################################
save_enrichmentList_xls <- function(enrichmentList,fileNameHead="deg_list_"){
	library(WriteXLS)
	fileNameTail=".xls"
	if(is.null(enrichmentList)){stop("Error! enrichmentList cannot be empty\n")}
	for(el in 1:length(enrichmentList)){
		file_name = gsub(pattern="\\s{1,}",replacement="\\_",x=names(enrichmentList[el]))
		saveName <- paste(fileNameHead,file_name,fileNameTail,sep="_")
		saveName <-  gsub(pattern="\\_{2,}",replacement="\\_",x=saveName)
		saveName <-  gsub(pattern="\\_{1,}\\.",replacement="\\.",x=saveName)
		saveTable <- enrichmentList[[el]]$enrichment
#		tableName <- deparse(substitute(saveTable))
		WriteXLS(x=list(quote(saveTable)),ExcelFileName=saveName,SheetNames = "Enrichment table",row.names=FALSE,col.names=TRUE)
		message(" saved ",file_name," enrichment list: ",saveName,"\n")
	}
} 

###############################################################################
#
# save the list of xgmml files generated by function generate_GO_graph_enrichmentList()
# input_params:
#              visualEnrichmentList : list of xgmml objects (generated by function generate_GO_graph_enrichmentList())
#                      fileNameHead : header for file names (for output file naming)
#                      fileNameTail : file name extension (for output files)
#
###############################################################################
save_enrichmentList_visual<- function(visualEnrichmentList,fileNameHead="deg_list_",fileNameTail=".xgmml"){
  library(XML)
  if(is.null(visualEnrichmentList)){stop("Error! visualEnrichmentList cannot be empty\n")}
  for(vel in 1:length(visualEnrichmentList)){
    file_name = gsub(pattern="\\s{1,}",replacement="\\_",x=names(visualEnrichmentList[vel]))
    saveName <- paste(fileNameHead,file_name,fileNameTail,sep="")
    message(" file name : ",saveName,"\n")
    saveXML(visualEnrichmentList[[vel]],file=saveName,encoding="UTF-8")
  }
}

