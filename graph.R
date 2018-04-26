#' ---
#' title: "Analysis of gene expression variance in schizophrenia using structural equation modeling"
#' purpose: "Construct and render graphs"
#' author: "Anna A. Igolkina (igolkinaanna11@gmail.com)"
#' date: "2017-2018"
#' ---


# Library to work with diagrams
library('devtools')
install("r-diagrammer-0.8.2-r3.3.1_0/lib/R/library/Diagrammer")
library('DiagrammeR')

get_graph_properties <- function(pathway, flag)
{
  genes_info = pathway$genes
  edges_info = pathway$edges
  genes_n_attr = 3;
  
  genes_mx_tmp = matrix(genes_info, ncol = length(genes_info)/genes_n_attr, nrow = genes_n_attr)
  
  genes_mx = list()
  genes_mx$name = genes_mx_tmp[1,]
  genes_mx$label = genes_mx_tmp[2,]
  genes_mx$type = genes_mx_tmp[3,]
  genes_mx$islatent = genes_mx_tmp[3,] == 'oval'
  genes_mx$n = length(genes_mx$name)
  
  genes_mx$variants = list()
  for (i in 1:length(genes_mx$name))
  {
    #if (genes_mx[3,i] == 'point')
    #  next
    s = strsplit(genes_mx$label[i], split = '\n');
    s = s[[1]]
    s_first <- s[1]
    colon <- ":"
    if (grepl(pattern = colon, x = s_first))
      s = s[-1]
    genes_mx$variants[[as.character(genes_mx$name[i])]] = s
  }
  

  edges_n_attr = 2;
  edges_mx_tmp = matrix(edges_info, ncol = length(edges_info)/edges_n_attr, nrow = edges_n_attr)
  edges_mx = list()
  edges_mx$from = edges_mx_tmp[1,]
  edges_mx$to = edges_mx_tmp[2,]
  edges_mx$label = ''
  edges_mx$signif = c()
  edges_mx$n = length(edges_mx$from)
  
  
 
  property_mx = list()
  property_mx$genes_mx = genes_mx
  property_mx$edges_mx = edges_mx
  return (property_mx)
}
#========================================================================
#
#========================================================================

construct_graph <- function(pathway, flag)
{  
  genes_mx = pathway$genes_mx
  edges_mx = pathway$edges_mx
  
  if (flag == FALSE) {
    egdes_label = ''
    edges_color = 'black'
  } else {
    egdes_label = edges_mx$label
    # according to the significance
    if (length(edges_mx$pval) == 0)
      edges_color = 'black'
    else
    {
      edge_col_range = c('SteelBlue', 'Goldenrod', 'OrangeRed', 'Crimson')
      thsholds = c(1, 0.05, 0.01, 0.001)
      
      edges_color = rep('', edges_mx$n)
      for (i in 1:length(thsholds))
        edges_color[edges_mx$pval < thsholds[i]] = edge_col_range[i]
    }
  }
  

  nodes = 
    create_nodes(n = length(genes_mx$name),
                 nodes = genes_mx$name,
                 label = genes_mx$label,
                 type = "lower",
                 style = "empty",
                 color = 'black',
                 shape = genes_mx$type
    );
  edges = 
    create_edges( from = edges_mx$from,
                  to = edges_mx$to,
                  label = egdes_label,
                  arrowhead = 'normal',
                  arrowhead = 'none',
                  color = edges_color);
  graph =
    create_graph(nodes = nodes,
                 edges = edges,
                 node_attrs = "fontname = Helvetica",
                 edge_attrs = c("color = black",
                                "arrowsize = 1",
                                "fontname = Helvetica",
                                "penwidth = 2"),
                 graph_attrs = "layout = dot")

  return (graph)
}





#========================================================================
#
#========================================================================


construct_graph_fit_single <- function(genes_mx, edges_mx, nodes, fit_single)
{

  
  params_sem = parameterEstimates(fit_single);
  regressions = params_sem[params_sem[,'op']=='~', c('lhs', 'rhs', 'est', 'pvalue')];
  
  genes = unique(c(regressions[,1], regressions[,2]));
  
  nodes_names = names(nodes)
  gene_names = rep('', length(nodes_names));
  for (i in 1:length(nodes))
  {
    intersection = intersect(nodes[[i]], genes)
    if (length(intersection) > 0)
    {
      gene_names[i] = intersection[1]
    }
  }
  
  for (i in 1:dim(regressions)[1])
  {
    i
    regressions[i,1] =  nodes_names[gene_names == regressions[i,1]]
    regressions[i,2] =  nodes_names[gene_names == regressions[i,2]]
  }
  
  charact_color = rep('',dim(edges_mx)[2])
  charact_label = rep('',dim(edges_mx)[2])
  for (i in 1:dim(edges_mx)[2])
  {
    idx = (regressions[,2] == edges_mx[1,i]) & (regressions[,1] == edges_mx[2,i])
    if (sum(idx) > 0)
    {
      charact_label[i] = round(regressions[idx, 3],digits=2);
    }
  }
    
  genes_mx_fit = genes_mx;
  for (i in 1:dim(genes_mx_fit)[2])
  {
    idx = (nodes_names == genes_mx_fit[1,i])

    if (sum(idx) > 0)
    {
      genes_mx_fit[2,i] = gene_names[idx]
      #genes_mx_fit[3,i] = 'rectangle'
      genes_mx_fit[6,i] = 1
    }
  }
  
  return (list(charact_label,genes_mx_fit))
}
