# =================== SOURCE ========================

# Own source code
# To use script in RStudio - set right working directory 
# or use the next code commented out:
path_wd = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path_wd)

func_graph = "graph.R"
func_models = "models.R"
func_pathways = "pathways.R"
source(func_graph)
source(func_models)
source(func_pathways)

# =================== DATA ===========================

dExpr = read.table('dExpr.txt')
dCon = dExpr[dExpr[,'disease'] == 0,]
dSch = dExpr[dExpr[,'disease'] == 1,]

# Covariates for the data, which are common for both disease groups
common_factors = c('Age', 'Sex', 'Pca1', 'Pca2', 'batch')

# =================== PIPLINE ========================

# Pathway to use
pathway_name = 'MAPK'
# pathway_name = 'FocalAdh'
# pathway_name = 'Neurotrophin'
# pathway_name = 'AKT'
# pathway_name = 'Serine'

pathway = get_pathway(pathway_name)

# choose the type of graph
flag = FALSE # inital graph

property_mx = get_graph_properties(pathway, flag)
pathway$genes_mx = property_mx$genes_mx
pathway$edges_mx = property_mx$edges_mx

graph = construct_graph(pathway, flag)

render_graph(graph)

#====================================================
#        Create all possible path models
#====================================================

path_mod = construct_path_mod(pathway, common_factors)

length(path_mod$group)


#====================================================
#             Estimate all models
#====================================================

est = estimate_path_mod(path_mod$group, dExpr, common_factors)

#====================================================
# Select the best model
#====================================================


# Failed models
idx_failed = est$indexes[,'df'] != ''

# Models with many non-significant interactions
min(est$pvals_c[idx_failed])
idx_signif = as.numeric(est$pvals_c) <= 4



idx = 1:length(path_mod$group)
idx = idx[idx_failed & idx_signif]
plot(est$indexes[idx,'cfi'],est$indexes[idx,'rmsea'])
# idx = idx[est$indexes[idx,'cfi'] == max(est$indexes[idx,'cfi'])]

idx = idx[est$indexes[idx,'cfi'] > 0.75]
idx = idx[est$indexes[idx,'rmsea'] < 0.2]
idx = idx[est$indexes[idx,'tli'] == min(est$indexes[idx,'tli'])]

mod_group = path_mod$group[[idx[1]]]

est$indexes[idx,'cfi']

0.750623685378339

0.76118576693227

#====================================================
#             WHEN MODEL IS PREPARED
#====================================================


fit_groups = sem(mod_group, data = dExpr, meanstructure = TRUE , fixed.x = FALSE, group = 'disease')

# Add label propery (estimation) to edges
pathway_fit = construct_pathway_fit(pathway, fit_groups, common_factors)

# Calculate Information criteria for separate models
pathway_fit = get_indexes(pathway_fit, mod_group, fit_groups, dCon, dSch, common_factors);

graph_fit = construct_graph(pathway_fit, TRUE)

for (i in 1:length(pathway_fit$indexes))
  graph_fit = add_node(graph_fit, label = pathway_fit$indexes[i], node = i)

graph_fit = set_node_attr(graph_fit, nodes = c(1:4), node_attr = 'color', values = 'white')

render_graph(graph_fit)






