# Library to work with SEM
library(lavaan)

#========================================================================
#
#========================================================================
modDouble2Single <-function(mod)
{
  idx_cut = c()
  for(i in 1:length(mod))
  {
    if (grepl(pattern = ':', x = mod[i]))
    {
      idx_cut = c(idx_cut, i);
    }else
    {
      mod[i] = gsub(pattern = 'c\\(\\w+, \\w+\\)\\*', replacement = '',x = mod[i])
    }
  }
  mod = as.matrix(mod[-idx_cut])
}
#========================================================================
#
#========================================================================
# split_mod_into_fit <-function(fit_groups)
# {
#   params_sem = parameterEstimates(fit_groups);
#   mod1 = c()
#   mod2 = c()
#   for (i in 1:dim(params_sem)[1])
#   {
#     if
#   }
# 
#   return(list(mod1, mod2))
# }

#========================================================================
#
#========================================================================
modDouble2Single_fit <-function(mod, fit_groups, common_factors)
{
  params_sem0 = parameterEstimates(fit_groups);
  params_sem = params_sem0;
  for (i in 1:length(common_factors))
  {
    params_sem = params_sem[params_sem[,1] != common_factors[i],]
    params_sem = params_sem[params_sem[,3] != common_factors[i],]
  }
  
  regr_est1 = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='1')&(params_sem[,'label']!=''), c('est')];
  regr_lab1 = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='1')&(params_sem[,'label']!=''), c('label')];
  regr_est2 = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='2')&(params_sem[,'label']!=''), c('est')];
  regr_lab2 = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='2')&(params_sem[,'label']!=''), c('label')];
  
  params_sem = params_sem0;
  cov_est1 = params_sem[(params_sem[,'op']=='~~') &(params_sem[,'group']=='1')&(params_sem[,'label']!=''), c('est')];
  cov_lab1 = params_sem[(params_sem[,'op']=='~~') &(params_sem[,'group']=='1')&(params_sem[,'label']!=''), c('label')];
  cov_est2 = params_sem[(params_sem[,'op']=='~~') &(params_sem[,'group']=='2')&(params_sem[,'label']!=''), c('est')];
  cov_lab2 = params_sem[(params_sem[,'op']=='~~') &(params_sem[,'group']=='2')&(params_sem[,'label']!=''), c('label')];
  
  cov_est1 = c(cov_est1, regr_est1)
  cov_est2 = c(cov_est2, regr_est2)
  cov_lab1 = c(cov_lab1, regr_lab1)
  cov_lab2 = c(cov_lab2, regr_lab2)
  

  mod1 = c();
  mod2 = c();
  for(i in 1:length(mod))
  {
    if (grepl(pattern = ':', x = mod[i]))  # do not use variables
    {
      next;
    }
    break_flag = FALSE
    for (j in 1:length(cov_lab1))
    {
      if (length(grep(cov_lab1[j], x =  mod[i])) != 0)
      {
        mod1 = c(mod1,gsub(pattern = 'c\\(\\w+, \\w+\\)', replacement = cov_est1[j],x = mod[i]))
        break_flag = TRUE
        break;
      }
    }
    for (j in 1:length(cov_lab2))
    {
      if (length(grep(cov_lab2[j], x =  mod[i])) != 0)
      {
        mod2 = c(mod2,gsub(pattern = 'c\\(\\w+, \\w+\\)', replacement = cov_est2[j],x = mod[i]))
        break_flag = TRUE
        break;
      }
    }
    if (break_flag)
      next;
    
    if (length(grep(' ~ ', x =  mod[i])) == 0)  # if the line is about covariance
    {
      mod1 = c(mod1,  mod[i]);
      mod2 = c(mod2,  mod[i]);
      next;
    }
    
  }
  mod1 = as.matrix(mod1);
  mod2 = as.matrix(mod2);
  
  return(list(mod1, mod2))
}

#========================================================================
#
#========================================================================
construct_path_mod <- function(pathway, common_factors)
{
  
  genes_mx = pathway$genes_mx
  edges_mx = pathway$edges_mx
  
  # construct the base model with common ID
  
  mod_base = c();
  mod_base_single = c();
  for (i in 1:edges_mx$n)
  {
    mod_base_single = rbind(mod_base_single,paste(c(edges_mx$to[i], ' ~ ', 
                                                    edges_mx$from[i]), collapse = ''))
    mod_base = rbind(mod_base,paste(c(edges_mx$to[i], ' ~ ', 
                                      'c(a',i,'_1, a',i,'_2)*',
                                      edges_mx$from[i]), collapse = ''))
    mod_base = rbind(mod_base, paste(c('d',i,' := a',i,'_1 - a',i,'_2'), collapse = ''))
  }
  genes_exo = setdiff(edges_mx$from, edges_mx$to);
  
  # Set covariances between exogen variables to zero
  if (length(genes_exo) > 1)
  {
    for (i in 1:(length(genes_exo)-1))
    {
      for (j in (i+1):length(genes_exo))
      {
        mod_base_single = rbind(mod_base_single,paste(c(genes_exo[i], ' ~~ 0*', genes_exo[j]), collapse = ''))
        mod_base = rbind(mod_base,paste(c(genes_exo[i], ' ~~ 0*', genes_exo[j]), collapse = ''))
      }
    }
  }
  
  # # Set variances ?
  # for (i in 1:genes_mx$n)
  # {
  #   if (!genes_mx$islatent[i])
  #   {
  # 
  #       mod_base = rbind(mod_base,paste(c(genes_mx$name[i], ' ~~ ',
  #                                         'c(v0',i,'_',j, ', v0',i,'_',j,')*',
  #                                         genes_mx$name[i]), collapse = ''))
  #       mod_base_single = rbind(mod_base_single,paste(c(genes_mx$name[i], ' ~~ ', genes_mx$name[i]), collapse = ''))
  #   }
  # }

  
  # Additional factors: common_factors
  for (i in 1:genes_mx$n)
  {
    if (!genes_mx$islatent[i])
    {
      for (j in 1:length(common_factors))
      {
        
        mod_base = rbind(mod_base,paste(c(genes_mx$name[i], ' ~ ',
                                          'c(a0',i,'_',j, ', a0',i,'_',j,')*',
                                          common_factors[j]), collapse = ''))
        mod_base_single = rbind(mod_base_single,paste(c(genes_mx$name[i], ' ~ ', common_factors[j]), collapse = ''))
      }
    }
  }
  
  
  # add latent expressions
  for (i in 1:genes_mx$n)
  {
    if (genes_mx$islatent[i])
    {
      mod_tmp_single = mod_base_single;
      mod_tmp = mod_base;
      mod_tmp_single = rbind(mod_tmp_single, paste(c(genes_mx$name[i], ' =~ NA*', paste(sort(nodes[[genes_mx$name[i]]]), collapse = ' + ')), collapse = '')) 
      mod_tmp_single = rbind(mod_tmp_single, paste(c(genes_mx$name[i], ' ~~ 1*', genes_mx$name[i]), collapse = ''))      
      mod_tmp = rbind(mod_tmp, paste(c(genes_mx$name[i], ' =~ NA*', paste(sort(nodes[[genes_mx$name[i]]]), collapse = ' + ')), collapse = '')) 
      mod_tmp = rbind(mod_tmp, paste(c(genes_mx$name[i], ' ~~ 1*', genes_mx$name[i]), collapse = ''))              
      
      mod_base_single = mod_tmp_single
      mod_base = mod_tmp
      next
    }
  }
  
  
  
  # contruct all possible models
  models_single = list()
  models = list()
  
  models_single[[1]] = mod_base_single
  models[[1]] = mod_base
  for (i in 1:genes_mx$n)
  {
    if (genes_mx$islatent[i])
    {
      next
    }
    n_models = length(models);
    for (j in 1:length(models))
    {
      for (gene_variant in genes_mx$variants[[genes_mx$name[i]]])
      {
        mod_tmp_single = models_single[[j]];
        mod_tmp = models[[j]];
        for (l in 1:length(mod_tmp_single))
        {
          mod_tmp_line_s = gsub(x = mod_tmp_single[l], pattern = genes_mx$name[i], gene_variant);
          mod_tmp_single[l] = mod_tmp_line_s
        }      
        for (l in 1:length(mod_tmp))
        {
          mod_tmp_line = gsub(x = mod_tmp[l], pattern = genes_mx$name[i], gene_variant);
          mod_tmp[l] = mod_tmp_line
        }
        models_single = c(models_single, list(mod_tmp_single))
        models = c(models, list(mod_tmp))
        
      }
    }

    # Remove models, where genes were not replaced
    models_single = models_single[-(1:n_models)]
    models = models[-(1:n_models)]

  }
  path_models = list()
  path_models$single = models_single
  path_models$group = models
  return (path_models)
}
#========================================================================
#
#========================================================================

estimate_path_mod <- function(models, dExpr, common_factors)
{
  # param "models" - group models
  # fing aic for models on both control and schuzophrenia data
  aic_res = c()
  fit_measures = c()
  pvals = c()
  pvals_s = c()
  pvals_c = c()
  for (i in 1:length(models))
  {
    print(i)
    mod = models[[i]];
    fit_tmp <- tryCatch(sem(mod, data = dExpr, meanstructure = TRUE , fixed.x = FALSE, group = 'disease'),
                        warning = function(w) {NaN;},
                        error = function(e) {NaN;})  
    if (class(fit_tmp) != "numeric")
    {
      #aic_res = rbind(aic_res, AIC(fit_tmp))
      aic_res = rbind(aic_res, -fitMeasures(fit_tmp, 'cfi'))
      fit_measures = rbind(fit_measures, fitMeasures(fit_tmp))
      
      params_sem = parameterEstimates(fit_tmp);
      for (i in 1:length(common_factors))
      {
        params_sem = params_sem[params_sem[,1] != common_factors[i],]
        params_sem = params_sem[params_sem[,3] != common_factors[i],]
      }
      regressions1 = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='1'), c('pvalue')];
      regressions2 = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='2'), c('pvalue')];
      pvals = rbind(pvals, sum((regressions1 > 0.01) & (regressions2 > 0.01)))
      pvals_c = rbind(pvals_c, sum(regressions1 > 0.01))
      pvals_s = rbind(pvals_s, sum(regressions2 > 0.01))
    }else
    {
      aic_res = rbind(aic_res, 1000)
      fit_measures = rbind(fit_measures, rep('',41))
      pvals = rbind(pvals, '')
      pvals_c = rbind(pvals_c, '')
      pvals_s = rbind(pvals_s, '')
    }
  }
  
  
  estimation = list()
  estimation$indexes = fit_measures
  estimation$pvals_g = pvals
  estimation$pvals_c = pvals_c
  estimation$pvals_s = pvals_s
  
  return (estimation)
}  

#========================================================================
#
#========================================================================


construct_pathway_fit <- function(pathway, fit_groups, common_factors)
{
  genes_mx = pathway$genes_mx
  edges_mx = pathway$edges_mx
  
  params_sem = parameterEstimates(fit_groups);
  for (i in 1:length(common_factors))
  {
    params_sem = params_sem[params_sem[,1] != common_factors[i],]
    params_sem = params_sem[params_sem[,3] != common_factors[i],]
  }
  
  edge_from_to = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='1'), c('lhs', 'rhs')];
  regressions_gr1 = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='1'), c( 'est', 'se','pvalue')];
  regressions_gr2 = params_sem[(params_sem[,'op']=='~') &(params_sem[,'group']=='2'), c( 'est', 'se','pvalue')];
  regressions_diff = params_sem[(params_sem[,'op']==':=') &(params_sem[,'group']=='0'), c('est', 'se','pvalue')];

  
  signif_diff = regressions_diff[,'pvalue']
  
  
  # Construct edge labels
  edge_labels = c()
  for (i in 1:length(signif_diff))
  {
    est_gr1 = sprintf('%.02f', as.numeric(regressions_gr1[i,'est']))
    est_gr2 = sprintf('%.02f', as.numeric(regressions_gr2[i,'est']))
    
    se_gr1 = sprintf('%.02f', as.numeric(regressions_gr1[i,'se']))
    se_gr2 = sprintf('%.02f', as.numeric(regressions_gr2[i,'se']))
    
    s_gr1 = paste(c(est_gr1, se_gr1), collapse = '±')
    s_gr2 = paste(c(est_gr2, se_gr2), collapse = '±')
    
    signif_tmp = 'n.s.'
    if (signif_diff[i] < 0.001) { signif_tmp = '***'
    } else if (signif_diff[i] < 0.01) {signif_tmp = '**'
    } else if (signif_diff[i] < 0.05) {signif_tmp = '*'}
    
    s_diff = paste(c('Δ', signif_tmp), collapse = ': ')
    
    edge_labels = c(edge_labels, paste(c(s_gr1, s_gr2,s_diff), collapse = '\n'))    
    
  }
  
  
  # Names of noves across variants
  genes = unique(c(edge_from_to[,1], edge_from_to[,2]));
  print(genes)
  nodes_name = genes_mx$name;
  for (i in 1:genes_mx$n)
  {
    intersection = intersect(genes_mx$variants[[i]], genes)
    if (length(intersection) > 0)
    {
      nodes_name[i] = intersection[1]
    }
  }
  
  pathway$genes_mx$label = nodes_name
  pathway$edges_mx$label = edge_labels
  pathway$edges_mx$pval = signif_diff
  
  # latent genename change
  # latent = params_sem[(params_sem[,'op']=='=~') &(params_sem[,'group']=='1'), c('lhs', 'rhs')];
  # latent_uniq = unique(latent[,1])
  # for (i in 1:length(gene_names))
  # {
  #   intersection = intersect(gene_names[i], latent_uniq)
  #   if (length(intersection) > 0)
  #   {
  #     gene_names[i] = paste(latent[latent[,1] == intersection,2], collapse = '\n')
  #   }
  # }
  
  
  
  return (pathway)
}
#========================================================================
#
#========================================================================

get_indexes <- function(pathway_fit, mod_group, fit_groups, dCon, dSch, common_factors)
{
  
  mod12 = modDouble2Single_fit(mod_group, fit_groups, common_factors)
  fit_c = sem(mod12[[1]], data = dCon, meanstructure = TRUE , fixed.x = FALSE)
  print(parameterestimates(fit_c))
  fit_s = sem(mod12[[2]], data = dSch, meanstructure = TRUE , fixed.x = FALSE)
  
  c <- paste(c('CFI:',round(fitMeasures(fit_c,'cfi'),digits=3), round(fitMeasures(fit_s,'cfi'),digits=3)), collapse='\n')
  d <- paste(c('RMSEA:',round(fitMeasures(fit_c,'rmsea'),digits=3), round(fitMeasures(fit_s,'rmsea'),digits=3)), collapse='\n')
  # e <- paste(c('PGFI:',round(fitMeasures(fit_c,'pgfi'),digits=3), round(fitMeasures(fit_s,'pgfi'),digits=3)), collapse='\n')
  a <- paste(c('AIC:',round(AIC(fit_c),digits=1), round(AIC(fit_s),digits=1)), collapse='\n')
  b <- paste(c('BIC:',round(BIC(fit_c),digits=1), round(BIC(fit_s),digits=1)), collapse='\n')
  
  pathway_fit$indexes = c(c,d,a,b)
  
  return (pathway_fit)
}

