#' ---
#' title: "Analysis of gene expression variance in schizophrenia using structural equation modeling"
#' purpose: "Topologies of pathways wto  analyse"
#' author: "Anna A. Igolkina (igolkinaanna11@gmail.com)"
#' date: "2017-2018"
#' ---



get_pathway <- function(pathway_name)
{
  
  # 'rectange' means observed variable
  # 'oval' means the latern variable
  
  
  #======================================================
  #                      MAPK
  #======================================================
  genes_info_mapk = c(
    'SOS','SOS1\nSOS2','rectangle',
    'Ras','HRAS\nKRAS\nNRAS\nMRAS','rectangle',
    'RASA2', 'RASA2','rectangle',
    'NF1', 'NF1','rectangle',
    'RAPGEF2', 'RAPGEF2', 'rectangle',
    'PRKC','PRKCA\nPRKCB\nPRKCG','rectangle',
    'Raf1', 'RAF1', 'rectangle',
    'BRaf', 'BRAF', 'rectangle',
    'MEK','MEK:\nMAP2K1\nMAP2K2','rectangle',
    'ERK','ERK:\nMAPK1\nMAPK3','rectangle',
    'Rap','RAP1A\nRAP1B','rectangle'
    
  );
  
  edges_info_mapk = c(
    'SOS', 'Ras',
    'RASA2', 'Ras',
    'NF1', 'Ras',
    'RAPGEF2', 'Ras',
    'PRKC', 'Ras',
    'PRKC', 'Raf1',
    'Ras', 'Raf1',
    'Raf1', 'MEK',
    'MEK', 'ERK',
    'ERK','Raf1',
    'Rap','BRaf',
    'BRaf', 'MEK'
    
  )
  
  
  #======================================================
  #                      FOCAL1
  #======================================================
  
  
  genes_info_focal1 = c(
    'PTK2','PTK2', 'rectangle',
    'BCAR1', 'BCAR1', 'rectangle',
    'Crk', 'CRKL\nCRK', 'rectangle',
    'DOCK1', 'DOCK1', 'rectangle',
    'Pik3', 'PIK3CA\nPIK3CB\nPIK3CD\nPIK3R1\nPIK3R2\nPIK3R3', 'rectangle',
    'Vav', 'VAV2\nVAV1\nVAV3', 'rectangle',
    'Rac', 'RAC1\nRAC2\nRAC3', 'rectangle',
    'Pak', 'PAK1\nPAK2\nPAK3\nPAK4\nPAK5\nPAK6', 'rectangle',
    'RAPGEF1','RAPGEF1', 'rectangle',
    'Rap', 'RAP1B\nRAP1A', 'rectangle',
    'PTEN','PTEN', 'rectangle'
    
  );
  
  
  edges_info_focal1 = c( 
    'PTK2','BCAR1',
    'PTK2','Pik3',
    'BCAR1', 'Crk',
    'Crk','DOCK1',
    'Pik3','Vav',
    'PTEN','Vav',
    'PTEN','Pik3',
    'Vav', 'Rac',
    'DOCK1', 'Rac',
    'Rac', 'Pak',
    'Crk','RAPGEF1',
    'RAPGEF1','Rap'
  );
  
  
  
  #======================================================
  #                      NEUROTROPHIN
  #======================================================
  genes_info_neuro = c(
    'NGFR','NGFR', 'rectangle',
    'TRAF6','TRAF6', 'rectangle',
    'RAC1','RAC1', 'rectangle',
    'MAP3K1','MAP3K1', 'rectangle',
    'MAP2K7','MAP2K7', 'rectangle',
    'Jnk','MAP8\nMAP9\nMAP10', 'rectangle',
    'TP73','TP73', 'rectangle',
    'TP53','TP53', 'rectangle',
    'JUN','JUN', 'rectangle'
    
  )
  
  
  edges_info_neuro = c( 
    'NGFR','TRAF6',
    'NGFR','RAC1',
    'RAC1','MAP3K1',
    'MAP3K1','MAP2K7',
    'MAP2K7','Jnk',
    'TRAF6','Jnk',
    'Jnk','TP53',
    'TP73','TP53',
    'Jnk','JUN'
    
  )
  
  
  #======================================================
  #                     PI3K-AKT
  #======================================================
  
  
  
  genes_info_akt = c(
    'IRS1','IRS1', 'rectangle',
    'Pik3', 'PIK3CA\nPIK3CB\nPIK3CG\nPIK3R1\nPIK3R2\nPIK3R3', 'rectangle',
    'PTK2','PTK2', 'rectangle',
    'PDPK1','PDPK1', 'rectangle',
    'PTEN','PTEN', 'rectangle',
    'Akt','AKT3', 'rectangle',
    'Mtor','CRTC2', 'rectangle',
    'Foxo','FOXO3', 'rectangle',
    'BAD','BAD', 'rectangle',
    'GSK3B','GSK3B', 'rectangle',
    'Tsc','TSC1\nTSC2', 'rectangle'
    
  )
  
  
  edges_info_akt = c(
    'IRS1','Pik3',
    'PTK2','Pik3',
    'PDPK1','Pik3',
    'Pik3','Akt',
    'PDPK1','Akt',
    'PTEN','Akt',
    'Mtor','Akt',
    'Akt','Foxo',
    'Akt','Tsc',
    'Akt','GSK3B',
    'Akt','BAD'
  )
  
  
  
  #======================================================	
  #          SERINE
  #======================================================
  
  
  
  genes_info_serine = c(
    'PHGDH', 'PHGDH', 'rectangle',
    'PSAT1', 'PSAT1', 'rectangle',
    'PSPH', 'PSPH', 'rectangle',
    'SHMT1','SHMT1', 'rectangle',
    'SHMT2','SHMT2', 'rectangle'
  )
  
  
  
  edges_info_serine = c('PHGDH', 'PSAT1',
                        'PSAT1', 'PSPH',
                        'PSPH','SHMT1',
                        'PSPH','SHMT2')
  
  
  if (pathway_name == 'MAPK'){
    genes_info = genes_info_mapk
    edges_info = edges_info_mapk
  } else if (pathway_name == 'FocalAdh'){
    genes_info = genes_info_focal1
    edges_info = edges_info_focal1
  } else if (pathway_name == 'Neurotrophin'){
    genes_info = genes_info_neuro
    edges_info = edges_info_neuro
  } else if (pathway_name == 'AKT'){
    genes_info = genes_info_akt
    edges_info = edges_info_akt
  } else if (pathway_name == 'Serine'){
    genes_info = genes_info_serine
    edges_info = edges_info_serine
  } else
    stop('Pathway not found')
  
  pathway = list()
  pathway$genes = genes_info
  pathway$edges = edges_info
  
  return(pathway)
  
}












