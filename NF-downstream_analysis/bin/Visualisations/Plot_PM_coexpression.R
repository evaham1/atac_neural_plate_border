#!/usr/bin/env Rscript

print("Plots PM feature and co-accessibility plots")

############################## Load libraries #######################################
library(optparse)
library(future)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)
library(data.table)
library(Seurat)
library(patchwork)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = FALSE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    # data paths for the different inputs
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/4_PM_GAMs/FullData/" # SEACells metadata + PM averages
    data_path = "./output/NF-downstream_analysis/Processing/ss4/SEACELLS_INTEGRATING_WF/Integrated_SEACells_label_transfer/rds_files/" # latent time on metacells metadata 
    # output paths:
    rds_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/5_PM_FeaturePlots/ss4/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/5_PM_FeaturePlots/ss4/plots/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


########################       FUNCTIONS    ########################################

## adapted from Alex's coexpression UMAP plot
plot_umap_pm_coaccessibility <- function(seurat_object, PM_avgs, pm_1, pm_2, col.threshold = 0, two.colors = c("red", "blue"), negative.color = 'gray80', limit = 0, 
                                        module_names = c('Gene module 1', 'Gene module 2'), show_legend = TRUE,
                                         axes_label_size = 12, axes_title_size = 10, axes_tick_size = 0.15){
  
  # set params
  start = 1
  end = 100
  width = end - start
  
  # extract PM average scores and scale them
  pm_1_scaled <- (PM_avgs[,pm_1] - min(PM_avgs[,pm_1]))/(max(PM_avgs[,pm_1]) - min(PM_avgs[,pm_1])) * width + start
  pm_2_scaled <- (PM_avgs[,pm_2] - min(PM_avgs[,pm_2]))/(max(PM_avgs[,pm_2]) - min(PM_avgs[,pm_2])) * width + start
  
  # extract colour matrix
  dat <- data.frame(pm_1_scaled, pm_2_scaled, row.names = rownames(PM_avgs))
  dat <-  round(dat, 0)
  col_mat = Seurat:::BlendMatrix(n = 100, col.threshold = col.threshold, two.colors =  two.colors, negative.color = negative.color)
  col_mat <- as.data.frame.table(col_mat, responseName = "value") %>% mutate_if(is.factor, as.integer)
  col_mat[!(col_mat$Var1 > limit*100 & col_mat$Var2 > limit*100), 'value'] <- negative.color
  colnames(col_mat) <- c('a', 'b', 'mix')
  
  # assign colours to cells
  cell_cols <- unlist(apply(dat, 1, function(x){filter(col_mat, a == x[[1]] & b == x[[2]])[[3]]}))
  col_mat[,1:2] <- col_mat[,1:2]/100
  key_plot <- ggplot(col_mat %>% filter(mix != !!negative.color), aes(x = a, y = b)) +
    xlab(module_names[1]) +
    ylab(module_names[2]) +
    geom_tile(aes(fill = mix)) +
    scale_fill_identity() +
    scale_x_continuous(breaks = c(limit, 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = c(limit, 1), expand = c(0.01, 0.01)) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(size = axes_title_size),
          axis.text.x = element_text(size = axes_label_size),
          axis.title.y = element_text(size = axes_title_size),
          axis.text.y = element_text(size = axes_label_size),
          axis.ticks.length=unit(axes_tick_size,"cm"))
  
  plot_data <- as.data.frame(seurat_object[["umap"]]@cell.embeddings)
  
  # add cell colours to plot_data
  plot_data <- merge(plot_data, as.data.frame(cell_cols), by=0) %>% column_to_rownames('Row.names')
  
  umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, colour = cell_cols)) +
    geom_point(colour = negative.color, size = 6) +
    geom_point(data = plot_data %>% filter(cell_cols != negative.color), size = 6) +
    scale_colour_manual(values=plot_data %>% filter(cell_cols != negative.color) %>% dplyr::pull(cell_cols))+
    theme_void() +
    NoLegend()
  
  if (show_legend == FALSE){
    print(umap_plot)
  } else {
    layout <- '
    BA
    B#
    '
    wrap_plots(A = key_plot, B = umap_plot, design = layout, widths = c(4,1), heights = c(1,3))
  }
}

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'Non-neural',
                              'PPR', 'aPPR', 'pPPR', 'Placodal',
                              'eNPB', 'NPB', 'aNPB', 'pNPB',
                              'NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'Neural',
                              'node', 'streak', 'PGC', 'BI', 'meso', 'endo', 'Contam',
                              'MIXED', 'Unmapped')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "#9792A3",
                                "#7C8483", "#EAEAEA")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam',
                                       'MIXED', 'Unmapped')
########################       STAGE COLOURS     ###########################################
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
names(stage_colours) <- stage_order
############################################################################################
lineage_colours = c('placodal' = '#3F918C', 'NC' = '#DE4D00', 'neural' = '#8000FF')

########################################################################################################
#                                 Read in data and clean up                               #
########################################################################################################

########## PEAK MODULE AVERAGE SCORES DF #############
PM_avg_scores <- read_csv(paste0(data_path, "rds_files/PM_avg_scores.csv"))
PM_avg_scores <- column_to_rownames(PM_avg_scores, "SEACell_ID")
print("PM average scores: ")
print(head(PM_avg_scores))

######### SEACell METADATA #############
metadata <- read_csv(paste0(data_path, "rds_files/Combined_SEACell_integrated_metadata_latent_time.csv"))
print("SEACell metadata:")
print(head(metadata))

########## ATAC METACELL SEURAT OBJECT ############# 
label <- setdiff(sub('_.*', '', list.files(data_path)), "seacells_seurat_integrated.RDS")
label <- label[label != 'rds']
print(label)
seurat <- readRDS(paste0(data_path, label, "_seacells_seurat_integrated.RDS"))

print("Data read in!")

########################################################################################################
#                         Add latent time + lineage probs to seurat object                             #
########################################################################################################

print("Adding latent time and lineage probs to seurat object...")

print(head(metadata))

## extract only the correct stage metacells
metadata <- metadata[grepl(unique(seurat@meta.data$stage), metadata$Rownames), ]

## reorg metadata so in order of cells as appear in the seurat object
metadata <- metadata %>%
  mutate(Rownames = gsub("_", "-", Rownames)) %>%
  mutate(Rownames = substr(Rownames, 1, nchar(Rownames) - 4))
metadata <- metadata %>% arrange(Rownames)
head(metadata)

## check that the seacells match and are in the same order
if (sum(rownames(seurat@meta.data) == metadata$Rownames) == nrow(seurat@meta.data)){
  "SEACell IDs match!"
} else {stop("ERROR! SEACell IDs dont match!")}

## add latent time and lineage probs as metadata to seurat object
seurat@meta.data[["rna_latent_time"]] <- metadata$rna_latent_time
seurat@meta.data[["rna_lineage_neural_probability"]] <- metadata$rna_lineage_neural_probability
seurat@meta.data[["rna_lineage_NC_probability"]] <- metadata$rna_lineage_NC_probability
seurat@meta.data[["rna_lineage_placodal_probability"]] <- metadata$rna_lineage_placodal_probability

head(seurat@meta.data)

## plot feature plots of these metadata
png(paste0(plot_path, 'rna_latent_time_feature_plot.png'), width = 15, height = 15, units='cm', res=200)
print(
  FeaturePlot(seurat, features = "rna_latent_time", pt.size = 6) +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.text = element_text(size=16),
          legend.key.size = unit(1, 'cm'))
) 
graphics.off()

png(paste0(plot_path, 'rna_lineage_neural_probability_feature_plot.png'), width = 15, height = 15, units='cm', res=200)
print(
  FeaturePlot(seurat, features = "rna_lineage_neural_probability", pt.size = 6) +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.text = element_text(size=16),
          legend.key.size = unit(1, 'cm'))
) 
graphics.off()

png(paste0(plot_path, 'rna_lineage_NC_probability_feature_plot.png'), width = 15, height = 15, units='cm', res=200)
print(
  FeaturePlot(seurat, features = "rna_lineage_NC_probability", pt.size = 2) +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.text = element_text(size=16),
          legend.key.size = unit(1, 'cm'))
) 
graphics.off()

png(paste0(plot_path, 'rna_lineage_placodal_probability_feature_plot.png'), width = 15, height = 15, units='cm', res=200)
print(
  FeaturePlot(seurat, features = "rna_lineage_placodal_probability", pt.size = 6) +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.text = element_text(size=16),
          legend.key.size = unit(1, 'cm'))
) 
graphics.off()

########################################################################################################
#                                 Add PMs as features to seurat object                               #
########################################################################################################

print("Adding average PM accessibility to seurat object...")

plot_path = "./plots/feature_plots/"
dir.create(plot_path, recursive = T)

## split up the PM df by stage
split_df <- split(PM_avg_scores, substring(rownames(PM_avg_scores), nchar(rownames(PM_avg_scores)) - 3))
names(split_df) <- substr(names(split_df), 2, 4)

## detect which stage is the seurat object and extract the correct df
df <- split_df[[unique(seurat@meta.data$stage)]]
head(df)

## reorder df so the metacells are in the same order as they appear in the seurat metadata
df <- rownames_to_column(df, "SEACell_ID")
df <- df %>%
  mutate(SEACell_ID = gsub("_", "-", SEACell_ID)) %>%
  mutate(SEACell_ID = substr(SEACell_ID, 1, nchar(SEACell_ID) - 4))
df <- df %>% arrange(SEACell_ID)
head(df)

## check that the seacells match and are in the same order
if (sum(rownames(seurat@meta.data) == df$SEACell_ID) == nrow(seurat@meta.data)){
  "SEACell IDs match!"
} else {stop("ERROR! SEACell IDs dont match!")}

## add each average PM score as metadata
SEACell_IDs <- df[,1]
df <- df[,-1]
for (module in colnames(df)){
  seurat@meta.data[[module]] <- df[[module]]
  
  # feature plot of that PM
  png(paste0(plot_path, module, '_feature_plot.png'), width = 15, height = 15, units='cm', res=200)
  print(
    FeaturePlot(seurat, features = module, pt.size = 6) +
      theme_void() +
      theme(plot.title = element_blank(),
            legend.text = element_text(size=16),
            legend.key.size = unit(1, 'cm'))
  ) 
  graphics.off()
}


########################################################################################################
#                                 Plot co-accessibility of PMs                       #
########################################################################################################

print("Plotting co-accessibility...")

rownames(df) <- SEACell_IDs
head(df)

plot_path = "./plots/coaccessibility_plots/limit_0.4/"
dir.create(plot_path, recursive = T)

limit = 0.4

# neural/NC with placodal PM1
PMA = "FullData_PM1"
PMB = "FullData_PM6"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                                     module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                                     limit = limit)
graphics.off()

PMB = "FullData_PM7"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM10"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM11"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM12"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMA = "FullData_PM13"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM14"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM15"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()


# neural/NC with placodal PM2
PMA = "FullData_PM2"
PMB = "FullData_PM6"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM7"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM10"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM11"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM12"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM13"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM14"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM15"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

# neural/NC with placodal PM3
PMA = "FullData_PM3"
PMB = "FullData_PM6"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM7"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM10"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM11"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM12"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM13"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM14"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM15"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

# neural/NC with placodal PM4
PMA = "FullData_PM4"
PMB = "FullData_PM6"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM7"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM10"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM11"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM12"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM13"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM14"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()

PMB = "FullData_PM15"
png(paste0(plot_path, 'Coaccessibility_plot_', substr(PMA, 10, 14), "-", substr(PMB, 10, 14), '.png'), width = 15, height = 15, units='cm', res=200)
plot_umap_pm_coaccessibility(seurat, df, PMA, PMB,  
                             module_names = c(substr(PMA, 10, 14), substr(PMB, 10, 14)),
                             limit = limit)
graphics.off()