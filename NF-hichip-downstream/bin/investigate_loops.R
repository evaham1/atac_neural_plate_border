#!/usr/bin/env Rscript

### script to investigate the loops made by HiCDCPlus
print("script to investigate the loops made by HiCDCPlus")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(HiCDCPlus)
library(gridExtra)
library(grid)
library(VennDiagram)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

############################## Set up s"cript options #######################################
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
    
    # for testing with one sample
    plot_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output_investigating/plots/"
    rds_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output_investigating/rds_files/"
    data_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output/" # for HiCDCPlus output of sample 1
    
    # for testing with stringent consensus file
    data_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_diff_interactions/rds_files/" # for consensus interactions across 3 NF samples - all 3
    plot_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_output_investigating/stringent/plots/"
    rds_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_output_investigating/stringent/rds_files/"
    
    # for testing with less stringent consensus file
    data_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_diff_interactions/rds_files/" # for consensus interactions across 3 NF samples - 2 of 3
    plot_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_output_investigating/less_stringent/plots/"
    rds_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_output_investigating/less_stringent/rds_files/"
    
    # for the other data inputs
    data_path = "./output/NF-hichip-downstream/bins/peaks_intersect/" # for filtering by peaks
    data_path = "./output/NF-hichip-downstream/bins/promoters_intersect/" # for filtering by promoters
    data_path = "./output/NF-hichip-downstream/bins/bed_files/" # for bins
    
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


######################################   Read in data   ####################################################

print("reading in data...")

# read in all bins
print("All bins:")
bins <- data.table::fread(paste0(data_path, "bins.bed"))
colnames(bins) <- c("chr", "start", "end", "bin_ID")
head(bins)

# read in peaks intersected with bins
print("Peaks intersected with bins:")
peaks_bins <- data.table::fread(paste0(data_path, "FullData_PeakSet_bins_intersected.bed"))
colnames(peaks_bins) <- c("bin_chr", "bin_start", "bin_end", "bin_ID", "peak_chr", "peak_start", "peak_end", "peak_ID", "dunno", "strand")
head(peaks_bins)
dim(peaks_bins)

# read in promoters intersected with bins
print("Promoters intersected with bins:")
promoters_bins <- data.table::fread(paste0(data_path, "tag_chroms_bins_intersected.bed"))
colnames(promoters_bins) <- c("bin_chr", "bin_start", "bin_end", "bin_ID", "promoter_chr", "promoter_start", "promoter_end", "dunno", "strand", "gene_ID", "gene_name")
head(promoters_bins)
dim(promoters_bins)

# read in HiCDCPlus output
print("HiCDCPlus output:")
# interactions <- data.table::fread(paste0(data_path, "rds_files/NF_HiChip_r1_HiCDC_output_filtered.txt"))
interactions <- data.table::fread(paste0(data_path, "Consensus_interactions.txt"))
interactions <- data.table::fread(paste0(data_path, "Consensus_interactions_stringent.txt"))
interactions <- interactions[,-1]
head(interactions)
dim(interactions)

print("data read in!")

##############  Looking at anchors and how many times they appear    ###########################################

print("Investigating interaction anchor frequencies...")

# how many significant interactions are there
how_many_interactions <- nrow(interactions)
print(paste0("Number of sig interactions: ", how_many_interactions))

# extract all anchors in the significant interactions
anchors_I <- interactions[, 1:3]
colnames(anchors_I) <- c("chr", "start", "end")
anchors_I <- anchors_I %>% mutate(bin_ID = paste0(chr, "-", start+1, "-", end))
anchors_J <- interactions[, 4:6]
colnames(anchors_J) <- c("chr", "start", "end")
anchors_J <- anchors_J %>% mutate(bin_ID = paste0(chr, "-", start+1, "-", end))

interactions <- interactions %>% 
  mutate(anchor_I = anchors_I$bin_ID) %>%
  mutate(anchor_J = anchors_J$bin_ID)

# rbind all these anchors together
all_anchors <- rbind(anchors_I, anchors_J)
print("All anchors:")
print(head(all_anchors))

# frequency at which these anchors appear
print("Frequency at which anchors occur in significant interactions: ")
print(table(table(all_anchors$bin_ID)))

png(paste0(plot_path, "freq_of_anchors_in_interactions.png"), width=60, height=30, units = 'cm', res = 200)
plot(table(table(all_anchors$bin_ID)))
graphics.off()

# pull out the anchors that appear more than 100 times - can vary this threshold
highly_interacting_anchors <- table(all_anchors$bin_ID)[table(all_anchors$bin_ID) > 100]
print("highly interacting anchors: ")
print(highly_interacting_anchors)

##############  Overlap of anchors with any promoter and any peak called across all scATAC data    ###########################################

## how many bins in each type of data
bins_numbers <- data.frame(
  Total_Bins = length(unique(bins$bin_ID)),
  Interacting_Bins = length(unique(all_anchors$bin_ID)),
  Promoter_Bins = length(unique(promoters_bins$bin_ID)),
  Peak_Bins = length(unique(peaks_bins$bin_ID))
)

png(paste0(plot_path, 'bin_counts_table.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(tableGrob(bins_numbers, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# sanity check that all anchors are in the bins
if (sum(unique(all_anchors$bin_ID) %in% bins$bin_ID)){
  print("All anchors are in the bins!")
}else{stop("PROBLEM: NOT ALL ANCHORS ARE IN THE BINS!")}

# Venn diagram of bins between interactions, peaks and genes
venn.diagram(
  x = list(unique(all_anchors$bin_ID), unique(promoters_bins$bin_ID), unique(peaks_bins$bin_ID)),
  category.names = c("Interactions", "Promoters", "Peaks"),
  filename = paste0(plot_path, 'bins_venn_diagram.png'),
  output=TRUE, disable.logging = TRUE,
  # Output features
  imagetype="png",
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
graphics.off()



##############  FILTER INTERACTIONS    ###########################################

## 1) Remove interactions in which the bins involved are NOT annotated as peaks or promoters

# all annotated bins
annotated_bins <- unique(c(peaks_bins$bin_ID, promoters_bins$bin_ID))
length(annotated_bins)

# filter interactions so only include ones where both bins are annotated
filtered_interactions <- interactions %>% filter(anchor_I %in% annotated_bins & anchor_J %in% annotated_bins)
dim(filtered_interactions)
dim(interactions)


## 2) Include only interactions where there is at least ONE promoter and ONE peak on either side of each interactions

# split interactions to 2 dataframes: one where both bins are associated with peaks + promoters, and one where at least one bin is associated to only peak or only promoters
filtered_interactions <- filtered_interactions %>% rownames_to_column('ID')
bins_in_both <- filtered_interactions %>% filter(anchor_I %in% peaks_bins$bin_ID & anchor_I %in% promoters_bins$bin_ID | anchor_J %in% peaks_bins$bin_ID & anchor_J %in% promoters_bins$bin_ID)
bins_in_one <- filtered_interactions[!filtered_interactions$ID %in% bins_in_both$ID]

# in the case where at least one bin is only associated to one, make sure one anochor is associated to peak and one to promoter
filtered_bins_in_one <- bins_in_one[bins_in_one[['anchor_I']] %in% peaks_bins$bin_ID & bins_in_one[['anchor_J']] %in% promoters_bins$bin_ID | 
  bins_in_one[['anchor_J']] %in% peaks_bins$bin_ID & bins_in_one[['anchor_I']] %in% promoters_bins$bin_ID,]

# rbind this filtered df with the bins_in_both df
filtered_interactions <- rbind(filtered_bins_in_one, bins_in_both)

nrow(filtered_interactions)

##############  Read in genes of interest    ###########################################
# from Alex's GMs analysis, GMs made on ss8 data, these are the ones used for coexpression and temporal analysis

# PPR specific GMs
# GM12; AKR1D1, ATP1A1, ATP1B1, ATP2B1, B3GNT7, CCDC25, CGNL1, DAG1, EMILIN2, ENSGALG00000004814, EPCAM, FAM184B, FAM89A, GATA2, GATA3, IRAK2, IRF6, MAP7, METRNL, MPZL3, NET1, PLEKHA5, POGLUT2, PPFIBP1, RGN, SLC16A10, SLC25A4, TSPAN13, TTC39A, UNC5B, Z-TJP2
# GM14; BRINP1, CITED4, DLX5, DLX6, ENSGALG00000023936, ENSGALG00000040010, FN1, HESX1, HIF1A, KCNAB1, LAMB1, LRP11, NFKB1, PAX6, PITX1, PITX2, SEM1, SFRP1, SH3D19, SHISA2, SIX3, SPON1, SST, WDR1, Z-HAPLN1
# GM13; AKAP12, ASS1, BASP1, CD99, ENSGALG00000011296, ENSGALG00000041054, ENSGALG00000042443, EYA2, FERMT2, LGMN, METTL24, NR2F2, NUCKS1, SIX1

PPR_genes <- c("ASS1", "EYA2", "SIX1")

PPR_genes <- c("AKR1D1", "ATP1A1", "ATP1B1", "ATP2B1", "B3GNT7", "CCDC25", "CGNL1", "DAG1", "EMILIN2", "ENSGALG00000004814", "EPCAM", "FAM184B", "FAM89A", "GATA2", "GATA3", "IRAK2", "IRF6", "MAP7", "METRNL", "MPZL3", "NET1", "PLEKHA5", "POGLUT2", "PPFIBP1", "RGN", "SLC16A10", "SLC25A4", "TSPAN13", "TTC39A", "UNC5B", "Z-TJP2",
               "BRINP1", "CITED4", "DLX5", "DLX6", "ENSGALG00000023936", "ENSGALG00000040010", "FN1", "HESX1", "HIF1A", "KCNAB1", "LAMB1", "LRP11", "NFKB1", "PAX6", "PITX1", "PITX2", "SEM1", "SFRP1", "SH3D19", "SHISA2", "SIX3", "SPON1", "SST", "WDR1", "Z-HAPLN1",
               "AKAP12", "ASS1", "BASP1", "CD99", "ENSGALG00000011296", "ENSGALG00000041054", "ENSGALG00000042443", "EYA2", "FERMT2", "LGMN", "METTL24", "NR2F2", "NUCKS1", "SIX1")

# NC specific GMs
# GM40; AGTRAP, BRINP2, CDH11, CMTM8, ENSGALG00000001136, FRZB, GADD45A, HUNK, LARP7, LMX1B, MRAS, MSX1, SFRP2, SOX11, SPSB4, Z-FST, ZEB2, ZIC1
# GM42; BMP5, CDH6, CSRNP1, DRAXIN, EN1, ENSGALG00000013505, FOXD3, NKD1, NRP2, OLFM1, PAX7, SNAI2, SOX9, TFAP2B, TMEM132C, TSPAN18, WLS, WNT6, Z-ENC1, ZFHX4, ZNF423
# GM43; COL9A3, ENSGALG00000037717, ENSGALG00000053185, ERMN, ETS1, LMO4, PPP1R1C, RASL11B, RFTN2, SOX10, SOX8, TNC, Z-MEF2C, Z-PLK2
# GM44; CXCR4, ENSGALG00000030512, ENSGALG00000031427, FABP7, FKBP11, GLIPR2, ID1, ID2, MEOX1, MYL4, OLFML3, SOX5, WNT1

NC_genes <- c("AGTRAP", "BRINP2", "CDH11", "CMTM8", "ENSGALG00000001136", "FRZB", "GADD45A", "HUNK", "LARP7", "LMX1B", "MRAS", "MSX1", "SFRP2", "SOX11", "SPSB4", "Z-FST", "ZEB2", "ZIC1",
              "BMP5", "CDH6", "CSRNP1", "DRAXIN", "EN1", "ENSGALG00000013505", "FOXD3", "NKD1", "NRP2", "OLFM1", "PAX7", "SNAI2", "SOX9", "TFAP2B", "TMEM132C", "TSPAN18", "WLS", "WNT6", "Z-ENC1", "ZFHX4", "ZNF423",
              "CDON", "COTL1", "ENSGALG00000048488", "MFAP2", "PRTG", "TUBB3",
              "COL9A3", "ENSGALG00000037717", "ENSGALG00000053185", "ERMN", "ETS1", "LMO4", "PPP1R1C", "RASL11B", "RFTN2", "SOX10", "SOX8", "TNC", "Z-MEF2C", "Z-PLK2",
              "CXCR4", "ENSGALG00000030512", "ENSGALG00000031427", "FABP7", "FKBP11", "GLIPR2", "ID1", "ID2", "MEOX1", "MYL4", "OLFML3", "SOX5", "WNT1")

##############  Read in peaks of interest    ###########################################

## PPR specific PMs: FullData_PM6, FullData_PM7, FullData_PM9
PPR_peaks <- c(
  "chr1-51885355-51885855", "chr1-56743328-56743828", "chr15-2282477-2282977", "chr17-8061757-8062257", "chr2-34364392-34364892", "chr2-72766649-72767149", "chr2-91522834-91523334", "chr25-3194857-3195357", "chr28-2880144-2880644", "chr3-100929473-100929973", "chr3-101800822-101801322", "chr3-18826788-18827288", "chr3-3153527-3154027", "chr3-37532963-37533463", "chr3-40860007-40860507", "chr3-42110813-42111313", "chr3-46021620-46022120", "chr3-9604325-9604825", "chr33-4973016-4973516", "chr4-15228534-15229034", "chr4-3545860-3546360", "chr4-39659304-39659804", "chr4-42955909-42956409", "chr4-61691248-61691748", "chr4-82791681-82792181", "chr5-13025286-13025786", "chr5-14360597-14361097", "chr5-2185340-2185840", "chr5-25039278-25039778", "chr5-26331436-26331936", "chr5-37774715-37775215", "chr6-13570150-13570650", "chr6-22584965-22585465", "chr6-3635106-3635606", "chr7-14213427-14213927", "chr7-26273335-26273835", "chr8-20658815-20659315", "chr8-20919981-20920481", "chr8-8130129-8130629", "chr9-1192259-1192759", "chr9-12061796-12062296",
  "chr1-105919837-105920337", "chr1-113306753-113307253", "chr1-15015808-15016308", "chr1-156920062-156920562", "chr1-157299139-157299639", "chr1-15892762-15893262", "chr1-15919402-15919902", "chr1-168641905-168642405", "chr1-168893830-168894330", "chr1-176740306-176740806", "chr1-194840235-194840735", "chr1-20566555-20567055", "chr1-33571364-33571864", "chr1-60416481-60416981", "chr1-62211336-62211836", "chr1-72723164-72723664", "chr1-8065715-8066215", "chr1-88390949-88391449", "chr1-93486724-93487224", "chr10-14491785-14492285", "chr10-14492580-14493080", "chr10-20340332-20340832", "chr10-20452537-20453037", "chr11-19273365-19273865", "chr11-1987100-1987600", "chr11-2257660-2258160", "chr12-18975182-18975682", "chr12-3599217-3599717", "chr12-6958745-6959245", "chr12-7090763-7091263", "chr12-7158977-7159477", "chr12-7163182-7163682", "chr12-7163824-7164324", "chr12-7292027-7292527", "chr12-8524299-8524799", "chr12-8995625-8996125", "chr13-12548521-12549021", "chr13-16386526-16387026", "chr13-1671556-1672056", "chr13-8986883-8987383", "chr13-9166470-9166970", "chr14-13306078-13306578", "chr14-5259731-5260231", "chr14-5273840-5274340", "chr14-7010284-7010784", "chr15-6045070-6045570", "chr15-767877-768377", "chr17-5823220-5823720", "chr18-3734720-3735220", "chr18-836196-836696", "chr18-9597171-9597671", "chr19-8931349-8931849", "chr2-111543253-111543753", "chr2-113095334-113095834", "chr2-126210787-126211287", "chr2-130112541-130113041", "chr2-147786771-147787271", "chr2-23194456-23194956", "chr2-24093951-24094451", "chr2-37262664-37263164", "chr2-40163142-40163642", "chr2-40666551-40667051", "chr2-40777614-40778114", "chr2-41076879-41077379", "chr2-41096965-41097465", "chr2-46579772-46580272", "chr2-46796604-46797104", "chr2-60768678-60769178", "chr2-90151761-90152261", "chr2-98858224-98858724", "chr20-10297374-10297874", "chr20-3417599-3418099", "chr20-5668390-5668890", "chr20-9882128-9882628", "chr20-9895339-9895839", "chr21-394147-394647", "chr23-2270591-2271091", "chr23-328199-328699", "chr23-4463337-4463837", "chr23-639645-640145", "chr25-2726380-2726880", "chr26-4786832-4787332", "chr26-4822131-4822631", "chr27-4828979-4829479", "chr27-6164738-6165238", "chr27-7146919-7147419", "chr28-2142363-2142863", "chr28-4636074-4636574", "chr3-100893435-100893935", "chr3-102692943-102693443", "chr3-103753008-103753508", "chr3-105402084-105402584", "chr3-108232358-108232858", "chr3-108928367-108928867", "chr3-13375460-13375960", "chr3-1449127-1449627", "chr3-16699856-16700356", "chr3-22710909-22711409", "chr3-32309110-32309610", "chr3-32517550-32518050", "chr3-39074445-39074945", "chr3-39370694-39371194", "chr3-64337111-64337611", "chr3-67485272-67485772", "chr3-71754992-71755492", "chr3-79368221-79368721", "chr3-8650978-8651478", "chr3-87433943-87434443", "chr3-97408796-97409296", "chr33-4665299-4665799", "chr4-1390646-1391146", "chr4-17556200-17556700", "chr4-17640299-17640799", "chr4-240453-240953", "chr4-33759038-33759538", "chr4-34644939-34645439", "chr4-354006-354506", "chr4-43388891-43389391", "chr4-83638368-83638868", "chr4-9476451-9476951", "chr5-1148018-1148518", "chr5-17414692-17415192", "chr5-24251959-24252459", "chr5-3948608-3949108", "chr5-40035538-40036038", "chr7-1009852-1010352", "chr7-33112267-33112767", "chr8-26298399-26298899", "chr9-18127591-18128091",
  "chr1-115742625-115743125", "chr1-118138264-118138764", "chr1-124572961-124573461", "chr1-131433117-131433617", "chr1-140680906-140681406", "chr1-141276478-141276978", "chr1-142414827-142415327", "chr1-1548751-1549251", "chr1-16017632-16018132", "chr1-16510397-16510897", "chr1-166577443-166577943", "chr1-169475933-169476433", "chr1-18484334-18484834", "chr1-186493813-186494313", "chr1-194742462-194742962", "chr1-20918719-20919219", "chr1-24222888-24223388", "chr1-25525043-25525543", "chr1-29787549-29788049", "chr1-29823089-29823589", "chr1-30789661-30790161", "chr1-34361501-34362001", "chr1-34362002-34362502", "chr1-38326892-38327392", "chr1-38366243-38366743", "chr1-39551257-39551757", "chr1-42364379-42364879", "chr1-43485243-43485743", "chr1-46478536-46479036", "chr1-60970364-60970864", "chr1-61696787-61697287", "chr1-85099460-85099960", "chr1-86749742-86750242", "chr1-8841329-8841829", "chr1-96390837-96391337", "chr1-98741966-98742466", "chr10-10272404-10272904", "chr10-12849026-12849526", "chr10-5692798-5693298", "chr10-7462858-7463358", "chr11-15521598-15522098", "chr11-15529701-15530201", "chr11-388210-388710", "chr12-11534333-11534833", "chr12-14776981-14777481", "chr12-17501359-17501859", "chr12-18364771-18365271", "chr12-18943301-18943801", "chr13-12527823-12528323", "chr13-13545910-13546410", "chr13-18191067-18191567", "chr13-18442010-18442510", "chr13-9207099-9207599", "chr14-5275234-5275734", "chr17-5108318-5108818", "chr17-8890725-8891225", "chr18-4077933-4078433", "chr19-9685399-9685899", "chr2-149081179-149081679", "chr2-19630486-19630986", "chr2-39091436-39091936", "chr2-50593634-50594134", "chr2-5316227-5316727", "chr2-63741706-63742206", "chr2-63931517-63932017", "chr20-1015766-1016266", "chr20-11539936-11540436", "chr20-12379867-12380367", "chr20-13165061-13165561", "chr21-1003730-1004230", "chr21-1362953-1363453", "chr21-2402710-2403210", "chr21-2705825-2706325", "chr22-1883318-1883818", "chr23-2249285-2249785", "chr23-4880807-4881307", "chr23-5189129-5189629", "chr26-1704254-1704754", "chr26-618206-618706", "chr28-416980-417480", "chr3-11684665-11685165", "chr3-17132759-17133259", "chr3-353593-354093", "chr3-48440133-48440633", "chr3-65064074-65064574", "chr4-23449835-23450335", "chr4-49589263-49589763", "chr4-50663077-50663577", "chr4-65285432-65285932", "chr4-77824880-77825380", "chr7-8114384-8114884"
)



##############  Look at which interactions include PPR genes    ###########################################


## should write a function which:
# inputs 1 or more features (genes or peaks)
# extracts the bins which overlaps with the input feature(s)
# extracts the interactions which overlap with the bins
# extracts the OTHER bins which these interactions connect with (ie not input bin)
# returns a table like this:
# input feature | input bin | interaction | other bin | other bin feature | other bin feature type


########################
## function to extract bins from features (gene_names, gene_IDs or peak_IDs)
extract_bins_from_features <- function(features, feature_type, bin_dictionary){
  if (!feature_type %in% c("gene_name", "gene_ID", "peak_ID")) {
    stop("feature_type must be 'gene_name', 'gene_ID' or 'peak_ID'!")
  } else {
    bins <- bin_dictionary %>% 
      dplyr::filter(get(feature_type) %in% features)
  }
  if (feature_type == "peak_ID"){
    bins <- bins %>% dplyr::select(peak_ID, bin_ID)
  } else {
    bins <- bins %>% dplyr::select(gene_ID, gene_name, bin_ID)
  }
  return(bins)
}
#############

##############
##### reads in input bins and finds interactions
## function to extract all bins that interact with input bins
# removes bins which are not interacting significantly
# annotated anchors of interactions with the information from the input bins
## input_bins is bins to search for, input_feature_type is the type of annotation of input bins
## interactions is the dictionary of interactions to search
extract_interacting_bins <- function(input_bins, input_feature_type, interactions){
  
  if (!input_feature_type %in% c("gene", "peak")) {
    stop("input_feature_type must be 'gene', 'peak'!")
  }
  
  # extract interactions with the input bins
  extracted_interactions <- interactions %>% 
    filter(anchor_I %in% input_bins$bin_ID | anchor_J %in% input_bins$bin_ID) %>%
    dplyr::select(anchor_I, anchor_J)
  
  # merge the input bins annotations with the anchors
  if(input_feature_type == "gene"){
    colnames(input_bins) <- c("anchor_I_gene_ID", "anchor_I_gene_name", "anchor_I")
    output_interactions <- merge(extracted_interactions, input_bins, by = "anchor_I", all.x = TRUE)
    colnames(input_bins) <- c("anchor_J_gene_ID", "anchor_J_gene_name", "anchor_J")
    output_interactions <- merge(output_interactions, input_bins, by = "anchor_J", all.x = TRUE)
  }
  if(input_feature_type == "peak"){
    colnames(input_bins) <- c("anchor_I_peak_ID", "anchor_I")
    output_interactions <- merge(extracted_interactions, input_bins, by = "anchor_I", all.x = TRUE)
    colnames(input_bins) <- c("anchor_J_peak_ID", "anchor_J")
    output_interactions <- merge(output_interactions, input_bins, by = "anchor_J", all.x = TRUE)
  }
  
  # identify whether the input anchors are in anchor I, anchor J or both
  output_interactions <- output_interactions %>%
    mutate(input_anchor = case_when(
      !is.na(anchor_I_peak_ID) & !is.na(anchor_J_peak_ID) ~ "both",
      !is.na(anchor_I_peak_ID) ~ "anchor_I",
      !is.na(anchor_J_peak_ID) ~ "anchor_J",
      TRUE ~ NA_character_
    ))
  # then swap around the anchors (they dont mean anything anyway) so all input anchors are anchor I
  df <- output_interactions %>%
    mutate(anchor_I_peak_ID = case_when(input_anchor == "anchor_J" ~ anchor_J_peak_ID, TRUE ~ anchor_I_peak_ID),
      anchor_J_peak_ID = case_when(input_anchor == "anchor_J" ~ anchor_I_peak_ID, TRUE ~ anchor_J_peak_ID),
      anchor_I = case_when(input_anchor == "anchor_J" ~ anchor_J, TRUE ~ anchor_I),
      anchor_J = case_when(input_anchor == "anchor_J" ~ anchor_I, TRUE ~ anchor_J)
    )
  
  
  return(output_interactions)
}
##############

########################
##### reads in interactions from gene bins and annotates to promoters
## function to extract features from binIDs
# output is collapsed so each row = binID, annotated with all interacting peaks and promoters
extract_peaks_from_interactions <- function(output_bins, peak_bin_dictionary){
  # extract all peaks in these bins and collapse
  peaks <- peak_bin_dictionary %>% 
    dplyr::filter(peak_bin_dictionary$bin_ID %in% output_bins$anchor_I | peak_bin_dictionary$bin_ID %in% output_bins$anchor_J) %>%
    dplyr::select(peak_ID, bin_ID) %>%
    group_by(bin_ID) %>%
    summarise(across(everything(), ~paste0(unique(.), collapse = ",")))
  # add peak annotations to the anchors in the interactions
  colnames(peaks) <- c("anchor_I", "anchor_I_peak_ID")
  output_interactions <- merge(output_bins, peaks, by = "anchor_I", all.x = TRUE)
  colnames(peaks) <- c("anchor_J", "anchor_J_peak_ID")
  output_interactions <- merge(output_interactions, peaks, by = "anchor_J", all.x = TRUE)
  return(output_interactions)
}

########################
## function to extract features from binIDs
# output is collapsed so each row = binID, annotated with all interacting peaks and promoters
extract_features_from_bins <- function(output_bins, peak_bin_dictionary, promoter_bin_dictionary){
  # extract all peaks in these bins and collapse
  peaks <- peak_bin_dictionary %>% 
    dplyr::filter(peak_bin_dictionary$bin_ID %in% output_bins$anchor_I | peak_bin_dictionary$bin_ID %in% output_bins$anchor_J) %>%
    dplyr::select(peak_ID, bin_ID) %>%
    group_by(bin_ID) %>%
    summarise(across(everything(), ~paste0(unique(.), collapse = ",")))
  # extract all gene promoters in these bins and collapse
  promoters <- promoter_bin_dictionary %>%
    dplyr::filter(promoter_bin_dictionary$bin_ID %in% output_bins$anchor_I | promoter_bin_dictionary$bin_ID %in% output_bins$anchor_J) %>%
    dplyr::select(gene_ID, gene_name, bin_ID) %>%
    group_by(bin_ID) %>%
    summarise(across(everything(), ~paste0(unique(.), collapse = ",")))
  # add peak annotations to the anchors in the interactions
  colnames(peaks) <- c("anchor_I", "anchor_I_peak_ID")
  output_interactions <- merge(output_bins, peaks, by = "anchor_I", all.x = TRUE)
  colnames(peaks) <- c("anchor_J", "anchor_J_peak_ID")
  output_interactions <- merge(output_interactions, peaks, by = "anchor_J", all.x = TRUE)
  # add gene annotations to the anchors in the interactions
  colnames(promoters) <- c("anchor_I", "anchor_I_gene_ID", "anchor_I_gene_name")
  output_interactions <- merge(output_interactions, promoters, by = "anchor_I", all.x = TRUE)
  colnames(promoters) <- c("anchor_J", "anchor_J_gene_ID", "anchor_J_gene_name")
  output_interactions <- merge(output_interactions, promoters, by = "anchor_J", all.x = TRUE)
  return(output_interactions)
}



####### Running through processes from genes to peaks:

# 1) Extract input_bins from PPR genes
input_bins <- extract_bins_from_features(features=PPR_genes, feature_type="gene_name", bin_dictionary=promoters_bins)

# 2) See which output_bins are interacting with these input_bins
output_bins <- extract_interacting_bins(input_bins = input_bins, interactions = filtered_interactions)

# 3) See which features are in these output_bins
output_features <- extract_peaks_from_interactions(output_bins, peaks_bins)

# 4) See which of the two anchors includes the input genes
output_features <- output_features %>% mutate(input_gene_id = coalesce(anchor_I_gene_ID, anchor_J_gene_ID)) %>%
  mutate(input_gene_name = coalesce(anchor_I_gene_name, anchor_J_gene_name)) %>% 
  select(anchor_I, anchor_J, input_gene_id, input_gene_name, anchor_I_peak_ID, anchor_J_peak_ID)

# Extract all peaks to plot in heatmap
peaks <- peaks_bins %>% dplyr::filter(bin_ID %in% output_bins$anchor_I | bin_ID %in% output_bins$anchor_J)
length(peaks$peak_ID)
write.csv(peaks$peak_ID, "./output/Rshiny_PPR_input_peaks.csv")

# 1) Extract input_bins from NC genes
input_bins <- extract_bins_from_features(features=NC_genes, feature_type="gene_name", bin_dictionary=promoters_bins)

# 2) See which output_bins are interacting with these input_bins (removes input bins)
output_bins <- extract_interacting_bins(input_bins = input_bins, interactions = interactions)

# 3) See which features are in these output_bins
output_features <- extract_features_from_bins(output_bins, peaks_bins, promoters_bins)
peaks <- peaks_bins %>% dplyr::filter(bin_ID %in% output_bins)
length(peaks$peak_ID)
write.csv(peaks$peak_ID, "./output/Rshiny_NC_input_peaks.csv")



# for each interaction, how many interactions are between input_bins and how many are input_bin to output_bin
investigate_interactions <- selected_interactions %>% 
  mutate(anchor_I_in_input = anchor_I %in% input_bins$bin_ID) %>% 
  mutate(anchor_J_in_input = anchor_J %in% input_bins$bin_ID) %>% 
  mutate(Number_input_bins = rowSums(across(where(is.logical))))
table(investigate_interactions$Number_input_bins)



####### Running through processes from peaks to genes:

# 1) Extract input_bins from PPR peaks
input_bins <- extract_bins_from_features(features=PPR_peaks, feature_type="peak_ID", bin_dictionary=peaks_bins)

# 2) See which output_bins are interacting with these input_bins
output_bins <- extract_interacting_bins(input_bins = input_bins, input_feature_type = "peak", interactions = filtered_interactions)

# 3) See which features are in these output_bins
output_features <- extract_peaks_from_interactions(output_bins, peaks_bins)

