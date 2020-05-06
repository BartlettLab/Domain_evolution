#200426 had to reinstall all packages in new lib loc.. something broke
#orig lib "/home/jm33a/R/x86_64-pc-linux-gnu-library/3.6"
.libPaths("/home/jm33a/R_2")
.libPaths()

#BiocManager::install("Biostrings")
require(Biostrings)

install.packages("tidyverse", lib = "/home/jm33a/R_2")
require(tidyverse)

#BiocManager::install("msa")
require(msa)

devtools::install_github("hrbrmstr/statebins")
require(statebins)

#install.packages("zoo")
require(zoo)

#install.packages("cowplot")
require(cowplot)

#install.packages("ggsignif")
require(ggsignif)

source("~/bin/R_rainclouds.R")

################################################
#Set up initial stuff

    setwd("/home/jm33a/domain_evolution/clade_trees/GSO/tree_from_hits/follow_up_analyses/site_rates/")
    anchor_gene <-"AT5G44700.1_GSO2"
    burnin <- 0.1
    lrr_range <- 94:839 
    tm_range <- 877:897
    rlk_range <- 948:1232
   
    setwd("/home/jm33a/domain_evolution/clade_trees/PEPR/tree_from_hits/follow_up_analyses/site_rates/")
    anchor_gene <-"AT1G73080.1_PEPR1"
    burnin <- 0.2
    lrr_range <- 31:721 
    tm_range <- 770:790
    rlk_range <- 827:1115
   
    setwd("/home/jm33a/domain_evolution/clade_trees/CLV1/tree_from_hits/follow_up_analyses/site_rates/")
    anchor_gene <-"AT1G75820.1_CLV1"
    burnin <- 0.6
    lrr_range <- 93:600 
    tm_range <- 639:659
    rlk_range <- 685:969
    
    setwd("/home/jm33a/domain_evolution/clade_trees/RGI5/tree_from_hits/follow_up_analyses/site_rates/")
    anchor_gene <-"AT1G34110.1_RGI5"
    burnin <- 0.1
    lrr_range <- 90:677 
    tm_range <- 707:727
    rlk_range <- 772:1067    
    
#    setwd("/home/jm33a/domain_evolution/clade_trees/PXY/tree_from_hits/follow_up_analyses/site_rates/")
#    anchor_gene <-"AT5G61480.1_PXY"
#    burnin <- 0.5
#    lrr_range <- 80:607 
#    tm_range <- 653:673
#    rlk_range <- 719:1001    
    
    setwd("/home/jm33a/domain_evolution/clade_trees/IKU2/tree_from_hits/follow_up_analyses/site_rates/")
    anchor_gene <-"AT3G19700.1_IKU2"
    burnin <- 0.1
    lrr_range <- 66:578 
    tm_range <- 617:637
    rlk_range <- 671:970    
    
    
    
    setwd("/home/jm33a/domain_evolution/clade_trees/FLS2/tree_from_hits/follow_up_analyses/site_rates/")
    anchor_gene <-"AT5G46330.1_FLS2"
    burnin <- 0.1
    lrr_range <- 97:769 
    tm_range <- 807:827
    rlk_range <- 870:1155      
    
    
    
####   IQtree site rate tests     #####    
  alignment_site_rates_table <- read.table(list.files(pattern = "*align.fasta.rate"),header=TRUE)  [, c("Site", "Rate")] #read in the site rates for the alignment, stripping categories data
  colnames(alignment_site_rates_table) <- c("MSA_position", "Rate")
  peptide_MSA <-  as(readAAMultipleAlignment(filepath = list.files(pattern = "*align.fasta$"), format = "fasta"), "BStringSet") #bring in the entire peptide MSA as biostrings set

#Build the site rates table based on the target gene
  peptide_anchor_sequence <- as.character(peptide_MSA[anchor_gene]) #get the sequence of gene in alignment, including dashes
  alignment_site_rates_table$peptide_anchor_sequence <- unlist(strsplit(peptide_anchor_sequence, "")) #show target sequence residues on site rates table
  anchor_gene_site_rates_table <- alignment_site_rates_table[alignment_site_rates_table$peptide_anchor_sequence != "-",] #rates from alignment table stripped of gaps in target gene
  anchor_gene_site_rates_table <- anchor_gene_site_rates_table %>% mutate(Gene_position = row_number()) #add the reference gene's residue numbering to table
   
#calculate sliding mean of site rates and add to site rate table
  sliding_window_size <- 50
  sliding_window_table <- anchor_gene_site_rates_table %>% 
    select(Gene_position, Rate) %>% 
    mutate(sliding_window = rollmean(Rate, k = sliding_window_size, fill = NA)) #calculate sliding mean rate
  anchor_gene_site_rates_table$sliding_window_rate <- sliding_window_table$sliding_window #add column with calculated sliding mean rate
  
  #transform table for visualization and plot
  vis_tab <- anchor_gene_site_rates_table[, c("Gene_position", "Rate", "sliding_window_rate") ] %>% gather(metric, evolution_rate, c(Rate, sliding_window_rate))

  
####   Selection tests ####


  
  
  setwd("../nucleotide_seqs/") #move to the same clade's nucleotide analysis location
   CDS_MSA <-  as(readDNAMultipleAlignment(filepath = list.files(pattern = "*translation.fasta"), format = "fasta"), "BStringSet") #bring in the entire CDS MSA as biostrings set
  CDS_anchor_sequence <- as.character(CDS_MSA[anchor_gene]) #get the sequence of gene in alignment, including dashes
  
    #### FEL test####
  fel_neg_pval <- 0.001
  fel_pos_pval <- 0.3
  
     fel_tbl <- read.csv(list.files(pattern = "FEL*"),header=TRUE)   #import FEL results output file 
    fel_tbl <- fel_tbl %>% #identify sites under selection according to manufacturer's instructions
    mutate(Selection_detected = case_when( 
                                  ( alpha > beta ) & ( p.value < fel_neg_pval ) ~ 'Negative',
                                  ( alpha < beta ) & ( p.value < fel_pos_pval ) ~ 'Positive',
                                  TRUE ~ '')) 
    sum(fel_tbl$Selection_detected == 'Positive') #check that results agree with web interface results
    sum(fel_tbl$Selection_detected == 'Negative')#check that results agree with web interface results
 
 #strip positions of FEL table MSA that are not present in anchor sequence
    fel_tbl$anchor_gene_seq <- unlist(strsplit(CDS_anchor_sequence, ""))[1:(str_length(CDS_anchor_sequence)/3) * 3] #print every third character of anchor gene in column, this is to check if anchor gene has sequence or gap at each position

 
####section to move selection hits to nearest anchor gene residue, probably don't use####    
##if selection detected was on position without anchor gene sequence, move it to nearest available anchor gene position
##this approach may be problematic because it is taking data from positions less likely to be represented by all sequences, ie. an insertion only found in one or a few genes. Therefore it is highly likely to have extreme signal. And it stacks up values, which for neg heatmap style overrepresents gappy locations.
#   fel_tbl$anchor_gene_residue_number <- NA
#   anchor_gene_position <- 1
#           for(msa_position in 1:nrow(fel_tbl)){ #add residue counting based on anchor gene sequence
#      if(fel_tbl[msa_position, 'anchor_gene_seq']!= "-") { #if the anchor gene has sequence there
#        fel_tbl[msa_position, 'anchor_gene_residue_number'] <- anchor_gene_position #give it a number
#        anchor_gene_position <- anchor_gene_position + 1 #and increment the number
#      }
#    }
# 
#   fel_tbl$anchor_gene_residue_number <- na.locf(fel_tbl$anchor_gene_residue_number) #replace NAs in gene residue number with nearest gene residue
#   colnames(fel_tbl)[1] <- 'Site' #rename column
   
####original version, loses all selection values not in anchor gene sequence ####   
    fel_tbl <- fel_tbl[fel_tbl$anchor_gene_seq != "-", c('Site', 'Selection_detected')] #remove all MSA positions that anchor gene didn't have sequence in, retain only site and selection test column
    rownames(fel_tbl) <- 1:(nrow(fel_tbl)) #fix row numbers
    fel_tbl$Site <- 1:(nrow(fel_tbl)) #remove alignment residue numbering, replace with gene residue
    fel_tbl <- fel_tbl[fel_tbl$Selection_detected != "",] #remove rows without selection detected, so table is now only a list of positions with selection, and which type.
    fel_tbl$Selection_detected <- as.factor(fel_tbl$Selection_detected) #store as factor
    
    
  #### SLAC test #####
    slac_neg_pval <- .001
    slac_pos_pval <- .2
    
            slac_tbl <- read.csv(list.files(pattern = "SLAC*"),header=TRUE)   #import FEL output file and identify sites that are under selection
    slac_tbl <- slac_tbl %>%
      mutate(Selection_detected = case_when( #call postive and negative selection for any given P value
        ( P..dN.dS...1..1 < slac_neg_pval ) ~ 'Negative',
        ( P..dN.dS...1. < slac_pos_pval ) ~ 'Positive', 
        TRUE ~ ''))
    
    #check if agrees with datamonkey call
    sum(slac_tbl$Selection_detected == "Negative")
    sum(slac_tbl$Selection_detected == "Positive")
    
    
    #strip positions of SLC table that are not present in anchor sequence
    slac_tbl$anchor_gene_seq <- unlist(strsplit(CDS_anchor_sequence, ""))[1:(str_length(CDS_anchor_sequence)/3) * 3] #print every third character of anchor gene in column
    slac_tbl <- slac_tbl[slac_tbl$anchor_gene_seq != "-", c('Site', 'Selection_detected')] #table with just the positions from anchor gene
    rownames(slac_tbl) <- 1:(nrow(slac_tbl)) #fix row numbers
    slac_tbl$Site <- 1:(nrow(slac_tbl)) #remove alignment residue numbering, replace with gene residue
    slac_tbl <- slac_tbl[slac_tbl$Selection_detected != "",] #remove rows without selection detected
    slac_tbl$Selection_detected <- as.factor(slac_tbl$Selection_detected) #store as factor
    
    
    
    ####FUBAR test####
    fubar_neg_prob <- 0.999
    fubar_pos_prob <- 0.7
    fubar_tbl <- read.csv(list.files(pattern = "FUBAR*"),header=TRUE)   #import FEL output file and identify sites that are under selection
    fubar_tbl <- fubar_tbl %>%
      mutate(Selection_detected = case_when( 
        (X.alpha. > X.beta.) & ( Prob..alpha...beta.. > fubar_neg_prob ) ~ 'Negative',
        (X.alpha. < X.beta.) & (Prob..alpha...beta...1 > fubar_pos_prob ) ~ 'Positive',
        TRUE ~ '')) #These are the criteria used to call postive and negative selection for any given prob value
    
    fubar_tbl$anchor_gene_seq <- unlist(strsplit(CDS_anchor_sequence, ""))[1:(str_length(CDS_anchor_sequence)/3) * 3] #print every third character of anchor gene in column
    fubar_tbl <- fubar_tbl[fubar_tbl$anchor_gene_seq != "-", c('Site', 'Selection_detected')] #table with just the positions from anchor gene
    rownames(fubar_tbl) <- 1:(nrow(fubar_tbl)) #fix row numbers
    fubar_tbl$Site <- 1:(nrow(fubar_tbl)) #remove alignment residue numbering, replace with gene residue
    fubar_tbl <- fubar_tbl[fubar_tbl$Selection_detected != "",] #remove rows without selection detected
    fubar_tbl$Selection_detected <- as.factor(fubar_tbl$Selection_detected) #store as factor
  
    sum(fubar_tbl$Selection_detected == 'Positive') #check that results agree with web interface results
    sum(fubar_tbl$Selection_detected == 'Negative')#check that results agree with web interface results
    
    
    ####MEME test####
    meme_pos_pval <- 0.1
    meme_tbl <- read.csv(list.files(pattern = "MEME*"),header=TRUE)   #import FEL output file and identify sites that are under selection
    meme_tbl <- meme_tbl %>%
      mutate(Selection_detected = case_when( 
        ( p.value <= meme_pos_pval ) ~ 'Positive',
        TRUE ~ '')) #These are the criteria used to call postive and negative selection for any given prob value
 
    #   sum(meme_tbl$Selection_detected == 'Positive') #check that results agree with web interface results
    
    meme_tbl$anchor_gene_seq <- unlist(strsplit(CDS_anchor_sequence, ""))[1:(str_length(CDS_anchor_sequence)/3) * 3] #print every third character of anchor gene in column
    meme_tbl <- meme_tbl[meme_tbl$anchor_gene_seq != "-", c('Site', 'Selection_detected')] #table with just the positions from anchor gene
    rownames(meme_tbl) <- 1:(nrow(meme_tbl)) #fix row numbers
    meme_tbl$Site <- 1:(nrow(meme_tbl)) #remove alignment residue numbering, replace with gene residue
    meme_tbl <- meme_tbl[meme_tbl$Selection_detected != "",] #remove rows without selection detected
    meme_tbl$Selection_detected <- as.factor(meme_tbl$Selection_detected) #store as factor
    


    
    
  ####   Viualize trace plot #######
    
    protein_y_pos <- 3 # set the position of the protein model over the chart
    protein_y_range <- .1 # set the thickness of the protein model over the chart

    trace_plot <- ggplot() +
      #Draw the evolution rate traces
      geom_line(data = vis_tab, aes(Gene_position, evolution_rate, color = factor(metric, labels = c("Site evolution rate", "50 residue mean")), alpha=metric)) + 
      scale_color_manual(values=c('gray','forestgreen')) + 
      scale_alpha_manual(values = c(.4,1),guide=F) + #set site rate to partially transparent
      #These draw lines at the mean evolution rate of each domain
      geom_segment(aes(x=min(lrr_range),xend=max(lrr_range),y=mean(na.omit(anchor_gene_site_rates_table[lrr_range,'Rate'])),yend=mean(na.omit(anchor_gene_site_rates_table[lrr_range,'Rate']))),linetype = "dotted", size = 1) + #LRR average
      annotate("text", hjust = 1, x=min(lrr_range)-10, y=mean(na.omit(anchor_gene_site_rates_table[lrr_range,'Rate'])), label = "LRR mean", size = 4) +
#      geom_segment(aes(x=min(tm_range),xend=max(tm_range),y=mean(na.omit(anchor_gene_site_rates_table[tm_range,'Rate'])),yend=mean(na.omit(anchor_gene_site_rates_table[tm_range,'Rate']))),linetype = "dotted", size = 1) + #TM average
#      annotate("text", hjust = 1, x=min(tm_range)-10, y=mean(na.omit(anchor_gene_site_rates_table[tm_range,'Rate'])), label = "TM mean", size = 4) +
      geom_segment(aes(x=min(rlk_range),xend=max(rlk_range),y=mean(na.omit(anchor_gene_site_rates_table[rlk_range,'Rate'])),yend=mean(na.omit(anchor_gene_site_rates_table[rlk_range,'Rate']))),linetype = "dotted", size = 1) + #RLK average
      annotate("text", hjust = 1, x=min(rlk_range)-10, y=mean(na.omit(anchor_gene_site_rates_table[rlk_range,'Rate'])), label = "RLK mean", size = 4) +
      # Visual language
      theme_bw() +
      theme(legend.position = c(.85, .5),legend.key.width = unit(.1, "npc") ,legend.background = element_rect(fill="white"), legend.box.background = element_rect(colour = "black", size = 2)) +
#      theme(legend.position = c(.15, .85),legend.key.width = unit(.1, "npc") ,legend.background = element_rect(fill="white"), legend.box.background = element_rect(colour = "black", size = 2)) +
      labs(title = paste(anchor_gene, "site rates and selection tests using", length(CDS_MSA), "genes", sep = " ") , x= "Position in anchor gene", y = "Rate of evolution") + 
      labs(color = "Evo Rate") +
      ylim(0, 4.15) +
      # These draw the protien domain shapes
      geom_rect(mapping=aes(xmin=0,xmax=nrow(anchor_gene_site_rates_table), ymin=protein_y_pos-.5*protein_y_range,  ymax=protein_y_pos+ .5*protein_y_range),colour = "darkgrey", fill = "grey") + 
      statebins:::geom_rrect(mapping=aes(xmin=min(lrr_range), xmax=max(lrr_range), ymin=protein_y_pos-protein_y_range,  ymax=protein_y_pos+protein_y_range), colour = "black", fill = "#1A96AA") +
      statebins:::geom_rrect(mapping=aes(xmin=min(tm_range), xmax=max(tm_range), ymin=protein_y_pos-protein_y_range,  ymax=protein_y_pos+protein_y_range), colour = "black", fill = "red") +
      statebins:::geom_rrect(mapping=aes(xmin=min(rlk_range), xmax=max(rlk_range), ymin=protein_y_pos-protein_y_range,  ymax=protein_y_pos+protein_y_range), colour = "black", fill = "#72E8D1") +
      #draw selection detected
      geom_rect(data = fel_tbl[fel_tbl$Selection_detected == "Positive",], aes(xmin = Site-1, xmax = Site +1, ymin = protein_y_pos + .2, ymax = protein_y_pos + .3), fill = "red", alpha=0.7) + #need to comment this out if no pos selected residues
      geom_rect(data = fel_tbl[fel_tbl$Selection_detected == "Negative",], aes(xmin = Site-2, xmax = Site +2, ymin = protein_y_pos + .1, ymax = protein_y_pos + .2), fill = "blue",  alpha=0.09) +
      annotate("text", hjust = 1, x=-5, y= protein_y_pos+.2, label = "FEL", size = 3) +
      
      geom_rect(data = fubar_tbl[fubar_tbl$Selection_detected == "Positive",], aes(xmin = Site-1, xmax = Site +1, ymin = protein_y_pos +.8, ymax = protein_y_pos +.9), fill = "red", alpha=0.7) + #need to comment this out if no pos selected residues
      geom_rect(data = fubar_tbl[fubar_tbl$Selection_detected == "Negative",], aes(xmin = Site-2, xmax = Site +2, ymin = protein_y_pos + .7, ymax = protein_y_pos +.8), fill = "blue",  alpha=0.09) +
      annotate("text", hjust = 1, x=-5, y= protein_y_pos+.8, label = "FUBAR", size = 3) +
      
      geom_rect(data = slac_tbl[slac_tbl$Selection_detected == "Positive",], aes(xmin = Site-1, xmax = Site +1, ymin = protein_y_pos+ .5, ymax = protein_y_pos +.6), fill = "red", alpha=0.7) + #need to comment this out if no pos selected residues
      geom_rect(data = slac_tbl[slac_tbl$Selection_detected == "Negative",], aes(xmin = Site-2, xmax = Site +2, ymin = protein_y_pos + .4, ymax = protein_y_pos +.5), fill = "blue",  alpha=0.09) +
      annotate("text", hjust = 1, x=-5, y= protein_y_pos+.5, label = "SLAC", size = 3) + 
      
      geom_rect(data = meme_tbl[meme_tbl$Selection_detected == "Positive",], aes(xmin = Site-1, xmax = Site +1, ymin = protein_y_pos+ 1, ymax = protein_y_pos +1.1), fill = "red", alpha=0.7) + #need to comment this out if no pos selected residues
      annotate("text", hjust = 1, x=-5, y= protein_y_pos+ 1.05, label = "MEME", size = 3)  
      
    
  trace_plot   
    
    #ggsave(plot = last_plot(), filename = paste(anchor_gene, "evolution_rate_and_selection_tests", "pdf", sep = "."), height = 9, width = 12)
 
  #check for zeros if above plot fails
  sum(fubar_tbl$Selection_detected == 'Positive')
  sum(fubar_tbl$Selection_detected == 'Negative')
  sum(slac_tbl$Selection_detected == "Negative")
  sum(slac_tbl$Selection_detected == "Positive")
  sum(fel_tbl$Selection_detected == 'Positive')
  sum(fel_tbl$Selection_detected == 'Negative')

    #find sites under pos selection
  which(slac_tbl$Selection_detected == "Positive")

  
  
 
  
  
   ##    BEAST section, for visualizing posterior probabilities of relative rates #########    
  #return to rate directory
  setwd("../site_rates/")

  #import log file from BEAST output to get the relative rates  
  beast_trace_table <- read.table(list.files(pattern = "*.log"), comment.char = "#", header=T, sep = "\t") #import beast output
  beast_trace_table <- beast_trace_table[round((burnin * nrow(beast_trace_table))):nrow(beast_trace_table), ] # remove burn in
  
  rate.mean <- mean(beast_trace_table$rate.mean) * mean(beast_trace_table$TreeHeight) #calculate rate.mean, the mean substitutaion rate over the entire tree
  beast_trace_table <- beast_trace_table[ ,c("Sample", "mutationRate.LRR","mutationRate.RLK")] # just retain relative rates
    #modify the relative mutation rates of both domains to be in units of substitutions per site per age of tree
    beast_trace_table$mutationRate.LRR <- beast_trace_table$mutationRate.LRR * rate.mean
    beast_trace_table$mutationRate.RLK <- beast_trace_table$mutationRate.RLK * rate.mean
  beast_trace_table_long <- beast_trace_table %>% pivot_longer(cols = c("mutationRate.LRR","mutationRate.RLK"), names_to = "Mutation_Rate_Estimate") #for visualiztion
                                     
 beast_plot <- ggplot(beast_trace_table_long, aes(`Mutation_Rate_Estimate`, value)) + 
    geom_boxplot(fill = c("#1A96AA", "#72E8D1"), outlier.alpha = 0) +
   labs(title = "BEAST posterior rates", x= "Domain") + 
    theme_bw() + theme(axis.title.y=element_blank()) +
    geom_signif(comparisons = list(c('mutationRate.LRR', 'mutationRate.RLK')), map_signif_level = T, test = 't.test') +
    expand_limits(y=0)
      

########visualize both together#######
  
 #plot_grid(trace_plot, beast_plot, labels = c('A', 'B'), rel_widths = c(4, 1)) #original version with regular box plot
 #   ggsave(plot = last_plot(), filename = paste(anchor_gene, "evolution_rate_and_selection_tests", "pdf", sep = "."), height = 6, width = 14)
  #niftier, with jittered data dots. half violin can be made if you want to fiddle with raincloud plot
      beast_plot <-      ggplot(beast_trace_table_long, aes(`Mutation_Rate_Estimate`, value)) + 
        geom_jitter(data = beast_trace_table_long, aes(x = Mutation_Rate_Estimate, y = value, color = Mutation_Rate_Estimate), size = .45, alpha = .15, width = .45) + 
        scale_color_manual(values = c("#1A96AA", "#72E8D1")) +
 #       geom_flat_violin(data = beast_trace_table_long, aes( x = Mutation_Rate_Estimate, y = value, fill = Mutation_Rate_Estimate), width = .8, position = position_nudge(x = .3, y =0)) +
        geom_boxplot(data = beast_trace_table_long, aes(x = Mutation_Rate_Estimate, y = value ,  fill = Mutation_Rate_Estimate), outlier.alpha = 0, width = .5, alpha = 0.4) +
            scale_fill_manual(values=c("#1A96AA", "#72E8D1")) +
        labs(title = "Posterior rate estimates", x= "Domain") + 
        theme_bw() + theme(axis.title.y=element_blank()) +
        geom_signif(comparisons = list(c('mutationRate.LRR', 'mutationRate.RLK')), map_signif_level = T, test = 't.test') +
        theme(legend.position = "none")
      plot_grid(trace_plot, beast_plot, labels = c('A', 'B'), rel_widths = c(4, 1))
      
#      ggsave(plot = last_plot(), filename = paste(anchor_gene, "evolution_rate_and_selection_tests", "pdf", sep = "."), height = 6, width = 14)
   
      
      
  #initital way to quantify selection: count number of positivley selected sites per sites in domain (rate)
      #set up array, 3 dimentions, 1 = test results (7 categories), 2 = gene (10 categories), 3 = domain (2 categories, [,,1] = LRR domain, [,,2] = RLK domain
      #selection_rate_tbl <- array(NA, c(7,11,2)) 
      #dimnames(selection_rate_tbl)[[1]] <- c('FUBAR_positive', 'FUBAR_negative', 'FEL_positive', 'FEL_negative', 'SLAC_positive', 'SLAC_negative', 'MEME_positive')
      #dimnames(selection_rate_tbl)[[2]] <- c('BAM1','CEPR2', 'CLV1', 'FLS2', 'GSO', 'HAE', 'HBD', 'IKU2', 'PEPR', 'PXY', 'RGI5')
      #dimnames(selection_rate_tbl)[[3]] <- c('LRR', 'RLK')
  
      gene_index <- 5

      fubar_lrr_tbl <- fubar_tbl[fubar_tbl$Site > min(lrr_range) & fubar_tbl$Site < max(lrr_range),] #just the region over LRR domain
      selection_rate_tbl['FUBAR_positive', gene_index, 1 ] <- sum(fubar_lrr_tbl$Selection_detected == "Positive") / length(lrr_range) #positively selected sites per residue in LRR domain
      selection_rate_tbl['FUBAR_negative', gene_index, 1 ] <- sum(fubar_lrr_tbl$Selection_detected == "Negative") / length(lrr_range) #negatively selected sites per residue in LRR domain
      
      fubar_rlk_tbl <- fubar_tbl[fubar_tbl$Site > min(rlk_range) & fubar_tbl$Site < max(rlk_range),] #just the region over LRR domain
      selection_rate_tbl['FUBAR_positive', gene_index, 2 ] <- sum(fubar_rlk_tbl$Selection_detected == "Positive") / length(rlk_range) #positively selected sites per residue in LRR domain
      selection_rate_tbl['FUBAR_negative', gene_index, 2 ] <- sum(fubar_rlk_tbl$Selection_detected == "Negative") / length(rlk_range) #negatively selected sites per residue in LRR domain
      
      fel_lrr_tbl <- fel_tbl[fel_tbl$Site > min(lrr_range) & fel_tbl$Site < max(lrr_range),] #just the region over LRR domain
      selection_rate_tbl['FEL_positive', gene_index, 1 ] <- sum(fel_lrr_tbl$Selection_detected == "Positive") / length(lrr_range) #positively selected sites per residue in LRR domain
      selection_rate_tbl['FEL_negative', gene_index, 1 ] <- sum(fel_lrr_tbl$Selection_detected == "Negative") / length(lrr_range) #negatively selected sites per residue in LRR domain
      
      fel_rlk_tbl <- fel_tbl[fel_tbl$Site > min(rlk_range) & fel_tbl$Site < max(rlk_range),] #just the region over LRR domain
      selection_rate_tbl['FEL_positive', gene_index, 2 ] <- sum(fel_rlk_tbl$Selection_detected == "Positive") / length(rlk_range) #positively selected sites per residue in LRR domain
      selection_rate_tbl['FEL_negative', gene_index, 2 ] <- sum(fel_rlk_tbl$Selection_detected == "Negative") / length(rlk_range) #negatively selected sites per residue in LRR domain
      
      slac_lrr_tbl <- slac_tbl[slac_tbl$Site > min(lrr_range) & slac_tbl$Site < max(lrr_range),] #just the region over LRR domain
      selection_rate_tbl['SLAC_positive', gene_index, 1 ] <- sum(slac_lrr_tbl$Selection_detected == "Positive") / length(lrr_range) #positively selected sites per residue in LRR domain
      selection_rate_tbl['SLAC_negative', gene_index, 1 ] <- sum(slac_lrr_tbl$Selection_detected == "Negative") / length(lrr_range) #negatively selected sites per residue in LRR domain
      
      slac_rlk_tbl <- slac_tbl[slac_tbl$Site > min(rlk_range) & slac_tbl$Site < max(rlk_range),] #just the region over LRR domain
      selection_rate_tbl['SLAC_positive', gene_index, 2 ] <- sum(slac_rlk_tbl$Selection_detected == "Positive") / length(rlk_range) #positively selected sites per residue in LRR domain
      selection_rate_tbl['SLAC_negative', gene_index, 2 ] <- sum(slac_rlk_tbl$Selection_detected == "Negative") / length(rlk_range) #negatively selected sites per residue in LRR domain
      
      meme_lrr_tbl <- meme_tbl[meme_tbl$Site > min(lrr_range) & meme_tbl$Site < max(lrr_range),] #just the region over LRR domain
      selection_rate_tbl['MEME_positive', gene_index, 1 ] <- sum(meme_lrr_tbl$Selection_detected == "Positive") / length(lrr_range) #positively selected sites per residue in LRR domain

      meme_rlk_tbl <- meme_tbl[meme_tbl$Site > min(rlk_range) & meme_tbl$Site < max(rlk_range),] #just the region over LRR domain
      selection_rate_tbl['MEME_positive', gene_index, 2 ] <- sum(meme_rlk_tbl$Selection_detected == "Positive") / length(rlk_range) #positively selected sites per residue in LRR domain

    
      
      
      
      
      #sum all before calculating rates
    sum(sum(fubar_lrr_tbl$Selection_detected == "Positive") + 
        sum(fel_lrr_tbl$Selection_detected == "Positive") + 
        sum(slac_lrr_tbl$Selection_detected == "Positive") + 
        sum(meme_lrr_tbl$Selection_detected == "Positive")) / 
        length(lrr_range) # LRR pos rate
      
    sum(sum(fubar_rlk_tbl$Selection_detected == "Positive") +
        sum(fel_rlk_tbl$Selection_detected == "Positive") +
        sum(slac_rlk_tbl$Selection_detected == "Positive") + 
        sum(meme_rlk_tbl$Selection_detected == "Positive")) /
        length(rlk_range) # RLK pos rate
        
    sum(sum(fubar_lrr_tbl$Selection_detected == "Negative") + 
          sum(fel_lrr_tbl$Selection_detected == "Negative") + 
          sum(slac_lrr_tbl$Selection_detected == "Negative") + 
          sum(meme_lrr_tbl$Selection_detected == "Negative")) / 
      length(lrr_range) # LRR neg rate
    
    sum(sum(fubar_rlk_tbl$Selection_detected == "Negative") +
          sum(fel_rlk_tbl$Selection_detected == "Negative") +
          sum(slac_rlk_tbl$Selection_detected == "Negative") + 
          sum(meme_rlk_tbl$Selection_detected == "Negative")) /
      length(rlk_range) # RLK neg rate
    
    
#this finds rate for each then averages the averages to weigh each test equally
        mean(c(sum(fubar_lrr_tbl$Selection_detected == "Positive") / length(lrr_range) ,
          sum(fel_lrr_tbl$Selection_detected == "Positive") / length(lrr_range)  ,
          sum(slac_lrr_tbl$Selection_detected == "Positive") / length(lrr_range)  ,
          sum(meme_lrr_tbl$Selection_detected == "Positive") / length(lrr_range))) #lrr pos rate

        mean(c(sum(fubar_rlk_tbl$Selection_detected == "Positive") / length(rlk_range) ,
               sum(fel_rlk_tbl$Selection_detected == "Positive") / length(rlk_range)  ,
               sum(slac_rlk_tbl$Selection_detected == "Positive") / length(rlk_range)  ,
               sum(meme_rlk_tbl$Selection_detected == "Positive") / length(rlk_range)))    # rlk pos rate
      
        mean(c(sum(fubar_lrr_tbl$Selection_detected == "Negative") / length(lrr_range) ,
               sum(fel_lrr_tbl$Selection_detected == "Negative") / length(lrr_range)  ,
               sum(slac_lrr_tbl$Selection_detected == "Negative") / length(lrr_range)  ,
               sum(meme_lrr_tbl$Selection_detected == "Negative") / length(lrr_range))) #lrr neg rate
        
        mean(c(sum(fubar_rlk_tbl$Selection_detected == "Negative") / length(rlk_range) ,
               sum(fel_rlk_tbl$Selection_detected == "Negative") / length(rlk_range)  ,
               sum(slac_rlk_tbl$Selection_detected == "Negative") / length(rlk_range)  ,
               sum(meme_rlk_tbl$Selection_detected == "Negative") / length(rlk_range)))         # rlk neg rate 
      
      
      
      
      
      
      
      ####experimental stuff        
      
      
   #line colored by domain   
   iqtree_tab_trace_table <-   vis_tab[vis_tab$metric == 'sliding_window_rate',	] #just sliding window values
#define positions in gene
      iqtree_tab_trace_table$domain <- 'none'
     iqtree_tab_trace_table[1:min(lrr_range), 'domain'] <- "N_term"
     iqtree_tab_trace_table[max(lrr_range):min(tm_range), 'domain'] <- "none1"
     iqtree_tab_trace_table[lrr_range, 'domain'] <- "lrr"
     iqtree_tab_trace_table[max(tm_range):min(rlk_range), 'domain'] <- "none2"
     iqtree_tab_trace_table[tm_range, 'domain'] <- "tm"
     iqtree_tab_trace_table[max(rlk_range):nrow(iqtree_tab_trace_table), 'domain'] <- "C_term"
     iqtree_tab_trace_table[rlk_range, 'domain'] <- "rlk"
   
   
            ggplot() +
        geom_path(data = iqtree_tab_trace_table, aes(Gene_position, evolution_rate, color = factor(domain), alpha=1)) 
            + scale_color_manual(values = c("N_term" = "gray","lrr"= "#1A96AA","none1"= "gray", "tm"= "41C1CE", "gray", "#72E8D1", "gray"))
            scale_color_manual(values = c("4" = "blue", "6" = "green", "8" = "purple", 
                                          "group1" = "black", 
                                          "group2" = "red"))
      
      
      
      
   length(CDS_MSA)
      
      ###old stuff
   
   
    geom_quasirandom(fill = "lightseagreen", colour = "lightseagreen", bandwidth=2,alpha=.1)
    theme_bw()

  +
    geom_quasirandom(fill = "lightseagreen", colour = "lightseagreen", bandwidth=2,alpha=.1)  +
    geom_jitter(width = 0.1, alpha = .1, shape = 1)
    
  
  
  beast_trace_table %>% 

  
  
  
  
  
  
  #this geom_ribbon fills in area between line and average. But still needs some work...
   #geom_ribbon(data = vis_tab, mapping = aes(x = Gene_position, y = evolution_rate, ymin = mean(anchor_gene_site_rates_table$Rate), ymax = evolution_rate, fill = metric), linetype = "dotted") 

#geom_point(data = contact_residues, aes(Site, evolution_rate), color = "red") + 
   # geom_rect(data=contact_residues, mapping=aes(xmin= Site -2, xmax=Site + 2, ymin=.5, ymax=1), fill = "orange", alpha = 1) 







ggplot() + geom_point(data = contact_residues, aes(Site, evolution_rate)) 





#hae_seq <- as.character(aa2$AT4G28490.1) #print sequence of HAE in alignment, including dashes
#hae_seq <- as.character(aa2$AT1G28440.1) #print sequence of HSL1 in alignment, including dashes
#hae_seq <- as.character(aa2$Solyc03g006300.1.1) #print sequence of SlHAE in alignment, including dashes
hae_seq <- as.character(aa2$LOC_Os01g13800.1) #print sequence of OsHAE in alignment, including dashes




#hae residues in contact with ligand  ####
contact_residues <- data.frame(c(409,407,385,361,364,364,383,339,337,337,315,313,290,290,268,268,313,290,266,266,264,264,240,242,218,196,196,196,196,172,148))
colnames(contact_residues) <- "Site"
df <- inner_join(contact_residues, tab4, by = "Site")
contact_residues <- df[df$metric == "Rate",c(1,3)]
typeof(contact_residues$evolution_rate)
contact_residues$evolution_rate <- as.numeric(contact_residues$evolution_rate)
contact_residues$Site <- as.numeric(contact_residues$Site)
######




