
pacman::p_load('dplyr', 'tidyr', 'htmltools', 'bbplot', 'scales', 'data.table', 'ggmsa',
               'ggplot2', 'rdrop2', 'shiny', 'BiocManager', 'DECIPHER', 'ComplexUpset',
               'dendextend', 'data.table', 'Biostrings', 'alakazam', "unikn", 'ggupset',
               'plotly', "jcolors", 'ggdendro', "RColorBrewer","kmer","heatmaply", "ggseqlogo", 
               "msa", "stringr", install = F)



bag_data <- fread("MUSA_with_rss_and_leader_and_asc_2025-02-05.csv")
setnames(bag_data, "allele", "old_name")
bag_data[, allele := new_tag]  # Assuming new names match old names

bag_data[, `:=`(
  gene = stringr::str_extract(allele, "[^*]+"),
  ASC = sub("\\*.*", "", allele),
  Family = sub("-.*", "", allele),
  chain = substr(allele, 1, 3)
)]

bag_data <- bag_data[!(bag_data$allele == ""),]


data_ <- fread("repertoire_genotype_filter_d_04_02_rep_data.csv")
setnames(data_, "allele", "old_allele")
# Update allele in data_ based on bag_data
data_[, allele := bag_data[match(old_allele, old_name), new_tag]]
data_ <- data_[!is.na(allele)]
data_[, sum_count := sum(count), by = sample]
data_[, frac := count / sum_count]
data_[, frac_allele := count / sum_count]

optimized_thresholds_f <- fread("04_02_thresholds_table.tsv")
setnames(optimized_thresholds_f, "allele", "old_allele")
optimized_thresholds_f[, allele := bag_data[match(old_allele, old_name), new_tag]]

bag_data[, main_allele := ifelse(
  allele == "" | is.na(allele), 
  allele, 
  sub("_.*", "", allele)
)]


# Step 3: Prepare a new `bag_data` table with only the `allele` and updated `threshold`
bag_threshold <- bag_data[, .(allele, main_allele)]

optimized_thresholds <- merge(
  bag_threshold,
  optimized_thresholds_f,
  by = "allele",
  all.x = TRUE
)

# Step 2: Assign default thresholds for rows with NA in `threshold`
optimized_thresholds[, threshold := fifelse(
  is.na(threshold), 
  fifelse(
    grepl("IGHV", main_allele), 1e-04,  # Default for IGHV
    1e-03                             # Default for other gene types
  ),
  threshold
)]


allele_appearance <- function(data_, imgt_genes,imgt_bag_data, chain = "IGH") {
  
  result <- imgt_bag_data %>%
    select(allele, samples_AIRRseq, samples_genomic) %>%
    group_by(allele)%>%
    mutate(
      samples_genomic = if_else(is.na(samples_genomic), list(""), strsplit(samples_genomic, ",")),
      samples_AIRRseq = if_else(is.na(samples_AIRRseq), list(""), strsplit(samples_AIRRseq, ","))
    ) %>%
    summarize(
      allele = allele,
      count_in_both = length(intersect(unlist(samples_genomic), unlist(samples_AIRRseq)) %>% 
                              .[. != "" & !is.na(.)]),
      count_only_genomic = length(setdiff(unlist(samples_genomic), unlist(samples_AIRRseq)) %>% 
                                    .[. != "" & !is.na(.)]),
      count_only_AIRRseq = length(setdiff(unlist(samples_AIRRseq), unlist(samples_genomic)) %>% 
                                    .[. != "" & !is.na(.)])
    )

  result <- unique(result)

  # Reshape for plotting
  result_long <- result %>%
    pivot_longer(
      cols = starts_with("count"),
      names_to = "category",
      values_to = "count"
    ) %>%
    mutate(category = recode(category, 
                            "count_in_both" = "Both",
                            "count_only_genomic" = "Genomic Only",
                            "count_only_AIRRseq" = "Repertoire Only",
                            .default = "Not in Database"))

  # Filter and adjust height
  data_ <- data_[grepl(imgt_genes, data_$allele),]
  alleles <- data_$allele
  data_$imgt_call <- data_$allele

  height = length(unique(data_$imgt_call))
  height = ifelse(height < 20, 25, height)
  height = height * 30

  # Define colors
  category_colors <- setNames(hcl.colors(4, palette = "Sunset"), 
                              c("Both", "Repertoire Only", "Genomic Only", "Not in Database"))

  # Plot the processed data
  p <- ggplot(result_long, aes(x = allele, y = count, fill = category)) + 
    geom_bar(stat = "identity") +  # Use stat = "identity" to plot the precomputed counts
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = category_colors) +
    labs(x = "Allele", y = "# Individuals", fill = "") + 
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(
        size = 10
      ), 
      axis.title = element_text(size = 14)
    )
    

  p1<-ggplotly(p, height = height-height/10, width = height)
  
  
  p <- ggplot(data_, 
       aes(x = allele, y = frac_allele,
           text = paste0("sample: ", sample))) +
    geom_boxplot() + # Add boxplot to show the distribution
    geom_point(data = subset(data_, frac_allele != 0), 
               aes(color = factor(frac_allele != 0)), # Add color aesthetic if needed
               position = position_jitter(width = 0.2), 
               alpha = 1, size = 1.5) +
    scale_color_manual(values = c("TRUE" = alpha("red", 1)),
                       guide = "none") + # Hide legend if not needed
    labs(x = "Allele", y = "Usage (%)") + 
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(
        size = 10
      ), 
      axis.title = element_text(size = 14)
    )



  p2<-ggplotly(p, height = height-height/10, width = height)
  
  sub_p1 <- subplot(p2, # boxplot and point - decrease size of points and remove legend
                    p1,
                    nrows = 2, margin = 0.007)
  
    return(list(p1 = ggplotly(p1, height = height-height/10, width = height-height/4),
              p2 = ggplotly(p2, height = height-height/10, width = height)))
  
}



heatmap_alleles <-function(data_, g_group, allele_db) {

    data_  <- data_[grepl(imgt_genes, data_$allele),]
    
    alleles <- data_$allele
    data_$imgt_call <- data_$allele
    data_upset <- data_ %>%
      select(sample, imgt_call, frac) %>%
      group_by(sample) %>%
      dplyr::mutate(frequency = sum(as.numeric(frac) / n(), na.rm = T),
                    val = 1, 
                    text = paste("</br>sample: ", sample,
                                 "</br>Alleles: ", paste0(sort(unique(imgt_call)), collapse = ","),
                                 "</br>Frequency: ", unique(frequency)),
                    text2 = paste("</br>Alleles: ", paste0(sort(unique(imgt_call)), collapse = ","))) %>%
      select(-frac) %>%
      pivot_wider(names_from = imgt_call, 
                  values_from = val, 
                  values_fill = list(val = 0))
    
    
    
    t<-setNames(data_$imgt_call,data_$imgt_call)
    missing_alleles <- names(t)[!(t) %in% colnames(data_upset[,5:ncol(data_upset)])]
    
    data_upset[,missing_alleles] <- 0
    
    p_upset <- upset(
      data_upset,
      unique(t),
      annotations = list(
        'Frequency'=ggplot(mapping = aes(x=intersection, y = frequency)) +
          smplot2::sm_boxplot() + upset_themes$default
      ),
      width_ratio=0.1, encode_sets = F
    )
    
    p1 <- ggplotly(p_upset[[2]]+ aes(text = text), tooltip = "text")
    
    p2 <- ggplotly(p_upset[[4]] + aes(text = text2), tooltip = "text")
    p3 <- ggplotly(p_upset[[5]], tooltip = "text")
    p_interaction <- ggplotly(p_upset, tooltip = "x")
    
    sub_p1 <- subplot(hide_legend(p1), # boxplot and point - decrease size of points and remove legend
                      p2, 
                      hide_legend(p_interaction), 
                      nrows = 3, margin = 0.0007)
    sub_p2 <- subplot(plotly_empty(), 
                      plotly_empty(), 
                      p3, 
                      nrows = 3, margin = 0.0007)
    sub_p3 <- subplot(sub_p2, sub_p1, nrows = 1, widths= c(0.4,0.6) ,margin = 0.11)
    
    height = (length(unique(data_$imgt_call)))
    height = ifelse(length(height)<20, 25, height)
    height = height*100
    d<-ggplotly(sub_p3, height = height*1.5, width = height*1.2)
    return(d)
  }


seq_align <- function(sequences, alleles, imgt_genes, chain = "IGH") { 

    temp_output <- file("nul", open = "w") # For Windows use "nul", for Linux/macOS use "/dev/null"
    sink(temp_output, type = "output")
    sink(temp_output, type = "message")

    # Perform the msa alignment
    alignment <- msa(sequences, type = "dna")

    # Restore normal output
    sink(type = "output")
    sink(type = "message")
    close(temp_output)

    # Convert msa alignment to DNAStringSet
    alignment <- DNAStringSet(alignment)

    mat_sub <-
      DistanceMatrix(
        alignment,
        includeTerminalGaps = FALSE,
        penalizeGapGapMatches = FALSE,
        penalizeGapLetterMatches = T,
        verbose = F
      )
    
    if(length(mat_sub)==1) return(NULL)
    colnames(mat_sub) <-  gsub("IGH", "", colnames(mat_sub))
    rownames(mat_sub) <-  gsub("IGH", "", rownames(mat_sub))
    
    matrix_sequences <-
      as.data.frame(sapply(as.character(alignment), seqinr::s2c), stringsAsFactors = F)

    nucs <-
       nrow(matrix_sequences) - sum(apply(matrix_sequences, 1, function(x)
         all(x == ".")))
    if (length(alleles) < 3) {
      hc <- hclust(as.dist(mat_sub))
    } else{
      hc <- ape::nj(as.dist(mat_sub))
      hc <- ape::ladderize(hc)
    }
    dend <- as.dendrogram(hc)
    labels(dend) <- alleles
    
    dend <- dendextend::set(dend, "labels_cex", 0.5)
    ggd1 <- as.ggdend(dend)
    ggd1$labels$y <- ggd1$labels$y-0.01

    p_dend <- ggplot(ggd1, horiz = T,  theme = NULL)  +
      theme(
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_line(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        panel.grid.minor = element_line(color = "gray"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"
      ) +
      scale_y_continuous(limits = c(-0.035, NA), sec.axis = sec_axis(~ . * nucs, name = "Mutations")) +
      ylab("distance(nt)")
    
    pg <- ggplotly(p_dend)%>% layout(margin = list(l = 75))
    return(pg)
  }



seq_align2 <-function(sequences, alleles, imgt_genes, chain = "IGH") {    

# Assuming you have an allele column in your data
    temp_output <- file("nul", open = "w") # For Windows use "nul", for Linux/macOS use "/dev/null"
    sink(temp_output, type = "output")
    sink(temp_output, type = "message")

    # Perform the msa alignment
    alignment <- msa(sequences, type = "dna")

    # Restore normal output
    sink(type = "output")
    sink(type = "message")
    close(temp_output)

    # Convert msa alignment to DNAStringSet
    alignment <- DNAStringSet(alignment)

  sequences <- as.character(alignment)  
  
  len <-nchar(sequences[1])

    matrix_sequences <-
      as.data.frame(sapply(sequences, seqinr::s2c), stringsAsFactors = F)


      
    names(matrix_sequences) = alleles

    nucs <-
      nrow(matrix_sequences) - sum(apply(matrix_sequences, 1, function(x)
        all(x == ".")))
    
    snps <-
      which(apply(matrix_sequences, 1, function(x)
        length(unique(x)) != 1))
    
    matrix_sequences$annot <-
      apply(matrix_sequences, 1, function(x)
        length(unique(x)) != 1)

    matrix_sequences$pos <- 1:len

    matrix_sequences_plot <-
      reshape2::melt(matrix_sequences, id.vars = c("pos", "annot"))
    matrix_sequences_plot$id <- matrix_sequences_plot$pos
    matrix_sequences_plot$allele <-
      gsub("IGH", "", matrix_sequences_plot$variable)
    matrix_sequences_plot$allele <-
      factor(matrix_sequences_plot$allele,
             levels = unique(matrix_sequences_plot$allele))
    
    # matrix_sequences_plot$value[matrix_sequences_plot$value == "."] <-
    #   NA
    matrix_sequences_plot$annot_text <-
      sapply(1:nrow(matrix_sequences_plot), function(i)
        ifelse(
          matrix_sequences_plot$annot[i],
          matrix_sequences_plot$value[i],
          ""
        ))
    
    hotspot <- c()
    if (length(snps) != 0) {
      for (s in snps) {
        ht_snp <- c()
        for (i in 1:(ncol(matrix_sequences) - 2)) {
          if (s > 3)
            ht_snp <- c(ht_snp,
                        grepl("G[TC][AT]",
                              paste0(matrix_sequences[(s - 1):(s + 2), i], collapse = "")) |
                          grepl("[AT][AG]C",
                                paste0(matrix_sequences[(s -
                                                           2):(s + 1), i], collapse = "")))
        }
        if (any(ht_snp))
          hotspot <- c(hotspot, s)
      }
    }

    plot_align <-
      function(dt,
               low_bound = 1,
               upper_boud = 80,
               hotspot) {
        if (length(hotspot) != 0)
          ht <-
            hotspot[which(as.numeric(hotspot) %in% low_bound:upper_boud)]
        else
          ht <- NULL
        p <- ggplot(dt[dt$id >= low_bound &
                           dt$id < upper_boud, ]) +
          geom_tile(aes(
            x = (pos),
            y = (allele),
            fill = value
          ), colour = "white") +
          geom_text(aes(
            x = (pos),
            y = (allele),
            label = annot_text
          ), color = "black") +
          #coord_equal(expand = F, xlim = c(low_bound, upper_boud), ratio = 9/5, clip = "off") +
          #bbplot::bbc_style() +
          scale_fill_manual(values = setNames(
            c(unname(jcolors("pal2")[c(1, 3, 4, 5)]), "gray50"),
            c("A","C","G","T",".")
            )) +
          theme_align # "#1380A1", "#FAAB18", "#990000", "#588300"
        
        if (!is.null(ht)) {
          for (h in ht) {
            df_rect <- data.frame(ymin = 1,
                                  ymax = length(unique(dt$allele)),
                                  x = h)
            
            
            p <- p + geom_rect(
              data = df_rect,
              size = 1,
              fill = NA,
              linejoin = "bevel",
              lty = 1,
              colour = jcolors("pal2")[2],
              aes(
                xmin = x - 0.5,
                xmax = x + 0.5,
                ymin = ymin - 0.5,
                ymax = ymax + 0.5
              )
            )
          }
        }
        
        p <- ggplotly(p, tooltip = "none")
        
        if (low_bound == 1)
          p <- p %>% layout(legend = list(orientation = "h", x = 0, y = 1.1))
        return(p)
      }
    
    theme_align <- theme(
      axis.line = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(
        size = 12
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    )
    
    if (len < 80) {
      # Handle the case where len is less than 80
      p_list <- list(
        plot_align(matrix_sequences_plot, 1, len, hotspot)
      )
    } else {
      # Define the chunk size
      chunk_size <- 80
      
      # Calculate the number of chunks
      num_chunks <- ceiling(len / chunk_size)
      
      # Adjust len to make it a multiple of chunk_size by padding
      padded_len <- num_chunks * chunk_size
      
      # Create bounds for all chunks with the same width
      bounds <- data.frame(
        low_bound = seq(1, padded_len, by = chunk_size),
        upper_bound = seq(chunk_size, padded_len, by = chunk_size)
      )
      
      # Apply plot_align to each range
      p_list <- apply(bounds, 1, function(x) {
        xmin <- x[1]  # Calculate xmin
        xmax <- x[2]  # Calculate xmax
        plot_align(matrix_sequences_plot, xmin, xmax, hotspot)
      })
    }

    p_list
    
}
