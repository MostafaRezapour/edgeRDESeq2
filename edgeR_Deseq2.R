main <- function() {
  tryCatch({
    print("Starting script execution.")
    cat("\014")
    # Clear the environment
    rm(list = ls())
    while (!is.null(dev.list())) dev.off()
    print("Script has started running.")
    library(edgeR)
    library(tcltk)
    library(rstudioapi)
    print("Guideline")
    intro_message <- function() {
      require(tcltk)
      
      # Create the main window
      tt <- tktoplevel()
      tkwm.title(tt, "RNA-Seq Analysis Software")
      tcl("wm", "attributes", tt, topmost = TRUE)
   
      msg <- "EdgeR: Please make sure your CSV file is formatted accordingly."
      tkgrid(tklabel(tt, text = msg, wraplength = 600, justify = "left"), padx = 20, pady = 20)
      user_choice <- tclVar("continue")
      button_frame <- tkframe(tt)
      continue_btn <- tkbutton(button_frame, text = "Continue", command = function() {
        tclvalue(user_choice) <- "continue"  # Set choice to continue
        tkdestroy(tt)  # Close the window
      })
      exit_btn <- tkbutton(button_frame, text = "Exit and Fix", command = function() {
        tclvalue(user_choice) <- "exit"  # Set choice to exit
        tkdestroy(tt)  # Close the window
      })
      tkgrid(continue_btn, padx = 10, pady = 10)
      tkgrid(exit_btn, padx = 10, pady = 10)
      tkgrid(button_frame)
      tkwait.window(tt)
      if (tclvalue(user_choice) == "exit") {
        cat("Please fix the input CSV file and rerun the script.\n")
        return(FALSE)     }
      
      return(TRUE)  
    }
    if (!intro_message()) {
      stop("Execution stopped by user. Fix the input CSV file and rerun the script.")
    }
    cat("Script continues...\n")
    load_data <- function(filepath) {
      data <- read.csv(filepath, header = TRUE, check.names = FALSE)
      gene_symbols <- data[, 1]  # First column as gene symbols
      counts <- data[, -1]       # All other columns as counts
      dge <- DGEList(counts = as.matrix(counts), genes = data.frame(GeneSymbol = gene_symbols))
      o <- order(rowSums(dge$counts), decreasing = TRUE)
      dge <- dge[o, ]
      d <- duplicated(dge$genes$GeneSymbol)
      dge <- dge[!d, ]
      rownames(dge$counts) <- dge$genes$GeneSymbol
      return(dge)    }
    detect_groups <- function(column_names) {
      groups <- sub(" \\(.*\\)", "", column_names)
      group_info <- table(groups)
      return(as.list(group_info))    }
    confirm_groups <- function(group_info) {
      tt <- tktoplevel()
      tkwm.title(tt, "Confirm Groups")
      tkgrid(tklabel(tt, text = "Detected Groups and Sample Counts:"))
      for (group_name in names(group_info)) {
        tkgrid(tklabel(tt, text = sprintf("%s: %d samples", group_name, group_info[[group_name]])))
      }
      user_decision <- tclVar("no")

      tkgrid(tkbutton(tt, text = "Confirm", command = function() {
        tclvalue(user_decision) <- "yes"
        tkdestroy(tt)
      }))
      tkgrid(tkbutton(tt, text = "Decline", command = function() {
        tclvalue(user_decision) <- "no"
        tkdestroy(tt)
      }))
      tkwait.window(tt)
      
      return(tclvalue(user_decision) == "yes")
    }

    run_analysis <- function() {
      filepath <- file.choose()
      dge <- load_data(filepath)
      column_names <- colnames(dge$counts)
      group_info <- detect_groups(column_names)
      if (confirm_groups(group_info)) {
        group_labels <- rep(NA, ncol(dge$counts))
        for (group_name in names(group_info)) {
          group_indices <- which(sub(" \\(.*\\)", "", column_names) == group_name)
          group_labels[group_indices] <- group_name
        }
        dge <- DGEList(counts = dge$counts, group = factor(group_labels))
        print(dge)
        
        return(dge)
      } else {
        tkmessageBox(title = "Error", message = "Please fix the CSV file and try again.")
        return(NULL)
      }
    }
    dge <- run_analysis()
    table(dge$samples$group)
    dge <- calcNormFactors(dge)
    norm_factors <- dge$samples$norm.factors 
    lib_sizes <- dge$samples$lib.size         

    effective_lib_sizes <- lib_sizes * norm_factors
    tmm_normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE)
    write.csv(tmm_normalized_counts, file = "TMM_normalized_counts.csv", row.names = TRUE)
    cat("TMM-normalized counts saved to 'TMM_normalized_counts.csv'.")
    
    select_groups <- function(dge) {
      require(tcltk) 
      group_names <- levels(dge$samples$group)
      
   
      tt <- tktoplevel()
      tkwm.title(tt, "Select Baseline and Treatment Groups")
      baseline <- tclVar(group_names[1])
      treatment <- tclVar(if (length(group_names) > 1) group_names[2] else group_names[1])
      tkgrid(tklabel(tt, text = "Select Baseline Group:"))
      baseline_menu <- ttkcombobox(tt, values = group_names, textvariable = baseline, state = "readonly")
      tkgrid(baseline_menu)
      tkgrid(tklabel(tt, text = "Select Treatment Group:"))
      treatment_menu <- ttkcombobox(tt, values = group_names, textvariable = treatment, state = "readonly")
      tkgrid(treatment_menu)
      confirm_btn <- tkbutton(tt, text = "Confirm", command = function() {
        if (tclvalue(baseline) == tclvalue(treatment)) {
          tkmessageBox(title = "Error", message = "Baseline and Treatment groups must be different.")
        } else {
          tkdestroy(tt)
        }
      })
      tkgrid(confirm_btn)
      tkwait.window(tt)
      baseline_group_size <- sum(dge$samples$group == tclvalue(baseline))
      treatment_group_size <- sum(dge$samples$group == tclvalue(treatment))
      print(paste("Size of Baseline group (", tclvalue(baseline), "): ", baseline_group_size, sep = ""))
      print(paste("Size of Treated group (", tclvalue(treatment), "): ", treatment_group_size, sep = ""))
      return(list(baseline = tclvalue(baseline), treated = tclvalue(treatment)))
    }
    perform_glm_qlf <- function(dge, selected_groups) {
      dge_subset <- dge[, dge$samples$group %in% c(selected_groups$baseline, selected_groups$treated)]
      dge_subset$samples$group <- relevel(dge_subset$samples$group, ref = selected_groups$baseline)
      design <- model.matrix(~ group, data = dge_subset$samples)
      dge_subset <- estimateDisp(dge_subset, design)
      fit <- glmQLFit(dge_subset, design)
      qlf <- glmQLFTest(fit, coef = 2)  # coef = 2 compares the treated group against the baseline group
      top_results <- topTags(qlf, n = Inf)
      print(head(top_results$table))
      top_results$table$Bonferroni_pvalue <- p.adjust(top_results$table$PValue, method = "bonferroni")
      write.csv(top_results$table, file = "differential_expression_results.csv", row.names = TRUE)
      plotMD(qlf, main = "MA Plot", xlab = "Average Log CPM", ylab = "Log Fold Change")
      abline(h = c(-1, 1), col = "blue")
      return(top_results)
    }
    ask_comparisons <- function() {
      require(tcltk)
      tt <- tktoplevel()
      tkwm.title(tt, "Number of Comparisons")
      tkgrid(tklabel(tt, text = "How many comparisons would you like to perform?"))
      num_comparisons <- tclVar("1")  # Default value
      entry <- tkentry(tt, textvariable = num_comparisons, width = 10)
      tkgrid(entry)
      confirm_btn <- tkbutton(tt, text = "Confirm", command = function() {
        tkdestroy(tt)
      })
      tkgrid(confirm_btn)
      tkwait.window(tt)
      
      return(as.integer(tclvalue(num_comparisons)))
    }

    num_comparisons <- ask_comparisons()
    for (i in 1:num_comparisons) {
      cat("\nPerforming comparison", i, "of", num_comparisons, "\n")
      
      selected_groups <- select_groups(dge)  # Select baseline and treated groups
      results <- perform_glm_qlf(dge, selected_groups)
      
      filename <- paste0("EdgeR_", selected_groups$baseline, "_vs_", selected_groups$treated, ".csv")
      results_table <- results$table
      results_table <- data.frame(GeneSymbol = rownames(results_table), results_table)
      write.csv(results_table, file = filename, row.names = FALSE)
      cat("Differential expression results saved to:", filename, "\n")
    }
    print("Script completed successfully.")
    
  }, error = function(e) {
    print(sprintf("An error occurred: %s", e$message))
    
    library(DESeq2)
    library(tcltk)
    
    intro_message <- function() {
      tt <- tktoplevel()
      tkwm.title(tt, "RNA-Seq Analysis Software")
      tcl("wm", "attributes", tt, topmost = TRUE)
      msg <- "DESeq2 Analysis. Please make sure your CSV file is formatted accordingly."
      tkgrid(tklabel(tt, text = msg, wraplength = 600, justify = "left"), padx = 20, pady = 20)
      
      user_choice <- tclVar("continue")
      button_frame <- tkframe(tt)
      tkgrid(tkbutton(button_frame, text = "Continue", command = function() {
        tclvalue(user_choice) <- "continue"
        tkdestroy(tt)
      }), padx = 10, pady = 10)
      tkgrid(tkbutton(button_frame, text = "Exit and Fix", command = function() {
        tclvalue(user_choice) <- "exit"
        tkdestroy(tt)
      }), padx = 10, pady = 10)
      tkgrid(button_frame)
      tkwait.window(tt)
      if (tclvalue(user_choice) == "exit") return(FALSE)
      return(TRUE)
    }
    if (!intro_message()) stop("Execution stopped by user.")
    
    load_data <- function(filepath) {
      data <- read.csv(filepath, header = TRUE, check.names = FALSE)
      gene_symbols <- data[, 1]
      counts <- data[, -1]
      row_sums <- rowSums(counts)
      data_ordered <- data[order(row_sums, decreasing = TRUE), ]
      data_unique <- data_ordered[!duplicated(data_ordered[, 1]), ]
      counts <- data_unique[, -1]
      gene_symbols <- data_unique[, 1]
      rownames(counts) <- gene_symbols
      return(counts)
    }
    detect_groups <- function(column_names) {
      groups <- sub(" \\(.+\\)", "", column_names)
      return(groups)
    }
    confirm_groups <- function(column_names) {
      groups <- sub(" \\(.+\\)", "", column_names)
      group_counts <- table(groups)
      
      tt <- tktoplevel()
      tkwm.title(tt, "Confirm Groups")
      tkgrid(tklabel(tt, text = "Detected Groups and Sample Counts:"))
      
      for (g in names(group_counts)) {
        label_text <- sprintf("%s: %d samples", g, group_counts[[g]])
        tkgrid(tklabel(tt, text = label_text, justify = "left"), padx = 10, sticky = "w")
      }
      
      decision <- tclVar("no")
      button_frame <- tkframe(tt)
      
      tkgrid(tkbutton(button_frame, text = "Confirm", command = function() {
        tclvalue(decision) <- "yes"
        tkdestroy(tt)
      }), padx = 10)
      
      tkgrid(tkbutton(button_frame, text = "Decline", command = function() {
        tclvalue(decision) <- "no"
        tkdestroy(tt)
      }), padx = 10)
      
      tkgrid(button_frame, pady = 10)
      tkwait.window(tt)
      
      return(tclvalue(decision) == "yes")
    }
    select_groups <- function(groups) {
      unique_groups <- unique(groups)
      tt <- tktoplevel()
      tkwm.title(tt, "Select Baseline and Treatment")
      
      baseline <- tclVar(unique_groups[1])
      treatment <- tclVar(unique_groups[2])
      tkgrid(tklabel(tt, text = "Select Baseline:"))
      tkgrid(ttkcombobox(tt, values = unique_groups, textvariable = baseline))
      tkgrid(tklabel(tt, text = "Select Treatment:"))
      tkgrid(ttkcombobox(tt, values = unique_groups, textvariable = treatment))
      tkgrid(tkbutton(tt, text = "Confirm", command = function() {
        if (tclvalue(baseline) == tclvalue(treatment)) {
          tkmessageBox(message = "Baseline and Treatment must differ.")
        } else {
          tkdestroy(tt)
        }
      }))
      tkwait.window(tt)
      return(list(baseline = tclvalue(baseline), treated = tclvalue(treatment)))
    }
    run_deseq2 <- function(counts, groups, baseline, treated) {
      group_factor <- factor(groups)
      design_df <- data.frame(row.names = colnames(counts), group = group_factor)
      dds <- DESeqDataSetFromMatrix(countData = counts, colData = design_df, design = ~group)
      dds$group <- relevel(dds$group, ref = baseline)
      dds <- DESeq(dds)
      res <- results(dds, contrast = c("group", treated, baseline))
      resOrdered <- res[order(res$pvalue), ]
      resDF <- as.data.frame(resOrdered)
      resDF$GeneSymbol <- rownames(resDF)
      colnames(resDF)[colnames(resDF) == "log2FoldChange"] <- "logFC"
      colnames(resDF)[colnames(resDF) == "pvalue"] <- "PValue"
      colnames(resDF)[colnames(resDF) == "padj"] <- "FDR"
      resDF$Bonferroni_pvalue <- p.adjust(resDF$PValue, method = "bonferroni")
      resDF <- resDF[, c("GeneSymbol", setdiff(names(resDF), "GeneSymbol"))]
      write.csv(resDF, file = paste0("DESeq2_", baseline, "_vs_", treated, "_with_FDR_and_Bonferroni.csv"), row.names = FALSE)
    }
    ask_comparisons <- function() {
      tt <- tktoplevel()
      tkwm.title(tt, "Number of Comparisons")
      tkgrid(tklabel(tt, text = "How many comparisons?"))
      num_comparisons <- tclVar("1")
      tkgrid(tkentry(tt, textvariable = num_comparisons))
      tkgrid(tkbutton(tt, text = "Confirm", command = function() tkdestroy(tt)))
      tkwait.window(tt)
      return(as.integer(tclvalue(num_comparisons)))
    }
    filepath <- file.choose()
    counts <- load_data(filepath)
    column_names <- colnames(counts)
    groups <- detect_groups(column_names)
    if (!confirm_groups(column_names)) stop("Group confirmation declined.")
    num_comp <- ask_comparisons()
    for (i in 1:num_comp) {
      selected <- select_groups(groups)
      run_deseq2(counts, groups, selected$baseline, selected$treated)
    }
    print("Script completed successfully.")
  }, error = function(e) {
    print(sprintf("An error occurred: %s", e$message))
  })

}








main()
