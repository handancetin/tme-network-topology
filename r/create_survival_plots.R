# Load libraries
library(GEOquery)
library(survival)
library(survminer)
library(Biobase) 
library(biomaRt)

# Create directory for saving data outputs for Python
if (!dir.exists("output_data")) {
  dir.create("output_data")
  cat("Created directory: output_data\n")
}

# Function to map ENSEMBL IDs to probe IDs
map_ensembl_to_probes <- function(ensembl_ids, probe_map_file) {
  # Read the probe to gene mapping file
  probe_gene_map <- read.delim(probe_map_file, header = FALSE, 
                               col.names = c("probe_id", "gene_symbol"),
                               stringsAsFactors = FALSE)
  
  # Connect to Ensembl BioMart
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Get gene symbols for the ENSEMBL IDs
  gene_data <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = ensembl_ids,
                     mart = ensembl)
  
  # Merge with probe mapping file
  results <- list()
  for (i in 1:nrow(gene_data)) {
    ensembl_id <- gene_data$ensembl_gene_id[i]
    gene_symbol <- gene_data$hgnc_symbol[i]
    
    # Find all probe IDs matching this gene symbol
    matching_probes <- probe_gene_map$probe_id[probe_gene_map$gene_symbol == gene_symbol]
    
    if (length(matching_probes) > 0) {
      results[[ensembl_id]] <- data.frame(
        ensembl_id = ensembl_id,
        gene_symbol = gene_symbol,
        probe_id = matching_probes,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all results
  if (length(results) > 0) {
    combined_results <- do.call(rbind, results)
    return(combined_results)
  } else {
    return(data.frame(ensembl_id = character(0), 
                      gene_symbol = character(0), 
                      probe_id = character(0)))
  }
}

# Convert ENSGs to Affy IDs
genes_data <- read.delim("genes.tsv", header = TRUE, stringsAsFactors = FALSE)
genes_list <- unique(sub("\\..*", "", genes_data$genes))
genes_list <- c(genes_list, "ENSG00000078098") # Add FAP for control
genes_information <- map_ensembl_to_probes(genes_list, "probe_sets_to_genes.txt")
genes_probe_ids <- unique(genes_information$probe_id)

# Save genes_information for Python
write.csv(genes_information, "output_data/genes_information.csv", row.names = FALSE)

# Create a data frame to store all results
all_results <- data.frame(
  ensembl_id = character(),
  gene_symbol = character(),
  probe_id = character(),
  os_cutpoint = numeric(),
  os_hr = numeric(),
  os_pval = numeric(),
  os_high_worse = logical(),
  dfs_cutpoint = numeric(),
  dfs_hr = numeric(),
  dfs_pval = numeric(),
  dfs_high_worse = logical(),
  stringsAsFactors = FALSE
)

# Download and load the GSE17536 dataset
cat("\nDownloading GSE17536 dataset. This may take some time...\n")
gset <- getGEO("GSE17536", GSEMatrix = TRUE, getGPL = FALSE)
gset <- gset[[1]]

# Get expression and phenotype data
exprs_data <- exprs(gset)
pheno_data <- pData(gset)

# Save expression data for Python
write.csv(exprs_data, "output_data/expression_data.csv")

# Extract survival information from characteristics columns
overall_survival <- as.numeric(pheno_data[, "overall survival follow-up time:ch1"])

# Convert text event values to numeric (1 = event/death, 0 = no event/alive)
event_values <- pheno_data[, "overall_event (death from any cause):ch1"]
overall_survival_event <- ifelse(event_values == "death", 1, 0)

# Extract disease-free survival information
dfs_time <- as.numeric(pheno_data[, "dfs_time:ch1"])
dfs_event_values <- pheno_data[, "dfs_event (disease free survival; cancer recurrence):ch1"]

# Convert event values to numeric (1 = event/recurrence, 0 = no event)
dfs_event <- ifelse(dfs_event_values == "recurrence" | dfs_event_values == "1", 1, 0)

# Save survival data for Python
survival_data <- data.frame(
  sample_id = rownames(pheno_data),
  os_time = overall_survival,
  os_event = overall_survival_event,
  dfs_time = dfs_time,
  dfs_event = dfs_event
)
write.csv(survival_data, "output_data/survival_data.csv", row.names = FALSE)

# Create directory for saving plots if it doesn't exist
if (!dir.exists("survival_plots")) {
  dir.create("survival_plots")
  cat("Created directory: survival_plots\n")
}

# Loop through each gene and its probe IDs
for (i in 1:nrow(genes_information)) {
  ensembl_id <- genes_information$ensembl_id[i]
  gene_symbol <- genes_information$gene_symbol[i]
  probe_id <- genes_information$probe_id[i]
  
  # Skip if the probe ID is not found in the expression data
  if (!(probe_id %in% rownames(exprs_data))) {
    cat(paste("Probe ID", probe_id, "not found in expression data. Skipping...\n"))
    next
  }
  
  # Extract expression values for the current probe
  gene_expr <- exprs_data[rownames(exprs_data) == probe_id, ]
  
  # Create a data frame for overall survival analysis
  surv_data <- data.frame(
    sample = colnames(exprs_data),
    time = overall_survival,
    event = overall_survival_event,
    gene_expr = gene_expr
  )
  
  # Remove rows with missing values
  surv_data <- surv_data[complete.cases(surv_data), ]
  
  # Print some summary statistics
  cat(paste("\nProcessing gene:", gene_symbol, "- ENSEMBL ID:", ensembl_id, "- Probe ID:", probe_id, "\n"))
  cat("Expression summary:\n")
  print(summary(surv_data$gene_expr))
  
  # Generate survival analysis and plots
  tryCatch({
    # Find optimal cutpoint for gene expression in OS analysis
    cutpoint <- surv_cutpoint(surv_data, time = "time", event = "event", variables = "gene_expr")
    cat("Cutpoint information:\n")
    print(cutpoint$cutpoint)
    
    # Check which samples have higher and lower expression
    high_count <- sum(surv_data$gene_expr > cutpoint$cutpoint$cutpoint)
    low_count <- sum(surv_data$gene_expr <= cutpoint$cutpoint$cutpoint)
    cat("Based on actual expression values:\n")
    cat("Samples with high expression (>", cutpoint$cutpoint$cutpoint, "):", high_count, "\n")
    cat("Samples with low expression (<=", cutpoint$cutpoint$cutpoint, "):", low_count, "\n")
    
    # Categorize samples based on the cutpoint
    surv_cat <- surv_categorize(cutpoint)
    surv_cat$gene_expr <- as.factor(surv_cat$gene_expr)
    
    # Check the counts after categorization
    cat("After surv_categorize:\n")
    print(table(surv_cat$gene_expr))
    
    # Create the survival fit
    fit <- survfit(Surv(time, event) ~ gene_expr, data = surv_cat)
    
    # Get summary of survival in each group to determine which performs worse
    fit_summary <- summary(fit)
    cat("Median survival times:\n")
    print(fit_summary$table)
    
    # Run a Cox model
    cox_model <- coxph(Surv(time, event) ~ gene_expr, data = surv_cat)
    
    # Get the coefficient to determine direction of effect
    cox_coef <- coef(cox_model)[1]
    hr <- exp(cox_coef)
    pval <- summary(cox_model)$waldtest["pvalue"]
    
    # Check which group is the reference level
    ref_levels <- levels(surv_cat$gene_expr)
    cat("Factor levels:", paste(ref_levels, collapse=", "), "\n")
    if (length(ref_levels) > 0) {
      ref_level <- ref_levels[1]
      cat("Reference level in Cox model:", ref_level, "\n")
    } else {
      cat("Warning: gene_expr has no factor levels\n")
      ref_level <- NULL
    }
    cat("Cox coefficient:", cox_coef, "\n")
    cat("Hazard ratio:", hr, "\n")
    cat("P-value:", pval, "\n")
    
    # Determine if high expression is worse based on median survival and HR
    # This is tricky due to how surv_categorize labels groups
    # We look at both the categorization and the actual survival outcomes
    
    # If high group has more events/shorter survival, then high is worse
    is_high_worse <- FALSE
    
    # Check raw data from fit summary to determine which group has worse survival
    high_group_median <- fit_summary$table["gene_expr=high", "median"]
    if (!is.na(high_group_median)) {
      # If high group has a defined median survival (some patients died)
      low_group_median <- fit_summary$table["gene_expr=low", "median"]
      if (is.na(low_group_median)) {
        # If low group doesn't have median (better survival), then high is worse
        is_high_worse <- TRUE
      } else {
        # Both groups have median survival, compare them
        is_high_worse <- high_group_median < low_group_median
      }
    } else {
      # High group doesn't have median survival, check if low group does
      low_group_median <- fit_summary$table["gene_expr=low", "median"]
      if (!is.na(low_group_median)) {
        # Low group has worse survival
        is_high_worse <- FALSE
      } else {
        # Neither has median, use hazard ratio direction
        # Check if "high" category has worse outcomes by comparing coefficient
        if (!is.null(ref_level)) {
          if (ref_level == "high") {
            # If high is reference, negative coef means other group (low) has better outcomes (high is worse)
            is_high_worse <- cox_coef < 0
          } else {
            # If low is reference, positive coef means high has worse outcomes
            is_high_worse <- cox_coef > 0
          }
        } else {
          # Default to using coefficient sign when reference level is unclear
          is_high_worse <- cox_coef > 0
        }
      }
    }
    
    cat("High expression is worse:", is_high_worse, "\n")
    
    # For DFS analysis
    dfs_data <- data.frame(
      sample = colnames(exprs_data),
      time = dfs_time,
      event = dfs_event,
      gene_expr = gene_expr
    )
    
    # Remove rows with missing values
    dfs_data <- dfs_data[complete.cases(dfs_data), ]
    
    # Find optimal cutpoint for DFS
    dfs_cutpoint <- surv_cutpoint(dfs_data, time = "time", event = "event", variables = "gene_expr")
    
    # Categorize DFS data
    dfs_cat <- surv_categorize(dfs_cutpoint)
    dfs_cat$gene_expr <- as.factor(dfs_cat$gene_expr)
    
    # Create DFS fit
    dfs_fit <- survfit(Surv(time, event) ~ gene_expr, data = dfs_cat)
    
    # Run Cox model for DFS
    dfs_cox_model <- coxph(Surv(time, event) ~ gene_expr, data = dfs_cat)
    dfs_coef <- coef(dfs_cox_model)[1]
    dfs_hr <- exp(dfs_coef)
    dfs_pval <- summary(dfs_cox_model)$waldtest["pvalue"]
    
    # Check DFS reference level
    dfs_ref_levels <- levels(dfs_cat$gene_expr)
    cat("DFS factor levels:", paste(dfs_ref_levels, collapse=", "), "\n")
    if (length(dfs_ref_levels) > 0) {
      dfs_ref_level <- dfs_ref_levels[1]
      cat("DFS reference level:", dfs_ref_level, "\n")
    } else {
      cat("Warning: DFS gene_expr has no factor levels\n")
      dfs_ref_level <- NULL
    }
    cat("DFS Cox coefficient:", dfs_coef, "\n")
    cat("DFS Hazard ratio:", dfs_hr, "\n")
    cat("DFS P-value:", dfs_pval, "\n")
    
    # Determine if high expression is worse for DFS
    dfs_is_high_worse <- FALSE
    
    # Check if "high" category has worse outcomes for DFS
    dfs_fit_summary <- summary(dfs_fit)
    
    # Similar logic as for OS, but for DFS
    if (!is.null(dfs_ref_level)) {
      if (dfs_ref_level == "high") {
        # If high is reference, negative coef means low has better outcomes (high is worse)
        dfs_is_high_worse <- dfs_coef < 0
      } else {
        # If low is reference, positive coef means high has worse outcomes
        dfs_is_high_worse <- dfs_coef > 0
      }
    } else {
      # Default to using coefficient sign
      dfs_is_high_worse <- dfs_coef > 0
    }
    
    cat("DFS: High expression is worse:", dfs_is_high_worse, "\n")
    
    # Store the results for this gene
    gene_result <- data.frame(
      ensembl_id = ensembl_id,
      gene_symbol = gene_symbol,
      probe_id = probe_id,
      os_cutpoint = cutpoint$cutpoint$cutpoint,
      os_hr = hr,
      os_pval = pval,
      os_high_worse = is_high_worse,
      dfs_cutpoint = dfs_cutpoint$cutpoint$cutpoint,
      dfs_hr = dfs_hr,
      dfs_pval = dfs_pval,
      dfs_high_worse = dfs_is_high_worse,
      stringsAsFactors = FALSE
    )
    
    all_results <- rbind(all_results, gene_result)
    
    # Save individual gene analysis data for Python plotting
    base_file_name <- paste0("output_data/", gene_symbol, "_", probe_id, "_", ensembl_id)
    
    # Save OS categorized data
    surv_cat$cutpoint <- cutpoint$cutpoint$cutpoint
    write.csv(surv_cat, paste0(base_file_name, "_os_data.csv"), row.names = FALSE)
    
    # Save DFS categorized data
    dfs_cat$cutpoint <- dfs_cutpoint$cutpoint$cutpoint
    write.csv(dfs_cat, paste0(base_file_name, "_dfs_data.csv"), row.names = FALSE)
    
    # Only save if significant and high expression is associated with worse OS
    if (pval < 0.005 && dfs_pval < 0.005 && is_high_worse) {
      cat("Saving plot: p-value < 0.005 and high expression is associated with worse survival\n")
      
      # For plotting, we need to define the color palette
      palette_colors <- c("#E7B800", "#2E9FDF")  # yellow, blue
      
      names(fit$strata) <- tolower(gsub("gene_expr=", "", names(fit$strata)))
      names(dfs_fit$strata) <- tolower(gsub("gene_expr=", "", names(dfs_fit$strata)))
      palette_colors <- c("high" = "#E7B800", "low" = "#2E9FDF")
      
      
      os_plot <- ggsurvplot(
        fit,
        data = surv_cat,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.height = 0.25,
        palette = palette_colors,
        title = paste("Overall Survival -", gene_symbol, "- HR:", round(hr, 2)),
        xlab = "Time (months)",
        ylab = "Overall Survival Probability",
        ggtheme = theme_bw(),    
        font.main = c(16, "bold"),
        font.x = c(14),
        font.y = c(14),
        font.tickslab = c(12)
      )
      
      # Plot DFS
      dfs_plot <- ggsurvplot(
        dfs_fit,
        data = dfs_cat,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.height = 0.25,
        palette = palette_colors,
        title = paste("Disease-Free Survival -", gene_symbol, "- HR:", round(dfs_hr, 2)),
        xlab = "Time (months)",
        ylab = "Disease-Free Survival Probability",
        surv.median.line = "hv",
        ggtheme = theme_bw(),     
        font.main = c(16, "bold"),
        font.x = c(14),
        font.y = c(14),
        font.tickslab = c(12)
      )
      
      # Arrange plots
      combined_plots <- arrange_ggsurvplots(
        list(OS = os_plot, DFS = dfs_plot),
        print = FALSE,
        ncol = 2, nrow = 1,
        title = paste("Survival Analysis by", gene_symbol, " (", probe_id, ") ", "Expression (HR:", round(hr, 2), "p =", format(pval, digits = 3), ")")
      )
      
      # Save plot
      file_name <- paste0("survival_plots/", gene_symbol, "_", "HR_", round(hr, 3), "_",  probe_id, "_", ensembl_id, ".png")
      pdf(file_name, width = 14, height = 8)
      print(combined_plots)
      dev.off()
      
      cat(paste("Saved plot to:", file_name, "\n"))
    } else {
      cat("Skipping plot: p-value or high expression criteria not met\n")
    }
    
  }, error = function(e) {
    cat(paste("Error processing gene:", gene_symbol, "- Probe ID:", probe_id, "-", e$message, "\n"))
  })
}

# Save all results to a CSV file
write.csv(all_results, "output_data/all_gene_results.csv", row.names = FALSE)

cat("\nAll genes processed. Results saved to output_data/all_gene_results.csv\n")