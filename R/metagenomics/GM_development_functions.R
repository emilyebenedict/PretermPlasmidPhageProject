#' PERMANOVA for repeated measures using adonis2
#'
#' @param D An N-by-N distance matrix (must be a \code{dist} object)
#' @param subject A per-sample (length-N) character vector containing the subject identifiers
#' @param subject_data Data frame with per-subject metadata. Must have rownames matching unique(subject)
#' @param sample_data Data frame with per-sample metadata. Must have rownames matching D
#' @param metadata_order Character vector specifying the variables to include in the model
#' @param permutations Number of permutations to perform
#' @param ncores Number of cores for parallelization
#' @param by_type Type of test for adonis2 ("terms" or "margin")
#' @return An object similar to adonis2 output with empirical p-values
PERMANOVA_repeat_measures <- function(
    D,
    subject, subject_data = NULL,
    sample_data = NULL,
    metadata_order = c(names(subject_data), names(sample_data)),
    permutations = 999, ncores = 1,
    by_type = "terms") {
  
  if (!inherits(D, "dist")) {
    stop("D must be a dist object")
  }
  
  if (!is.character(subject)) stop("subject must be a character vector")
  if (length(subject) != nrow(as.matrix(D))) {
    stop("Length of subject must match number of samples in D")
  }
  if (!identical(names(subject), rownames(as.matrix(D)))) {
    stop("Names of subject vector must match sample names in D")
  }
  
  if (is.null(subject_data) & is.null(sample_data)) {
    stop("At least one of subject_data or sample_data must be provided")
  }
  
  if (is.null(subject_data)) {
    subject_data <- data.frame(.placeholder = 1, row.names = unique(subject))
  }
  
  if (length(unique(subject)) != nrow(subject_data)) {
    stop("Number of unique subjects must match rows in subject_data")
  }
  if (!setequal(unique(subject), rownames(subject_data))) {
    stop("Row names of subject_data must match unique values in subject")
  }
  
  if (is.null(sample_data)) {
    sample_data <- data.frame(.placeholder = 1, row.names = rownames(as.matrix(D)))
  }
  if (nrow(as.matrix(D)) != nrow(sample_data)) {
    stop("sample_data must have a row for each sample in D")
  }
  if (!identical(rownames(as.matrix(D)), rownames(sample_data))) {
    stop("Row names of sample_data must match row names of D")
  }
  
  if (length(intersect(names(subject_data), names(sample_data))) > 0) {
    stop("Column names must not overlap between subject_data and sample_data")
  }
  
  if (!all(metadata_order %in% c(names(subject_data), names(sample_data)))) {
    stop("metadata_order must only include columns from subject_data or sample_data")
  }
  
  mtdat <- cbind(subject_data[subject, , drop = FALSE], sample_data)
  
  if (any(is.na(mtdat[, metadata_order, drop = FALSE]))) {
    stop("Missing values found in model metadata")
  }
  
  # Construct formula explicitly
  formula <- as.formula(paste("D ~", paste(metadata_order, collapse = " + ")))
  ad <- vegan::adonis2(formula, permutations = 0, data = mtdat[, metadata_order, drop = FALSE], by = by_type)
  
  R2 <- ad$R2
  names(R2) <- rownames(ad)
  
  library(permute)
  doParallel::registerDoParallel(ncores)
  
  nullsamples <- foreach::`%dopar%`(
    foreach::foreach(i = seq_len(permutations), .combine = cbind),
    {
      subject.i <- sample(unique(subject))
      sample.i <- shuffle(nrow(sample_data), control = how(blocks = subject))
      
      i.subject_data <- subject_data[subject.i, , drop = FALSE]
      rownames(i.subject_data) <- rownames(subject_data)
      
      mtdat_perm <- cbind(i.subject_data[subject, , drop = FALSE],
                          sample_data[sample.i, , drop = FALSE])
      perm_ad <- vegan::adonis2(formula, permutations = 0, data = mtdat_perm[, metadata_order, drop = FALSE], by = by_type)
      r2_perm <- perm_ad$R2
      matrix(r2_perm, ncol = 1, dimnames = list(names(r2_perm), NULL))
    })
  
  doParallel::stopImplicitCluster()
  
  n <- length(R2)
  R2[n - 1] <- 1 - R2[n - 1]
  nullsamples[n - 1, ] <- 1 - nullsamples[n - 1, ]
  
  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)
  P[n] <- NA
  
  ad$`Pr(>F)` <- P
  return(ad)
}

PERMANOVA_repeat_measures <- function(
    D,
    subject, subject_data = NULL,
    sample_data = NULL,
    metadata_order = c(names(subject_data), names(sample_data)),
    permutations = 999, ncores = 1,
    by_type = "terms") {
  
  if (!inherits(D, "dist")) stop("D must be a dist object")
  if (!is.character(subject)) stop("subject must be a character vector")
  if (length(subject) != nrow(as.matrix(D))) {
    stop("Length of subject must match number of samples in D")
  }
  if (!identical(names(subject), rownames(as.matrix(D)))) {
    stop("Names of subject vector must match sample names in D")
  }

  if (is.null(subject_data) & is.null(sample_data)) {
    stop("At least one of subject_data or sample_data must be provided")
  }
  
  if (is.null(subject_data)) {
    subject_data <- data.frame(.placeholder = 1, row.names = unique(subject))
  }
  
  if (length(unique(subject)) != nrow(subject_data)) {
    stop("Number of unique subjects must match rows in subject_data")
  }
  if (!setequal(unique(subject), rownames(subject_data))) {
    stop("Row names of subject_data must match unique values in subject")
  }
  
  if (is.null(sample_data)) {
    sample_data <- data.frame(.placeholder = 1, row.names = rownames(as.matrix(D)))
  }
  if (nrow(as.matrix(D)) != nrow(sample_data)) {
    stop("sample_data must have a row for each sample in D")
  }
  if (!identical(rownames(as.matrix(D)), rownames(sample_data))) {
    stop("Row names of sample_data must match row names of D")
  }
  
  if (length(intersect(names(subject_data), names(sample_data))) > 0) {
    stop("Column names must not overlap between subject_data and sample_data")
  }
  
  if (!all(metadata_order %in% c(names(subject_data), names(sample_data)))) {
    stop("metadata_order must only include columns from subject_data or sample_data")
  }
  
  mtdat <- cbind(subject_data[subject, , drop = FALSE], sample_data)
  
  if (any(is.na(mtdat[, metadata_order, drop = FALSE]))) {
    stop("Missing values found in model metadata")
  }
  
  formula <- as.formula(paste("D ~", paste(metadata_order, collapse = " + ")))
  ad <- vegan::adonis2(formula, permutations = 0, data = mtdat[, metadata_order, drop = FALSE], by = by_type)
  
  R2 <- ad$R2
  names(R2) <- rownames(ad)
  
  library(permute)
  doParallel::registerDoParallel(ncores)

  nullsamples <- foreach::`%dopar%`(
    foreach::foreach(i = seq_len(permutations), .combine = cbind),
    {
      # Permute subjects and samples
      subject.i <- sample(unique(subject))
      sample.i <- shuffle(nrow(sample_data), control = how(blocks = subject))
      
      # Get permuted subject data
      i.subject_data <- subject_data[subject.i, , drop = FALSE]
      rownames(i.subject_data) <- rownames(subject_data)  # restore original rownames

      # Reconstruct metadata
      mtdat_perm <- tryCatch({
        cbind(i.subject_data[subject, , drop = FALSE],
              sample_data[sample.i, , drop = FALSE])
      }, error = function(e) {
        message("Error during cbind in permutation ", i, ": ", conditionMessage(e))
        return(NULL)  # Skip this iteration
      })

      if (is.null(mtdat_perm)) return(matrix(NA, nrow = length(R2), ncol = 1))

      # Run adonis2 on permuted metadata
      perm_ad <- vegan::adonis2(formula, permutations = 0, data = mtdat_perm[, metadata_order, drop = FALSE], by = by_type)
      r2_perm <- perm_ad$R2
      matrix(r2_perm, ncol = 1, dimnames = list(names(r2_perm), NULL))
    })

  doParallel::stopImplicitCluster()
  
  # Post-processing: residual adjustment
  n <- length(R2)
  R2[n - 1] <- 1 - R2[n - 1]
  nullsamples[n - 1, ] <- 1 - nullsamples[n - 1, ]

  # Handle possible NA rows (e.g. from skipped permutations)
  nullsamples <- nullsamples[, colSums(is.na(nullsamples)) == 0, drop = FALSE]

  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (ncol(nullsamples) + 1)
  P[n] <- NA  # No p-value for "Total"

  ad$`Pr(>F)` <- P
  return(ad)
}

# Function to extract and predict from a new phyloseq object
predict_microbiota_age <- function(ps_obj, model, top_taxa) {
  # Extract OTU table and normalize to relative abundance
  relabund <- otu_table(ps_obj)
  if (taxa_are_rows(ps_obj)) {
    relabund <- t(relabund)
  }
  #relabund <- sweep(relabund, 1, rowSums(relabund), "/")
  
  colnames(relabund) <- sapply(str_split(colnames(relabund), "\\|"), function(x) if(length(x) >= 7) x[7] else NA)
  
  # Convert to data frame and filter for top taxa
  df <- as.data.frame(relabund)
  df <- df[, intersect(colnames(df), top_taxa), drop = FALSE]
  
  # Add missing taxa as zeros
  for (taxon in setdiff(top_taxa, colnames(df))) {
    df[[taxon]] <- 0
  }
  df <- df[, top_taxa]  # ensure column order matches
  
  # Predict microbiota age
  predicted <- predict(model, newdata = df)
  
  # Combine with metadata
  metadata <- sample_data(ps_obj) %>% as.data.frame()
  metadata$Predicted_Microbiota_Age <- predicted
  
  return(metadata)
}

perform_wilcox_test <- function(df, facet_var, x_var, y_var, 
                                p.adjust.method = "BH",
                                label_column = c("p.adj.signif", "p.adj"),
                                y_pos_multiplier = 1.1, y_spacing = 0.1) {
  
  label_column <- match.arg(label_column)
  
  # Convert to factors for consistent ordering
  df[[facet_var]] <- as.factor(df[[facet_var]])
  df[[x_var]] <- as.factor(df[[x_var]])
  
  # Convert variable names to symbols
  facet_sym <- rlang::sym(facet_var)
  x_sym <- rlang::sym(x_var)
  y_sym <- rlang::sym(y_var)
  
  # ---- Run Wilcoxon tests per facet ----
  stat.test <- df %>%
    group_by(!!facet_sym) %>%
    rstatix::pairwise_wilcox_test(
      formula = as.formula(paste(y_var, "~", x_var)),
      p.adjust.method = p.adjust.method
    ) %>%
    ungroup()
  
  # ---- Add y positions ----
  y_max <- df %>%
    group_by(!!facet_sym) %>%
    summarise(y.max = max(!!y_sym, na.rm = TRUE), .groups = "drop")
  
  stat.test <- stat.test %>%
    left_join(y_max, by = facet_var) %>%
    group_by(!!facet_sym) %>%
    mutate(y.position = y.max * (y_pos_multiplier + y_spacing * row_number())) %>%
    ungroup()
  
  # ---- Add dynamic column matching x_var name ----
  stat.test <- stat.test %>%
    mutate(!!x_sym := group1)  # dynamically create a column named after x_var
  
  # ---- Reorder columns ----
  stat.test <- stat.test %>%
    relocate(!!x_sym, .after = group2) %>%
    relocate(y.position, .after = p.adj.signif)
  
  return(stat.test)
}

create_pairwise_beta_div_matrix <- function(ps_obj, distance_metric, save_matrix = FALSE, distance_metric_string, output_name){
  meta <- data.frame(phyloseq::sample_data(ps_obj))
  # Create distance matrix table
  ps.dist <- phyloseq::distance(ps_obj, method=distance_metric)
  ###Looking at pairwise  distances
  ps.dist.pairwise <- as.matrix(ps.dist)
  ps.dist.pairwise2 <- melt(ps.dist.pairwise)
  colnames(ps.dist.pairwise2)[1] <- "Sample_ID1"
  colnames(ps.dist.pairwise2)[2] <- "Sample_ID2"
  colnames(ps.dist.pairwise2)[3] <- distance_metric
  
  # Removing duplicated rows
  ps.dist.pairwise_uniq <- ps.dist.pairwise2[!duplicated(apply(ps.dist.pairwise2, 1, function(x) paste(sort(x), collapse = ""))),]
  # Removing self comparisons with values of zero
  ps.dist.pairwise_uniq <- ps.dist.pairwise_uniq[which(ps.dist.pairwise_uniq$Sample_ID1 != ps.dist.pairwise_uniq$Sample_ID2),] 
  # Merge final distance matrix with metadata
  ps.dist.pairwise_uniq <- merge(ps.dist.pairwise_uniq, meta, by.x = "Sample_ID1", by.y = "Sample_ID_R_final")
  ps.dist.pairwise_uniq <- merge(ps.dist.pairwise_uniq, meta, by.x = "Sample_ID2", by.y = "Sample_ID_R_final")
  
  ### Comparing within and between participant distances 
  # Adding comparison type column 
  ps.dist.pairwise_uniq$Comparison[ps.dist.pairwise_uniq$ID.x == ps.dist.pairwise_uniq$ID.y] <- "Within participant"
  ps.dist.pairwise_uniq$Comparison[ps.dist.pairwise_uniq$ID.x != ps.dist.pairwise_uniq$ID.y] <- "Between participant"
  
  if(save_matrix == TRUE) {
    write.csv(ps.dist.pairwise_uniq, paste0(output_path, deparse(substitute(ps_obj)), "_", distance_metric_string, "_",
                                            output_name, ".csv"), row.names = F)
  }
  
  ## Prepare for Wilcoxon ranksum tests between groups and participant pair groups
  ps.dist.pairwise_uniq$Comparison <- factor(ps.dist.pairwise_uniq$Comparison, levels = c("Between participant","Within participant"))
  
  return(ps.dist.pairwise_uniq)
  
}

perform_emmeans_lm <- function(df, facet_var, x_var, y_var,
                               fixed_effects = NULL,   # covariates in the model
                               p.adjust.method = "BH",
                               label_column = c("p.adj.signif", "p.adj"),
                               y_pos_multiplier = 1.1,
                               y_spacing = 0.1) {
  label_column <- match.arg(label_column)
  
  # Ensure factors for consistent ordering
  df[[facet_var]] <- as.factor(df[[facet_var]])
  df[[x_var]]     <- as.factor(df[[x_var]])
  
  # Symbols for tidy eval
  facet_sym <- rlang::sym(facet_var)
  x_sym     <- rlang::sym(x_var)
  y_sym     <- rlang::sym(y_var)
  
  # ---------- Build model formula y ~ x_var (+ fixed effects) ----------
  rhs_terms <- c(x_var)
  if (!is.null(fixed_effects)) {
    # fixed_effects can be a vector or a single string "cov1 + cov2"
    if (length(fixed_effects) == 1) {
      extra <- strsplit(fixed_effects, "[+,]")[[1]]
      extra <- trimws(extra)
      extra <- extra[extra != ""]
      rhs_terms <- c(rhs_terms, extra)
    } else {
      rhs_terms <- c(rhs_terms, fixed_effects)
    }
  }
  rhs_terms <- unique(rhs_terms)
  model_formula <- as.formula(paste(y_var, "~", paste(rhs_terms, collapse = " + ")))
  
  # ---------- Loop over facets and run lm + emmeans ----------
  facet_levels <- levels(df[[facet_var]])
  res_list <- vector("list", length(facet_levels))
  
  for (i in seq_along(facet_levels)) {
    lev <- facet_levels[i]
    sub_df <- df[df[[facet_var]] == lev, , drop = FALSE]
    
    # Need at least 2 levels of x_var within this facet
    if (dplyr::n_distinct(sub_df[[x_var]]) < 2) {
      res_list[[i]] <- NULL
      next
    }
    
    # Fit linear model
    fit <- stats::lm(model_formula, data = sub_df)
    
    # Estimated marginal means for x_var
    emm <- emmeans::emmeans(fit, specs = x_var)
    
    # Pairwise contrasts with p-value adjustment
    contr <- emmeans::contrast(emm, method = "pairwise")
    contr_sum <- summary(contr, adjust = p.adjust.method)
    
    contr_df <- as.data.frame(contr_sum)
    
    # contrast column is like "A - B"; split into group1, group2
    pairs <- strsplit(as.character(contr_df$contrast), " - ")
    contr_df$group1 <- vapply(pairs, `[`, character(1), 1)
    contr_df$group2 <- vapply(pairs, `[`, character(1), 2)
    
    # ---- NEW: strip outer parentheses, e.g. "(None Near-term)" -> "None Near-term" ----
    contr_df$group1 <- gsub("^\\((.*)\\)$", "\\1", contr_df$group1)
    contr_df$group2 <- gsub("^\\((.*)\\)$", "\\1", contr_df$group2)
    
    # Rename p.value -> p and duplicate to p.adj (already adjusted)
    contr_df <- contr_df %>%
      dplyr::rename(p = p.value) %>%
      dplyr::mutate(p.adj = p)
    
    # Add significance labels based on adjusted p
    contr_df <- contr_df %>%
      rstatix::add_significance("p.adj") %>%
      dplyr::rename(p.adj.signif = p.adj.signif)
    
    # Add .y. and facet column
    contr_df$.y. <- y_var
    contr_df[[facet_var]] <- lev
    
    res_list[[i]] <- contr_df
  }
  
  # Bind all facets
  stat.test <- dplyr::bind_rows(res_list)
  
  if (nrow(stat.test) == 0L) {
    warning("No facets had at least two levels of ", x_var,
            " — no pairwise emmeans contrasts were computed.")
    return(stat.test)
  }
  
  # ---- (Optional but helpful) drop rows with NA p.adj ----
  # stat.test <- stat.test %>%
  #   dplyr::filter(!is.na(p.adj))
  
  # ---------- Compute y positions for plotting ----------
  y_max <- df %>%
    dplyr::group_by(!!facet_sym) %>%
    dplyr::summarise(y.max = max(!!y_sym, na.rm = TRUE), .groups = "drop")
  
  stat.test <- stat.test %>%
    dplyr::left_join(y_max, by = facet_var) %>%
    dplyr::group_by(!!facet_sym) %>%
    dplyr::mutate(
      y.position = y.max * (y_pos_multiplier + y_spacing * dplyr::row_number())
    ) %>%
    dplyr::ungroup()
  
  # ---------- Dynamic x column matching x_var name ----------
  stat.test <- stat.test %>%
    dplyr::mutate(!!x_sym := group1)
  
  # Make sure this column uses the same levels as the original factor
  stat.test[[x_var]] <- factor(stat.test[[x_var]], levels = levels(df[[x_var]]))
  
  # ---------- Tidy column order ----------
  stat.test <- stat.test %>%
    dplyr::relocate(!!x_sym, .after = group2) %>%
    dplyr::relocate(y.position, .after = p.adj.signif)
  
  return(stat.test)
}

perform_emmeans_lmer <- function(df, facet_var, x_var, y_var,
                                 fixed_effects = NULL,          # covariates in the model
                                 random_effects = "(1|Patient)",# random effects (string or vector)
                                 p.adjust.method = "BH",
                                 label_column = c("p.adj.signif", "p.adj"),
                                 y_pos_multiplier = 1.1,
                                 y_spacing = 0.1,
                                 lmer_method = c("REML", "ML"),
                                 df_method = c("kenward-roger", "satterthwaite", "asymptotic")) {
  label_column <- match.arg(label_column)
  lmer_method  <- match.arg(lmer_method)
  df_method    <- match.arg(df_method)
  
  # Ensure factors for consistent ordering
  df[[facet_var]] <- as.factor(df[[facet_var]])
  df[[x_var]]     <- as.factor(df[[x_var]])
  
  # Symbols for tidy eval
  facet_sym <- rlang::sym(facet_var)
  x_sym     <- rlang::sym(x_var)
  y_sym     <- rlang::sym(y_var)
  
  # ---------- Build fixed-effect RHS: x_var (+ fixed effects) ----------
  rhs_terms <- c(x_var)
  if (!is.null(fixed_effects)) {
    if (length(fixed_effects) == 1) {
      extra <- strsplit(fixed_effects, "[+,]")[[1]]
      extra <- trimws(extra)
      extra <- extra[extra != ""]
      rhs_terms <- c(rhs_terms, extra)
    } else {
      rhs_terms <- c(rhs_terms, fixed_effects)
    }
  }
  rhs_terms <- unique(rhs_terms)
  
  # ---------- Build random effects string ----------
  if (is.null(random_effects) || length(random_effects) == 0) {
    stop("random_effects must be provided for lmer(), e.g. '(1|Patient)'.")
  }
  rand_str <- if (length(random_effects) == 1) {
    random_effects
  } else {
    paste(random_effects, collapse = " + ")
  }
  rand_str <- trimws(rand_str)
  
  # ---------- Full model formula: y ~ fixed + random ----------
  model_formula <- stats::as.formula(
    paste(y_var, "~", paste(rhs_terms, collapse = " + "), "+", rand_str)
  )
  
  # ---------- Loop over facets and run lmer + emmeans ----------
  facet_levels <- levels(df[[facet_var]])
  res_list <- vector("list", length(facet_levels))
  
  for (i in seq_along(facet_levels)) {
    lev <- facet_levels[i]
    sub_df <- df[df[[facet_var]] == lev, , drop = FALSE]
    
    # Need at least 2 levels of x_var within this facet
    if (dplyr::n_distinct(sub_df[[x_var]]) < 2) {
      res_list[[i]] <- NULL
      next
    }
    
    # Fit mixed model (guard against failures)
    fit <- tryCatch(
      lme4::lmer(model_formula, data = sub_df, REML = (lmer_method == "REML")),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      res_list[[i]] <- NULL
      next
    }
    
    # Estimated marginal means for x_var
    # NOTE: df_method requires lmerTest and/or pbkrtest depending on choice.
    emm <- emmeans::emmeans(fit, specs = x_var, lmer.df = df_method)
    
    # Pairwise contrasts with p-value adjustment
    contr <- emmeans::contrast(emm, method = "pairwise")
    contr_sum <- summary(contr, adjust = p.adjust.method)
    
    contr_df <- as.data.frame(contr_sum)
    
    # contrast column like "A - B"; split into group1, group2
    pairs <- strsplit(as.character(contr_df$contrast), " - ")
    contr_df$group1 <- vapply(pairs, `[`, character(1), 1)
    contr_df$group2 <- vapply(pairs, `[`, character(1), 2)
    
    # strip outer parentheses, e.g. "(None Near-term)" -> "None Near-term"
    contr_df$group1 <- gsub("^\\((.*)\\)$", "\\1", contr_df$group1)
    contr_df$group2 <- gsub("^\\((.*)\\)$", "\\1", contr_df$group2)
    
    # Rename p.value -> p and duplicate to p.adj (already adjusted)
    contr_df <- contr_df %>%
      dplyr::rename(p = p.value) %>%
      dplyr::mutate(p.adj = p)
    
    # Add significance labels based on adjusted p
    contr_df <- contr_df %>%
      rstatix::add_significance("p.adj") %>%
      dplyr::rename(p.adj.signif = p.adj.signif)
    
    # Add .y. and facet column
    contr_df$.y. <- y_var
    contr_df[[facet_var]] <- lev
    
    res_list[[i]] <- contr_df
  }
  
  # Bind all facets
  stat.test <- dplyr::bind_rows(res_list)
  
  if (nrow(stat.test) == 0L) {
    warning("No facets had at least two levels of ", x_var,
            " (or model fits failed) — no pairwise emmeans contrasts were computed.")
    return(stat.test)
  }
  
  # ---------- Compute y positions for plotting ----------
  y_max <- df %>%
    dplyr::group_by(!!facet_sym) %>%
    dplyr::summarise(y.max = max(!!y_sym, na.rm = TRUE), .groups = "drop")
  
  stat.test <- stat.test %>%
    dplyr::left_join(y_max, by = facet_var) %>%
    dplyr::group_by(!!facet_sym) %>%
    dplyr::mutate(
      y.position = y.max * (y_pos_multiplier + y_spacing * dplyr::row_number())
    ) %>%
    dplyr::ungroup()
  
  # ---------- Dynamic x column matching x_var name ----------
  stat.test <- stat.test %>%
    dplyr::mutate(!!x_sym := group1)
  
  # Make sure this column uses the same levels as the original factor
  stat.test[[x_var]] <- factor(stat.test[[x_var]], levels = levels(df[[x_var]]))
  
  # ---------- Tidy column order ----------
  stat.test <- stat.test %>%
    dplyr::relocate(!!x_sym, .after = group2) %>%
    dplyr::relocate(y.position, .after = p.adj.signif)
  
  return(stat.test)
}

plot_pcoa <- function(ps_obj, group_var = "abx_group", custom_colors = NULL, adonisformula = "ps.dist ~ Gestational_age_weeks + abx_group_names", 
                      output_dir = "DOL_bin_PCoA_plots") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Extract metadata and group levels
  meta <- data.frame(sample_data(ps_obj))
  
  # Perform ordination
  ord <- ordinate(ps_obj, method = "PCoA", distance = "bray")
  
  # Plot ordination
  p <- plot_ordination(ps_obj, ord, type = "samples", color = group_var) +
    geom_point(shape = 19, size = 3) +
    stat_ellipse() +
    scale_color_manual(values = custom_colors) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "left",
      plot.title = element_text(size = 16, hjust = 0.5),
      panel.border = element_rect(fill = NA, colour = "black", size = 1)
    ) +
    labs(title = paste(levels(ps_obj@sam_data$DOL_bin)[1]))
  
  # Add marginal boxplots
  p_with_marginals <- ggMarginal(p, type = "boxplot", groupFill = TRUE, size = 10)
  
  # Save the plot
  file_name <- paste0("PCoA_", levels(ps_obj@sam_data$DOL_bin)[1], ".pdf")
  ggsave(file.path(output_dir, file_name), plot = p_with_marginals, bg = "transparent", width = 9, height = 6)
  message("Saved: ", file_name)
  
  ## PERMANOVA
  ps.dist <- phyloseq::distance(ps_obj, method="bray") 
  
  ps.adonis <- vegan::adonis2(formula(adonisformula),  
                              data=meta,
                              na=na.omit,
                              permutations = 10000,
                              by = "terms",
                              subset=complete.cases(ord))  
  adonis_name <- paste0("PCoA_", levels(ps_obj@sam_data$DOL_bin)[1], "PERMANOVA.csv")
  write.csv(ps.adonis, file.path(output_dir, adonis_name))


}