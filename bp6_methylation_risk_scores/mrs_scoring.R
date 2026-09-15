# Read in beta matrix, coverage matrix, cell fractions, and chronologic age. 
beta_mat <- as.matrix(read.table("/path/to/Illumina/ref/beta.tsv", sep = "\t", header = TRUE, row.names = 1))
coverage_mat <- read.table("/path/to/Illumina/ref/coverage.tsv", sep = "\t", header = TRUE, row.names = 1)
cell_fractions <- read.table("/path/to/cell/fractions.tsv", sep = "\t", header = TRUE, row.names = 1)
chrono_age <- readRDS("/path/to/chronologic/age.rds")

# EWAS summary statistics, must label the columns the following: probe, weight, trait, threshold 
score_weights <- fread("/path/to/EWAS/summary/stats.csv")

# Placeholder — fill in real trait/threshold definitions before running.
score_defs <- data.table(
  trait = "X",
  threshold = "X"
)

# Defining helper functions
beta_to_m <- function(x, epsilon = 1e-6) {
  x <- pmin(pmax(x, epsilon), 1 - epsilon)
  log2(x / (1 - x))
}

m_to_beta <- function(x) {
  2^x / (1 + 2^x)
}

# Regress out the given covariates (e.g. age, cell fractions) from each CpG's methylation values, returning covariate-adjusted beta values.
adjust_methylation <- function(beta_mat, covariates) {
  stopifnot(
    identical(colnames(beta_mat), rownames(covariates))
  )
  m_mat <- beta_to_m(beta_mat)

  adjusted <- t(
    apply(m_mat, 1, function(y) {
      complete <- complete.cases(y, covariates)
      result <- rep(NA_real_, length(y))

      # Require enough complete samples to fit the model
      if (sum(complete) > ncol(covariates) + 2) {
        fit <- lm(
          y[complete] ~ .,
          data = covariates[complete, , drop = FALSE]
        )

        # Residual + original mean preserves the CpG methylation scale
        result[complete] <- residuals(fit) + mean(y[complete])
      }

      result
    })
  )

  dimnames(adjusted) <- dimnames(beta_mat)
  m_to_beta(adjusted)
}

# Score qc function
score_qc_function <- function(observed,
                              expected,
                              min_prop_large,
                              require_all_if_n_le) {
  if (expected <= require_all_if_n_le) {
    pass_qc <- observed == expected
  } else {
    pass_qc <- observed / expected >= min_prop_large
  }
  pass_qc
}

# Calculate weighted score 
calculate_weighted_scores <- function(beta_mat,
                                      score_weights,
                                      score_defs,
                                      min_prop_large,
                                      require_all_if_n_le,
                                      coverage_mat = NULL,
                                      qc_function = NULL) {
  results <- data.table(
    sample = colnames(beta_mat)
  )

  qc_records <- list()

  for (i in seq_len(nrow(score_defs))) {
    current_trait <- score_defs$trait[i]
    current_definition <- score_defs$threshold[i]

    weights <- unique(
      score_weights[
        trait == current_trait &
          threshold == current_definition,
        .(probe, weight)
      ],
      by = "probe"
    )

    available_probes <- intersect(
      weights$probe,
      rownames(beta_mat)
    )

    score_name <- paste(
      current_trait,
      current_definition,
      sep = "_"
    )

# No score probes available
    if (length(available_probes) == 0) {
      results[, (score_name) := NA_real_]

      qc_records[[score_name]] <- data.table(
        sample = colnames(beta_mat),
        trait = current_trait,
        definition = current_definition,
        n_expected = nrow(weights),
        n_observed = 0L,
        pct_observed = 0,
        pass_qc = FALSE
      )

      next
    }

  # Extract methylation values
    methylation <- beta_mat[
      available_probes,
      drop = FALSE
    ]

    weights <- weights[
      match(available_probes, probe)
    ]

# Count CpGs that are actually observed
    observed_probes <- colSums(
      !is.na(methylation)
    )

    expected_probes <- length(
      score_weights[
        trait == current_trait &
          threshold == current_definition,
        unique(probe)
      ]
    )

    pct_observed <- observed_probes / expected_probes

# Score level QC
    pass_qc <- if (is.null(qc_function)) {
      rep(TRUE, ncol(methylation))
    } else {
      qc_function(
        observed = observed_probes,
        expected = expected_probes,
        min_prop_large = min_prop_large,
        require_all_if_n_le = require_all_if_n_le
      )
    }

# Calculate weighted score
    score <- as.numeric(
      crossprod(
        weights$weight,
        methylation
      )
    )

    # Scores failing QC are not retained
    score[!pass_qc] <- NA_real_

    results[, (score_name) := score]

# Record QC
    qc_records[[score_name]] <- data.table(
      sample = colnames(beta_mat),
      trait = current_trait,
      definition = current_definition,
      n_expected = expected_probes,
      n_observed = as.integer(observed_probes),
      pct_observed = pct_observed,
      pass_qc = pass_qc
    )
  }

  list(
    scores = results,
    qc = rbindlist(
      qc_records,
      fill = TRUE
    )
  )
}

# Scoring
score_probes <- intersect(
  rownames(beta_mat),
  unique(score_weights$probe)
)

score_beta <- beta_mat[
  score_probes,
  ,
  drop = FALSE
]

score_coverage <- if (is.null(coverage_mat)) {
  NULL
} else {
  coverage_mat[
    score_probes,
    ,
    drop = FALSE
  ]
}

# Ensure age covariate order matches beta matrix order
age_covariates <- data.frame(
  Chrono_age = chrono_age[colnames(beta_mat)],
  row.names = colnames(beta_mat)
)
head(age_covariates)

# Remove one cell-fraction column because estimated fractions sum to 1. Reorder cell_fractions to match beta_mat's sample order by name first.
cell_covariates <- cell_fractions[
  colnames(beta_mat),
  -ncol(cell_fractions),
  drop = FALSE
]

age_cell_covariates <- cbind(
  age_covariates,
  cell_covariates
)

# Generate three methylation matrices
beta_unadjusted <- score_beta

beta_age_adjusted <- adjust_methylation(
  score_beta,
  age_covariates
)

beta_age_cell_adjusted <- adjust_methylation(
  score_beta,
  age_cell_covariates
)

# Calculating MRS
min_prop_large <- 0.2       # require ≥80% of expected CpGs for larger scores
require_all_if_n_le <- 5    # scores with ≤5 CpGs must have all of them observed

raw_scores <- calculate_weighted_scores(
  beta_mat = beta_unadjusted,
  coverage_mat = score_coverage,
  score_weights = score_weights,
  score_defs = score_defs,
  qc_function = score_qc_function,
  min_prop_large = min_prop_large,
  require_all_if_n_le = require_all_if_n_le
)
 
 
age_adjusted_scores <- calculate_weighted_scores(
  beta_mat = beta_age_adjusted,
  coverage_mat = score_coverage,
  score_weights = score_weights,
  score_defs = score_defs,
  qc_function = score_qc_function,
  min_prop_large = min_prop_large,
  require_all_if_n_le = require_all_if_n_le
)
 
 
age_cell_adjusted_scores <- calculate_weighted_scores(
  beta_mat = beta_age_cell_adjusted,
  coverage_mat = score_coverage,
  score_weights = score_weights,
  score_defs = score_defs,
  qc_function = score_qc_function,
  min_prop_large = min_prop_large,
  require_all_if_n_le = require_all_if_n_le
)

# Merging final socres 
final_scores <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "sample",
      all = TRUE
    )
  },
  list(
    data.table(
      sample = colnames(beta_mat)
    ),
    data.table(
      sample = names(chrono_age),
      Chrono_age = chrono_age
    ),
    data.table(
      sample = rownames(cell_fractions),
      cell_fractions
    ),
    raw_scores$scores,
    age_adjusted_scores$scores,
    age_cell_adjusted_scores$scores
  )
)
 
# Identify score columns
score_columns <- setdiff(
  names(final_scores),
  c(
    "sample",
    "Chrono_age",
    colnames(cell_fractions)
  )
)

# Standardize MRS values
final_scores[
  ,
  paste0(score_columns, "_z") :=
    lapply(
      .SD,
      function(x) {
        as.numeric(scale(x))
      }
    ),
  .SDcols = score_columns
]

# Exporting results
# Saving final scores
write(
  final_scores,
  "methylation_risk_scores.tsv",
  sep = "\t"
)
# Saving 
score_qc <- rbindlist(
  list(
 
    raw_scores$qc[
      ,
      adjustment := "unadjusted"
    ],
 
    age_adjusted_scores$qc[
      ,
      adjustment := "age_adjusted"
    ],
 
    age_cell_adjusted_scores$qc[
      ,
      adjustment := "age_cell_adjusted"
    ]
  ),
  fill = TRUE
)
fwrite(
  score_qc,
  "methylation_risk_score_QC.tsv",
  sep = "\t"
)
