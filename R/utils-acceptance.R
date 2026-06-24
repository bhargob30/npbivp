format_acceptance <- function(acc_raw, m_fit) {
  acc_raw <- as.numeric(acc_raw)

  raw_names <- c(
    "lambda0_intercept_precision",
    paste0("lambda", seq_len(m_fit), "_common_scale"),
    paste0("weight", seq_len(m_fit), "_transformed"),
    "error_precision"
  )

  names(acc_raw) <- raw_names

  active_names <- c(
    "lambda0_intercept_precision",
    "lambda_common_scale",
    if (m_fit > 1) paste0("weight", seq_len(m_fit - 1), "_transformed"),
    "error_precision"
  )

  active_values <- c(
    acc_raw["lambda0_intercept_precision"],
    acc_raw["lambda1_common_scale"],
    if (m_fit > 1) acc_raw[paste0("weight", seq_len(m_fit - 1), "_transformed")],
    acc_raw["error_precision"]
  )

  data.frame(
    parameter = active_names,
    acceptance_rate = as.numeric(active_values),
    row.names = NULL
  ) |>
    list(active = _, raw = acc_raw)
}
