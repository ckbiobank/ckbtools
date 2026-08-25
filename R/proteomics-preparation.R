#' Check CKB data release
#'
#' Checks the data release number is suitable.
#' Use at the start of proteomics data preparation functions.
#'
#' @keywords internal
#' @noRd
#'
check_data_release <- function() {
  data_release <- getOption("ckb.data.release")
  if (is.null(data_release)) {
    cli::cli_abort(c("CKB data release not set.",
                     "i" = "Use {.fn options} to set {.field ckb.data.release} to your data release number."),
                   call = rlang::caller_env())
  }

  if (!data_release %in% c("19.02", "19.03")) {
    cli::cli_abort(c("x" = "This function is not validated for CKB data release {data_release}."),
                   call = rlang::caller_env())
  }

  return(invisible(data_release))
}


#' Preparation of Olink meta data
#'
#' @param data Raw CKB Olink meta data (imported from olink_explore_meta.rds or
#' olink_explore_meta.csv)
#' @param raw_olink_protein_data Raw Olink protein data (imported from
#' data_baseline_olink_explore.rds or data_baseline_olink_explore.csv)
#' @param names Naming scheme to use for proteins: "protein" (default) to use
#' assay names or "id" to use Olink IDs.
#' @param columns Optional columns to include in the meta data.
#' (Default: c("uniprot", "olinkid"))
#'
#' @export
prepare_olink_meta_data <- function(data,
                                    raw_olink_protein_data,
                                    names = c("protein", "id"),
                                    columns = c("uniprot",
                                                "olinkid")) {

  # check data release
  check_data_release()

  # check arugments
  names <- rlang::arg_match(names)

  # Obtain Olink ID from the raw protein data set
  olinkid_data <- raw_olink_protein_data %>%
    dplyr::distinct(.data$assay, .data$panel_full, .data$olinkid) %>%
    dplyr::mutate(panel = factor(dplyr::case_when(
      .data$panel_full == "Cardiometabolic" ~ "Cardiometabolic_I",
      .data$panel_full == "Inflammation" ~ "Inflammation_I",
      .data$panel_full == "Neurology" ~ "Neurology_I",
      .data$panel_full == "Oncology" ~ "Oncology_I",
      TRUE ~ .data$panel_full
    )))

  # Obtain Olink normalization type from raw protein data set
  normalization_data <- raw_olink_protein_data %>%
    dplyr::select("csid", "olinkid", "normalization") %>%
    dplyr::group_by(.data$olinkid) %>%
    dplyr::summarise(normalization = factor(paste(unique(.data$normalization), collapse = "")))

  ## Expect that each protein should only have one of two possible
  ## normalization types
  if (!all(levels(normalization_data$normalization) %in% c("Intensity", "Plate control"))) {
    cli::cli_abort("A protein assay appears to have unexpected or multiple normalization types. There is a problem with the raw data.")
  }


  # Add batch identifier for explore panel to match expand panel
  data <- data %>%
    dplyr::mutate(panel = factor(dplyr::case_when(
      .data$panel == "Cardiometabolic" ~ "Cardiometabolic_I",
      .data$panel == "Inflammation" ~ "Inflammation_I",
      .data$panel == "Neurology" ~ "Neurology_I",
      .data$panel == "Oncology" ~ "Oncology_I",
      .data$panel == "Cardiometabolic II" ~ "Cardiometabolic_II",
      .data$panel == "Inflammation II" ~ "Inflammation_II",
      .data$panel == "Neurology II" ~ "Neurology_II",
      .data$panel == "Oncology II" ~ "Oncology_II",
      TRUE ~ .data$panel
    )))

  # Fix some assay names to match protein data
  data <- data %>%
    dplyr::mutate(assay = dplyr::case_when(
      .data$assay == "BOLA2" ~ "BOLA2_BOLA2B",
      .data$assay == "MICA_MICB" ~ "MICB_MICA",
      .data$assay == "CERT1" ~ "CERT",
      .data$assay == "WARS1" ~ "WARS",
      .data$assay == "DEFB103A" ~ "DEFB103A_DEFB103B",
      .data$assay == "CGB3" ~ "CGB3_CGB5_CGB8",
      .data$assay == "AMY1A" ~ "AMY1A_AMY1B_AMY1C",
      .data$assay == "CTAG1A" ~ "CTAG1A_CTAG1B",
      .data$assay == "SPACA5" ~ "SPACA5_SPACA5B",
      .data$assay == "NT-proBNP" ~ "NTproBNP",
      TRUE ~ .data$assay
    )) %>%
    dplyr::filter(!.data$assay %in% c("TNFSF9", "TOM1L2"))

  # Join the Olink ID and normalization data
  data <- data %>%
    dplyr::left_join(olinkid_data, by = c("assay", "panel")) %>%
    dplyr::left_join(normalization_data, by = "olinkid")

  # Create protein name
  if (names == "protein") {
    data <- data %>%
      dplyr::add_count(.data$assay) %>%
      dplyr::mutate(name = dplyr::if_else(.data$n > 1,
                                          paste0(.data$assay,
                                                 "_",
                                                 .data$panel),
                                          .data$assay))
  } else if (names == "id") {
    data$name <- data$olinkid
  }

  # Fix uniprot format when there are multiple
  if ("uniprot" %in% columns) {
    data$uniprot <- gsub("_", "|", data$uniprot)
  }

  # Tidy data frame and select columns
  data <- data %>%
    dplyr::relocate(tidyselect::all_of("name"), .before = 1) %>%
    dplyr::select("name",
                  "assay",
                  "panel",
                  "normalization",
                  tidyselect::any_of(columns))

  return(data)
}


#' Preparation of Olink protein data
#'
#' @param data Raw Olink protein data (imported from
#' data_baseline_olink_explore.rds or data_baseline_olink_explore.csv)
#' @param meta_data Olink meta data produced by [prepare_olink_meta_data()].
#' @param transform Transformation to apply to protein values.
#' Currently only "none" is available.
#' @param qc
#' QC steps to apply:
#' \itemize{
#'    \item "default" (default): replace protein values with a QC warning with
#'    a missing value
#'    \item "none"
#'    }
#' @param columns Other columns in raw data to include in output. For example,
#' set to "lod" to include the limit of detection for each protein.
#'
#' @export
prepare_olink_data <- function(data,
                               meta_data,
                               transform = c("none"),
                               qc = c("default", "none"),
                               columns = NULL) {
  # argument verification
  transform <- rlang::arg_match(transform)
  qc <- rlang::arg_match(qc)

  # check data release
  check_data_release()

  # Add batch identifier for explore panel to match expand panel
  data <- data %>%
    dplyr::mutate(panel_full = dplyr::case_when(
      .data$panel_full == "Cardiometabolic" ~ "Cardiometabolic_I",
      .data$panel_full == "Inflammation" ~ "Inflammation_I",
      .data$panel_full == "Neurology" ~ "Neurology_I",
      .data$panel_full == "Oncology" ~ "Oncology_I",
      TRUE ~ .data$panel_full
    ))

  # Default QC: Exclude a protein value if `qc_warning == 1`
  if (qc == "default") {
    data$npx[data$qc_warning == 1] <- NA_real_
    cli::cli_alert_success("Protein values with QC warnings have been replaced with missing values. {.emph [Default]}",
                           wrap = TRUE)
  } else if (qc == "none") {
    cli::cli_alert_success("Protein values with QC warnings have {.strong not} been removed.",
                           wrap = TRUE)
  }


  data$name <- meta_data$name[match(data$olinkid, meta_data$olinkid)]

  # Convert to wide format
  wide_data <- data %>%
    dplyr::filter(!is.na(.data$name)) %>%
    dplyr::select("csid", "name", "npx") %>%
    tidyr::pivot_wider(
      names_from = "name",
      values_from = "npx"
    )

  other_wide_data <- NULL
  if (!is.null(columns)) {
    other_wide_data <- data %>%
      dplyr::select("csid", "name", tidyselect::all_of(columns)) %>%
      tidyr::pivot_wider(
        names_from = "name",
        names_glue = "{.value}_{name}",
        values_from = tidyselect::all_of(columns)
      )
  }

  # Assign plate id for each participant
  data_long_plateid_I <- data %>%
    dplyr::select("csid", "batch", "plateid") %>%
    dplyr::filter(.data$batch == "explore") %>%
    dplyr::mutate(plateid_I = .data$plateid)

  data_plateid_I <- data_long_plateid_I[!duplicated(data_long_plateid_I[c("csid")]), ]
  data_plateid_I <- data_plateid_I[, c("csid", "plateid_I")]

  data_long_plateid_II <- data %>%
    dplyr::select("csid", "batch", "plateid") %>%
    dplyr::filter(.data$batch == "expand") %>%
    dplyr::mutate(plateid_II = .data$plateid)

  data_plateid_II <- data_long_plateid_II[!duplicated(data_long_plateid_II[c("csid")]), ]
  data_plateid_II <- data_plateid_II[, c("csid", "plateid_II")]

  plates <- data %>%
    dplyr::group_by(.data$csid) %>%
    dplyr::summarise(plateid = paste(sort(unique(.data$plateid)),
                                     collapse = " ")) %>%
    dplyr::mutate(plateid = factor(.data$plateid))

  data <- list(wide_data,
               other_wide_data,
               data_plateid_I,
               data_plateid_II,
               plates) %>%
    purrr::discard(is.null) %>%
    purrr::reduce(dplyr::inner_join, by = "csid") %>%
    dplyr::relocate("plateid_I", "plateid_II", "plateid",
                    .after = "csid")

  expected_plates <- getOption("olink_expected_plates") %||% 46
  if (length(levels(plates$plateid)) != expected_plates) {
    cli::cli_warn("We expect {expected_plates} plates, but this data includes {length(levels(plates$plateid))} plates. Please check you are using correct data.")
  }

  data <- data %>%
    dplyr::rename_with(function(x) gsub("^npx_", "", x))


  return(data)
}




#' Preparation of SomaScan meta data
#'
#' @param data Raw CKB SomaScan meta data (imported from somalogic_meta.rds or
#' somalogic_meta.csv)
#' @param names Naming scheme to use for proteins: "protein" (default) to use
#' target name or "id" to use SomaLogic sequence IDs.
#' @param columns Optional columns to include in the meta data.
#' (Default: c("seqid", "somaid", "target_full_name", "target", "uniprot",
#' "organism", "type", "dilution"))
#' @param filter Filter to apply to rows:
#' \itemize{
#'  \item "HumanProtein" (default): Only return meta data for aptamers where
#'  organism is "Human" and type is "Protein".
#'  \item "All": Return all rows of meta data.
#' }
#'
#' @export
prepare_somascan_meta_data <- function(data,
                                       names = c("protein", "id"),
                                       columns = c("seqid", "somaid",
                                                   "target_full_name", "target",
                                                   "uniprot",
                                                   "organism", "type",
                                                   "dilution"),
                                       filter = "HumanProtein") {

  names <- rlang::arg_match(names)
  filter <- rlang::arg_match(filter, values = c("HumanProtein", "All"))

  # check data release
  check_data_release()

  # Create protein name
  if (names == "protein") {
    data <- data %>%
      dplyr::add_count(.data$target) %>%
      dplyr::mutate(name = dplyr::if_else(.data$n > 1,
                                          paste0(.data$target,
                                                 "_",
                                                 .data$seqid),
                                          .data$target)) %>%
      dplyr::mutate(name = gsub("\\.|\\-|\\s|\\:", "_", .data$name))

  } else if (names == "id") {
    data$name <- data$aptname
  }

  # filter
  if (filter == "HumanProtein") {
    data <- dplyr::filter(data,
                          .data$organism == "Human",
                          .data$type == "Protein")
  }

  #
  data <- dplyr::select(data, "name", "aptname", tidyselect::any_of(columns))

  return(data)
}




#' Prepare SomaScan protein data
#'
#' This functions prepares a data set of SomaScan protein data
#' given the raw CKB SomaScan data.
#'
#' @details
#' # QC acceptance criteria
#' Hybridization control scale factor (`hybcontrolnormscale`) and median signal
#' normalization scale factors (`normscale_20`, `normscale_0_5`,
#' `normscale_0_005`) have an acceptance range of 0.4-2.5.
#'
#' Percentage of measurements used for ANML scale factor calculation
#' (`anmlfractionused_20`, `anmlfractionused_0_5`, `anmlfractionused_0_005`) has
#' a criteria of > 30% (0.3).
#'
#' Samples that meet these criteria are given a PASS for normalization
#' acceptance criteria for all row scale factors (`rowcheck`), while samples
#' that are outside the criteria are given a FLAG.
#'
#' Hybridization control scale factor out of range suggests a technical issue.
#' Set `qc = "default"` to exclude participant samples with hybridization
#' control scale factor (`hybcontrolnormscale`) out of range (i.e. less than 0.4
#' or greater than 2.5).
#'
#' Other criteria out of range suggests that the sample is fairly different from
#' a typical population.
#' Set `qc = "flagged"` to exclude participant samples with a FLAG for
#' normalization acceptance criteria for all row scale factors (`rowcheck`).
#'
#' @param data Raw CKB SomaScan protein data (imported from
#' data_baseline_somalogic.rds or data_baseline_somalogic.csv)
#' @param sample_data Raw CKB SomaScan protein sample information data (imported
#' from data_baseline_somalogic_samples.rds or
#' data_baseline_somalogic_samples.csv)
#' @param meta_data SomaScan meta data produced by
#' [prepare_somascan_meta_data()].
#' @param transform Transformation to apply to protein values:
#' \itemize{
#'  \item "log2" (default): log base 2
#'  \item "log10": log base 10
#'  \item "none"
#' }
#' @param qc
#' QC steps to apply.
#' \itemize{
#'  \item "default" (default): Exclude participant samples with
#'  hybridization control scale factor out of range.
#'  \item "flagged": Exclude participants samples flagged for any normalization
#'  acceptance criteria.
#'  \item "none": Do not exclude any participant samples.
#' }
#' @param repeat_samples Default "remove" will remove data for the second sample
#' for participants who had a sample run twice. Use "retain" to keep this data;
#' the second samples will have sample_source=1.
#' @param columns Names of any additional columns from sample_data to include in
#' the output.
#'
#' @export
prepare_somascan_data <- function(data,
                                  sample_data,
                                  meta_data,
                                  transform = c("log2", "log", "log10", "none"),
                                  qc = c("default", "flagged", "none"),
                                  repeat_samples = c("remove", "retain"),
                                  columns = NULL) {

  # argument verification
  transform <- rlang::arg_match(transform)
  qc <- rlang::arg_match(qc)
  repeat_samples <- rlang::arg_match(repeat_samples)

  # check data release
  check_data_release()

  data <- data[, c("csid", "sample_source", "variable", "value")]

  # Transformations
  if (transform == "log2") {
    data$value <- log2(data$value)
    cli::cli_alert_success("RFU values have been log transformed (base 2). {.emph [Default]}",
                           wrap = TRUE)
  } else if (transform == "log") {
    data$value <- log(data$value)
    cli::cli_alert_success("RFU values have been (natural) log transformed.",
                           wrap = TRUE)
  } else if (transform == "log10") {
    data$value <- log10(data$value)
    cli::cli_alert_success("RFU values have been log transformed (base 10).",
                           wrap = TRUE)
  }

  data$name <- meta_data$name[match(data$variable, meta_data$aptname)]
  # Pivot to wide data
  wide_data <- data %>%
    dplyr::filter(!is.na(.data$name)) %>%
    dplyr::select("csid", "sample_source", "name", "value") %>%
    tidyr::pivot_wider(
      names_from = "name",
      values_from = "value"
    )

  # QC
  keep_samples <- sample_data %>%
    ## hybcontrolnormscale is a character vector in the raw data
    ## convert to numeric
    dplyr::mutate(hybcontrolnormscale = as.numeric(.data$hybcontrolnormscale))

  if (qc == "default") {
    ## exclude participant samples with hybcontrolnormscale
    ## out of range (i.e. less than 0.4 or greater than 2.5)
    keep_samples <- keep_samples %>%
      dplyr::filter(.data$hybcontrolnormscale >= 0.4,
                    .data$hybcontrolnormscale <= 2.5)
    cli::cli_alert_success("Participant samples with hybcontrolnormscale out of range (i.e. less than 0.4 or greater than 2.5) have been removed. {.emph [Default]}",
                           wrap = TRUE)
  } else if (qc == "flagged") {
    keep_samples <- keep_samples %>%
      dplyr::filter(.data$rowcheck == "PASS")
    cli::cli_alert_success("Participant samples flagged for any normalization acceptance criteria have been removed.",
                           wrap = TRUE)
  } else if (qc == "none") {
    cli::cli_alert_success("Participant samples with QC flags have {.strong not} been removed.",
                           wrap = TRUE)
  }

  # Select columns and merge data frames
  data <- keep_samples %>%
    dplyr::select("csid",
                  "sample_source",
                  "plateid",
                  tidyselect::any_of(columns)) %>%
    dplyr::mutate(plateid = factor(.data$plateid)) %>%
    dplyr::left_join(wide_data, by = c("csid", "sample_source"))

  # Handle repeat samples
  if (repeat_samples == "remove") {
    data <- data %>%
      dplyr::filter(.data$sample_source == 0) %>%
      dplyr::select(-"sample_source")

    cli::cli_alert_success("Data for the second sample (for participants who had a sample run twice) have been removed. {.emph [Default]}",
                           wrap = TRUE)
  }

  return(data)
}


#' Preparation of a proteomics analysis data set
#'
#'
#' @param participant_data A data frame of participant data
#' @param protein_data Protein data produced by [prepare_olink_data()] or
#' [prepare_somascan_data()].
#' @param protein_meta_data Corresponding protein meta data produed by
#' [prepare_olink_meta_data()] or [prepare_somascan_meta_data()].
#' @param plate_correction Type of plate correction to apply:
#' \itemize{
#'  \item "none" (default)
#'  \item "median-subcohort": For each protein, subtract the plate median (of
#'  subcohort values) and add the overall median (of subcohort values).
#'  \item "mean-subcohort": For each protein, subtract the plate mean (of
#'  subcohort values) and add the overall mean (of subcohort values).
#' }
#' @param scaling Scaling to apply to values for each protein:
#' \itemize{
#'  \item "none" (default)
#'  \item "center-subcohort": Subtract the overall mean of subcohort values.
#'  \item "scale-subcohort": Subtract the overall mean and divide by the sample
#'  standard deviation of subcohort values.
#'  \item "center": Subtract the overall mean.
#'  \item "scale": Subtract the overall mean and divide by the sample standard
#'  deviation.
#'  \item "rint": Rank inverse normal transformation.
#' }
#' @param id Name of column used to join participant_data and protein_data.
#' (Default: "csid")
#'
#' @export
prepare_proteomics_analysis_dataset <- function(
  participant_data,
  protein_data,
  protein_meta_data,
  plate_correction = c("none",
                       "median-subcohort",
                       "mean-subcohort"),
  scaling = c("none",
              "center-subcohort",
              "scale-subcohort",
              "center",
              "scale",
              "rint"),
  id = "csid"
) {
  plate_correction <- rlang::arg_match(plate_correction)
  scaling <- rlang::arg_match(scaling)

  if (any(grepl("subcohort", c(plate_correction, scaling))) &&
      !"subcohort" %in% names(participant_data)) {
    cli::cli_abort("There must be a {.field subcohort} column in the {.arg participant_data} data frame.")
  }

  protein_names <- protein_meta_data$name

  data <- participant_data %>%
    dplyr::inner_join(protein_data, by = id)

  # Plate correction
  if (plate_correction == "median-subcohort") {

    proteins_for_plate_correction <- protein_names
    proteins_no_plate_correction <- NULL

    ## For Olink, only apply this correction to assays
    ## that had "Intensity" normalization applied
    if ("normalization" %in% names(protein_meta_data)) {
      proteins_for_plate_correction <- protein_names[protein_meta_data$normalization == "Intensity"]
      proteins_no_plate_correction <- setdiff(protein_names, proteins_for_plate_correction)
    }

    overall_medians <- data %>%
      dplyr::summarise(dplyr::across(
        tidyselect::all_of(proteins_for_plate_correction),
        function(x) stats::median(x[.data$subcohort == 1], na.rm = TRUE)
      )) %>%
      as.list()

    data <- data %>%
      dplyr::group_by(.data$plateid) %>%
      dplyr::mutate(dplyr::across(
        tidyselect::all_of(proteins_for_plate_correction),
        function(x) x - stats::median(x[.data$subcohort == 1], na.rm = TRUE)
      )) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(
        tidyselect::all_of(proteins_for_plate_correction),
        function(x) x + overall_medians[[dplyr::cur_column()]]
      ))

    cli::cli_alert_success("Plate correction 'median-subcohort' has been applied to protein values.",
                           wrap = TRUE)
    if (length(proteins_no_plate_correction) > 0) {
      cli::cli_alert_info("{proteins_no_plate_correction} had no plate correction applied because they were not originally 'Intensity' normalised.",
                          wrap = TRUE)
    }


  } else if (plate_correction == "mean-subcohort") {
    proteins_for_plate_correction <- protein_names

    overall_means <- data %>%
      dplyr::summarise(dplyr::across(
        tidyselect::all_of(proteins_for_plate_correction),
        function(x) mean(x[.data$subcohort == 1], na.rm = TRUE)
      )) %>%
      as.list()

    data <- data %>%
      dplyr::group_by(.data$plateid) %>%
      dplyr::mutate(dplyr::across(
        tidyselect::all_of(proteins_for_plate_correction),
        function(x) x - mean(x[.data$subcohort == 1], na.rm = TRUE)
      )) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(
        tidyselect::all_of(proteins_for_plate_correction),
        function(x) x + overall_means[[dplyr::cur_column()]]
      ))

    cli::cli_alert_success("Plate correction 'mean-subcohort' has been applied to protein values.",
                           wrap = TRUE)
  } else if (plate_correction == "none") {
    cli::cli_alert_success("No plate correction has been applied to protein values. {.emph [Default]}",
                           wrap = TRUE)
  }

  # Scaling
  if (scaling == "none") {
    cli::cli_alert_success("No scaling of protein values has been applied. {.emph [Default]}",
                           wrap = TRUE)
    return(data)
  } else if (scaling %in% c("rint", "center", "scale")) {
    scaling_func <- switch(
      scaling,
      "rint" = function(x) rint(x),
      "center" = function(x) as.numeric(scale(x, center = TRUE, scale = FALSE)),
      "scale" = function(x) as.numeric(scale(x, center = TRUE, scale = TRUE))
    )
    data <- data %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(protein_names), scaling_func))

    switch(
      scaling,
      "rint" = cli::cli_alert_success("Rank inverse normal transformation has been applied to protein values.", wrap = TRUE),
      "center" = cli::cli_alert_success("Protein values have been centred using overall means.", wrap = TRUE),
      "scale" = cli::cli_alert_success("Protein values have been centred and scaled using overall means and standard deviations.", wrap = TRUE)
    )

  } else if (scaling == "center-subcohort") {
    data <- data %>%
      dplyr::mutate(dplyr::across(
        tidyselect::all_of(protein_names),
        function(x) (x - mean(x[.data$subcohort == 1], na.rm = TRUE))
      ))
    cli::cli_alert_success("Protein values have been centred using subcohort means.", wrap = TRUE)
  } else if (scaling == "scale-subcohort") {
    data <- data %>%
      dplyr::mutate(dplyr::across(
        tidyselect::all_of(protein_names),
        function(x) (x - mean(x[.data$subcohort == 1], na.rm = TRUE)) / stats::sd(x[.data$subcohort == 1], na.rm = TRUE)
      ))
    cli::cli_alert_success("Protein values have been centred and scaled using subcohort means and standard deviations.",
                           wrap = TRUE)
  }

  return(data)
}
