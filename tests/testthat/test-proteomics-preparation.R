
withr::local_options(list(ckb.data.release = "19.02",
                          olink_expected_plates = 26))

# --- check_data_release -------------------------------------------------------
test_that("check_data_release enforces presence and specific value", {
  withr::with_options(
    list(ckb.data.release = NULL),
    expect_error(check_data_release(), "ckb.data.release")
  )
  withr::with_options(
    list(ckb.data.release = "99.99"),
    expect_error(check_data_release(), "not validated")
  )
  withr::with_options(
    list(ckb.data.release = "19.02"),
    expect_invisible(expect_equal(check_data_release(), "19.02"))
  )
})




# --- prepare_olink_meta_data --------------------------------------------------
test_that("prepare_olink_meta_data: joins olinkid + normalization and builds names", {
  olink_meta_data <- prepare_olink_meta_data(
    data = ckbtools_raw_olink_meta_data,
    raw_olink_protein_data = ckbtools_raw_olink_protein_data,
    names = "protein"
  )

  expect_s3_class(olink_meta_data, "data.frame")
  expect_true(all(c("name", "assay", "panel", "normalization") %in% names(olink_meta_data)))
  expect_false(any(is.na(olink_meta_data$name)))
  expect_false(any(duplicated(olink_meta_data$name)))
  expect_true(all(levels(olink_meta_data$normalization) %in% c("Intensity", "Plate control")))
  expect_true(all(olink_meta_data$panel %in% c("Cardiometabolic_I",
                                               "Cardiometabolic_II",
                                               "Inflammation_I",
                                               "Inflammation_II",
                                               "Neurology_I",
                                               "Neurology_II",
                                               "Oncology_I",
                                               "Oncology_II")))

  # names="id" path uses olinkid as name
  olink_meta_data_id <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data, ckbtools_raw_olink_protein_data, names = "id")
  expect_true(all(olink_meta_data_id$name == olink_meta_data_id$olinkid))
  expect_true(all(grepl("^OID", olink_meta_data_id$name)))
})

test_that("prepare_olink_meta_data errors on inconsistent normalization", {
  bad_raw_olink_protein_data <- ckbtools_raw_olink_protein_data
  bad_raw_olink_protein_data$normalization[1] <- "Plate control"
  expect_error(prepare_olink_meta_data(ckbtools_raw_olink_meta_data,
                                       bad_raw_olink_protein_data,
                                       names = "id"),
               "A protein assay appears to have unexpected or multiple normalization types. There is a problem with the raw data.")
})

test_that("prepare_olink_data QC='default' increases NAs vs QC='none'", {
  olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data, ckbtools_raw_olink_protein_data)

  out_none <- suppressMessages(prepare_olink_data(ckbtools_raw_olink_protein_data, olink_meta_data, qc = "none"))
  out_def  <- suppressMessages(prepare_olink_data(ckbtools_raw_olink_protein_data, olink_meta_data, qc = "default"))

  expect_gt(
    sum(is.na(out_def[olink_meta_data$name])),
    sum(is.na(out_none[olink_meta_data$name]))
  )
})

test_that("prepare_olink_data includes requested extra columns when present", {
  olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data, ckbtools_raw_olink_protein_data)
  out <- suppressMessages(
    prepare_olink_data(ckbtools_raw_olink_protein_data,
                       olink_meta_data,
                       qc = "none",
                       columns = "lod")
  )
  expect_all_true(paste0("lod_", olink_meta_data$name) %in% names(out))
})


# --- prepare_olink_data -------------------------------------------------------
test_that("prepare_olink_data: QC, widening, and plate columns", {
  olink_meta_data <- prepare_olink_meta_data(
    data = ckbtools_raw_olink_meta_data,
    raw_olink_protein_data = ckbtools_raw_olink_protein_data
  )

  suppressMessages(
    olink_data <- prepare_olink_data(
      data = ckbtools_raw_olink_protein_data,
      meta_data = olink_meta_data,
      qc = "default"  # should set NPX to NA where qc_warning==1
    )
  )

  # structure
  expect_true(all(c("csid", "plateid", "plateid_I", "plateid_II") %in% names(olink_data)))
  expect_gt(ncol(olink_data), 10)  # wide format with protein columns
  expect_true(any(olink_meta_data$name %in% names(olink_data)))

  # spot-check QC: if any qc_warning==1 for a protein, those values become NA
  one <- ckbtools_raw_olink_protein_data$olinkid[1]
  nm  <- olink_meta_data$name[match(one, olink_meta_data$olinkid)]
  flagged_ids <- ckbtools_raw_olink_protein_data$csid[ckbtools_raw_olink_protein_data$olinkid == one & ckbtools_raw_olink_protein_data$qc_warning == 1]
  if (length(flagged_ids) > 0 && nm %in% names(olink_data)) {
    expect_true(all(is.na(olink_data[olink_data$csid %in% flagged_ids, nm, drop = TRUE])))
  }
})




# --- prepare_somascan_meta_data ----------------------------------------------
test_that("prepare_somascan_meta_data: naming + optional filtering", {
  # names=protein + default filter
  somascan_meta_data <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data,
                                                   names = "protein",
                                                   filter = "HumanProtein")
  expect_true(all(c("name", "aptname") %in% names(somascan_meta_data)))
  expect_true(all(somascan_meta_data$aptname %in% ckbtools_raw_somascan_meta_data$aptname))
  expect_true(all(somascan_meta_data$aptname %in% ckbtools_raw_somascan_meta_data$aptname[ckbtools_raw_somascan_meta_data$organism == "Human" & ckbtools_raw_somascan_meta_data$type == "Protein"]))

  # names=id + include all
  somascan_meta_data_id <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data, names = "id", filter = "All")
  expect_true(all(somascan_meta_data_id$name == somascan_meta_data_id$aptname))
})

test_that("prepare_somascan_meta_data filtering All vs HumanProtein", {
  somascan_meta_data_all <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data, names = "protein", filter = "All")
  somascan_meta_data_hp  <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data, names = "protein", filter = "HumanProtein")

  expect_gte(nrow(somascan_meta_data_all), nrow(somascan_meta_data_hp))
  expect_all_true(somascan_meta_data_hp$aptname %in% somascan_meta_data_all$aptname)
})


# --- prepare_somascan_data ----------------------------------------------------
test_that("prepare_somascan_data: transforms, QC, and repeat handling", {
  somascan_meta_data  <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data, names = "protein", filter = "HumanProtein")

  # default: log2 + QC(default) + remove repeats
  suppressMessages(
    somascan_data <- prepare_somascan_data(
      data = ckbtools_raw_somascan_protein_data,
      sample_data = ckbtools_raw_somascan_sample_data,
      meta_data = somascan_meta_data,
      trans = "log2",
      qc = "default",
      repeat_samples = "remove"
    )
  )

  expect_true(all(c("csid", "plateid") %in% names(somascan_data)))
  expect_false("sample_source" %in% names(somascan_data))   # removed when repeat_samples="remove"
  prot_cols <- setdiff(names(somascan_data), c("csid", "plateid"))
  expect_true(any(prot_cols %in% somascan_meta_data$name))

  # keep repeats + no transform + no QC removal
  suppressMessages(
    somascan_data_keep <- prepare_somascan_data(
      data = ckbtools_raw_somascan_protein_data,
      sample_data = ckbtools_raw_somascan_sample_data,
      meta_data = somascan_meta_data,
      trans = "none",
      qc = "none",
      repeat_samples = "retain"
    )
  )
  expect_true("sample_source" %in% names(somascan_data_keep))
  expect_true(all(somascan_data_keep$sample_source %in% c(0, 1)))
})



test_that("prepare_somascan_data transformations correct", {
  somascan_meta_data  <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data,
                                                    names = "protein",
                                                    filter = "HumanProtein")
  suppressMessages(
    somascan_data_none <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                                ckbtools_raw_somascan_sample_data,
                                                somascan_meta_data,
                                                trans = "none",
                                                qc = "none",
                                                repeat_samples = "retain")
  )
  suppressMessages(
    somascan_data_log2 <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                                ckbtools_raw_somascan_sample_data,
                                                somascan_meta_data,
                                                trans = "log2",
                                                qc = "none",
                                                repeat_samples = "retain")
  )
  suppressMessages(
    somascan_data_log <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                               ckbtools_raw_somascan_sample_data,
                                               somascan_meta_data,
                                               trans = "log",
                                               qc = "none",
                                               repeat_samples = "retain")
  )
  suppressMessages(
    somascan_data_log10 <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                                 ckbtools_raw_somascan_sample_data,
                                                 somascan_meta_data,
                                                 trans = "log10",
                                                 qc = "none",
                                                 repeat_samples = "retain")
  )

  protein <- somascan_meta_data$name[[1]]

  expect_equal(
    log2(somascan_data_none[[protein]]), somascan_data_log2[[protein]]
  )
  expect_equal(
    log(somascan_data_none[[protein]]), somascan_data_log[[protein]]
  )
  expect_equal(
    log10(somascan_data_none[[protein]]), somascan_data_log10[[protein]]
  )
})


test_that("prepare_somascan_data default QC removes out-of-range hybcontrolnormscale", {
  somascan_meta_data  <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data,
                                                    names = "protein",
                                                    filter = "HumanProtein")

  # Force one sample out of acceptable hybcontrolnormscale range
  raw_somascan_sample_data_bad <- ckbtools_raw_somascan_sample_data
  raw_somascan_sample_data_bad$hybcontrolnormscale[1] <- "100"  # coerces to numeric, clearly > 2.5

  suppressMessages(
    keep <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                  ckbtools_raw_somascan_sample_data,
                                  somascan_meta_data,
                                  qc = "none")
  )
  suppressMessages(
    drop <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                  raw_somascan_sample_data_bad,
                                  somascan_meta_data,
                                  qc = "default")
  )

  csid_bad <- raw_somascan_sample_data_bad$csid[1]
  expect_true(csid_bad %in% keep$csid)
  expect_false(csid_bad %in% drop$csid)
})



# --- prepare_proteomics_analysis_dataset -------------------------------------
test_that("prepare_proteomics_analysis_dataset: join, plate-correction, scaling", {

  csid_in_subcohort <- ckbtools_ascertainments$csid[ckbtools_ascertainments$olinkexpexpan_chd_b1_subcohort == 1]
  ckbtools_participant_data$subcohort <- 0
  ckbtools_participant_data$subcohort[ckbtools_participant_data$csid %in% csid_in_subcohort] <- 1

  # Olink
  olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data, ckbtools_raw_olink_protein_data)
  suppressMessages(
    olink_data <- prepare_olink_data(ckbtools_raw_olink_protein_data, olink_meta_data, qc = "default")
  )

  # No plate correction, no scaling
  suppressMessages(
    olink_analysis_data_none <- prepare_proteomics_analysis_dataset(
      participant_data = ckbtools_participant_data,
      protein_data = olink_data,
      protein_meta_data = olink_meta_data,
      plate_correction = "none",
      scaling = "none",
      id = "csid"
    )
  )
  expect_true("plateid" %in% names(olink_analysis_data_none))
  expect_true(all(olink_meta_data$name %in% names(olink_analysis_data_none)))

  # Plate correction median-subcohort + rint scaling
  suppressMessages(
    olink_analysis_data_pc <- prepare_proteomics_analysis_dataset(
      participant_data = ckbtools_participant_data,
      protein_data = olink_data,
      protein_meta_data = olink_meta_data,
      plate_correction = "median-subcohort",
      id = "csid"
    )
  )

  suppressMessages(
    olink_analysis_data_pc_rint <- prepare_proteomics_analysis_dataset(
      participant_data = ckbtools_participant_data,
      protein_data = olink_data,
      protein_meta_data = olink_meta_data,
      plate_correction = "median-subcohort",
      scaling = "rint",
      id = "csid"
    )
  )

  # After RINT, columns should be roughly centred around 0 (not exact due to ties/missing)
  means <- vapply(olink_analysis_data_pc_rint[olink_meta_data$name],
                  function(x) mean(x, na.rm = TRUE),
                  numeric(1))
  expect_all_true(abs(means) < 0.01)


  # SomaScan
  somascan_meta_data <- prepare_somascan_meta_data(ckbtools_raw_somascan_meta_data)
  suppressMessages(
    somascan_data <- prepare_somascan_data(ckbtools_raw_somascan_protein_data,
                                           ckbtools_raw_somascan_sample_data,
                                           somascan_meta_data,
                                           qc = "default")
  )

  somascan_analysis_data <- suppressMessages(
    prepare_proteomics_analysis_dataset(
      participant_data = ckbtools_participant_data,
      protein_data = somascan_data,
      protein_meta_data = somascan_meta_data,
      plate_correction = "none",
      scaling = "scale",
      id = "csid"
    )
  )

  expect_true(all(somascan_meta_data$name %in% names(somascan_analysis_data)))
})

testthat::test_that("prepare_proteomics_analysis_dataset errors when subcohort is required but missing", {
  olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data,
                                             ckbtools_raw_olink_protein_data)
  suppressMessages(
    olink_data <- prepare_olink_data(ckbtools_raw_olink_protein_data,
                                     olink_meta_data,
                                     qc = "none")
  )

  testthat::expect_error(
    prepare_proteomics_analysis_dataset(
      participant_data = ckbtools_participant_data,
      protein_data = olink_data,
      protein_meta_data = olink_meta_data,
      plate_correction = "median-subcohort",
      scaling = "none",
      id = "csid"
    ),
    "There must be a subcohort column"
  )
})


testthat::test_that("prepare_proteomics_analysis_dataset tolerates missing protein values (no crash)", {

  csid_in_subcohort <- ckbtools_ascertainments$csid[ckbtools_ascertainments$olinkexpexpan_chd_b1_subcohort == 1]
  ckbtools_participant_data$subcohort <- 0
  ckbtools_participant_data$subcohort[ckbtools_participant_data$csid %in% csid_in_subcohort] <- 1

  olink_meta_data <- prepare_olink_meta_data(ckbtools_raw_olink_meta_data, ckbtools_raw_olink_protein_data)
  suppressMessages(
    olink_data <- prepare_olink_data(ckbtools_raw_olink_protein_data, olink_meta_data, qc = "default")
  )

  # Expect no errors (functions emit cli messages and plate-count warnings in synthetic data)
  testthat::expect_no_error(suppressMessages(
    prepare_proteomics_analysis_dataset(ckbtools_participant_data,
                                        olink_data,
                                        olink_meta_data,
                                        "none",
                                        "none")
  ))
  testthat::expect_no_error(suppressMessages(
    prepare_proteomics_analysis_dataset(ckbtools_participant_data,
                                        olink_data,
                                        olink_meta_data,
                                        "none",
                                        "center")
  ))
  testthat::expect_no_error(suppressMessages(
    prepare_proteomics_analysis_dataset(ckbtools_participant_data,
                                        olink_data,
                                        olink_meta_data,
                                        "none",
                                        "scale")
  ))
  testthat::expect_no_error(suppressMessages(
    prepare_proteomics_analysis_dataset(ckbtools_participant_data,
                                        olink_data,
                                        olink_meta_data,
                                        "none",
                                        "rint")
  ))
})
