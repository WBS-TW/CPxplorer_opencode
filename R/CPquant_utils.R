


# Various utilities and helper functions for CPquant

#-------------------------------
parse_transition_rel_ab <- function(transition_note) {
    matches <- stringr::str_match(transition_note, "\\{[^}]*\\}\\{([^}]*)\\}")
    as.numeric(matches[, 2])
}

#-------------------------------
normalize_cpquant_skyline <- function(df, standard_types = "Group Mixtures") {
    df <- df |>
        #make sure to update any_of() if Skyline changes the variable names in future versions
        dplyr::rename(Replicate_Name = tidyr::any_of(c("Replicate Name", "ReplicateName", "Replicate"))) |>
        dplyr::rename(Sample_Type = tidyr::any_of(c("Sample Type", "SampleType"))) |>
        dplyr::rename(Molecule_List = tidyr::any_of(c("Molecule List", "MoleculeList"))) |>
        dplyr::rename(Mass_Error_PPM = tidyr::any_of(c("Mass Error PPM", "MassErrorPPM"))) |>
        dplyr::rename(Isotope_Label_Type = tidyr::any_of(c("Isotope Label Type", "IsotopeLabelType"))) |>
        dplyr::rename(Chromatogram_Precursor_MZ = tidyr::any_of(c("Chromatogram Precursor M/Z", "ChromatogramPrecursorMz"))) |>
        dplyr::rename(Analyte_Concentration = tidyr::any_of(c("Analyte Concentration", "AnalyteConcentration"))) |>
        dplyr::rename(Batch_Name = tidyr::any_of(c("Batch Name", "BatchName"))) |>
        dplyr::rename(Transition_Note = tidyr::any_of(c("Transition Note", "TransitionNote"))) |>
        dplyr::rename(Sample_Dilution_Factor = tidyr::any_of(c("Sample Dilution Factor", "SampleDilutionFactor"))) |>
        dplyr::mutate(
            Analyte_Concentration = as.numeric(Analyte_Concentration),
            Area = as.numeric(Area),
            Area = tidyr::replace_na(Area, 0),
            C_homologue = stringr::str_extract(Molecule, "C\\d+"),
            Cl_homologue = stringr::str_extract(Molecule, "Cl\\d+"),
            C_number = as.numeric(stringr::str_extract(C_homologue, "\\d+")),
            Cl_number = as.numeric(stringr::str_extract(Cl_homologue, "\\d+")),
            PCA = stringr::str_c(C_homologue, Cl_homologue, sep = ""),
            Rel_Ab = parse_transition_rel_ab(Transition_Note)
        )

    if (identical(standard_types, "Group Mixtures")) {
        df <- df |>
            dplyr::mutate(Quantification_Group = stringr::str_extract(Batch_Name, "^[^_]+"))
    }

    df
}

#-------------------------------
read_cpquant_skyline <- function(path, standard_types = "Group Mixtures") {
    readxl::read_excel(
        path,
        guess_max = 5000,
        na = c("", "NA", "#N/A", "N/A")
    ) |>
        normalize_cpquant_skyline(standard_types = standard_types)
}

#-------------------------------
apply_cpquant_preprocessing <- function(
    df,
    remove_samples = NULL,
    quan_sum = "Quan only",
    correct_with_rs = "No",
    choose_rs = NULL,
    blank_subtraction = "No",
    standard_types = "Group Mixtures"
) {
    if (!is.null(remove_samples) && length(remove_samples) > 0) {
        df <- df |>
            dplyr::filter(!Replicate_Name %in% remove_samples)
    }

    if (identical(quan_sum, "Sum Quan+Qual")) {
        df <- df |>
            dplyr::group_by(Replicate_Name, Molecule_List, Molecule) |>
            dplyr::mutate(Area = sum(Area)) |>
            dplyr::ungroup()
    }

    if (identical(correct_with_rs, "Yes") && any(df$Molecule_List == "RS")) {
        if (!is.null(choose_rs) && nzchar(choose_rs)) {
            df <- df |>
                dplyr::group_by(Replicate_Name) |>
                dplyr::mutate(
                    .rs_val = dplyr::first(Area[
                        Molecule == choose_rs &
                            Molecule_List == "RS" &
                            Isotope_Label_Type == "Quan"
                    ]),
                    Area = dplyr::if_else(!is.na(.rs_val) & .rs_val > 0, Area / .rs_val, NA_real_)
                ) |>
                dplyr::ungroup() |>
                dplyr::select(-.rs_val)
        }
    }

    if (identical(blank_subtraction, "Yes, by avg area of blanks")) {
        df_blank <- df |>
            dplyr::filter(Sample_Type == "Blank") |>
            dplyr::group_by(Molecule, Molecule_List, Isotope_Label_Type) |>
            dplyr::summarise(AverageBlank = mean(Area, na.rm = TRUE), .groups = "drop") |>
            dplyr::filter(!Molecule_List %in% c("IS", "RS", "VS"))

        df <- df |>
            dplyr::left_join(
                df_blank,
                by = c("Molecule", "Molecule_List", "Isotope_Label_Type")
            ) |>
            dplyr::mutate(AverageBlank = tidyr::replace_na(AverageBlank, 0)) |>
            dplyr::mutate(
                Area = dplyr::case_when(
                    Sample_Type == "Unknown" ~ Area - AverageBlank,
                    TRUE ~ Area
                ),
                Area = ifelse(Area < 0, 0, Area)
            )
    }

    if (identical(standard_types, "Group Mixtures")) {
        df <- df |>
            dplyr::mutate(
                Quantification_Group = stringr::str_extract(Batch_Name, "^[^_]+")
            )
    }

    df
}

#-------------------------------
prepare_cpquant_standards <- function(df, min_rsquared = 0.9) {
    df |>
        dplyr::filter(
            Sample_Type == "Standard",
            !Molecule_List %in% c("IS", "RS"),
            Isotope_Label_Type == "Quan",
            !is.na(Batch_Name),
            Batch_Name != "NA"
        ) |>
        dplyr::mutate(
            C_range = stringr::str_extract_all(Quantification_Group, "\\d+"),
            C_min = as.numeric(purrr::map_chr(C_range, ~ .x[1])),
            C_max = as.numeric(purrr::map_chr(C_range, function(x) {
                if (length(x) > 1) x[2] else x[1]
            }))
        ) |>
        dplyr::select(-C_range) |>
        dplyr::group_by(
            Batch_Name, Sample_Type, Molecule, Molecule_List,
            C_number, Cl_number, PCA, Quantification_Group, C_min, C_max
        ) |>
        tidyr::nest() |>
        dplyr::filter(C_number >= C_min & C_number <= C_max) |>
        dplyr::mutate(models = purrr::map(data, ~ stats::lm(Area ~ Analyte_Concentration, data = .x))) |>
        dplyr::mutate(coef = purrr::map(models, stats::coef)) |>
        dplyr::mutate(RF_raw = purrr::map_dbl(models, ~ as.numeric(stats::coef(.x)["Analyte_Concentration"]))) |>
        dplyr::mutate(intercept = purrr::map_dbl(coef, ~ as.numeric(purrr::pluck(.x, "(Intercept)", .default = NA_real_)))) |>
        dplyr::mutate(cal_rsquared = purrr::map(models, summary)) |>
        dplyr::mutate(cal_rsquared = purrr::map_dbl(cal_rsquared, ~ as.numeric(purrr::pluck(.x, "r.squared", .default = 0)))) |>
        dplyr::select(-coef) |>
        dplyr::mutate(
            cal_rsquared = ifelse(is.nan(cal_rsquared) | is.na(cal_rsquared), 0, cal_rsquared),
            RF = dplyr::if_else(is.na(RF_raw) | RF_raw < 0 | cal_rsquared < min_rsquared, 0, RF_raw)
        ) |>
        dplyr::ungroup() |>
        dplyr::group_by(Batch_Name) |>
        dplyr::mutate(Sum_RF_group = sum(RF, na.rm = TRUE)) |>
        dplyr::ungroup()
}

#-------------------------------
format_calibration_export_table <- function(standards_tbl) {
    if (!"RF_raw" %in% names(standards_tbl)) {
        standards_tbl <- dplyr::mutate(standards_tbl, RF_raw = RF)
    }

    standards_tbl |>
        dplyr::ungroup() |>
        dplyr::select(
            Batch_Name,
            Quantification_Group,
            Molecule_List,
            Molecule,
            C_number,
            Cl_number,
            PCA,
            C_min,
            C_max,
            RF_raw,
            RF,
            intercept,
            cal_rsquared,
            Sum_RF_group
        ) |>
        dplyr::mutate(
            C_number = as.numeric(C_number),
            Cl_number = as.numeric(Cl_number),
            C_min = as.numeric(C_min),
            C_max = as.numeric(C_max),
            RF_raw = as.numeric(RF_raw),
            RF = as.numeric(RF),
            intercept = as.numeric(intercept),
            cal_rsquared = as.numeric(cal_rsquared),
            Sum_RF_group = as.numeric(Sum_RF_group)
        )
}

#-------------------------------
write_cpquant_calibration_xlsx <- function(
    path,
    standards_tbl,
    raw_standards_df = NULL,
    metadata = list()
) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Calibration")
    openxlsx::addWorksheet(wb, "Standards_Data")
    openxlsx::addWorksheet(wb, "Metadata")

    cal_export <- format_calibration_export_table(standards_tbl)
    openxlsx::writeData(wb, "Calibration", cal_export)

    if (is.null(raw_standards_df)) {
        raw_standards_df <- tibble::tibble()
    } else {
        raw_standards_df <- raw_standards_df |>
            dplyr::filter(Sample_Type == "Standard")
    }
    openxlsx::writeData(wb, "Standards_Data", raw_standards_df)

    default_meta <- list(
        created_date = as.character(Sys.Date()),
        standardTypes = "Group Mixtures",
        quanSum = "Quan only",
        blankSubtraction = "No",
        correctWithRS = "No",
        chooseRS = "",
        removeRsquared = "0.9",
        quanUnit = ""
    )
    merged_meta <- utils::modifyList(default_meta, lapply(metadata, function(x) {
        if (is.null(x) || length(x) == 0 || all(is.na(x))) "" else as.character(x[[1]])
    }))
    meta_tbl <- tibble::tibble(
        Setting = names(merged_meta),
        Value = unlist(merged_meta, use.names = FALSE)
    )
    openxlsx::writeData(wb, "Metadata", meta_tbl)

    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    invisible(path)
}

#-------------------------------
read_cpquant_calibration_xlsx <- function(
    path,
    min_rsquared = 0.9,
    quan_sum = NULL,
    correct_with_rs = NULL,
    choose_rs = NULL,
    standard_types = "Group Mixtures"
) {
    sheets <- readxl::excel_sheets(path)
    if (!"Calibration" %in% sheets) {
        stop("Invalid calibration file: missing 'Calibration' worksheet.")
    }

    cal_df <- readxl::read_excel(path, sheet = "Calibration", guess_max = 5000)
    required_cols <- c(
        "Batch_Name", "Quantification_Group", "Molecule_List", "Molecule",
        "C_number", "Cl_number", "PCA", "RF", "intercept", "cal_rsquared"
    )
    missing_cols <- setdiff(required_cols, names(cal_df))
    if (length(missing_cols) > 0) {
        stop(sprintf(
            "Invalid calibration file: missing required column(s) in 'Calibration' sheet: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    metadata <- list()
    if ("Metadata" %in% sheets) {
        meta_df <- readxl::read_excel(path, sheet = "Metadata")
        if (all(c("Setting", "Value") %in% names(meta_df))) {
            meta_vals <- ifelse(is.na(meta_df$Value), "", as.character(meta_df$Value))
            metadata <- stats::setNames(as.list(meta_vals), as.character(meta_df$Setting))
        }
    }

    raw_standards <- tibble::tibble()
    if ("Standards_Data" %in% sheets) {
        raw_in <- readxl::read_excel(
            path,
            sheet = "Standards_Data",
            guess_max = 5000,
            na = c("", "NA", "#N/A", "N/A")
        )
        if (nrow(raw_in) > 0) {
            raw_standards <- normalize_cpquant_skyline(raw_in, standard_types = standard_types)
        }
    }

    settings_changed <- FALSE
    if (nrow(raw_standards) > 0 && length(metadata) > 0) {
        if (!is.null(quan_sum) && !is.null(metadata$quanSum) && nzchar(metadata$quanSum) &&
            !identical(as.character(quan_sum), as.character(metadata$quanSum))) {
            settings_changed <- TRUE
        }
        if (!is.null(correct_with_rs) && !is.null(metadata$correctWithRS) && nzchar(metadata$correctWithRS) &&
            !identical(as.character(correct_with_rs), as.character(metadata$correctWithRS))) {
            settings_changed <- TRUE
        }
        if (identical(correct_with_rs, "Yes") && !is.null(choose_rs) && nzchar(choose_rs) &&
            !is.null(metadata$chooseRS) && !identical(as.character(choose_rs), as.character(metadata$chooseRS))) {
            settings_changed <- TRUE
        }
    }

    if (settings_changed && nrow(raw_standards) > 0) {
        processed_std <- apply_cpquant_preprocessing(
            raw_standards,
            remove_samples = NULL,
            quan_sum = if (!is.null(quan_sum)) quan_sum else "Quan only",
            correct_with_rs = if (!is.null(correct_with_rs)) correct_with_rs else "No",
            choose_rs = choose_rs,
            blank_subtraction = "No",
            standard_types = standard_types
        )
        standards <- prepare_cpquant_standards(processed_std, min_rsquared = min_rsquared)
    } else {
        if (!"RF_raw" %in% names(cal_df)) {
            cal_df <- dplyr::mutate(cal_df, RF_raw = as.numeric(RF))
        }
        if (!all(c("C_min", "C_max") %in% names(cal_df))) {
            cal_df <- cal_df |>
                dplyr::mutate(
                    C_range = stringr::str_extract_all(Quantification_Group, "\\d+"),
                    C_min = as.numeric(purrr::map_chr(C_range, ~ .x[1])),
                    C_max = as.numeric(purrr::map_chr(C_range, function(x) {
                        if (length(x) > 1) x[2] else x[1]
                    }))
                ) |>
                dplyr::select(-C_range)
        }
        standards <- cal_df |>
            dplyr::mutate(
                Sample_Type = "Standard",
                C_number = as.numeric(C_number),
                Cl_number = as.numeric(Cl_number),
                C_min = as.numeric(C_min),
                C_max = as.numeric(C_max),
                RF_raw = as.numeric(RF_raw),
                intercept = as.numeric(intercept),
                cal_rsquared = as.numeric(cal_rsquared),
                cal_rsquared = ifelse(is.nan(cal_rsquared) | is.na(cal_rsquared), 0, cal_rsquared),
                RF = dplyr::if_else(is.na(RF_raw) | RF_raw < 0 | cal_rsquared < min_rsquared, 0, RF_raw)
            ) |>
            dplyr::group_by(Batch_Name) |>
            dplyr::mutate(Sum_RF_group = sum(RF, na.rm = TRUE)) |>
            dplyr::ungroup()
    }

    list(
        standards = standards,
        raw_standards = raw_standards,
        metadata = metadata
    )
}

#-------------------------------
compare_calibration_rf_qc <- function(
    saved_standards,
    current_skyline_df,
    pass_pct = 20,
    review_pct = 30
) {
    empty_batch <- tibble::tibble(
        Batch_Name = character(),
        Quantification_Group = character(),
        n_homologues = integer(),
        n_standards = integer(),
        Sum_RF_saved = numeric(),
        Sum_RF_current = numeric(),
        RF_ratio = numeric(),
        pct_diff = numeric(),
        qc_flag = character()
    )
    empty_homologue <- tibble::tibble(
        Batch_Name = character(),
        Quantification_Group = character(),
        Molecule_List = character(),
        Molecule = character(),
        C_number = numeric(),
        Cl_number = numeric(),
        PCA = character(),
        n_standards = integer(),
        RF_saved = numeric(),
        RF_current = numeric(),
        RF_ratio = numeric(),
        pct_diff = numeric(),
        qc_flag = character()
    )

    if (is.null(saved_standards) || nrow(saved_standards) == 0 ||
        is.null(current_skyline_df) || nrow(current_skyline_df) == 0) {
        return(list(batch_qc = empty_batch, homologue_qc = empty_homologue))
    }

    saved_active <- saved_standards |>
        dplyr::ungroup() |>
        dplyr::filter(RF > 0) |>
        dplyr::select(
            Batch_Name, Quantification_Group, Molecule_List, Molecule,
            C_number, Cl_number, PCA, RF_saved = RF
        )

    if (nrow(saved_active) == 0) {
        return(list(batch_qc = empty_batch, homologue_qc = empty_homologue))
    }

    check_points <- current_skyline_df |>
        dplyr::filter(
            Sample_Type == "Standard",
            !Molecule_List %in% c("IS", "RS", "VS"),
            Isotope_Label_Type == "Quan",
            !is.na(Batch_Name),
            Batch_Name != "NA",
            !is.na(Analyte_Concentration),
            Analyte_Concentration > 0
        ) |>
        dplyr::mutate(RF_point = Area / Analyte_Concentration)

    if (nrow(check_points) == 0) {
        return(list(batch_qc = empty_batch, homologue_qc = empty_homologue))
    }

    joined <- saved_active |>
        dplyr::inner_join(
            check_points |>
                dplyr::select(Batch_Name, Molecule, Replicate_Name, Analyte_Concentration, Area, RF_point),
            by = c("Batch_Name", "Molecule")
        )

    if (nrow(joined) == 0) {
        return(list(batch_qc = empty_batch, homologue_qc = empty_homologue))
    }

    homologue_qc <- joined |>
        dplyr::group_by(
            Batch_Name, Quantification_Group, Molecule_List, Molecule,
            C_number, Cl_number, PCA, RF_saved
        ) |>
        dplyr::summarise(
            n_standards = as.integer(dplyr::n()),
            RF_current = mean(RF_point, na.rm = TRUE),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            RF_ratio = dplyr::if_else(RF_saved > 0, RF_current / RF_saved, NA_real_),
            pct_diff = dplyr::if_else(RF_saved > 0, (RF_current - RF_saved) / RF_saved * 100, NA_real_),
            qc_flag = dplyr::case_when(
                is.na(pct_diff) | !is.finite(pct_diff) ~ "No signal",
                abs(pct_diff) <= pass_pct ~ "Pass",
                abs(pct_diff) <= review_pct ~ "Review",
                TRUE ~ "Fail"
            )
        ) |>
        dplyr::arrange(Batch_Name, C_number, Cl_number)

    batch_qc <- homologue_qc |>
        dplyr::group_by(Batch_Name, Quantification_Group) |>
        dplyr::summarise(
            n_homologues = as.integer(dplyr::n()),
            n_standards = as.integer(max(n_standards, na.rm = TRUE)),
            Sum_RF_saved = sum(RF_saved, na.rm = TRUE),
            Sum_RF_current = sum(RF_current, na.rm = TRUE),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            RF_ratio = dplyr::if_else(Sum_RF_saved > 0, Sum_RF_current / Sum_RF_saved, NA_real_),
            pct_diff = dplyr::if_else(Sum_RF_saved > 0, (Sum_RF_current - Sum_RF_saved) / Sum_RF_saved * 100, NA_real_),
            qc_flag = dplyr::case_when(
                is.na(pct_diff) | !is.finite(pct_diff) ~ "No signal",
                abs(pct_diff) <= pass_pct ~ "Pass",
                abs(pct_diff) <= review_pct ~ "Review",
                TRUE ~ "Fail"
            )
        ) |>
        dplyr::arrange(Batch_Name)

    list(batch_qc = batch_qc, homologue_qc = homologue_qc)
}

# =========================================================
# Helper UI builders
# =========================================================


#-------------------------------
defineVariablesUI <- function(df_for_choices) {
    rep_names <- df_for_choices$Replicate_Name
    rep_names <- unique(stats::na.omit(rep_names))
    rep_names <- rep_names[order(rep_names)]

    shiny::fluidRow(
        shiny::column(
            6,
            shiny::selectInput(
                inputId = "removeSamples", # select samples to remove from quantification
                label   = "Remove samples from quantification?",
                choices = rep_names,
                selected = NULL,
                multiple = TRUE
            )
        ),

        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),

        shiny::column(
            6,
            shiny::sliderInput(
                inputId = "removeRsquared", # keep only Molecule from standard calibration curves above this R^2
                label   = "Keep the calibration curves above this R-squared (0 means keep everything)",
                min     = 0,
                max     = 1,
                value   = 0.90,
                step    = 0.05
            )
            # shiny::checkboxInput(
            #   inputId = "zerointercept", # force y-intercept through zero
            #   label   = "Set intercept to zero",
            #   value   = FALSE
            # )
        )
    )
}

#-------------------------------
defineCorrectionUI <- function(df_for_choices) {
    # Extract RS molecules for selection; keep it stable (from base data)
    rs_choices <- df_for_choices |>
        dplyr::filter(.data$Molecule_List == "RS") |>
        dplyr::pull(.data$Molecule) |>
        unique() |>
        stats::na.omit() |>
        as.character()

    rs_choices <- rs_choices[order(rs_choices)]

    shiny::fluidRow(
        shiny::column(
            6,
            shiny::selectInput(
                inputId = "chooseRS",      # select which will be the RS
                label   = "Choose RS for correction",
                choices = rs_choices,
                selected = if (length(rs_choices)) rs_choices[1] else NULL,
                multiple = FALSE
            )
        )
    )
}

# =========================================================
# Helper plotting functions
# =========================================================

#-------------------------------
plot_skyline_output <- function(Skyline_output){

    Skyline_output |>
        dplyr::filter(Isotope_Label_Type == "Quan") |>
        dplyr::mutate(OrderedMolecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)]))) |>
        plotly::plot_ly(
            x = ~ OrderedMolecule,
            y = ~ Area,
            color = ~ Sample_Type,
            type = "box",
            text = ~paste(
                "Homologue: ", PCA,
                "<br>Sample: ", Replicate_Name,
                "<br>Area:", round(Area, 2)
            ),
            hoverinfo = "text"
        ) |>
        plotly::layout(xaxis = list(title = 'Molecule'),
                       yaxis = list(title = 'Area'))
}

#-------------------------------
make_calibration_table <- function(CPs_standards_input) {

    CPs_standards_input %>%
        dplyr::ungroup() %>%
        dplyr::select(
            Batch_Name,
            Quantification_Group,
            Molecule_List,
            Molecule,
            C_number,
            Cl_number,
            PCA,
            RF,
            intercept,
            cal_rsquared
        ) %>%
        dplyr::mutate(
            RF = as.numeric(RF),
            intercept = as.numeric(intercept),
            cal_rsquared = as.numeric(cal_rsquared)
        ) %>%
        dplyr::arrange(
            Batch_Name, Quantification_Group,
            Molecule_List, Molecule
        )
}


#-------------------------------
##currently not used##
plot_calibration_curves <- function(CPs_standards, quantUnit) {
    # Unnest the data first
    CPs_standards_unnested <- CPs_standards |>
        dplyr::filter(RF > 0) |>
        tidyr::unnest(data) |>
        dplyr::mutate(Molecule = factor(Molecule,
                                        levels = unique(Molecule[order(C_number, Cl_number)])))

    # Get unique Quantification Groups
    groups <- unique(CPs_standards_unnested$Quantification_Group)

    # Create a list to store individual plots
    plot_list <- list()

    # Create individual plots for each group
    for(i in seq_along(groups)) {
        group_data <- CPs_standards_unnested |>
            dplyr::filter(Quantification_Group == groups[i])

        plot_list[[i]] <- plotly::plot_ly() |>
            plotly::add_trace(
                data = group_data,
                x = ~Analyte_Concentration,
                y = ~Area,
                color = ~PCA,
                type = 'scatter',
                mode = 'markers',
                legendgroup = ~Molecule_List,
                legendgrouptitle = list(text = ~Molecule_List),
                showlegend = FALSE,
                name = ~paste(PCA, "(points)"),
                text = ~paste(
                    "Molecule:", Molecule,
                    "<br>Area:", round(Area, 2),
                    "<br>Concentration:", round(Analyte_Concentration, 2),
                    "<br>R2:", round(cal_rsquared, 3)
                ),
                hoverinfo = 'text'
            ) |>
            plotly::add_trace(
                data = group_data,
                x = ~Analyte_Concentration,
                y = ~RF * Analyte_Concentration + intercept,
                color = ~PCA,
                type = 'scatter',
                mode = 'lines',
                legendgroup = ~Molecule_List,
                name = ~PCA,
                hoverinfo = 'none'
            )
    }

    # Calculate layout
    subplot_cols <- min(length(groups), 2)  # Maximum 2 columns
    subplot_rows <- ceiling(length(groups) / subplot_cols)

    # Create annotations for titles
    annotations <- list()
    for(i in seq_along(groups)) {
        row <- ceiling(i/subplot_cols)
        col <- if(i %% subplot_cols == 0) subplot_cols else i %% subplot_cols

        annotations[[i]] <- list(
            text = groups[i],
            font = list(size = 14),
            xref = "paper",
            yref = "paper",
            x = (col - 0.5)/subplot_cols,
            y = 1 - (row - 1)/subplot_rows,
            xanchor = "center",
            yanchor = "bottom",
            showarrow = FALSE
        )
    }

    # Combine plots using subplot
    final_plot <- plotly::subplot(
        plot_list,
        nrows = subplot_rows,
        shareX = FALSE,
        shareY = FALSE,
        margin = 0.1
    ) |>
        plotly::layout(
            height = 400 * subplot_rows,
            showlegend = TRUE,
            annotations = annotations,
            margin = list(t = 50, b = 50, l = 50, r = 50),
            legend = list(
                groupclick = "togglegroup",
                tracegroupgap = 10,
                itemsizing = "constant"
            ),
            xaxis = list(title = paste0("Analyte Concentration/Amount (", quantUnit, ")"))
        )

    return(final_plot)
}

#-------------------------------
plot_quanqualratio <- function(Skyline_output_filt) {

    Skyline_output_filt |>
        dplyr::group_by(Replicate_Name, Molecule) |>
        dplyr::mutate(Quan_Area = ifelse(Isotope_Label_Type == "Quan", Area, NA)) |>
        tidyr::fill(Quan_Area, .direction = "downup") |>
        dplyr::mutate(QuanMZ = ifelse(Isotope_Label_Type == "Quan", Chromatogram_Precursor_MZ, NA)) |>
        tidyr::fill(QuanMZ, .direction = "downup") |>
        dplyr::mutate(QuanQualRatio = ifelse(Isotope_Label_Type == "Qual", Quan_Area/Area, 1)) |>
        dplyr::mutate(QuanQualRatio = dplyr::na_if(QuanQualRatio, Inf)) |>
        dplyr::mutate(QuanQualRatio = dplyr::na_if(QuanQualRatio, NaN)) |>
        tidyr::replace_na(list(QuanQualRatio = 0)) |>
        dplyr::mutate(QuanQualMZ = paste0(QuanMZ,"/",Chromatogram_Precursor_MZ)) |>
        dplyr::ungroup() |>
        dplyr::select(Replicate_Name, Sample_Type, Molecule_List, Molecule, QuanQualMZ, QuanQualRatio) |>
        plotly::plot_ly(x = ~Replicate_Name, y = ~QuanQualRatio, type = 'violin', color = ~Sample_Type,
                        text = ~paste("Sample: ", Replicate_Name,
                                      "<br>Molecule List: ", Molecule_List,
                                      "<br>Molecule: ", Molecule,
                                      "<br>Quan/Qual MZ: ", QuanQualMZ,
                                      "<br>Ratio: ", round(QuanQualRatio, 2)),
                        hoverinfo = "text") |>
        plotly::layout(title = 'Quan-to-Qual Ratio',
                       xaxis = list(title = 'Replicate Name'),
                       yaxis = list(title = 'Quan-to-Qual Ratio'))
}

#-------------------------------
plot_meas_vs_theor_ratio <- function(Skyline_output_filt) {

    Skyline_output_filt |>
        dplyr::group_by(Replicate_Name, Molecule) |>
        dplyr::mutate(Quan_Area = ifelse(Isotope_Label_Type == "Quan", Area, NA)) |>
        tidyr::fill(Quan_Area, .direction = "downup") |>
        dplyr::mutate(QuanMZ = ifelse(Isotope_Label_Type == "Quan", Chromatogram_Precursor_MZ, NA)) |>
        tidyr::fill(QuanMZ, .direction = "downup") |>
        dplyr::mutate(QuanQualRatio = ifelse(Isotope_Label_Type == "Qual", Quan_Area/Area, 1)) |>
        dplyr::mutate(QuanQualRatio = dplyr::na_if(QuanQualRatio, Inf)) |>
        dplyr::mutate(QuanQualRatio = dplyr::na_if(QuanQualRatio, NaN)) |>
        tidyr::replace_na(list(QuanQualRatio = 0)) |>
        dplyr::mutate(QuanQualMZ = paste0(QuanMZ,"/",Chromatogram_Precursor_MZ)) |>

        dplyr::mutate(Quan_Rel_Ab = ifelse(Isotope_Label_Type == "Quan", Rel_Ab, NA)) |>
        tidyr::fill(Quan_Rel_Ab, .direction = "downup") |>
        dplyr::mutate(QuanQual_Rel_Ab_Ratio = ifelse(Isotope_Label_Type == "Qual", Quan_Rel_Ab/Rel_Ab, 1)) |>
        dplyr::mutate(QuanQual_Rel_Ab_Ratio = dplyr::na_if(QuanQual_Rel_Ab_Ratio, Inf)) |>
        dplyr::mutate(QuanQual_Rel_Ab_Ratio = dplyr::na_if(QuanQual_Rel_Ab_Ratio, NaN)) |>
        tidyr::replace_na(list(QuanQual_Rel_Ab_Ratio = 0)) |>
        dplyr::ungroup() |>
        dplyr::mutate(MeasVSTheo = QuanQualRatio/QuanQual_Rel_Ab_Ratio) |>
        dplyr::mutate(MeasVSTheo = dplyr::na_if(MeasVSTheo, Inf)) |>
        dplyr::mutate(MeasVSTheo = dplyr::na_if(MeasVSTheo, NaN)) |>
        dplyr::mutate(Is_Outlier = MeasVSTheo > 3 | MeasVSTheo < 0.3) |>
        dplyr::mutate(Is_Outlier = factor(Is_Outlier, levels = c(FALSE, TRUE), labels = c("Within Limit", "Outlier"))) |>
        dplyr::select(Replicate_Name, Sample_Type, Molecule_List, Molecule, QuanQualMZ, QuanQualRatio, QuanQual_Rel_Ab_Ratio, MeasVSTheo, Is_Outlier) |>
        plotly::plot_ly(x = ~Replicate_Name, y = ~MeasVSTheo,
                        type = 'scatter', mode = 'markers',
                        color = ~Is_Outlier,
                        colors = c('blue', 'red'),
                        text = ~paste("Replicate:", Replicate_Name,
                                      "<br>Homologue Group: ", Molecule,
                                      "<br>Measured against Theoretical Ratio:", round(MeasVSTheo, 1)),
                        marker = list(size = 10)) |>
        layout(title = "Measured/Theoretical ratio >3 or <0.3 are marked in red (ratio of 1 means perfect match",
               xaxis = list(title = "Sample Name"),
               yaxis = list(title = "MeasVSTheo"))

}

#-------------------------------
prepare_isotope_pattern_qc <- function(Skyline_output_filt) {

    ion_qc <- Skyline_output_filt |>
        dplyr::filter(!Molecule_List %in% c("IS", "RS", "VS")) |>
        dplyr::group_by(Replicate_Name, Sample_Type, Molecule) |>
        dplyr::mutate(
            Observed_Pct = if (sum(Area, na.rm = TRUE) > 0) Area / sum(Area, na.rm = TRUE) * 100 else NA_real_,
            Theoretical_Pct = if (sum(Rel_Ab, na.rm = TRUE) > 0) Rel_Ab / sum(Rel_Ab, na.rm = TRUE) * 100 else NA_real_,
            Measured_Theoretical_Ratio = Observed_Pct / Theoretical_Pct,
            Ion_Label = paste(Isotope_Label_Type, round(Chromatogram_Precursor_MZ, 6), sep = " ")
        ) |>
        dplyr::mutate(Measured_Theoretical_Ratio = dplyr::na_if(Measured_Theoretical_Ratio, Inf)) |>
        dplyr::mutate(Measured_Theoretical_Ratio = dplyr::na_if(Measured_Theoretical_Ratio, NaN)) |>
        dplyr::ungroup()

    summary_qc <- ion_qc |>
        dplyr::group_by(Replicate_Name, Sample_Type, Molecule) |>
        dplyr::summarise(
            n_ions_observed = sum(!is.na(Observed_Pct) & Area > 0),
            n_ions_expected = sum(!is.na(Theoretical_Pct) & Theoretical_Pct > 0),
            total_area = sum(Area, na.rm = TRUE),
            cosine_similarity = {
                observed <- Observed_Pct
                theoretical <- Theoretical_Pct
                denom <- sqrt(sum(observed^2, na.rm = TRUE)) * sqrt(sum(theoretical^2, na.rm = TRUE))
                ifelse(denom > 0, sum(observed * theoretical, na.rm = TRUE) / denom, NA_real_)
            },
            pearson_r = ifelse(
                sum(!is.na(Observed_Pct) & !is.na(Theoretical_Pct)) >= 2,
                stats::cor(Observed_Pct, Theoretical_Pct, use = "complete.obs", method = "pearson"),
                NA_real_
            ),
            weighted_abs_percent_error = {
                theoretical_sum <- sum(Theoretical_Pct, na.rm = TRUE)
                if (theoretical_sum == 0) NA_real_ else {
                    weights <- Theoretical_Pct / theoretical_sum
                    sum(abs(Observed_Pct - Theoretical_Pct) * weights, na.rm = TRUE)
                }
            },
            max_ion_ratio_error = {
                ratio_error <- abs(log2(Measured_Theoretical_Ratio))
                ratio_error <- dplyr::na_if(ratio_error, NaN)
                if (any(!is.na(ratio_error) & is.finite(ratio_error))) max(ratio_error, na.rm = TRUE) else NA_real_
            },
            .groups = "drop"
        ) |>
        dplyr::mutate(
            max_ion_ratio_error = dplyr::if_else(is.infinite(max_ion_ratio_error) | is.nan(max_ion_ratio_error), NA_real_, max_ion_ratio_error),
            qc_flag = dplyr::case_when(
                is.na(cosine_similarity) ~ "No signal",
                cosine_similarity >= 0.90 ~ "Pass",
                cosine_similarity >= 0.75 ~ "Review",
                TRUE ~ "Fail"
            )
        )

    list(ion_qc = ion_qc, summary_qc = summary_qc)
}

#-------------------------------
plot_isotope_pattern_overlay <- function(ion_qc, selected_replicate, selected_molecule) {

    plot_data <- ion_qc |>
        dplyr::filter(Replicate_Name == selected_replicate,
                      Molecule == selected_molecule) |>
        dplyr::arrange(Chromatogram_Precursor_MZ)

    shiny::req(nrow(plot_data) > 0)

    plotly::plot_ly(plot_data, x = ~Ion_Label) |>
        plotly::add_bars(
            y = ~Observed_Pct,
            name = "Observed",
            text = ~paste(
                "Sample:", Replicate_Name,
                "<br>Molecule:", Molecule,
                "<br>m/z:", round(Chromatogram_Precursor_MZ, 6),
                "<br>Observed (%):", round(Observed_Pct, 2),
                "<br>Theoretical (%):", round(Theoretical_Pct, 2)
            ),
            hoverinfo = "text"
        ) |>
        plotly::add_trace(
            y = ~Theoretical_Pct,
            name = "Theoretical",
            type = "scatter",
            mode = "markers",
            marker = list(color = "black", size = 8)
        ) |>
        plotly::layout(
            title = paste("Observed vs theoretical isotope pattern:", selected_molecule),
            xaxis = list(title = "Transition", tickangle = 45),
            yaxis = list(title = "Relative abundance (%)"),
            barmode = "group"
        )
}

#-------------------------------
plot_isotope_similarity_heatmap <- function(summary_qc) {

    plot_data <- summary_qc |>
        dplyr::mutate(cosine_similarity = round(cosine_similarity, 3))

    plotly::plot_ly(
        plot_data,
        x = ~Replicate_Name,
        y = ~Molecule,
        z = ~cosine_similarity,
        type = "heatmap",
        colorscale = "Viridis",
        zmin = 0,
        zmax = 1,
        text = ~paste(
            "Sample:", Replicate_Name,
            "<br>Molecule:", Molecule,
            "<br>Cosine similarity:", cosine_similarity,
            "<br>QC flag:", qc_flag,
            "<br>Total area:", round(total_area, 0)
        ),
        hoverinfo = "text"
    ) |>
        plotly::layout(
            title = "Isotope-pattern similarity to theoretical pattern",
            xaxis = list(title = "Sample"),
            yaxis = list(title = "Molecule")
        )
}

#-------------------------------
plot_isotope_ratio_residuals <- function(ion_qc) {

    plot_data <- ion_qc |>
        dplyr::filter(!is.na(Measured_Theoretical_Ratio), Measured_Theoretical_Ratio > 0) |>
        dplyr::mutate(
            Is_Outlier = Measured_Theoretical_Ratio > 3 | Measured_Theoretical_Ratio < 0.3,
            Is_Outlier = factor(Is_Outlier, levels = c(FALSE, TRUE), labels = c("Within Limit", "Outlier"))
        )

    plotly::plot_ly(
        plot_data,
        x = ~Replicate_Name,
        y = ~Measured_Theoretical_Ratio,
        type = "scatter",
        mode = "markers",
        color = ~Is_Outlier,
        colors = c("blue", "red"),
        text = ~paste(
            "Sample:", Replicate_Name,
            "<br>Molecule:", Molecule,
            "<br>Transition:", Ion_Label,
            "<br>Observed/Theoretical:", round(Measured_Theoretical_Ratio, 2),
            "<br>Observed (%):", round(Observed_Pct, 2),
            "<br>Theoretical (%):", round(Theoretical_Pct, 2)
        ),
        hoverinfo = "text"
    ) |>
        plotly::layout(
            title = "Observed/theoretical isotope ion ratios (>3 or <0.3 marked red)",
            xaxis = list(title = "Sample"),
            yaxis = list(title = "Observed/Theoretical ratio", type = "log")
        )
}

#-------------------------------
plot_sample_contribution <- function(deconvolution) {

    # How much contribution of each sample to the final deconvoluted homologue group pattern
    plot_data <- deconvolution |>
        tidyr::unnest(deconv_coef) |>
        tidyr::unnest_longer(c(deconv_coef, Batch_Name)) |>
        dplyr::select(Replicate_Name, Batch_Name, deconv_coef)

    # Create the plotly stacked bar plot
    plotly::plot_ly(plot_data,
                    x = ~Replicate_Name,
                    y = ~deconv_coef,
                    type = "bar",
                    color = ~Batch_Name,
                    colors = "Spectral") |>
        plotly::layout(
            title = list(
                text = "Contributions from standards to deconvoluted homologue pattern",
                x = 0.5,
                y = 0.95
            ),
            barmode = "stack",
            xaxis = list(title = "Replicate Name"),
            yaxis = list(title = "Relative Contribution",
                         tickformat = ".2%"),
            showlegend = TRUE,
            legend = list(title = list(text = "Batch Name"))
        )
}

#-------------------------------
plot_homologue_group_pattern_comparison <- function(Sample_distribution, input_selectedSamples){

    # Filter data for selected samples and reshape data
    selected_samples <- Sample_distribution |>
        dplyr::filter(Replicate_Name %in% input_selectedSamples) |>
        dplyr::mutate(Molecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)])))

    # Get unique homologue groups for consistent coloring
    homologue_groups <- unique(selected_samples$C_homologue)

    # Create a list of plots, one for each Replicate_Name
    plot_list <- selected_samples |>
        split(selected_samples$Replicate_Name) |>
        purrr::map(function(df) {
            # Create base plot
            p <- plotly::plot_ly()

            # Add bars for each homologue group
            for(hg in homologue_groups) {
                df_filtered <- df[df$C_homologue == hg,]
                p <- p |>
                    plotly::add_trace(
                        data = df_filtered,
                        x = ~Molecule,
                        y = ~Relative_Area,
                        name = hg,
                        legendgroup = hg,
                        showlegend = (df_filtered$Replicate_Name[1] == input_selectedSamples[1]),
                        type = 'bar',
                        opacity = 1
                    )
            }

            # Add the black line for resolved_distribution
            p <- p |>
                plotly::add_trace(
                    data = df,
                    x = ~Molecule,
                    y = ~resolved_distribution,
                    name = "Deconvoluted Distribution",
                    legendgroup = "DeconvDistr",
                    showlegend = (df$Replicate_Name[1] == input_selectedSamples[1]),
                    type = 'scatter',
                    mode = 'lines+markers',
                    line = list(color = 'black'),
                    marker = list(color = 'black', size = 6),
                    opacity = 0.7
                )

            # Add layout
            p <- p |>
                plotly::layout(
                    xaxis = list(
                        title = "Homologue",
                        tickangle = 45
                    ),
                    yaxis = list(title = "Value"),
                    barmode = 'group',
                    annotations = list(
                        x = 0.5,
                        y = 1.1,
                        text = unique(df$Replicate_Name),
                        xref = 'paper',
                        yref = 'paper',
                        showarrow = FALSE
                    )
                )

            return(p)
        })

    # Combine the plots using subplot
    plotly::subplot(plot_list,
                    nrows = ceiling(length(plot_list)/2),
                    shareX = TRUE,
                    shareY = TRUE) |>
        plotly::layout(
            showlegend = TRUE,
            hovermode = 'closest',
            hoverlabel = list(bgcolor = "white"),
            barmode = 'group'
        ) |>
        plotly::config(displayModeBar = TRUE) |>
        htmlwidgets::onRender("
      function(el) {
        var plotDiv = document.getElementById(el.id);
        plotDiv.on('plotly_legendclick', function(data) {
          Plotly.restyle(plotDiv, {
            visible: data.data[data.curveNumber].visible === 'legendonly' ? true : 'legendonly'
          }, data.fullData.map((trace, i) => i).filter(i =>
            data.fullData[i].legendgroup === data.fullData[data.curveNumber].legendgroup
          ));
          return false;
        });
      }
    ")
}

# =========================================================
# Deconvolution function
# =========================================================

perform_deconvolution <- function(df, combined_standard, CPs_standards_sum_RF,
                                  sample_name = NA_character_) {

    df_matrix <- as.matrix(df)

    message(paste("df_matrix dimensions:", paste(dim(df_matrix), collapse = " x ")))
    message(paste("combined_standard dimensions:", paste(dim(combined_standard), collapse = " x ")))

    # For messages
    sample_label <- if (is.na(sample_name)) "<unknown sample>" else sample_name

    # Build df_vector
    if (ncol(df_matrix) == 1) {
        df_vector <- as.vector(df_matrix)

    } else {

        if (!"Relative_Area" %in% colnames(df_matrix)) {
            stop(sprintf(
                paste0(
                    "Deconvolution input error for sample '%s':\n",
                    "Column 'Relative_Area' not found in data passed to perform_deconvolution().\n",
                    "This usually means that the Relative_Area column was not created correctly ",
                    "before deconvolution (check that Isotope_Label_Type == 'Quan' exists and Area is numeric)."
                ),
                sample_label
            ))
        }

        df_vector <- df_matrix[, "Relative_Area", drop = TRUE]
    }

    # Dimension compatibility: rows (molecules) of combined_standard must match length of df_vector
    if (nrow(combined_standard) != length(df_vector)) {
        stop(sprintf(
            paste0(
                "Deconvolution failed for sample '%s': incompatible pattern dimensions.\n\n",
                "  • Length of sample relative pattern = %d\n",
                "  • Number of molecules in standards matrix = %d\n\n",
                "This usually means that the set of homologues (Molecule) in the samples\n",
                "does not match the homologues in the standards.\n",
                "Check that:\n",
                "  - Your Skyline export contains the same homologue groups in standards and samples, and\n",
                "  - You might have removed Molecules in Skyline that were in the original transition list."
            ),
            sample_label, length(df_vector), nrow(combined_standard)
        ))
    }

    # Validate inputs for nnls
    if (any(!is.finite(df_vector))) {
        stop(sprintf(
            "Deconvolution input error for sample '%s': df_vector contains NA / NaN / Inf values.",
            sample_label
        ))
    }
    if (any(!is.finite(combined_standard))) {
        stop(
            "Deconvolution input error: 'combined_standard' matrix contains NA / NaN / Inf values."
        )
    }

    # ---- existing nnls and output code stays the same ----
    deconv <- nnls::nnls(combined_standard, df_vector)

    deconv_coef <- deconv$x
    if (sum(deconv_coef) > 0) {
        deconv_coef <- deconv_coef / sum(deconv_coef)
    } else {
        deconv_coef[] <- 0
    }

    deconv_resolved <- as.vector(combined_standard %*% deconv_coef)
    names(deconv_resolved) <- rownames(combined_standard)

    sum_deconv_RF <- as.numeric(as.matrix(CPs_standards_sum_RF) %*% deconv_coef)

    sst <- sum((df_vector - mean(df_vector))^2)

    pred_norm_den <- sum(deconv_resolved)
    pred_norm <- if (is.finite(pred_norm_den) && pred_norm_den > 0) {
        deconv_resolved / pred_norm_den
    } else {
        rep(0, length(deconv_resolved))
    }

    ssr <- sum((df_vector - pred_norm)^2)
    deconv_rsquared <- if (sst > 0) 1 - (ssr / sst) else NA_real_

    names(deconv_coef) <- colnames(combined_standard)

    list(
        sum_deconv_RF   = sum_deconv_RF,
        deconv_coef     = deconv_coef,
        deconv_resolved = deconv_resolved,
        deconv_rsquared = deconv_rsquared
    )
}
