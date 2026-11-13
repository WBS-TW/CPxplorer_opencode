#' CPquant: Shiny CPs Quantification for Skyline Output
#' @param ...
#'
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @import htmlwidgets
#' @rawNamespace import(ggplot2, except=c(last_plot))
#' @import readxl
#' @import nnls
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import DT
#' @import plotly
#' @import purrr
#' @import markdown
#' @import openxlsx
#' @import stringr
#' @importFrom stats lm
#' @export


CPquant <- function(...){


    # =========================================================
    # UI
    # =========================================================

    options(shiny.maxRequestSize = 500 * 1024^2)

    ui <- shiny::navbarPage(
        "CPquant",
        shiny::tabPanel(
            "Quantification Inputs",
            shiny::fluidPage(
                shiny::sidebarLayout(
                    shiny::sidebarPanel(
                        width = 3,
                        shiny::fileInput("fileInput", "Import excel file from Skyline",
                                         accept = c('xlsx')),
                        shiny::textInput("quanUnit", "(Optional) Concentration unit:"),
                        shiny::radioButtons("quanSum",
                                            label = "Choose ions for quantification",
                                            choices = c("Quan only", "Sum Quan+Qual"),
                                            selected = "Quan only"),
                        shiny::radioButtons("blankSubtraction",
                                            label = "Subtraction with blank?",
                                            choices = c("Yes, by avg area of blanks", "No"),
                                            selected = "No"),
                        shiny::radioButtons("correctWithRS",
                                            label = "Correct with RS area?",
                                            choices = c("Yes", "No"),
                                            selected = "No"),
                        shiny::uiOutput("correctionUI"), # render UI if correctwithRS == "Yes"
                        shiny::radioButtons("calculateRecovery",
                                            label = "Calculate recovery? (req QC samples)",
                                            choices = c("Yes", "No"),
                                            selected = "No"),
                        shiny::radioButtons("calculateMDL",
                                            label = "Calculate MDL? (req blank samples)",
                                            choices = c("Yes", "No"),
                                            selected = "No"),
                        shiny::radioButtons("standardTypes",
                                            label = "Types of standards",
                                            choices = c("Group Mixtures"),
                                            selected = "Group Mixtures"),
                        shiny::radioButtons("calibrationSource",
                                            label = "Calibration curve source",
                                            choices = c("From Skyline file",
                                                        "Upload saved calibration (.xlsx)"),
                                            selected = "From Skyline file"),
                        shiny::conditionalPanel(
                            condition = "input.calibrationSource == 'Upload saved calibration (.xlsx)'",
                            shiny::fileInput("calibrationFile",
                                             "Import saved calibration (.xlsx)",
                                             accept = c("xlsx", ".xlsx"))
                        ),
                        shiny::tags$div(
                            title = "Wait until import file is fully loaded before pressing!",
                            shiny::actionButton('go', 'Proceed', width = "100%")
                        ),
                        shiny::uiOutput("defineVariables")
                    ),
                    shiny::mainPanel(
                        width = 9,
                        plotly::plotlyOutput("plot_Skyline_output"),
                        DT::DTOutput("table_Skyline_output")
                    )
                )
            )
        ),
        shiny::tabPanel(
            "Input Summary",
            shiny::sidebarPanel(
                width = 2,
                shiny::radioButtons("navSummary", "Choose tab:",
                                    choices = c("Included Standards",
                                                "Removed from Calibration",
                                                "Quan to Qual ratio",
                                                "Measured vs Theor Quan/Qual ratio"),
                                    selected = "Included Standards"),
                shiny::tags$hr(),
                shiny::tags$p("'Quan to Qual ratio' and 'Measured vs Theor Quan/Qual ratio' only works for Quan only option",
                              style = "font-style: italic;")
            ),
            shiny::mainPanel(
                width = 10,
                shiny::conditionalPanel(
                    condition = "input.navSummary == 'Included Standards'",
                    shiny::tags$h3("Standards included in calibration curves"),
                    shiny::downloadButton("downloadCalibration", "Save Calibration Curve (.xlsx)"),
                    shiny::tags$br(), shiny::tags$br(),
                    DT::DTOutput("CalibrationIncluded")
                ),
                shiny::conditionalPanel(
                    condition = "input.navSummary == 'Removed from Calibration'",
                    shiny::tags$h3("Calibration series removed from quantification"),
                    DT::DTOutput("CalibrationRemoved")
                ),
                shiny::conditionalPanel(
                    condition = "input.navSummary == 'Quan to Qual ratio'",
                    shiny::tags$h3("Violin plots of Quant/Qual ions"),
                    plotly::plotlyOutput("RatioQuantToQual", height = "80vh", width = "100%")
                ),
                shiny::conditionalPanel(
                    condition = "input.navSummary == 'Measured vs Theor Quan/Qual ratio'",
                    shiny::tags$h3("Measured divided by Theoretical Quant/Qual ratios"),
                    plotly::plotlyOutput("MeasVSTheor", height = "80vh", width = "100%")
                )
            )
        ),
        shiny::tabPanel(
            "Quantification Summary",
            shiny::fluidPage(
                downloadButton("downloadResults", "Export all results to Excel"),
                shiny::tags$br(), shiny::tags$br(),
                DT::DTOutput("quantTable")
            )
        ),
        shiny::tabPanel(
            "Standard contributions",
            shiny::fluidPage(
                plotly::plotlyOutput("sampleContributionPlot")
            )
        ),
        shiny::tabPanel(
            "Homologue Group Patterns",
            shiny::fluidRow(
                shiny::column(
                    width = 2,
                    shiny::radioButtons("plotHomologueGroups", "Choose tab:",
                                        choices = c("All Samples Overview", "Samples Overlay", "Samples Panels"),
                                        selected = "All Samples Overview"),
                    shiny::conditionalPanel(
                        condition = "input.plotHomologueGroups == 'Samples Overlay'",
                        shiny::uiOutput("sampleSelectionUIOverlay")
                    ),
                    shiny::conditionalPanel(
                        condition = "input.plotHomologueGroups == 'Samples Panels'",
                        shiny::uiOutput("sampleSelectionUIComparisons")
                    ),
                    shiny::tags$div(
                        title = "WAIT after pressing..Might take some time before plot shows!",
                        shiny::actionButton('go2', 'Plot', width = "100%")
                    )
                ),
                shiny::column(
                    width = 10,
                    shiny::tags$div(
                        style = "margin-top: 20px;",
                        shiny::conditionalPanel(
                            condition = "input.plotHomologueGroups == 'All Samples Overview'",
                            shiny::plotOutput("plotHomologuePatternStatic", height = "80vh", width = "100%")
                        ),
                        shiny::conditionalPanel(
                            condition = "input.plotHomologueGroups == 'Samples Overlay'",
                            plotly::plotlyOutput("plotHomologuePatternOverlay", height = "80vh", width = "100%")
                        ),
                        shiny::conditionalPanel(
                            condition = "input.plotHomologueGroups == 'Samples Panels'",
                            plotly::plotlyOutput("plotHomologuePatternComparisons", height = "80vh", width = "100%")
                        )
                    )
                )
            )
        ),
        shiny::tabPanel(
            "QA/QC",
            shiny::sidebarLayout(
                shiny::sidebarPanel(
                    width = 2,
                    shiny::radioButtons(
                        "navQAQC",
                        "Choose tab:",
                        choices = c(
                            "Recovery and MDL",
                            "Calibration RF Check",
                            "Isotope Similarity Heatmap",
                            "Isotope QC Summary",
                            "Isotope Pattern Overlay"

                        ),
                        selected = "Recovery and MDL"
                    ),
                    shiny::conditionalPanel(
                        condition = "input.navQAQC == 'Isotope Pattern Overlay'",
                        shiny::uiOutput("isotopeQCSampleUI"),
                        shiny::uiOutput("isotopeQCMoleculeUI")
                    )
                ),
                shiny::mainPanel(
                    width = 10,
                    shiny::conditionalPanel(
                        condition = "input.navQAQC == 'Recovery and MDL'",
                        DT::DTOutput("table_recovery"),
                        shiny::br(),
                        DT::DTOutput("table_MDL")
                    ),
                    shiny::conditionalPanel(
                        condition = "input.navQAQC == 'Calibration RF Check'",
                        shiny::uiOutput("calibrationRFCheckStatus"),
                        shiny::tags$h4("Batch-level Response Factor Comparison (Sum RF)"),
                        DT::DTOutput("table_cal_rf_batch_qc"),
                        shiny::br(),
                        shiny::tags$h4("Homologue-level Response Factor Comparison"),
                        DT::DTOutput("table_cal_rf_homologue_qc")
                    ),
                    shiny::conditionalPanel(
                        condition = "input.navQAQC == 'Isotope QC Summary'",
                        DT::DTOutput("table_isotope_qc")
                    ),
                    shiny::conditionalPanel(
                        condition = "input.navQAQC == 'Isotope Pattern Overlay'",
                        plotly::plotlyOutput("plotIsotopePatternOverlay", height = "80vh", width = "100%")
                    ),
                    shiny::conditionalPanel(
                        condition = "input.navQAQC == 'Isotope Similarity Heatmap'",
                        plotly::plotlyOutput("plotIsotopeSimilarityHeatmap", height = "80vh", width = "100%")
                    ),
                )
            )
        ),
        shiny::tabPanel(
            "Instructions",
            shiny::sidebarLayout(
                shiny::sidebarPanel(shiny::h3("Manual"), width = 3),
                shiny::mainPanel(
                    shiny::includeMarkdown(system.file("instructions_CPquant.md", package = "CPxplorer"))
                )
            )
        ),
        shiny::tabPanel(
            "Log",
            shiny::fluidPage(
                shiny::verbatimTextOutput("console_log"),
                shiny::tags$style("#console_log {max-height: 70vh;
                                    overflow-y: auto;
                                    background: #111827; /* dark */
                                    color: #e5e7eb;      /* light */
                                    font-family: monospace;
                                    font-size: 12px;
                                    padding: 12px;
                                    border-radius: 6px;}")
            )
        )
    )

    # =========================================================
    # Server
    # =========================================================

    server <- function(input, output, session) {



        # ---------------
        # In-memory console sink setup
        # ---------------
        log_buf <- character()
        if (!inherits(session, "MockShinySession")) {
            tc <- textConnection("log_buf", open = "w", local = TRUE)
            sink(tc, append = TRUE)                     # stdout
            sink(tc, append = TRUE, type = "message")   # stderr (conditions)
        }
        max_lines <- 600
        output$console_log <- shiny::renderText({
            shiny::invalidateLater(500, session)
            if (length(log_buf) > max_lines) {
                log_buf <<- tail(log_buf, max_lines)
            }
            paste(rev(log_buf), collapse = "\n")
        })
        # ---------------

        # ---- Quant unit passthrough ----
        quantUnit <- reactive({ input$quanUnit })

        # ---- 1) Load + minimal tidy (independent of removeSamples) ----
        base_df <- reactive({
            req(input$fileInput)

            df <- read_cpquant_skyline(input$fileInput$datapath, input$standardTypes)

            df
        })

        # ---- 2) Dynamic UI from base_df() ----
        output$defineVariables <- shiny::renderUI({
            req(base_df())
            defineVariablesUI(base_df())
        })

        output$correctionUI <- shiny::renderUI({
            req(base_df())
            if (input$correctWithRS == "Yes") {
                defineCorrectionUI(base_df())
            }
        })

        # ---- 3) Full processing with removal applied EARLY ----
        Skyline_output <- reactive({
            req(base_df())

            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "WAIT! Loading data...", value = 0)

            progress$set(value = 0.3, detail = "Starting from base data")
            df <- base_df()

            # Apply removal of samples BEFORE any downstream stats (RS, blanks, etc.)
            rmv <- input$removeSamples
            if (!is.null(rmv) && length(rmv) > 0) {
                df <- df |>
                    dplyr::filter(!Replicate_Name %in% rmv)
            }

            progress$set(value = 0.5, detail = "Applying quantification settings")

            df <- apply_cpquant_preprocessing(
                df,
                remove_samples = input$removeSamples,
                quan_sum = input$quanSum,
                correct_with_rs = input$correctWithRS,
                choose_rs = input$chooseRS,
                blank_subtraction = input$blankSubtraction,
                standard_types = input$standardTypes
            )

            progress$set(value = 1, detail = "Complete")
            print("Upload successful!")
            df
        })

        # ---- Event-gated input (for R^2) ----
        removeRsquared <- shiny::eventReactive(input$go, { as.numeric(input$removeRsquared) })
        Samples_Concentration <- reactiveVal() # store deconvolution for other tabs
        Isotope_QC <- reactiveVal(NULL)
        Calibration_RF_QC <- reactiveVal(NULL)

        # ---- Raw table & plot (reflect auto removal) ----
        output$table_Skyline_output <- DT::renderDT({
            DT::datatable(Skyline_output(),
                          options = list(paging = TRUE, pageLength = 50))
        })

        output$plot_Skyline_output <- plotly::renderPlotly({
            plot_skyline_output(Skyline_output())
        })

        # =========================================================
        # Deconvolution pipeline — runs on Proceed button
        # =========================================================
        shiny::observeEvent(input$go, {

            progress <- shiny::Progress$new()
            on.exit(progress$close())

            progress$set(message = "WAIT! Processing data...", value = 0)

            # Skyline_output already reflects sample removals
            Skyline_output_filt <- Skyline_output()


            #######################################
            ##### PREPARE FOR DECONVOLUTION #######
            #######################################
            progress$set(value = 0.2, detail = "Preparing standards data")

            if (input$standardTypes == "Group Mixtures") {
                using_saved_cal <- identical(input$calibrationSource, "Upload saved calibration (.xlsx)")

                if (using_saved_cal) {
                    if (is.null(input$calibrationFile) || is.null(input$calibrationFile$datapath)) {
                        showModal(
                            modalDialog(
                                title = "Missing saved calibration file",
                                "Please upload a saved calibration (.xlsx) file or switch 'Calibration curve source' to 'From Skyline file'.",
                                easyClose = TRUE,
                                footer = modalButton("OK")
                            )
                        )
                        return()
                    }

                    loaded_cal <- tryCatch({
                        read_cpquant_calibration_xlsx(
                            path = input$calibrationFile$datapath,
                            min_rsquared = removeRsquared(),
                            quan_sum = input$quanSum,
                            correct_with_rs = input$correctWithRS,
                            choose_rs = input$chooseRS,
                            standard_types = input$standardTypes
                        )
                    }, error = function(e) {
                        showModal(
                            modalDialog(
                                title = "Error loading saved calibration file",
                                tagList(
                                    p("Could not load the uploaded calibration (.xlsx) file."),
                                    p("Details:"),
                                    tags$pre(conditionMessage(e))
                                ),
                                easyClose = TRUE,
                                footer = modalButton("OK")
                            )
                        )
                        NULL
                    })

                    if (is.null(loaded_cal)) {
                        return()
                    }

                    CPs_standards <- loaded_cal$standards
                    raw_standards_for_export <- if (nrow(loaded_cal$raw_standards) > 0) {
                        loaded_cal$raw_standards
                    } else {
                        base_df() |>
                            dplyr::filter(
                                Sample_Type == "Standard",
                                !Replicate_Name %in% input$removeSamples
                            )
                    }
                } else {
                    # All Molecule that does not meet RF, cal_rsquared criteria will be mutated to RF=0 (not included into calibration)
                    CPs_standards <- prepare_cpquant_standards(
                        Skyline_output_filt,
                        min_rsquared = removeRsquared()
                    )
                    raw_standards_for_export <- base_df() |>
                        dplyr::filter(
                            Sample_Type == "Standard",
                            !Replicate_Name %in% input$removeSamples
                        )
                }

                cal_rf_qc <- compare_calibration_rf_qc(CPs_standards, Skyline_output_filt)
                Calibration_RF_QC(cal_rf_qc)

                output$downloadCalibration <- shiny::downloadHandler(
                    filename = function() {
                        paste("CPquant_Calibration_", Sys.Date(), ".xlsx", sep = "")
                    },
                    content = function(file) {
                        write_cpquant_calibration_xlsx(
                            path = file,
                            standards_tbl = CPs_standards,
                            raw_standards_df = raw_standards_for_export,
                            metadata = list(
                                created_date = as.character(Sys.Date()),
                                standardTypes = input$standardTypes,
                                quanSum = input$quanSum,
                                blankSubtraction = input$blankSubtraction,
                                correctWithRS = input$correctWithRS,
                                chooseRS = if (identical(input$correctWithRS, "Yes") && !is.null(input$chooseRS)) input$chooseRS else "",
                                removeRsquared = removeRsquared(),
                                quanUnit = quantUnit()
                            )
                        )
                    }
                )

                progress$set(value = 0.6, detail = "Preparing sample data")

                CPs_samples <- Skyline_output_filt |>
                    dplyr::filter(
                        Sample_Type %in% c("Unknown", "Blank"),
                        !Molecule_List %in% c("IS", "RS", "VS"),
                        Isotope_Label_Type == "Quan"
                    ) |>
                    dplyr::group_by(Replicate_Name) |>
                    dplyr::mutate(Relative_Area = Area / sum(Area, na.rm = TRUE)) |>
                    dplyr::ungroup() |>
                    dplyr::select(-Mass_Error_PPM, -Isotope_Label_Type, -Chromatogram_Precursor_MZ, -Analyte_Concentration, -Batch_Name) |>
                    dplyr::mutate(dplyr::across(Relative_Area, ~ replace(., is.nan(.), 0)))

                # Combined sample (nested) for per-sample deconvolution
                combined_sample <- CPs_samples |>
                    dplyr::group_by(Replicate_Name, Sample_Type) |>
                    tidyr::nest() |>
                    dplyr::ungroup()

                ###### Table for Calibration Curves ######
                output$CalibrationIncluded <- DT::renderDT({
                    req(CPs_standards)

                    # Fully flatten CPs_standards (removes list-cols, ensures every column is atomic)
                    CPs_standards |>
                        dplyr::ungroup() |>
                        filter(RF > 0) |>
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
                        ) |>
                        dplyr::mutate(
                            RF = round(as.numeric(RF), 7),
                            intercept = round(as.numeric(intercept), 3),
                            cal_rsquared = round(as.numeric(cal_rsquared), 3)
                        ) |>
                        DT::datatable(options = list(pageLength = 40))

                })


                output$CalibrationRemoved <- DT::renderDT({
                    req(CPs_standards)
                    CPs_standards |>
                        dplyr::ungroup() |>
                        dplyr::filter(RF <= 0) |>
                        dplyr::mutate(RF_display = if ("RF_raw" %in% names(CPs_standards)) as.numeric(RF_raw) else as.numeric(RF)) |>
                        dplyr::select(Batch_Name, Molecule, Quantification_Group, RF = RF_display, intercept, cal_rsquared) |>
                        dplyr::mutate(across(where(is.numeric), ~ signif(.x, digits = 4))) |>
                        DT::datatable(options = list(pageLength = 40))
                })

                ###### Calibration RF Check QA/QC outputs ######
                output$calibrationRFCheckStatus <- shiny::renderUI({
                    qc <- Calibration_RF_QC()
                    if (!using_saved_cal) {
                        shiny::tags$div(
                            style = "padding: 10px 14px; margin-bottom: 12px; background-color: #eef2ff; border-left: 4px solid #4f46e5; border-radius: 4px;",
                            shiny::tags$strong("Note: "),
                            "Calibration curves were built from the current Skyline file. When 'Upload saved calibration (.xlsx)' is selected and the Skyline file includes verification standards with matching Batch Name, this section compares saved vs. current Response Factors (RF = Area / Concentration; Pass <= 20%, Review <= 30%, Fail > 30%)."
                        )
                    } else if (is.null(qc) || nrow(qc$batch_qc) == 0) {
                        shiny::tags$div(
                            style = "padding: 10px 14px; margin-bottom: 12px; background-color: #fff3cd; border-left: 4px solid #ffc107; border-radius: 4px;",
                            shiny::tags$strong("No matching check standards found: "),
                            "Using saved calibration file (",
                            shiny::tags$code(input$calibrationFile$name),
                            "), but the current Skyline file does not contain any Sample_Type == 'Standard' rows with a matching Batch Name and Analyte Concentration > 0."
                        )
                    } else {
                        shiny::tags$div(
                            style = "padding: 10px 14px; margin-bottom: 12px; background-color: #d1e7dd; border-left: 4px solid #198754; border-radius: 4px;",
                            shiny::tags$strong("Saved Calibration QA/QC: "),
                            "Comparing saved calibration (",
                            shiny::tags$code(input$calibrationFile$name),
                            ") against verification standards in the current Skyline run (RF = Area / Concentration; Pass: |% diff| <= 20%, Review: <= 30%, Fail: > 30%)."
                        )
                    }
                })

                output$table_cal_rf_batch_qc <- DT::renderDT({
                    shiny::req(Calibration_RF_QC())
                    Calibration_RF_QC()$batch_qc |>
                        dplyr::mutate(
                            Sum_RF_saved = signif(Sum_RF_saved, 5),
                            Sum_RF_current = signif(Sum_RF_current, 5),
                            RF_ratio = round(RF_ratio, 3),
                            pct_diff = round(pct_diff, 1)
                        ) |>
                        DT::datatable(
                            filter = "top",
                            extensions = c("Buttons", "Scroller"),
                            options = list(
                                scrollY = 300, scrollX = 500, deferRender = TRUE, scroller = TRUE,
                                buttons = list(
                                    list(extend = "excel", filename = "Calibration_RF_Batch_QC", title = NULL,
                                         exportOptions = list(modifier = list(page = "all"))),
                                    list(extend = "csv", filename = "Calibration_RF_Batch_QC", title = NULL,
                                         exportOptions = list(modifier = list(page = "all")))
                                ),
                                dom = "lBfrtip"
                            ),
                            rownames = FALSE
                        ) |>
                        DT::formatStyle(
                            "qc_flag",
                            backgroundColor = DT::styleEqual(
                                c("Pass", "Review", "Fail", "No signal"),
                                c("#d1e7dd", "#fff3cd", "#f8d7da", "#e2e3e5")
                            ),
                            fontWeight = "600"
                        )
                })

                output$table_cal_rf_homologue_qc <- DT::renderDT({
                    shiny::req(Calibration_RF_QC())
                    Calibration_RF_QC()$homologue_qc |>
                        dplyr::mutate(
                            RF_saved = signif(RF_saved, 5),
                            RF_current = signif(RF_current, 5),
                            RF_ratio = round(RF_ratio, 3),
                            pct_diff = round(pct_diff, 1)
                        ) |>
                        DT::datatable(
                            filter = "top",
                            extensions = c("Buttons", "Scroller"),
                            options = list(
                                scrollY = 500, scrollX = 500, deferRender = TRUE, scroller = TRUE,
                                buttons = list(
                                    list(extend = "excel", filename = "Calibration_RF_Homologue_QC", title = NULL,
                                         exportOptions = list(modifier = list(page = "all"))),
                                    list(extend = "csv", filename = "Calibration_RF_Homologue_QC", title = NULL,
                                         exportOptions = list(modifier = list(page = "all")))
                                ),
                                dom = "lBfrtip"
                            ),
                            rownames = FALSE
                        ) |>
                        DT::formatStyle(
                            "qc_flag",
                            backgroundColor = DT::styleEqual(
                                c("Pass", "Review", "Fail", "No signal"),
                                c("#d1e7dd", "#fff3cd", "#f8d7da", "#e2e3e5")
                            ),
                            fontWeight = "600"
                        )
                })

                ###### Plots: Quan/Qual ratios ######
                output$RatioQuantToQual <- plotly::renderPlotly({
                    plot_quanqualratio(Skyline_output_filt)
                })
                output$MeasVSTheor <- plotly::renderPlotly({
                    plot_meas_vs_theor_ratio(Skyline_output_filt)
                })

                isotope_qc_input <- Skyline_output_filt
                rmv_qc <- input$removeSamples
                if (!is.null(rmv_qc) && length(rmv_qc) > 0) {
                    isotope_qc_input <- isotope_qc_input |>
                        dplyr::filter(!Replicate_Name %in% rmv_qc)
                }

                isotope_qc <- prepare_isotope_pattern_qc(isotope_qc_input)
                Isotope_QC(isotope_qc)

                output$table_isotope_qc <- DT::renderDT({
                    shiny::req(Isotope_QC())
                    Isotope_QC()$summary_qc |>
                        dplyr::mutate(
                            cosine_similarity = round(cosine_similarity, 3),
                            pearson_r = round(pearson_r, 3),
                            weighted_abs_percent_error = round(weighted_abs_percent_error, 2),
                            max_ion_ratio_error = round(max_ion_ratio_error, 2),
                            total_area = round(total_area, 0)
                        ) |>
                        DT::datatable(
                            filter = "top",
                            extensions = c("Buttons", "Scroller"),
                            options = list(
                                scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE,
                                buttons = list(
                                    list(extend = "excel", filename = "Isotope_pattern_QC", title = NULL,
                                         exportOptions = list(modifier = list(page = "all"))),
                                    list(extend = "csv", filename = "Isotope_pattern_QC", title = NULL,
                                         exportOptions = list(modifier = list(page = "all"))),
                                    list(extend = "colvis", targets = 0, visible = FALSE)
                                ),
                                dom = "lBfrtip",
                                fixedColumns = TRUE
                            ),
                            rownames = FALSE
                        ) |>
                        DT::formatStyle(
                            "qc_flag",
                            backgroundColor = DT::styleEqual(c("Pass", "Review", "Fail", "No signal"),
                                                             c("#d1e7dd", "#fff3cd", "#f8d7da", "#e2e3e5")),
                            fontWeight = "600"
                        )
                })

                output$isotopeQCSampleUI <- shiny::renderUI({
                    req(Isotope_QC())
                    samples <- unique(Isotope_QC()$ion_qc$Replicate_Name)
                    req(length(samples) > 0)
                    shiny::selectInput("selectedIsotopeQCSample", "Sample:", choices = samples, selected = samples[[1]])
                })

                output$isotopeQCMoleculeUI <- shiny::renderUI({
                    req(Isotope_QC(), input$selectedIsotopeQCSample)
                    molecules <- Isotope_QC()$ion_qc |>
                        dplyr::filter(Replicate_Name == input$selectedIsotopeQCSample) |>
                        dplyr::pull(Molecule) |>
                        unique()
                    req(length(molecules) > 0)
                    shiny::selectInput("selectedIsotopeQCMolecule", "Molecule:", choices = molecules, selected = molecules[[1]])
                })

                output$plotIsotopePatternOverlay <- plotly::renderPlotly({
                    req(Isotope_QC(), input$selectedIsotopeQCSample, input$selectedIsotopeQCMolecule)
                    plot_isotope_pattern_overlay(
                        Isotope_QC()$ion_qc,
                        input$selectedIsotopeQCSample,
                        input$selectedIsotopeQCMolecule
                    )
                })

                output$plotIsotopeSimilarityHeatmap <- plotly::renderPlotly({
                    req(Isotope_QC())
                    plot_isotope_similarity_heatmap(Isotope_QC()$summary_qc)
                })

                # ---- DECONVOLUTION ----
                progress$set(value = 0.8, detail = "Performing deconvolution")

                CPs_standards_input <- CPs_standards |>
                    dplyr::select(Molecule, Batch_Name, RF) |>
                    tidyr::pivot_wider(names_from = Batch_Name, values_from = "RF") |>
                    dplyr::mutate(across(everything(), ~ tidyr::replace_na(., 0)))

                CPs_standards_sum_RF <- CPs_standards |>
                    dplyr::select(Batch_Name, Sum_RF_group) |>
                    dplyr::distinct() |>
                    dplyr::ungroup() |>
                    tidyr::pivot_wider(names_from = Batch_Name, values_from = "Sum_RF_group") |>
                    dplyr::mutate(across(everything(), ~ tidyr::replace_na(., 0)))

                # Matrix with molecules as rows, batches as columns
                combined_standard <- CPs_standards_input |>
                    tibble::column_to_rownames(var = "Molecule") |>
                    as.matrix()

                std_molecules <- rownames(combined_standard)
                combined_sample <- combined_sample |>
                    dplyr::mutate(
                        data = purrr::map(data, function(sample_df) {
                            if (nrow(sample_df) == length(std_molecules) &&
                                setequal(sample_df$Molecule, std_molecules)) {
                                sample_df[match(std_molecules, sample_df$Molecule), , drop = FALSE]
                            } else {
                                sample_df
                            }
                        })
                    )

                # ---- DECONVOLUTION with UI error popup ----
                deconvolution <- tryCatch({

                    combined_sample |>
                        dplyr::mutate(result = purrr::map(
                            data,
                            ~ perform_deconvolution(
                                dplyr::select(.x, Relative_Area),
                                combined_standard,
                                CPs_standards_sum_RF
                            )
                        )) |>
                        dplyr::mutate(sum_Area = purrr::map_dbl(data, ~ sum(.x$Area))) |>
                        dplyr::mutate(sum_deconv_RF = as.numeric(purrr::map(result, purrr::pluck("sum_deconv_RF")))) |>
                        dplyr::mutate(Sample_Dilution_Factor = purrr::map_dbl(data, ~ dplyr::first(.x$Sample_Dilution_Factor))) |>
                        dplyr::mutate(
                            Concentration = dplyr::if_else(sum_deconv_RF > 0,
                                                           sum_Area / sum_deconv_RF * Sample_Dilution_Factor,
                                                           NA_real_)
                        ) |>
                        dplyr::mutate(Unit = quantUnit()) |>
                        dplyr::mutate(
                            deconv_coef = purrr::map(result, ~ tibble::tibble(
                                Batch_Name  = names(.x$deconv_coef),
                                deconv_coef = as.numeric(.x$deconv_coef)
                            )),
                            deconv_rsquared = as.numeric(purrr::map(result, purrr::pluck("deconv_rsquared"))),
                            deconv_resolved = purrr::map(result, ~ tibble::tibble(
                                Molecule         = names(.x$deconv_resolved),
                                deconv_resolved  = as.numeric(.x$deconv_resolved)
                            ))
                        ) |>
                        dplyr::select(-result)

                }, error = function(e) {

                    # Show the error as a popup in the Shiny UI
                    showModal(
                        modalDialog(
                            title = "Error during quantification / deconvolution",
                            tagList(
                                p("The quantification stopped due to an error in the deconvolution step."),
                                p("Details:"),
                                tags$pre(conditionMessage(e))
                            ),
                            easyClose = TRUE,
                            footer = modalButton("OK")
                        )
                    )

                    # Return NULL so we can handle it gracefully outside
                    return(NULL)
                })

                # If deconvolution failed, stop the rest of this observer
                if (is.null(deconvolution)) {
                    return()
                }

                progress$set(value = 0.9, detail = "Calculating final results")

                # Store for downstream tabs
                Samples_Concentration(deconvolution)

                progress$set(value = 1, detail = "Complete")

                # ---- Download handler (Excel) ----
                output$downloadResults <- shiny::downloadHandler(
                    filename = function() {
                        paste("CPquant_Results_", Sys.Date(), ".xlsx", sep = "")
                    },
                    content = function(file) {
                        wb <- openxlsx::createWorkbook()

                        # Sheets
                        openxlsx::addWorksheet(wb, "Quantification")
                        openxlsx::addWorksheet(wb, "StandardsContribution")
                        openxlsx::addWorksheet(wb, "HomologueDistribution")

                        # Quantification
                        openxlsx::writeData(
                            wb, "Quantification",
                            deconvolution |>
                                dplyr::select(Replicate_Name, Sample_Type, Concentration, Unit, deconv_rsquared) |>
                                dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3))
                        )

                        # StandardsContribution
                        openxlsx::writeData(
                            wb, "StandardsContribution",
                            deconvolution |>
                                dplyr::select(Replicate_Name, deconv_coef) |>
                                tidyr::unnest(deconv_coef) |>
                                dplyr::mutate(deconv_coef = deconv_coef * 100) |>
                                tidyr::pivot_wider(names_from = Batch_Name, values_from = deconv_coef)
                        )

                        # HomologueDistribution
                        openxlsx::writeData(
                            wb, "HomologueDistribution",
                            deconvolution |>
                                dplyr::mutate(data = purrr::map2(data, deconv_resolved, ~ dplyr::inner_join(.x, .y, by = "Molecule"))) |>
                                dplyr::mutate(data = purrr::map(data, ~ .x |>
                                                                    dplyr::mutate(resolved_distribution = deconv_resolved / sum(deconv_resolved)))) |>
                                dplyr::select(-deconv_resolved, -Sample_Dilution_Factor) |>
                                tidyr::unnest(data) |>
                                dplyr::mutate(Deconvoluted_Distribution = as.numeric(resolved_distribution)) |>
                                dplyr::rename(Relative_Distribution = Relative_Area) |>
                                dplyr::mutate(Molecule_Concentration = Deconvoluted_Distribution * Concentration) |>
                                dplyr::select(
                                    Replicate_Name, Sample_Type, Molecule_List, Molecule, C_homologue, Cl_homologue, PCA,
                                    Relative_Distribution, Deconvoluted_Distribution, Molecule_Concentration, Unit
                                )
                        )

                        openxlsx::saveWorkbook(wb, file)
                    }
                )

                # ---- Quant table ----
                output$quantTable <- DT::renderDT({
                    deconvolution |>
                        dplyr::select(Replicate_Name, Sample_Type, Concentration, Unit, deconv_rsquared) |>
                        dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3)) |>
                        DT::datatable(
                            filter = "top", extensions = c("Buttons", "Scroller"),
                            options = list(
                                scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE,
                                buttons = list(
                                    list(extend = "excel", filename = "Samples_concentration", title = NULL,
                                         exportOptions = list(modifier = list(page = "all"))),
                                    list(extend = "csv", filename = "Samples_concentration", title = NULL,
                                         exportOptions = list(modifier = list(page = "all"))),
                                    list(extend = "colvis", targets = 0, visible = FALSE)
                                ),
                                dom = "lBfrtip",
                                fixedColumns = TRUE
                            ),
                            rownames = FALSE
                        )
                })

                # ---- Contribution plot ----
                output$sampleContributionPlot <- plotly::renderPlotly({
                    plot_sample_contribution(deconvolution)
                })

                # ---------------- QA/QC ----------------
                QAQC <- deconvolution |>
                    dplyr::select(Replicate_Name, Sample_Type, deconv_rsquared) |>
                    dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3))

                if (input$calculateRecovery == "Yes") {
                    recovery_data <- Skyline_output_filt |>
                        dplyr::filter(Isotope_Label_Type == "Quan",
                                      Molecule_List %in% c("IS", "RS"))

                    # Split IS and RS tables (keeping key context columns)
                    IS_tbl <- recovery_data %>%
                        dplyr::filter(Molecule_List == "IS") %>%
                        dplyr::select(Replicate_Name, Sample_Type, Molecule_IS = Molecule, Area_IS = Area)

                    RS_tbl <- recovery_data %>%
                        dplyr::filter(Molecule_List == "RS") %>%
                        dplyr::select(Replicate_Name, Sample_Type, Molecule_RS = Molecule, Area_RS = Area)

                    # Cross-join IS and RS within each replicate
                    ratio_all_pairs <- IS_tbl %>%
                        dplyr::inner_join(RS_tbl, by = c("Replicate_Name", "Sample_Type")) %>%
                        dplyr::mutate(
                            RatioISRS = dplyr::if_else(Area_RS > 0, Area_IS / Area_RS, NA_real_)
                        )

                    # Calculate QC ratio
                    qc_ratio <- ratio_all_pairs |>
                        dplyr::filter(Sample_Type == "Quality Control") |>
                        dplyr::group_by(Molecule_IS, Molecule_RS) |>
                        dplyr::summarize(AverageRatio = mean(RatioISRS, na.rm = TRUE), .groups = "drop")

                    ratio_with_sample <- ratio_all_pairs %>%
                        dplyr::filter(Sample_Type %in% c("Unknown", "Blank")) %>%
                        dplyr::mutate(
                            Molecule_IS = as.character(Molecule_IS),
                            Molecule_RS = as.character(Molecule_RS)
                        ) %>%
                        dplyr::left_join(
                            qc_ratio %>%
                                dplyr::mutate(
                                    Molecule_IS = as.character(Molecule_IS),
                                    Molecule_RS = as.character(Molecule_RS)
                                ),
                            by = c("Molecule_IS", "Molecule_RS")
                        ) %>%
                        dplyr::mutate(
                            Recovery = dplyr::if_else(!is.na(AverageRatio) & AverageRatio != 0,
                                                      round(RatioISRS / AverageRatio * 100, 1),
                                                      NA_real_)
                        ) %>%
                        dplyr::select(Replicate_Name, Sample_Type, Molecule_IS, Molecule_RS, Recovery)

                    QAQC <- QAQC |>
                        dplyr::left_join(ratio_with_sample, by = c("Replicate_Name", "Sample_Type"))
                }

                # Recovery table
                output$table_recovery <- DT::renderDT({
                    DT::datatable(
                        QAQC,
                        filter = "top",
                        extensions = c("Buttons", "Scroller"),
                        options = list(
                            scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE,
                            buttons = list(
                                list(extend = "excel", filename = "Samples_recovery", title = NULL,
                                     exportOptions = list(modifier = list(page = "all"))),
                                list(extend = "csv", filename = "Samples_recovery", title = NULL,
                                     exportOptions = list(modifier = list(page = "all"))),
                                list(extend = "colvis", targets = 0, visible = FALSE)
                            ),
                            dom = "lBfrtip",
                            fixedColumns = TRUE
                        ),
                        rownames = FALSE
                    )
                })

                # MDL
                if (input$calculateMDL == "Yes") {
                    if (input$blankSubtraction == "No") {
                        MDL_data <- deconvolution |>
                            dplyr::filter(Sample_Type == "Blank") |>
                            dplyr::summarize(
                                MDL_sumPCA = mean(Concentration) + 3 * stats::sd(Concentration, na.rm = TRUE),
                                number_of_blanks = dplyr::n_distinct(Replicate_Name)
                            )
                    } else if (input$blankSubtraction == "Yes, by avg area of blanks") {
                        MDL_data <- deconvolution |>
                            dplyr::filter(Sample_Type == "Blank") |>
                            dplyr::summarize(
                                MDL_sumPCA = 3 * stats::sd(Concentration, na.rm = TRUE),
                                number_of_blanks = dplyr::n_distinct(Replicate_Name)
                            )
                    }
                }

                # MDL table
                output$table_MDL <- DT::renderDT({
                    if (exists("MDL_data")) {
                        DT::datatable(MDL_data, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
                    } else {
                        NULL
                    }
                })


                # --- success popup when everything is done ---
                showModal(
                    modalDialog(
                        title = "Quantification successful",
                        "Quantification successful.",
                        easyClose = TRUE,
                        footer = modalButton("OK")
                    )
                )


            } # end if (input$standardTypes == "Group Mixtures")
        })  # end observeEvent(input$go)

        # ---- Homologue Group Patterns tab ----
        output$sampleSelectionUIOverlay <- renderUI({
            req(Samples_Concentration())
            sample_names <- unique(Samples_Concentration()$Replicate_Name)
            selectInput("selectedSamples", "Select samples to compare:",
                        choices = sample_names,
                        multiple = TRUE,
                        selected = NULL)
        })

        output$sampleSelectionUIComparisons <- renderUI({
            req(Samples_Concentration())
            sample_names <- unique(Samples_Concentration()$Replicate_Name)
            selectInput("selectedSamples", "Select samples to compare:",
                        choices = sample_names,
                        multiple = TRUE,
                        selected = NULL)
        })

        shiny::observeEvent(input$go2, {
            req(Samples_Concentration())

            withProgress(message = 'Generating plot...', value = 0, {
                Sample_distribution <- Samples_Concentration() |>
                    dplyr::mutate(data = purrr::map2(data, deconv_resolved, ~ dplyr::inner_join(.x, .y, by = "Molecule"))) |>
                    dplyr::mutate(data = purrr::map(data, ~ .x |>
                                                        dplyr::mutate(resolved_distribution = deconv_resolved / sum(deconv_resolved)))) |>
                    dplyr::select(-deconv_resolved, -Sample_Dilution_Factor) |>
                    tidyr::unnest(data) |>
                    dplyr::mutate(resolved_distribution = as.numeric(resolved_distribution)) |>
                    dplyr::mutate(deconv_resolved = as.numeric(deconv_resolved)) |>
                    dplyr::mutate(Molecule_concentration = resolved_distribution * Concentration)

                incProgress(0.5)

                if (input$plotHomologueGroups == "All Samples Overview") {
                    output$plotHomologuePatternStatic <- shiny::renderPlot({
                        ggplot2::ggplot(Sample_distribution, aes(x = Molecule, y = Relative_Area, fill = C_homologue)) +
                            ggplot2::geom_col() +
                            ggplot2::facet_wrap(~ Replicate_Name) +
                            ggplot2::theme_minimal() +
                            ggplot2::theme(axis.text.x = element_blank()) +
                            ggplot2::labs(
                                title = "Relative Distribution (before deconvolution)",
                                x = "Homologue",
                                y = "Relative Distribution"
                            )
                    })
                } else if (input$plotHomologueGroups == "Samples Overlay") {
                    output$plotHomologuePatternOverlay <- plotly::renderPlotly({
                        req(input$selectedSamples)
                        req(Sample_distribution)

                        selected_samples <- Sample_distribution |>
                            dplyr::filter(Replicate_Name %in% input$selectedSamples) |>
                            dplyr::mutate(Molecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)])))

                        req(nrow(selected_samples) > 0)

                        plotly::plot_ly(
                            data = selected_samples,
                            x = ~ Molecule,
                            y = ~ Relative_Area,   # switch to resolved_distribution if you prefer
                            color = ~ Replicate_Name,
                            type = "bar",
                            text = ~ paste(
                                "Sample:", Replicate_Name,
                                "<br>Homologue:", Molecule,
                                "<br>Distribution:", round(Relative_Area, 3),
                                "<br>C-atoms:", C_homologue
                            ),
                            hoverinfo = "text"
                        ) |>
                            layout(
                                xaxis = list(title = "Homologue", tickangle = 45),
                                yaxis = list(title = "Distribution"),
                                barmode = 'group',
                                showlegend = TRUE,
                                height = 600,
                                margin = list(b = 100)
                            )
                    })
                } else if (input$plotHomologueGroups == "Samples Panels") {
                    output$plotHomologuePatternComparisons <- plotly::renderPlotly({
                        req(input$selectedSamples)
                        plot_homologue_group_pattern_comparison(Sample_distribution, input$selectedSamples)
                    })
                }

                incProgress(1)
            })
        })

        # ---- Close the app when the session ends (non-interactive) ----
        if (!interactive() && !inherits(session, "MockShinySession")) {
            session$onSessionEnded(function() {
                stopApp()
                q("no")
            })
        }
    }

    # =========================================================
    # Run the application
    # =========================================================
    shinyApp(ui = ui, server = server)




}
