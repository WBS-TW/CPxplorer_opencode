test_that("CPions Skyline fixture matches expected output", {
    actual <- skyline_fixture()

    expect_identical(
        names(actual),
        c(
            "Molecule List Name",
            "Molecule Name",
            "Parent_Formula",
            "Adduct",
            "Precursor Charge",
            "Label Type",
            "Precursor m/z",
            "Explicit Retention Time",
            "Explicit Retention Time Window",
            "Ion Selection Rank",
            "Original Rel_ab Rank",
            "Selected Reason",
            "Interfering Candidate Count",
            "Closest Interference m/z",
            "Closest Interference Molecule",
            "Resolution Needed",
            "Interference from other filtered ions?",
            "Note"
        )
    )
    expect_true(all(actual$`Label Type` %in% c("Quan", "Qual")))
    expect_equal(sum(actual$`Label Type` == "Quan"), 2L)
    expect_equal(sum(actual$`Label Type` == "Qual"), 4L)
    expect_true(all(actual$`Interference from other filtered ions?` == "NO" | is.na(actual$`Interference from other filtered ions?`)))
    expect_false(any(grepl("[INTERFERENCE]", actual$Note, fixed = TRUE)))
})

test_that("Skyline Note warns when final exported transition interferes", {
    source_data <- dplyr::bind_rows(
        CPxplorer:::getAdduct_normal("[PCA+Cl]-", 10:10, 3:3, 3L, 5L),
        CPxplorer:::getAdduct_normal("[PCO+Cl]-", 10:10, 3:3, 3L, 5L)
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Most intense",
        ms_resolution = 20000L,
        strategy = "balanced",
        preferred_qual_n = 2L
    )

    final_interference <- CPxplorer:::compute_transition_interference(actual, 20000L)
    warned_rows <- grepl("[INTERFERENCE]", actual$Note, fixed = TRUE)

    expect_identical(actual$`Interference from other filtered ions?`, final_interference$interference)
    expect_identical(warned_rows, actual$`Interference from other filtered ions?` == "YES")
    expect_false(any(warned_rows))
})

test_that("Filtered Skyline selection avoids interfered Quan ions", {
    actual <- skyline_dense_fixture()

    quan_rows <- dplyr::filter(actual, `Label Type` == "Quan")

    expect_gt(nrow(quan_rows), 0)
    expect_false(any(grepl("[INTERFERENCE]", quan_rows$Note, fixed = TRUE)))
    expect_true(all(table(actual$`Molecule Name`, actual$`Label Type`)[, "Quan"] == 1))
})

test_that("Filtered Skyline selection yields a final non-conflicting transition set", {
    actual <- skyline_dense_fixture()

    final_interference <- actual |>
        CPxplorer:::compute_transition_interference(ms_resolution = 20000L)

    non_interfered <- dplyr::filter(final_interference, interference == "NO")

    expect_equal(nrow(non_interfered), nrow(actual))
    expect_true(all(actual$`Label Type` %in% c("Quan", "Qual")))
    expect_true(all(actual$`Interference from other filtered ions?` == "NO"))
})

test_that("Filtered Skyline selection limits Qual ions per molecule", {
    actual <- skyline_dense_fixture()

    qual_counts <- actual |>
        dplyr::filter(`Label Type` == "Qual") |>
        dplyr::count(`Molecule Name`, name = "n_qual")

    expect_true(all(qual_counts$n_qual <= 2L))
})

test_that("Abundance Skyline strategy selects highest abundance Quan ions", {
    source_data <- tibble::tibble(
        Molecule_Formula = c("C10H17Cl3", "C10H17Cl3", "C11H19Cl3", "C11H19Cl3"),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100, 101, 200, 201),
        Rel_ab = c(40, 80, 60, 20)
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "abundance",
        preferred_qual_n = 0L
    )

    quan_rows <- dplyr::filter(actual, `Label Type` == "Quan")

    expect_equal(nrow(quan_rows), 2L)
    expect_true(all(quan_rows$`Original Rel_ab Rank` == 1L))
    expect_setequal(quan_rows$Note, c("{[PCA+Cl]-}{80}", "{[PCA+Cl]-}{60}"))
})

test_that("Interference Skyline strategy prefers non-interfering Quan ions", {
    source_data <- tibble::tibble(
        Molecule_Formula = c("C10H17Cl3", "C10H17Cl3", "C11H19Cl3", "C11H19Cl3"),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100, 101, 200, 201),
        Rel_ab = c(90, 80, 70, 60),
        interference = c("YES", "NO", "YES", "NO")
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "interference",
        preferred_qual_n = 0L
    )

    quan_rows <- dplyr::filter(actual, `Label Type` == "Quan")

    expect_equal(nrow(quan_rows), 2L)
    expect_true(all(quan_rows$`Original Rel_ab Rank` == 2L))
    expect_true(all(quan_rows$`Selected Reason` == "Least interference: higher-abundance ion interfered"))
})

test_that("Interference Skyline strategy falls back to least-interfering Quan ions", {
    source_data <- tibble::tibble(
        Molecule_Formula = c("C10H17Cl3", "C10H17Cl3", "C11H19Cl3"),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100.000, 300.000, 100.001),
        Rel_ab = c(90, 80, 70),
        interference = "YES"
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "interference",
        preferred_qual_n = 0L
    )

    target_quan <- actual |>
        dplyr::filter(`Molecule Name` == "C10H17Cl3", `Label Type` == "Quan")

    expect_equal(target_quan$`Precursor m/z`, 300)
    expect_equal(target_quan$`Original Rel_ab Rank`, 2L)
    expect_equal(target_quan$`Selected Reason`, "Least interference: all ions interfere")
})

test_that("Skyline export from unit-mass interference table uses integer precursor m/z", {
    source_data <- dplyr::bind_rows(
        CPxplorer:::getAdduct_normal("[PCA+Cl]-", 10:10, 3:3, 3L, 5L),
        CPxplorer:::getAdduct_normal("[PCO+Cl]-", 10:10, 3:3, 3L, 5L)
    )

    unit_mass_ions <- source_data |>
        CPxplorer:::apply_unit_mass_mz() |>
        CPxplorer:::compute_interference(ms_resolution = .Machine$integer.max)

    actual <- CPxplorer:::build_skyline_transition_list(
        unit_mass_ions,
        mode = "normal",
        quant_ion = "Most intense",
        ms_resolution = 20000L,
        strategy = "balanced",
        preferred_qual_n = 2L
    )

    expect_true(nrow(actual) > 0)
    expect_true(all(actual$`Precursor m/z` == round(actual$`Precursor m/z`, 0)))
})

test_that("rank_skyline_qual_candidates orders leftovers by strategy", {
    leftovers <- tibble::tibble(
        `Precursor m/z` = c(100, 200, 300, 400),
        Rel_ab = c(40, 90, 70, 80),
        .original_interference = c("YES", "YES", "NO", "NO"),
        `Interfering Candidate Count` = c(3L, 1L, 0L, 2L),
        `Resolution Needed` = c(50000, 30000, NA_real_, 25000)
    )

    abundance <- CPxplorer:::rank_skyline_qual_candidates(leftovers, "abundance")
    expect_equal(abundance$Rel_ab, c(90, 80, 70, 40))

    balanced <- CPxplorer:::rank_skyline_qual_candidates(leftovers, "balanced")
    expect_equal(balanced$Rel_ab, c(80, 70, 90, 40))

    interference <- CPxplorer:::rank_skyline_qual_candidates(leftovers, "interference")
    expect_equal(interference$Rel_ab, c(70, 80, 90, 40))
})

test_that("Abundance Skyline strategy keeps highest-abundance Qual leftovers including interfering ions", {
    source_data <- tibble::tibble(
        Molecule_Formula = rep("C10H17Cl3", 3),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100, 200, 300),
        Rel_ab = c(90, 80, 70),
        interference = c("NO", "YES", "NO")
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "abundance",
        preferred_qual_n = 2L
    )

    quan <- dplyr::filter(actual, `Label Type` == "Quan")
    quals <- dplyr::filter(actual, `Label Type` == "Qual") |>
        dplyr::arrange(`Original Rel_ab Rank`)

    expect_equal(quan$`Precursor m/z`, 100)
    expect_equal(quals$`Precursor m/z`, c(200, 300))
    expect_equal(quals$`Selected Reason`, c("Qual: fallback interfering", "Qual: non-interfering"))
    expect_true(grepl("[INTERFERENCE]", quals$Note[[1]], fixed = TRUE))
    expect_false(grepl("[INTERFERENCE]", quals$Note[[2]], fixed = TRUE))
})

test_that("Balanced Skyline strategy fills remaining Qual slots with interfering leftovers", {
    source_data <- tibble::tibble(
        Molecule_Formula = rep("C10H17Cl3", 4),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100, 200, 300, 400),
        Rel_ab = c(90, 80, 70, 60),
        interference = c("YES", "NO", "YES", "NO")
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "balanced",
        preferred_qual_n = 2L
    )

    quan <- dplyr::filter(actual, `Label Type` == "Quan")
    quals <- dplyr::filter(actual, `Label Type` == "Qual") |>
        dplyr::arrange(
            match(`Selected Reason`, c("Qual: non-interfering", "Qual: fallback interfering")),
            `Original Rel_ab Rank`
        )

    expect_equal(quan$`Precursor m/z`, 200)
    expect_equal(quals$`Precursor m/z`, c(400, 100))
    expect_equal(quals$`Selected Reason`, c("Qual: non-interfering", "Qual: fallback interfering"))
    expect_true(grepl("[INTERFERENCE]", quals$Note[quals$`Precursor m/z` == 100], fixed = TRUE))
})

test_that("Interference Skyline strategy ranks Qual leftovers by interference then abundance", {
    source_data <- tibble::tibble(
        Molecule_Formula = rep("C10H17Cl3", 5),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100, 250, 300, 400, 500),
        Rel_ab = c(90, 50, 80, 40, 70),
        interference = c("YES", "NO", "YES", "NO", "YES")
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "interference",
        preferred_qual_n = 3L
    )

    quan <- dplyr::filter(actual, `Label Type` == "Quan")
    quals <- dplyr::filter(actual, `Label Type` == "Qual")
    qual_order <- quals |>
        dplyr::arrange(match(`Selected Reason`, c("Qual: non-interfering", "Qual: fallback interfering")), `Original Rel_ab Rank`)

    expect_equal(quan$`Precursor m/z`, 250)
    expect_equal(qual_order$`Precursor m/z`, c(400, 100, 300))
    expect_equal(
        qual_order$`Selected Reason`,
        c("Qual: non-interfering", "Qual: fallback interfering", "Qual: fallback interfering")
    )
})

test_that("Skyline Qual selection skips leftovers that collide with another molecule", {
    source_data <- tibble::tibble(
        Molecule_Formula = c("C10H17Cl3", "C10H17Cl3", "C11H19Cl3", "C11H19Cl3", "C11H19Cl3"),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100.000, 200.000, 300.000, 100.001, 400.000),
        Rel_ab = c(90, 10, 90, 80, 20),
        interference = "NO"
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "abundance",
        preferred_qual_n = 2L
    )

    c11 <- dplyr::filter(actual, `Molecule Name` == "C11H19Cl3")
    c11_quals <- dplyr::filter(c11, `Label Type` == "Qual")

    expect_false(any(abs(c11$`Precursor m/z` - 100.001) < 1e-8))
    expect_equal(c11_quals$`Precursor m/z`, 400)
})
