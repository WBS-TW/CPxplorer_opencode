test_that("CPions advanced settings fixture matches expected output", {
    actual <- advanced_fixture()

    expect_gt(nrow(actual), 0)
    expect_true(all(c(
        "Molecule_Formula", "Halo_perc", "Compound_Class", "TP",
        "Charge", "Adduct", "Adduct_Isotopologue",
        "Adduct_Formula", "Isotope_Formula", "m/z", "Rel_ab"
    ) %in% names(actual)))

    expect_snapshot_value(actual, style = "json2")
})

expect_delta <- function(notation, ...) {
    expected <- c(C = 0L, H = 0L, Cl = 0L, Br = 0L, O = 0L, S = 0L, F = 0L)
    updates <- c(...)
    expected[names(updates)] <- as.integer(updates)
    expect_identical(CPxplorer:::parse_tp_notation(notation), expected)
}

test_that("parse_tp_notation maps catalog strings to element deltas", {
    expect_delta("None")
    expect_delta("-H+OH", O = 1L)
    expect_delta("-Cl+OH", H = 1L, Cl = -1L, O = 1L)
    expect_delta("-2Cl+2OH", H = 2L, Cl = -2L, O = 2L)
    expect_delta("-2H+2OH", O = 2L)
    expect_delta("-2H+O", H = -2L, O = 1L)
    expect_delta("-H+SO4H", O = 4L, S = 1L)
    expect_delta("-H+C6H10O7", C = 6L, H = 9L, O = 7L)
    expect_delta("-2H+2O", H = -2L, O = 2L)
    expect_delta("-Br+OH", H = 1L, Br = -1L, O = 1L)
    expect_delta("-2Br+2OH", H = 2L, Br = -2L, O = 2L)
    expect_delta("-H+OCH3", C = 1L, H = 2L, O = 1L)
    expect_delta("-Cl+OCH3", C = 1L, H = 3L, Cl = -1L, O = 1L)
    expect_delta("-4H+2O", H = -4L, O = 2L)
})

test_that("parse_tp_list splits custom TP text", {
    expect_identical(CPxplorer:::parse_tp_list(""), "None")
    expect_identical(CPxplorer:::parse_tp_list("None"), "None")
    expect_identical(
        CPxplorer:::parse_tp_list(" -H+OH ; -Cl+OH "),
        c("-H+OH", "-Cl+OH")
    )
})

test_that("parse_tp_notation rejects bad grammar", {
    expect_error(CPxplorer:::parse_tp_notation("-H+COOH"), "COOH")
    expect_error(CPxplorer:::parse_tp_notation("H+OH"), "sign")
    expect_error(CPxplorer:::parse_tp_notation(""), "empty")
})

test_that("tp_catalog lists predefined notations with None first", {
    notations <- CPxplorer:::tp_catalog_notations()
    expect_identical(notations[1], "None")
    expect_true(all(c(
        "-Cl+OH", "-H+OH", "-2Cl+2OH", "-2H+2OH", "-2H+O", "-H+SO4H",
        "-H+C6H10O7", "-2H+2O", "-Br+OH", "-2Br+2OH", "-H+OCH3",
        "-Cl+OCH3", "-4H+2O"
    ) %in% notations))
    for (notation in setdiff(notations, "None")) {
        expect_silent(CPxplorer:::parse_tp_notation(notation))
    }
})

test_that("parent_h uses class-specific hydrogen counts", {
    expect_identical(CPxplorer:::parent_h("PCA", 10L, 6L, 0L), 16L)
    expect_identical(CPxplorer:::parent_h("PCO", 10L, 6L, 0L), 14L)
    expect_identical(CPxplorer:::parent_h("PCdiO", 10L, 6L, 0L), 12L)
    expect_identical(CPxplorer:::parent_h("PCtriO", 10L, 6L, 0L), 10L)
    expect_identical(CPxplorer:::parent_h("BCA", 10L, 4L, 2L), 16L)
})

test_that("parent_class_row errors on unknown class", {
    expect_error(CPxplorer:::parent_h("UNKNOWN", 10L, 6L, 0L), "unknown compound class")
})
