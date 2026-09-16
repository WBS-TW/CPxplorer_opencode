test_that("CPions interfering ions fixture matches expected output", {
    actual <- interference_fixture()

    expect_true(all(c(
        "difflag", "difflead", "reslag", "reslead", "interference"
    ) %in% names(actual)))
    expect_type(actual$interference, "character")

    expect_snapshot_value(
        dplyr::select(
            actual,
            Molecule_Formula,
            Adduct,
            `m/z`,
            Rel_ab,
            difflag,
            difflead,
            reslag,
            reslead,
            interference
        ),
        style = "json2"
    )
})

test_that("apply_unit_mass_mz rounds m/z to nearest integer", {
    input <- tibble::tibble(`m/z` = c(376.8912, 377.4, 378.6))
    actual <- CPxplorer:::apply_unit_mass_mz(input)
    expect_identical(actual$`m/z`, c(377, 377, 379))
})

test_that("unit mass resolution flags only shared integer m/z as interfering", {
    ions <- tibble::tibble(
        Molecule_Formula = c("A", "B", "C"),
        `m/z` = c(376.4, 375.6, 378.1)
    )

    actual <- ions |>
        CPxplorer:::apply_unit_mass_mz() |>
        CPxplorer:::compute_interference(ms_resolution = .Machine$integer.max)

    expect_identical(actual$`m/z`, c(376, 376, 378))
    expect_identical(actual$interference, c("YES", "YES", "NO"))
})
