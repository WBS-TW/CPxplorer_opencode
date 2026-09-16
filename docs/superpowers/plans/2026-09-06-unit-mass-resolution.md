# Unit Mass Resolution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an unchecked **Unit mass resolution** checkbox on the CPions Interfering ions tab that hides MS Resolution, rounds m/z to integers, flags same-integer masses as interfering, and lets Skyline export those integer masses from the interference table.

**Architecture:** Add `apply_unit_mass_mz()` in `R/CPions_utils.R`. On Calculate with the checkbox on, round `m/z` then call existing `compute_interference()` with `.Machine$integer.max` so only exact-zero deltas are `YES`. Store the result in `CP_allions_compl2`. Skyline uses that table when it exists. No Skyline checkbox. No change to Normal/Advanced ion generation.

**Tech Stack:** R package CPxplorer, Shiny, dplyr, testthat edition 3.

## Global Constraints

- Checkbox label is exactly `Unit mass resolution`; input id is `unit_mass_resolution`; default `FALSE`.
- When checked, hide `MS Resolution`; do not add a Skyline checkbox.
- Round with `round(x, 0)` (R default, including half-to-even).
- Unit-mass interference: same integer m/z → `YES`; different integers → `NO`. Use `ms_resolution = .Machine$integer.max`, never `1L`.
- Do not add a `unit_mass` argument to `compute_interference()`.
- Do not change Normal/Advanced ion generation; those tables stay high-res.
- Do not add a separate nominal-mass column; `m/z` itself becomes integer in the interference table.
- If `CP_allions_compl2()` exists, Skyline `mz` export (Most intense and Interference-filtered) uses that table. If NULL, Most intense keeps Normal/Advanced tables; Interference-filtered still requires Calculate.
- Existing high-res interference snapshot must stay unchanged.
- Do not hand-edit `NAMESPACE` or `man/*.Rd`.
- Do not edit the embedded `isotopes` table in `R/CPions_utils.R`.
- Update `inst/instructions_CPions.md`. There is no `inst/instructions_CPions.html` in this repo; skip HTML.
- TDD: no production code without a failing test first.
- Do not commit unless the user explicitly asks; skip commit steps.
- Verification: `devtools::test(filter = "cpions-interference")`, `devtools::test(filter = "cpions-skyline")`, `devtools::test(filter = "cpions-app-smoke")`.

## File map

- `R/CPions_utils.R` — `apply_unit_mass_mz()`
- `R/CPions.R` — checkbox, hide MS Resolution, Calculate rounding, Skyline source table
- `inst/instructions_CPions.md` — Interfering ions docs
- `tests/testthat/test-cpions-interference.R` — rounding and unit-mass interference
- `tests/testthat/test-cpions-skyline.R` — integer Precursor m/z from unit-mass table
- `tests/testthat/test-cpions-app-smoke.R` — checkbox default plus unit-mass Calculate

---

### Task 1: Round m/z helper

**Files:**
- Modify: `R/CPions_utils.R` (insert `apply_unit_mass_mz()` immediately before `compute_interference()` around line 981)
- Test: `tests/testthat/test-cpions-interference.R`

**Interfaces:**
- Consumes: a data frame/tibble with an `m/z` column
- Produces: `apply_unit_mass_mz(CP_allions)` → same object with `m/z` replaced by `round(`m/z`, 0)`

- [ ] **Step 1: Write the failing rounding test**

Append to `tests/testthat/test-cpions-interference.R` (keep the existing snapshot test):

```r
test_that("apply_unit_mass_mz rounds m/z to nearest integer", {
    input <- tibble::tibble(`m/z` = c(376.8912, 377.4, 378.6))
    actual <- CPxplorer:::apply_unit_mass_mz(input)
    expect_identical(actual$`m/z`, c(377, 377, 379))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "cpions-interference")`

Expected: FAIL because `apply_unit_mass_mz` is not found.

- [ ] **Step 3: Write minimal implementation**

Insert immediately before `compute_interference()` in `R/CPions_utils.R`:

```r
apply_unit_mass_mz <- function(CP_allions) {
    dplyr::mutate(CP_allions, `m/z` = round(`m/z`, 0))
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `devtools::test(filter = "cpions-interference")`

Expected: PASS, including the existing high-res snapshot.

- [ ] **Step 5: Commit**

Skip unless the user asks.

---

### Task 2: Unit-mass interference flags

**Files:**
- Test: `tests/testthat/test-cpions-interference.R`
- Production: none if Task 1 helper plus existing `compute_interference()` already satisfy the assertions

**Interfaces:**
- Consumes: `apply_unit_mass_mz()`, `compute_interference(CP_allions, ms_resolution)`
- Produces: unit-mass table where shared integer m/z is `YES` and different integers are `NO`

- [ ] **Step 1: Write the failing interference test**

Append to `tests/testthat/test-cpions-interference.R`:

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "cpions-interference")`

Expected: FAIL on the interference assertion if current default-resolution logic were used, or PASS if `.Machine$integer.max` already yields this result. If it PASSES immediately, that is acceptable: the test locks the spec rule and is not testing a missing function. Do not change `compute_interference()`.

If it fails because first/last `case_when` branches disagree with the expected vector, fix the test data or expected vector to match the spec rule (same integer → YES, different integer → NO) without changing `compute_interference()`.

- [ ] **Step 3: Implementation**

No production change unless Step 2 fails for a reason other than assertion mismatch. Do not add a `unit_mass` argument.

- [ ] **Step 4: Run test to verify it passes**

Run: `devtools::test(filter = "cpions-interference")`

Expected: PASS. Existing snapshot unchanged.

- [ ] **Step 5: Commit**

Skip unless the user asks.

---

### Task 3: Skyline integer precursor m/z

**Files:**
- Test: `tests/testthat/test-cpions-skyline.R`
- Production: none in this task (export already copies `m/z` to `Precursor m/z`)

**Interfaces:**
- Consumes: `apply_unit_mass_mz()`, `compute_interference()`, `build_skyline_transition_list()`
- Produces: Skyline table whose `Precursor m/z` values are integers

- [ ] **Step 1: Write the failing Skyline test**

Append to `tests/testthat/test-cpions-skyline.R`:

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "cpions-skyline")`

Expected: PASS immediately if `Precursor m/z` already copies `m/z`. That locks the export contract. If it fails because `apply_unit_mass_mz` is missing, Task 1 is not done.

- [ ] **Step 3: Implementation**

No production change unless the test fails because Skyline re-computes high-res m/z. If that happens, stop and fix `build_skyline_transition_list()` so `Precursor m/z` comes from the input `m/z` column.

- [ ] **Step 4: Run test to verify it passes**

Run: `devtools::test(filter = "cpions-skyline")`

Expected: PASS. Existing Skyline tests unchanged.

- [ ] **Step 5: Commit**

Skip unless the user asks.

---

### Task 4: Interfering ions UI and Calculate path

**Files:**
- Modify: `R/CPions.R` (Interfering ions sidebar around lines 147-150; `observeEvent(input$go2)` around lines 475-488; Skyline `observeEvent(input$go3)` around lines 632-687)
- Test: `tests/testthat/test-cpions-app-smoke.R`

**Interfaces:**
- Consumes: `apply_unit_mass_mz()`, `compute_interference()`, `CP_allions_compl2`
- Produces: Shiny input `unit_mass_resolution`; Calculate writes integer-m/z interference table when checked; Skyline reads `CP_allions_compl2()` when non-NULL

- [ ] **Step 1: Write the failing app-smoke assertions**

In `tests/testthat/test-cpions-app-smoke.R`, change the existing go2 block to set the checkbox default, then add a unit-mass Calculate after the high-res check:

Replace:

```r
        session$setInputs(MSresolution = 20000, interfere_mode = "normal", go2 = 1)
        session$flushReact()
        expect_true(!is.null(output$Plotly))
        expect_true(!is.null(output$Plotly2))
        expect_true(!is.null(output$Table2))
```

with:

```r
        session$setInputs(
            unit_mass_resolution = FALSE,
            MSresolution = 20000,
            interfere_mode = "normal",
            go2 = 1
        )
        session$flushReact()
        expect_true(!is.null(output$Plotly))
        expect_true(!is.null(output$Plotly2))
        expect_true(!is.null(output$Table2))

        session$setInputs(unit_mass_resolution = TRUE, go2 = 2)
        session$flushReact()
        expect_true(!is.null(output$Table2))
```

Keep the existing Skyline block after this.

- [ ] **Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "cpions-app-smoke")`

Expected: FAIL or warn on unknown input `unit_mass_resolution` until the checkbox exists. If shiny testServer silently accepts unknown inputs, the test may still pass; continue to Step 3 because the UI is still required by the spec.

- [ ] **Step 3: Implement UI, Calculate, and Skyline source**

In `R/CPions.R` Interfering ions sidebar, replace the numeric input with:

```r
                            shiny::checkboxInput("unit_mass_resolution", "Unit mass resolution", value = FALSE),
                            shiny::conditionalPanel(
                                condition = "!input.unit_mass_resolution",
                                shiny::numericInput("MSresolution", "MS Resolution", value = 20000, min = 100, max = 5000000)
                            ),
```

In `observeEvent(input$go2)`, replace the `compute_interference(...)` call with:

```r
        if (isTRUE(input$unit_mass_resolution)) {
            CP_allions_interfere <- apply_unit_mass_mz(CP_allions_interfere)
            CP_allions_interfere <- compute_interference(CP_allions_interfere, .Machine$integer.max)
        } else {
            CP_allions_interfere <- compute_interference(CP_allions_interfere, MSresolution())
        }
```

In `observeEvent(input$go3)`, replace the two Most intense branches with one branch that prefers `CP_allions_compl2()`:

```r
    if (input$QuantIon == "Most intense") {
        source_ions <- if (!is.null(CP_allions_compl2())) {
            CP_allions_compl2()
        } else if (input$skyline_mode == "advanced") {
            CP_allions_glob_adv()
        } else {
            CP_allions_glob()
        }
        CP_allions_skyline <- build_skyline_transition_list(
            CP_allions = source_ions,
            mode = input$skyline_mode,
            quant_ion = "Most intense",
            ms_resolution = input$MSresolution,
            strategy = input$skyline_strategy,
            preferred_qual_n = as.integer(input$skyline_qual_n)
        )
    } else if (input$QuantIon == "Interference-filtered" & input$skyline_mode == "advanced") {
```

Leave the two Interference-filtered branches unchanged (they already use `CP_allions_compl2()`).

- [ ] **Step 4: Run tests to verify they pass**

Run: `devtools::test(filter = "cpions-app-smoke")`

Expected: PASS.

Also run: `devtools::test(filter = "cpions-interference")` and `devtools::test(filter = "cpions-skyline")`

Expected: PASS.

- [ ] **Step 5: Commit**

Skip unless the user asks.

---

### Task 5: Instructions

**Files:**
- Modify: `inst/instructions_CPions.md` (Introduction paragraph around line 13 and Interfering ions tab around lines 85-93)

**Interfaces:**
- Consumes: UI label `Unit mass resolution`
- Produces: user-facing docs for the checkbox, hidden MS Resolution, integer m/z, same-nominal interference, and Skyline integer export

- [ ] **Step 1: Update introduction**

Replace the Interfering ions sentence in the Instructions section:

```markdown
The _Interfering ions_ tab can be used to check for ions that interfer with each other at the estimated resolution of the mass spectrometer. Default is set to R=20,000. 
```

with:

```markdown
The _Interfering ions_ tab can be used to check for ions that interfer with each other at the estimated resolution of the mass spectrometer. Default is set to R=20,000. Check _Unit mass resolution_ (unchecked by default) to hide MS Resolution, round m/z values to integers for a low-resolution quadrupole, and treat ions with the same nominal m/z as interfering. After Calculate, Skyline export uses those integer masses.
```

- [ ] **Step 2: Update Interfering ions tab section**

After the `__interference__` paragraph, add:

```markdown
__Unit mass resolution__: unchecked by default. When checked, the MS Resolution input is hidden and all m/z values are rounded to integers to reflect a low-resolution (quadrupole) mass spectrometer. Ions that share the same nominal m/z are flagged as interfering. Skyline `mz` export uses these integer masses after Calculate.
```

- [ ] **Step 3: No automated test**

Instructions are rendered via `includeMarkdown`. No test change.

- [ ] **Step 4: Commit**

Skip unless the user asks.

---

## Spec coverage

- Checkbox + hide MS Resolution → Task 4
- Round m/z on Calculate → Tasks 1 and 4
- Same integer m/z → YES → Task 2
- Skyline uses interference table integers → Tasks 3 and 4
- Instructions → Task 5
- High-res snapshot unchanged → Tasks 1–2
- No Skyline checkbox, no CPquant, no isotopes table, no NAMESPACE/Rd edits → Global Constraints
