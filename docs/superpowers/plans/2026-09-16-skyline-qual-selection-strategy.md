# Skyline Qual Ion Selection by Strategy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make CPions Skyline interference-filtered **Selection strategy** apply to Qual ions as well as Quan, including fill-to-N with interfering leftovers and a skip for collisions with other molecules.

**Architecture:** Keep Quan selection unchanged. Add `rank_skyline_qual_candidates()` and `select_qual_by_strategy()` in `R/CPions_utils.R`. After all Quans are chosen, process each molecule’s leftovers in `unique(Molecule Name)` order, rank them with the same strategy, and keep up to `preferred_qual_n` Quals. Skip leftovers that collide with other molecules’ already-selected ions. Update Skyline UI copy and `inst/instructions_CPions.md`.

**Tech Stack:** R package CPxplorer, Shiny, dplyr, testthat edition 3.

## Global Constraints

- Strategy values stay `abundance`, `balanced`, `interference`; input id stays `skyline_strategy`.
- Do not change **Most intense** Quant Ion mode.
- Do not change Skyline export column names.
- Do not change interference calculation, Normal/Advanced ion generation, or CPquant.
- Do not pick Quan and Qual from one combined ranked list. All Quans first via `select_quan_by_strategy()`, then Quals from leftovers.
- Other-molecule collisions are never kept as fill-to-N Quals. Test `has_ms_interference()` against other molecules’ selected m/z only.
- Qual `Selected Reason`: `Qual: non-interfering` or `Qual: fallback interfering`.
- Qual `Ion Selection Rank` stays `Original Rel_ab Rank`.
- `preferred_qual_n = 0` still means no Quals.
- Do not hand-edit `NAMESPACE` or `man/*.Rd`.
- Do not edit the embedded `isotopes` table in `R/CPions_utils.R`.
- Update `inst/instructions_CPions.md`. There is no `inst/instructions_CPions.html` in this repo; skip HTML.
- TDD: no production code without a failing test first.
- Do not commit unless the user explicitly asks; skip commit steps.
- Verification: `devtools::test(filter = "cpions-skyline")`, and `devtools::test(filter = "cpions-app-smoke")` after UI copy changes.

## File map

- `R/CPions_utils.R` — `rank_skyline_qual_candidates()`, `select_qual_by_strategy()`, rewrite Qual loop in `select_skyline_filtered_candidates()`
- `R/CPions.R` — Skyline help text and strategy labels
- `inst/instructions_CPions.md` — Skyline strategy and Qual-count docs
- `tests/testthat/test-cpions-skyline.R` — ranking helper tests and Qual strategy tests

---

### Task 1: Rank leftover Qual candidates by strategy

**Files:**
- Modify: `R/CPions_utils.R` (insert `rank_skyline_qual_candidates()` immediately before `select_skyline_filtered_candidates()` around line 1215)
- Test: `tests/testthat/test-cpions-skyline.R`

**Interfaces:**
- Consumes: a tibble of leftover candidates with `Precursor m/z`, `Rel_ab`, `.original_interference`, `Interfering Candidate Count`, `Resolution Needed`
- Produces: `rank_skyline_qual_candidates(candidates, strategy)` → same rows reordered
  - `abundance`: descending `Rel_ab`, then `Precursor m/z`
  - `balanced`: `.original_interference` (`NO` before `YES`), then descending `Rel_ab`, then `Precursor m/z`
  - `interference`: `.original_interference`, then `Interfering Candidate Count`, then `dplyr::coalesce(Resolution Needed, 0)`, then descending `Rel_ab`, then `Precursor m/z`

- [ ] **Step 1: Write the failing ranking test**

Append to `tests/testthat/test-cpions-skyline.R`:

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "cpions-skyline")`

Expected: FAIL because `rank_skyline_qual_candidates` is not found.

- [ ] **Step 3: Write minimal implementation**

Insert immediately before `select_skyline_filtered_candidates()` in `R/CPions_utils.R`:

```r
rank_skyline_qual_candidates <- function(candidates, strategy) {
    switch(
        strategy,
        abundance = candidates |>
            dplyr::arrange(dplyr::desc(Rel_ab), `Precursor m/z`),
        balanced = candidates |>
            dplyr::arrange(
                .original_interference,
                dplyr::desc(Rel_ab),
                `Precursor m/z`
            ),
        interference = candidates |>
            dplyr::arrange(
                .original_interference,
                `Interfering Candidate Count`,
                dplyr::coalesce(`Resolution Needed`, 0),
                dplyr::desc(Rel_ab),
                `Precursor m/z`
            )
    )
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `devtools::test(filter = "cpions-skyline")`

Expected: PASS for the new ranking test. Existing Skyline tests still PASS.

---

### Task 2: Select Qual ions with the same strategy as Quan

**Files:**
- Modify: `R/CPions_utils.R` (`select_qual_by_strategy()` and Qual loop inside `select_skyline_filtered_candidates()`)
- Test: `tests/testthat/test-cpions-skyline.R`

**Interfaces:**
- Consumes: `rank_skyline_qual_candidates()`, `select_quan_by_strategy()`, `has_ms_interference()`, `prepare_skyline_candidates()`
- Produces: `select_qual_by_strategy(leftovers, strategy, preferred_qual_n, other_mz, own_mz, ms_resolution)` → zero or more Qual rows with `Label Type`, `Ion Selection Rank`, `Selected Reason`
- `select_skyline_filtered_candidates()` still returns the Quan+Qual table with Note flags

Qual acceptance:

1. All Quans are selected first. Then process molecules in `unique(candidates$Molecule Name)` order.
2. Skip a leftover if it collides with another molecule’s already-selected m/z.
3. `abundance`: keep remaining leftovers in ranked order up to N, including original/same-molecule interference.
4. `balanced` / `interference`: first keep leftovers that are not original-YES and do not collide with this molecule’s selected ions; if still short of N, fill with remaining leftovers that do not collide with other molecules.
5. `Selected Reason` is `Qual: non-interfering` or `Qual: fallback interfering`.
6. `Ion Selection Rank` is `Original Rel_ab Rank`.

- [ ] **Step 1: Write the four Qual selection tests**

Append to `tests/testthat/test-cpions-skyline.R`:

```r
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
        dplyr::arrange(dplyr::desc(Rel_ab))

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
            dplyr::desc(Rel_ab)
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
        dplyr::arrange(match(`Selected Reason`, c("Qual: non-interfering", "Qual: fallback interfering")), dplyr::desc(Rel_ab))

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
```

- [ ] **Step 2: Run tests to verify the first three fail**

Run: `devtools::test(filter = "cpions-skyline")`

Expected:

- Abundance Qual test FAIL: current Qual loop skips `interference == "YES"` leftovers, so m/z 200 is missing.
- Balanced fill-to-N FAIL: current loop does not fill remaining Qual slots with interfering leftovers, so m/z 100 is missing.
- Interference Qual order FAIL: current loop does not fill remaining slots, so only the non-interfering Qual is kept.
- Cross-molecule test may already PASS on current code (it already skips collisions with the growing selected set). Keep it as a regression guard.

- [ ] **Step 3: Implement Qual selection**

Insert `select_qual_by_strategy()` immediately after `rank_skyline_qual_candidates()` in `R/CPions_utils.R`:

```r
select_qual_by_strategy <- function(leftovers, strategy, preferred_qual_n, other_mz, own_mz, ms_resolution) {
    if (preferred_qual_n <= 0L || nrow(leftovers) == 0L) {
        return(leftovers[0, , drop = FALSE])
    }

    ranked <- rank_skyline_qual_candidates(leftovers, strategy)
    selected <- list()
    selected_own_mz <- own_mz

    accept_candidate <- function(candidate, reason) {
        candidate$`Label Type` <- "Qual"
        candidate$`Ion Selection Rank` <- candidate$`Original Rel_ab Rank`
        candidate$`Selected Reason` <- reason
        selected[[length(selected) + 1]] <<- candidate
        selected_own_mz <<- c(selected_own_mz, candidate$`Precursor m/z`[[1]])
    }

    candidate_is_interfering <- function(candidate) {
        mz <- candidate$`Precursor m/z`[[1]]
        original_yes <- identical(candidate$.original_interference[[1]], "YES")
        collides_own <- has_ms_interference(mz, selected_own_mz, ms_resolution)
        original_yes || collides_own
    }

    for (i in seq_len(nrow(ranked))) {
        if (length(selected) >= preferred_qual_n) {
            break
        }
        candidate <- ranked[i, , drop = FALSE]
        mz <- candidate$`Precursor m/z`[[1]]
        if (has_ms_interference(mz, other_mz, ms_resolution)) {
            next
        }
        is_interfering <- candidate_is_interfering(candidate)
        if (strategy == "abundance") {
            accept_candidate(
                candidate,
                if (is_interfering) "Qual: fallback interfering" else "Qual: non-interfering"
            )
            next
        }
        if (!is_interfering) {
            accept_candidate(candidate, "Qual: non-interfering")
        }
    }

    if (strategy != "abundance" && length(selected) < preferred_qual_n) {
        kept_ids <- vapply(selected, function(row) row$.candidate_id[[1]], numeric(1))
        for (i in seq_len(nrow(ranked))) {
            if (length(selected) >= preferred_qual_n) {
                break
            }
            candidate <- ranked[i, , drop = FALSE]
            if (candidate$.candidate_id[[1]] %in% kept_ids) {
                next
            }
            mz <- candidate$`Precursor m/z`[[1]]
            if (has_ms_interference(mz, other_mz, ms_resolution)) {
                next
            }
            accept_candidate(candidate, "Qual: fallback interfering")
        }
    }

    if (length(selected) == 0L) {
        return(ranked[0, , drop = FALSE])
    }
    dplyr::bind_rows(selected)
}
```

Replace the Qual loop in `select_skyline_filtered_candidates()` (from `selected_rows <- ...` through the `for` loop) with:

```r
    selected_rows <- if (nrow(quan_candidates) > 0) list(quan_candidates) else list()
    selected_ids <- quan_candidates$.candidate_id
    molecule_names <- unique(candidates$`Molecule Name`)
    selected_mz_by_molecule <- split(quan_candidates$`Precursor m/z`, quan_candidates$`Molecule Name`)
    for (molecule_name in molecule_names) {
        if (is.null(selected_mz_by_molecule[[molecule_name]])) {
            selected_mz_by_molecule[[molecule_name]] <- numeric(0)
        }
    }

    for (molecule_name in molecule_names) {
        leftovers <- candidates |>
            dplyr::filter(`Molecule Name` == molecule_name, !(.candidate_id %in% selected_ids))
        other_names <- setdiff(names(selected_mz_by_molecule), molecule_name)
        other_mz <- unlist(selected_mz_by_molecule[other_names], use.names = FALSE)
        own_mz <- selected_mz_by_molecule[[molecule_name]]
        quals <- select_qual_by_strategy(
            leftovers = leftovers,
            strategy = strategy,
            preferred_qual_n = preferred_qual_n,
            other_mz = other_mz,
            own_mz = own_mz,
            ms_resolution = ms_resolution
        )
        if (nrow(quals) > 0) {
            selected_rows[[length(selected_rows) + 1]] <- quals
            selected_ids <- c(selected_ids, quals$.candidate_id)
            selected_mz_by_molecule[[molecule_name]] <- c(own_mz, quals$`Precursor m/z`)
        }
    }
```

Keep the existing `dplyr::bind_rows(selected_rows) |> compute_transition_interference(...)` Note-flag tail unchanged.

Do not walk leftovers globally by `Original Rel_ab Rank`.

- [ ] **Step 4: Run Skyline tests**

Run: `devtools::test(filter = "cpions-skyline")`

Expected: the four new Qual tests PASS. Existing `preferred_qual_n = 0L` Quan tests still PASS.

If `Filtered Skyline selection yields a final non-conflicting transition set` or `CPions Skyline fixture matches expected output` fails because balanced fill-to-N now keeps an interfering Qual in the dense/small fixtures, update only those assertions to allow flagged interfering Quals. Do not weaken them if Qual picks did not actually change.

- [ ] **Step 5: Review any snapshot diffs**

If `_snaps/cpions-skyline.md` changes, confirm the change is from intentional Qual picks. Do not accept snapshot updates that only silence a failure.

---

### Task 3: Skyline UI copy and instructions

**Files:**
- Modify: `R/CPions.R` (Skyline sidebar around lines 174–186)
- Modify: `inst/instructions_CPions.md` (Skyline strategy section around lines 120–128)
- Test: `tests/testthat/test-cpions-app-smoke.R` (run existing smoke; no new assertion required)

**Interfaces:**
- Consumes: existing `skyline_strategy` input id and values `abundance` / `balanced` / `interference`
- Produces: labels without Quan-only wording; help text and instructions that cover Qual

- [ ] **Step 1: Update Skyline sidebar copy**

In `R/CPions.R`, replace the interference-filtered help text and strategy choices with:

```r
                                 shiny::helpText("When using interference filtering, choose a strategy to select Quan and Qual ions and set the preferred number of Qual ions."),
                                 shiny::selectInput(
                                     "skyline_strategy",
                                     label = "Selection strategy",
                                     choices = c(
                                        "Highest abundance" = "abundance",
                                        "Balanced" = "balanced",
                                        "Least interference" = "interference"
                                     ),
                                     selected = "balanced"
                                 ),
```

Keep `skyline_qual_n` unchanged.

- [ ] **Step 2: Update instructions**

In `inst/instructions_CPions.md`, replace the Selection strategy and Preferred number of Qual ions paragraphs with:

```markdown
__Selection strategy__ (only applies when _"Interference-filtered"_ is selected):

  _Highest abundance_: selects the most abundant ion as Quan regardless of interference. Qual ions are the next-most-abundant leftovers, including interfering ions. Interfering Quan or Qual ions are flagged [INTERFERENCE] in the Note column.

  _Balanced_: selects the highest-abundance ion that does not have mass interference at the set MS resolution as Quan. Qual ions prefer the highest-abundance non-interfering leftovers. If all ions for a _Molecule Name_ interfere, the most abundant ion is used as Quan and flagged [INTERFERENCE]. If fewer than the requested Qual ions are non-interfering, remaining Qual slots are filled with the next-best interfering leftovers and flagged [INTERFERENCE]. This balances signal strength with interference avoidance.

  _Least interference_: prioritizes minimizing interference over abundance for both Quan and Qual. Ions are sorted first by interference status and count, then by abundance. This may select a lower-abundance Quan or Qual ion to avoid interference. If fewer than the requested Qual ions are non-interfering, remaining Qual slots are filled with the next-best interfering leftovers and flagged [INTERFERENCE].

__Preferred number of Qual ions__ (only applies when _"Interference-filtered"_ is selected): sets the target number of Qual ions kept per _Molecule Name_ in the strategy-based selection. The actual count may be lower if too few leftovers remain after skipping collisions with other molecules' selected ions.
```

Do not create `inst/instructions_CPions.html`.

- [ ] **Step 3: Run Skyline and smoke tests**

Run: `devtools::test(filter = "cpions-skyline")`

Expected: PASS.

Run: `devtools::test(filter = "cpions-app-smoke")`

Expected: PASS. Existing smoke still sets `skyline_strategy = "Balanced"`; `tolower()` in `build_skyline_transition_list()` maps that to `balanced`.
