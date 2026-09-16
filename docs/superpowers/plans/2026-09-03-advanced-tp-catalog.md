# Advanced TP Catalog and Custom Input Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace hardcoded Advanced transformation-product branches with a notation catalog plus parser, add custom TP text input that blocks on invalid/infeasible chemistry, add extra TPs, and add Advanced-only PCdiO/PCtriO parent classes.

**Architecture:** Keep logic in `R/CPions_utils.R`. A `tp_catalog` (notation + name) and `parent_class_catalog` (H formula + has_Br) feed one parser (`parse_tp_notation()`) and one applier used by both the dropdown and custom text. `getAdduct_advanced()` builds parent atoms, applies TP deltas, then existing adduct/isotope code. Excel is regenerated from the catalog. Advanced UI in `R/CPions.R` gains PCdiO/PCtriO and a custom-TP checkbox.

**Tech Stack:** R package CPxplorer, Shiny, dplyr/tidyr/stringr, enviPat, testthat edition 3, openxlsx/readxl.

## Global Constraints

- Do not change Normal settings adducts or add `[PCdiO+Cl]-` / `[PCtriO+Cl]-` there.
- Excel is not the runtime TP source; `tp_catalog` is.
- Invalid grammar or infeasible TPs block the whole Advanced run; no partial table.
- If some homologues in a C/Cl range are valid and some are not, drop impossible rows and keep the rest. Block only when a selected class+TP leaves zero valid rows.
- Mixed PCA+BCA with `-Br+OH` blocks (any selected class without Br plus a Br-loss TP).
- Skyline `Molecule List Name` uses **parent** carbon count, not post-TP `Molecule_Formula`.
- Existing advanced snapshot (`TP = "None"`, PCA C10 Cl3) must stay unchanged.
- Do not hand-edit `NAMESPACE` or `man/*.Rd`.
- Do not edit the embedded `isotopes` table in `R/CPions_utils.R`.
- Update `inst/instructions_CPions.md` when Advanced UI changes. There is no `inst/instructions_CPions.html` in this repo; skip HTML.
- TDD: no production code without a failing test first.
- Do not commit unless the user explicitly asks; skip commit steps until then.
- Verification: `devtools::test(filter = "cpions-advanced")` for generation/parser; `devtools::test(filter = "cpions-app-smoke")` if Advanced inputs change.

## File map

- `R/CPions_utils.R` — catalogs, parser, parent atoms, apply-TP, feasibility, refactor `getAdduct_advanced()` / `generateInput_Envipat_advanced()`, Skyline class naming, Excel rebuild helper
- `R/CPions.R` — Compound Class choices, custom TP checkbox/text, submit validation
- `inst/CPions_TP_formula.xlsx` — regenerated from catalog
- `inst/instructions_CPions.md` — Advanced docs
- `tests/testthat/test-cpions-advanced.R` — parser, apply-TP, feasibility, extra TPs, PCdiO/PCtriO, Excel examples
- `tests/testthat/test-cpions-app-smoke.R` — checkbox default + new input ids

---

### Task 1: TP notation parser

**Files:**
- Modify: `R/CPions_utils.R` (insert new helpers after `create_formula()` around line 38)
- Test: `tests/testthat/test-cpions-advanced.R`

**Interfaces:**
- Consumes: none
- Produces:
  - `empty_tp_delta()` → named integer vector `c(C=0L, H=0L, Cl=0L, Br=0L, O=0L, S=0L, F=0L)`
  - `parse_tp_notation(notation)` → that named integer vector; errors with a message containing the notation
  - `parse_tp_list(text)` → character vector of trimmed notations; empty/`None` → `"None"`

- [ ] **Step 1: Write the failing parser tests**

Append to `tests/testthat/test-cpions-advanced.R` (keep the existing snapshot test):

```r
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: FAIL with `parse_tp_notation` not found (or equivalent).

- [ ] **Step 3: Implement parser**

Insert in `R/CPions_utils.R` after `create_formula()`:

```r
empty_tp_delta <- function() {
    c(C = 0L, H = 0L, Cl = 0L, Br = 0L, O = 0L, S = 0L, F = 0L)
}

.tp_group_atoms <- list(
    H = c(C = 0L, H = 1L, Cl = 0L, Br = 0L, O = 0L, S = 0L, F = 0L),
    Cl = c(C = 0L, H = 0L, Cl = 1L, Br = 0L, O = 0L, S = 0L, F = 0L),
    Br = c(C = 0L, H = 0L, Cl = 0L, Br = 1L, O = 0L, S = 0L, F = 0L),
    OH = c(C = 0L, H = 1L, Cl = 0L, Br = 0L, O = 1L, S = 0L, F = 0L),
    O = c(C = 0L, H = 0L, Cl = 0L, Br = 0L, O = 1L, S = 0L, F = 0L),
    SO4H = c(C = 0L, H = 1L, Cl = 0L, Br = 0L, O = 4L, S = 1L, F = 0L),
    OCH3 = c(C = 1L, H = 3L, Cl = 0L, Br = 0L, O = 1L, S = 0L, F = 0L),
    C6H10O7 = c(C = 6L, H = 10L, Cl = 0L, Br = 0L, O = 7L, S = 0L, F = 0L)
)

parse_tp_notation <- function(notation) {
    notation <- stringr::str_trim(as.character(notation))
    if (length(notation) != 1L || is.na(notation) || notation == "") {
        stop("empty TP notation", call. = FALSE)
    }
    if (identical(notation, "None")) {
        return(empty_tp_delta())
    }
    token_re <- "([+-])(\\d+)?(C6H10O7|SO4H|OCH3|OH|Cl|Br|O|H)"
    matches <- gregexpr(token_re, notation, perl = TRUE)[[1]]
    if (matches[1] == -1L) {
        stop(sprintf("invalid TP notation '%s': missing sign", notation), call. = FALSE)
    }
    consumed <- sum(attr(matches, "match.length"))
    if (consumed != nchar(notation)) {
        stop(sprintf("invalid TP notation '%s'", notation), call. = FALSE)
    }
    delta <- empty_tp_delta()
    starts <- as.integer(matches)
    lens <- attr(matches, "match.length")
    for (i in seq_along(starts)) {
        token <- substr(notation, starts[i], starts[i] + lens[i] - 1L)
        sign <- if (starts_with_minus <- substr(token, 1, 1) == "-") -1L else 1L
        rest <- substr(token, 2L, nchar(token))
        coef <- 1L
        coef_match <- regexpr("^\\d+", rest, perl = TRUE)
        if (coef_match == 1L) {
            coef <- as.integer(regmatches(rest, coef_match))
            rest <- substr(rest, attr(coef_match, "match.length") + 1L, nchar(rest))
        }
        group <- .tp_group_atoms[[rest]]
        if (is.null(group)) {
            stop(sprintf("unknown TP group '%s' in '%s'", rest, notation), call. = FALSE)
        }
        delta <- delta + sign * coef * group
    }
    delta
}

parse_tp_list <- function(text) {
    text <- if (length(text) == 0L || is.na(text)) "" else as.character(text)
    parts <- stringr::str_trim(unlist(strsplit(text, ";", fixed = TRUE)))
    parts <- parts[parts != ""]
    if (length(parts) == 0L) {
        return("None")
    }
    parts
}
```

Fix the sign-error path so `H+OH` (no leading `+/-`) hits the missing-sign message: if `substr(notation, 1, 1)` is not `+` or `-` and notation is not `None`, `stop(..., "missing sign")`. Keep the unknown-group path so `-H+COOH` mentions `COOH`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: new parser tests PASS; existing snapshot still PASS.

- [ ] **Step 5: Commit (only if the user asked)**

```bash
git add R/CPions_utils.R tests/testthat/test-cpions-advanced.R
git commit -m "feat: parse transformation-product notation into element deltas"
```

---

### Task 2: Catalogs

**Files:**
- Modify: `R/CPions_utils.R` (next to parser)
- Test: `tests/testthat/test-cpions-advanced.R`

**Interfaces:**
- Consumes: `parse_tp_notation()`
- Produces:
  - `tp_catalog` tibble columns `notation`, `name`
  - `parent_class_catalog` tibble columns `class`, `h_offset`, `has_br`, `advanced_only`
  - `tp_catalog_notations()` → character vector, `None` first
  - `parent_h(class, C, Cl, Br)` → integer H count

Parent H formulas (spec):

| class  | H                         | has_br | advanced_only |
|--------|---------------------------|--------|---------------|
| PCA    | `2*C + 2 - Cl`            | FALSE  | FALSE         |
| PCO    | `2*C - Cl`                | FALSE  | FALSE         |
| PCdiO  | `2*C - 2 - Cl`            | FALSE  | TRUE          |
| PCtriO | `2*C - 4 - Cl`            | FALSE  | TRUE          |
| BCA    | `2*C + 2 - Cl - Br`       | TRUE   | FALSE         |

Store `h_offset` as the saturated-alkane offset before removing Cl/Br: PCA/BCA `+2`, PCO `0`, PCdiO `-2`, PCtriO `-4`. Then `H = 2*C + h_offset - Cl - ifelse(has_br, Br, 0)`.

- [ ] **Step 1: Write failing catalog tests**

```r
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: FAIL, `tp_catalog_notations` not found.

- [ ] **Step 3: Implement catalogs**

```r
tp_catalog <- tibble::tribble(
    ~notation,       ~name,
    "None",          "None",
    "-Cl+OH",        "Cl-hydroxylation",
    "-H+OH",         "H-hydroxylation",
    "-2Cl+2OH",      "Double Cl-hydroxylation",
    "-2H+2OH",       "Double H-hydroxylation",
    "-2H+O",         "Oxidation (ketone/aldehyde)",
    "-H+SO4H",       "Sulfonation",
    "-H+C6H10O7",    "Glucuronidation",
    "-2H+2O",        "Carboxylic acid / omega-oxidation",
    "-Br+OH",        "Br-hydroxylation",
    "-2Br+2OH",      "Double Br-hydroxylation",
    "-H+OCH3",       "Methoxylation",
    "-Cl+OCH3",      "Cl-methoxylation",
    "-4H+2O",        "Diketone"
)

tp_catalog_notations <- function() {
    tp_catalog$notation
}

parent_class_catalog <- tibble::tribble(
    ~class,    ~h_offset, ~has_br, ~advanced_only,
    "PCA",      2L,        FALSE,   FALSE,
    "PCO",      0L,        FALSE,   FALSE,
    "PCdiO",   -2L,        FALSE,   TRUE,
    "PCtriO",  -4L,        FALSE,   TRUE,
    "BCA",      2L,        TRUE,    FALSE
)
```

Implement:

```r
parent_class_row <- function(class) {
    row <- parent_class_catalog[parent_class_catalog$class == class, , drop = FALSE]
    if (nrow(row) != 1L) {
        stop(sprintf("unknown compound class '%s'", class), call. = FALSE)
    }
    row
}

parent_h <- function(class, C, Cl, Br = 0L) {
    row <- parent_class_row(class)
    as.integer(2L * C + row$h_offset - Cl - if (row$has_br) Br else 0L)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: PASS.

- [ ] **Step 5: Commit (only if the user asked)**

```bash
git add R/CPions_utils.R tests/testthat/test-cpions-advanced.R
git commit -m "feat: add TP and parent-class catalogs"
```

---

### Task 3: Apply TP deltas and feasibility

**Files:**
- Modify: `R/CPions_utils.R`
- Test: `tests/testthat/test-cpions-advanced.R`

**Interfaces:**
- Consumes: `parse_tp_notation()`, `parent_h()`, `parent_class_row()`, `create_formula()`
- Produces:
  - `build_parent_grid(class, C, Cl, Clmax, Br, Brmax)` → tibble with `C,H,Cl,Br,S,O,F,Parent_Formula` (S/O/F = 0; Br = 0 unless `has_br`)
  - `apply_tp_to_parents(parents, notation)` → same columns plus `Molecule_Formula`, `TP`; rows with any negative atom count dropped
  - `tp_feasibility_error(classes, notations, C, Cl, Clmax, Br, Brmax)` → `NULL` if OK, otherwise a single error string naming the TP and reason

Feasibility rules (spec):

1. Parse error → message includes notation.
2. Br-loss (`delta["Br"] < 0`) if **any** selected class has `has_br == FALSE` → class mismatch, block even if BCA is also selected.
3. After apply, **zero** remaining rows for a class+notation over the selected range → infeasible, block.
4. Partial range (some Cl values work) → not an error; dropped rows only.

- [ ] **Step 1: Write failing feasibility tests**

```r
test_that("apply_tp_to_parents yields Excel example formulas for PCA C10H16Cl6", {
    parents <- CPxplorer:::build_parent_grid("PCA", 10L, 6L, 6L, 0L, 0L)
    expect_identical(parents$Parent_Formula, "C10H16Cl6")
    mol <- function(tp) {
        CPxplorer:::apply_tp_to_parents(parents, tp)$Molecule_Formula
    }
    expect_identical(mol("None"), "C10H16Cl6")
    expect_identical(mol("-Cl+OH"), "C10H17Cl5O")
    expect_identical(mol("-H+OH"), "C10H16Cl6O")
    expect_identical(mol("-2Cl+2OH"), "C10H18Cl4O2")
    expect_identical(mol("-2H+2OH"), "C10H16Cl6O2")
    expect_identical(mol("-2H+O"), "C10H14Cl6O")
    expect_identical(mol("-H+SO4H"), "C10H16Cl6O4S")
    expect_identical(mol("-H+C6H10O7"), "C16H25Cl6O7")
    expect_identical(mol("-2H+2O"), "C10H14Cl6O2")
    expect_identical(mol("-H+OCH3"), "C11H18Cl6O")
    expect_identical(mol("-Cl+OCH3"), "C11H19Cl5O")
    expect_identical(mol("-4H+2O"), "C10H12Cl6O2")
})

test_that("Br TPs apply on BCA and are rejected on PCA", {
    bca <- CPxplorer:::build_parent_grid("BCA", 10L, 4L, 4L, 2L, 2L)
    expect_identical(bca$Parent_Formula, "C10H16Cl4Br2")
    expect_identical(
        CPxplorer:::apply_tp_to_parents(bca, "-Br+OH")$Molecule_Formula,
        "C10H17Cl4BrO"
    )
    expect_match(
        CPxplorer:::tp_feasibility_error("PCA", "-Br+OH", 10L, 6L, 6L, 0L, 0L),
        "Br"
    )
    expect_match(
        CPxplorer:::tp_feasibility_error(c("PCA", "BCA"), "-Br+OH", 10L, 4L, 4L, 2L, 2L),
        "Br"
    )
    expect_null(
        CPxplorer:::tp_feasibility_error("BCA", "-Br+OH", 10L, 4L, 4L, 2L, 2L)
    )
})

test_that("infeasible Cl range blocks and partial range keeps valid rows", {
    expect_match(
        CPxplorer:::tp_feasibility_error("PCA", "-2Cl+2OH", 10L, 1L, 1L, 0L, 0L),
        "-2Cl\\+2OH"
    )
    expect_null(
        CPxplorer:::tp_feasibility_error("PCA", "-2Cl+2OH", 10L, 1L, 4L, 0L, 0L)
    )
    parents <- CPxplorer:::build_parent_grid("PCA", 10L, 1:4, 4L, 0L, 0L)
    out <- CPxplorer:::apply_tp_to_parents(parents, "-2Cl+2OH")
    expect_true(all(out$Cl >= 0L))
    expect_equal(nrow(out), 3L)
})
```

`build_parent_grid()` for PCA should use the same filters as today: `C >= Cl`, `Cl <= Clmax`. For BCA also `Br <= Brmax` and `Br + Cl <= C`. When `Cl` is an integer sequence, pass it as the vector `Cl` argument (same as `getAdduct_advanced`).

Signature to implement:

```r
build_parent_grid <- function(class, C, Cl, Clmax, Br, Brmax)
apply_tp_to_parents <- function(parents, notation)
tp_feasibility_error <- function(classes, notations, C, Cl, Clmax, Br, Brmax)
```

`Cl` / `C` / `Br` are integer vectors as in `getAdduct_advanced`. For the Cl=1..4 test, call `build_parent_grid("PCA", 10L, 1:4, 4L, 0L, 0L)`.

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: FAIL, functions not found.

- [ ] **Step 3: Implement apply/feasibility**

`build_parent_grid`: look up class; `tidyr::crossing` C/Cl (and Br if `has_br`); apply current filters; `H = parent_h(...)`; `S=O=F=0`; `Br=0` if not `has_br`; `Parent_Formula = create_formula(C,H,Cl,Br,S,O,F)`.

`apply_tp_to_parents`: `delta <- parse_tp_notation(notation)`; add delta to C/H/Cl/Br/O/S/F; drop rows with any of those `< 0`; `TP <- notation`; `Molecule_Formula <- create_formula(...)`.

`tp_feasibility_error`: for each notation, `tryCatch(parse_tp_notation)` and return parse error text; if `delta["Br"] < 0` and any class `!has_br`, return a message naming the TP and class; for each class+notation, `apply_tp_to_parents(build_parent_grid(...), notation)` and if `nrow == 0` return infeasible message. Return `NULL` if all OK.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: PASS. If `create_formula` omits count `1` (`Cl` not `Cl1`) that is correct; expected strings above match that.

- [ ] **Step 5: Commit (only if the user asked)**

```bash
git add R/CPions_utils.R tests/testthat/test-cpions-advanced.R
git commit -m "feat: apply TP deltas and reject infeasible chemistry"
```

---

### Task 4: Refactor `getAdduct_advanced()` onto the catalog

**Files:**
- Modify: `R/CPions_utils.R` `getAdduct_advanced()` (~560–645) and `generateInput_Envipat_advanced()` (~211–265)
- Test: `tests/testthat/test-cpions-advanced.R` (existing snapshot + new BCA Cl-decrement test)

**Interfaces:**
- Consumes: `build_parent_grid()`, `apply_tp_to_parents()`, `tp_feasibility_error()`
- Produces: same `getAdduct_advanced(Class, Adduct_Ion, TP, Charge, C, Cl, Clmax, Br, Brmax, threshold)` return shape as today

- [ ] **Step 1: Write failing BCA Cl-decrement test**

```r
test_that("BCA -Cl+OH decrements chlorine on the molecule formula", {
    actual <- CPxplorer:::getAdduct_advanced(
        Class = "BCA",
        Adduct_Ion = "+Cl",
        TP = "-Cl+OH",
        Charge = "-",
        C = 10:10,
        Cl = 4:4,
        Clmax = 4L,
        Br = 2:2,
        Brmax = 2L,
        threshold = 5L
    )
    expect_true(all(actual$Molecule_Formula == "C10H17Cl3Br2O"))
    expect_true(all(actual$Parent_Formula == "C10H16Cl4Br2"))
    expect_true(all(actual$TP == "-Cl+OH"))
})
```

Current code does **not** decrement BCA Cl, so this test should fail until the refactor.

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: FAIL, `Molecule_Formula` still has Cl4 (or similar), not Cl3.

- [ ] **Step 3: Refactor generation**

Replace the PCA/PCO/BCA `if` / `case_when(TP == ...)` block in `getAdduct_advanced()` with:

```r
err <- tp_feasibility_error(Class, TP, C, Cl, Clmax, Br, Brmax)
if (!is.null(err)) {
    stop(err, call. = FALSE)
}
data <- apply_tp_to_parents(
    build_parent_grid(Class, C, Cl, Clmax, Br, Brmax),
    TP
)
```

Then keep the existing Charge mutate, `generateInput_Envipat_advanced()`, `filter(Cl > 0)`, and isotope loop.

If `nrow(data) == 0` after `filter(Cl > 0)`, `stop()` with a message that the TP plus adduct leaves no Cl (spec rule 4).

In `generateInput_Envipat_advanced()`, **delete** the `O` and `S` `case_when(TP == ...)` blocks. Keep adduct Cl/H/Br/F mutations. Ensure incoming `O`, `S`, `F` columns exist:

```r
if (!"O" %in% names(data)) data$O <- 0
if (!"S" %in% names(data)) data$S <- 0
if (!"F" %in% names(data)) data$`F` <- 0
```

Change `Br = ifelse(Compound_Class == "BCA", Br, 0)` to `Br = ifelse(is.na(Br), 0, Br)` so catalog classes with Br still work and non-Br classes stay 0 from the parent grid.

Do not change adduct string construction: `None` → `[Class Adduct_Ion]Charge`, else `[Class TP Adduct_Ion]Charge`.

- [ ] **Step 4: Run tests**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: BCA test PASS; existing `TP = "None"` snapshot PASS (same columns/values). If the snapshot fails, inspect the diff; `None` output must not change. Do not accept a snapshot change unless it is an accidental column-order-only issue — then fix `select()` to preserve the old order.

- [ ] **Step 5: Commit (only if the user asked)**

```bash
git add R/CPions_utils.R tests/testthat/test-cpions-advanced.R
git commit -m "refactor: generate advanced TPs from catalog deltas"
```

---

### Task 5: Extra TPs and PCdiO/PCtriO through `getAdduct_advanced()`

**Files:**
- Test: `tests/testthat/test-cpions-advanced.R`
- Modify: only if Task 4 missed a class in `build_parent_grid` / Skyline is later

**Interfaces:**
- Consumes: `getAdduct_advanced()`
- Produces: molecule formulas for new TPs and new classes

- [ ] **Step 1: Write failing generation tests**

```r
test_that("getAdduct_advanced supports PCdiO, PCtriO, and new TPs", {
    pca_gluc <- CPxplorer:::getAdduct_advanced(
        "PCA", "+Cl", "-H+C6H10O7", "-", 10:10, 6:6, 6L, 0:0, 0L, 5L
    )
    expect_true(all(pca_gluc$Molecule_Formula == "C16H25Cl6O7"))
    expect_true(all(pca_gluc$Parent_Formula == "C10H16Cl6"))

    pcdio <- CPxplorer:::getAdduct_advanced(
        "PCdiO", "+Cl", "None", "-", 10:10, 6:6, 6L, 0:0, 0L, 5L
    )
    expect_true(all(pcdio$Molecule_Formula == "C10H12Cl6"))
    expect_true(all(pcdio$Compound_Class == "PCdiO"))

    pctrio <- CPxplorer:::getAdduct_advanced(
        "PCtriO", "+Cl", "-H+OH", "-", 10:10, 6:6, 6L, 0:0, 0L, 5L
    )
    expect_true(all(pctrio$Molecule_Formula == "C10H10Cl6O"))
    expect_true(all(pctrio$Compound_Class == "PCtriO"))
})
```

If Task 4 already implemented parent classes, this may pass immediately. If it passes on first run, keep the test (it locks the spec). If `getAdduct_advanced` still `stop`s on unknown class, that is the intended RED.

- [ ] **Step 2: Run tests**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: PASS if Task 4 used `parent_class_catalog`; otherwise implement the missing class path (no extra `if (Class == "PCdiO")` formula tables — grid/catalog only).

- [ ] **Step 3: Commit (only if the user asked)**

```bash
git add tests/testthat/test-cpions-advanced.R R/CPions_utils.R
git commit -m "test: cover PCdiO, PCtriO, and glucuronidation generation"
```

---

### Task 6: Regenerate `inst/CPions_TP_formula.xlsx`

**Files:**
- Modify: `R/CPions_utils.R` (helper)
- Modify: `inst/CPions_TP_formula.xlsx`
- Test: `tests/testthat/test-cpions-advanced.R`

**Interfaces:**
- Consumes: `tp_catalog`, `build_parent_grid()`, `apply_tp_to_parents()`, enviPat `isotopes`
- Produces: `build_tp_formula_table()` tibble with the current Excel columns; `write_tp_formula_xlsx(path)` writes it

Excel columns (keep names exactly):

1. `Name of TP`
2. `Transformation product`
3. `General formula`
4. `Example parent formula`
5. `Example molecule formula`
6. `Example adduct ion`
7. `Example adduct ion mz (calculated exact mass of the monoisotopic adduct ion of the transformation product Use: www.envipat.eawag.ch`
8. `Note`

Example parent: PCA `C10H16Cl6` except Br TPs, which use BCA `C10H16Cl4Br2`. Example adduct is `[M-H]-`. m/z is monoisotopic (highest 12C/35Cl/79Br, i.e. enviPat first/lowest-mass ion or the ion with Rel_ab of the monoisotopic peak). Use `getAdduct_advanced(..., Adduct_Ion = "-H", Charge = "-", threshold = 99)` and take the row with the lowest `m/z` among Rel_ab == max for that formula, or simpler: take the isotopologue `""` (monoisotopic) `m/z`.

General formula column (PCA parent `CxH2x+2-yCly`, Br TPs use BCA `CxH2x+2-y-zClyBrz`):

- `None`: `CxH2x+2-yCly`
- `-Cl+OH`: `CxH2x+2-y+1Cly-1O`
- `-H+OH`: `CxH2x+2-yClyO`
- `-2Cl+2OH`: `CxH2x+2-y+2Cly-2O2`
- `-2H+2OH`: `CxH2x+2-yClyO2`
- `-2H+O`: `CxH2x+2-y-2ClyO`
- `-H+SO4H`: `CxH2x+2-yClyO4S`
- `-H+C6H10O7`: `Cx+6H2x+2-y+9ClyO7`
- `-2H+2O`: `CxH2x+2-y-2ClyO2`
- `-Br+OH`: `CxH2x+2-y-z+1ClyBrz-1O`
- `-2Br+2OH`: `CxH2x+2-y-z+2ClyBrz-2O2`
- `-H+OCH3`: `Cx+1H2x+2-y+2ClyO`
- `-Cl+OCH3`: `Cx+1H2x+2-y+3Cly-1O`
- `-4H+2O`: `CxH2x+2-y-4ClyO2`

Note column: for `-Cl+OH` keep “Give exact same chemical formula as -H+OH with one less Cl …”; for `-H+OH` the inverse; others `NA`.

- [ ] **Step 1: Write failing Excel-shape test**

```r
test_that("TP formula table covers the catalog with Excel example formulas", {
    tbl <- CPxplorer:::build_tp_formula_table()
    expect_identical(tbl[["Transformation product"]], CPxplorer:::tp_catalog_notations())
    pca_none <- tbl[tbl[["Transformation product"]] == "None", ]
    expect_identical(pca_none[["Example parent formula"]], "C10H16Cl6")
    expect_identical(pca_none[["Example molecule formula"]], "C10H16Cl6")
    gluc <- tbl[tbl[["Transformation product"]] == "-H+C6H10O7", ]
    expect_identical(gluc[["Example molecule formula"]], "C16H25Cl6O7")
    br <- tbl[tbl[["Transformation product"]] == "-Br+OH", ]
    expect_identical(br[["Example parent formula"]], "C10H16Cl4Br2")
    expect_identical(br[["Example molecule formula"]], "C10H17Cl4BrO")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: FAIL, `build_tp_formula_table` not found.

- [ ] **Step 3: Implement helper and write the xlsx**

Implement `build_tp_formula_table()` looping `tp_catalog`. For each notation, choose class PCA unless `parse_tp_notation(notation)["Br"] < 0` then BCA. Build parent/molecule via `apply_tp_to_parents`. For m/z, call `getAdduct_advanced` with `-H` / `-` / threshold 5 or 99 and take `Isotopologue == ""`.

```r
write_tp_formula_xlsx <- function(path = "inst/CPions_TP_formula.xlsx") {
    openxlsx::write.xlsx(build_tp_formula_table(), path, overwrite = TRUE)
}
```

From the package root, after `devtools::load_all()`, run `CPxplorer:::write_tp_formula_xlsx()`.

- [ ] **Step 4: Run tests**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: PASS. Spot-check that old `[M-H]-` masses still match the previous file for `None` (344.9310) within rounding to 4 decimals if you round in the table; if current code uses 6 decimals, keep 4 in the Excel column to match the old file (`round(mz, 4)`).

- [ ] **Step 5: Commit (only if the user asked)**

```bash
git add R/CPions_utils.R inst/CPions_TP_formula.xlsx tests/testthat/test-cpions-advanced.R
git commit -m "feat: regenerate CPions TP formula Excel from catalog"
```

---

### Task 7: Skyline parent carbon in `Molecule List Name`

**Files:**
- Modify: `R/CPions_utils.R` `build_skyline_transition_list()` advanced branch (~1057–1068)
- Test: `tests/testthat/test-cpions-advanced.R`

**Interfaces:**
- Consumes: advanced ion table with `Compound_Class`, `TP`, `Parent_Formula`, `Molecule_Formula`
- Produces: `{Class}-C{parent n}` or `{Class}-C{parent n}_{TP}`

- [ ] **Step 1: Write failing Skyline naming test**

```r
test_that("Skyline list name uses parent carbon for TPs that add carbon", {
    ions <- CPxplorer:::getAdduct_advanced(
        "PCA", "+Cl", "-H+C6H10O7", "-", 10:10, 6:6, 6L, 0:0, 0L, 5L
    )
    sky <- CPxplorer:::build_skyline_transition_list(
        ions, mode = "advanced", quant_ion = "Most intense"
    )
    expect_true(all(sky[["Molecule List Name"]] == "PCA-C10_-H+C6H10O7"))
    expect_false(any(grepl("PCA-C16", sky[["Molecule List Name"]])))
})

test_that("Skyline list name supports PCdiO", {
    ions <- CPxplorer:::getAdduct_advanced(
        "PCdiO", "+Cl", "None", "-", 10:10, 6:6, 6L, 0:0, 0L, 5L
    )
    sky <- CPxplorer:::build_skyline_transition_list(
        ions, mode = "advanced", quant_ion = "Most intense"
    )
    expect_true(all(sky[["Molecule List Name"]] == "PCdiO-C10"))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: FAIL, list name `PCA-C16_...` and/or PCdiO not handled (`NA`).

- [ ] **Step 3: Replace class case_when**

In the advanced branch of `build_skyline_transition_list()`, replace the PCA/PCO/BCA `case_when` with:

```r
parent_c <- stringr::str_extract(Parent_Formula, "(?<=C)\\d+(?=H)")
dplyr::mutate(
    `Molecule List Name` = dplyr::case_when(
        stringr::str_detect(Compound_Class, "^IS$") ~ Compound_Class,
        stringr::str_detect(Compound_Class, "^RS$") ~ Compound_Class,
        TP == "None" ~ paste0(Compound_Class, "-C", parent_c),
        TRUE ~ paste0(Compound_Class, "-C", parent_c, "_", TP)
    )
)
```

Keep the rest of the function unchanged.

- [ ] **Step 4: Run tests**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: PASS. Also run `Rscript -e "devtools::test(filter = 'cpions-skyline')"` if Skyline snapshots exist for advanced mode; current skyline tests use normal mode and should be unaffected.

- [ ] **Step 5: Commit (only if the user asked)**

```bash
git add R/CPions_utils.R tests/testthat/test-cpions-advanced.R
git commit -m "fix: Skyline molecule list uses parent carbon count"
```

---

### Task 8: Advanced UI — custom TP checkbox, PCdiO/PCtriO, block submit

**Files:**
- Modify: `R/CPions.R` UI Compound Class (~91–118) and server Advanced reactives (~284–405)
- Test: `tests/testthat/test-cpions-app-smoke.R`

**Interfaces:**
- Consumes: `tp_catalog_notations()`, `parse_tp_list()`, `tp_feasibility_error()`, `getAdduct_advanced()`
- Produces: `input$TP_custom_adv` checkbox (default `FALSE`); `input$TP_text_adv` text; Compound Class includes `PCdiO`, `PCtriO`

- [ ] **Step 1: Extend app smoke inputs so missing new ids fail closed**

In `tests/testthat/test-cpions-app-smoke.R`, add to the Advanced `setInputs`:

```r
TP_custom_adv = FALSE,
TP_text_adv = "",
```

Keep `TP_adv = "None"`. After UI exists, add a second smoke that custom mode with `-H+COOH` errors:

```r
test_that("custom TP invalid grammar blocks advanced submit", {
    skip_if_not_installed("shiny")
    shiny::testServer(CPxplorer:::CPions_server, {
        session$setInputs(
            Cmin_adv = 10, Cmax_adv = 10,
            Clmin_adv = 3, Clmax_adv = 3,
            Brmin_adv = 1, Brmax_adv = 1,
            Compclass_adv = "PCA",
            Adducts_adv = "+Cl",
            Charge_adv = "-",
            TP_adv = "None",
            TP_custom_adv = TRUE,
            TP_text_adv = "-H+COOH",
            threshold_adv = 5,
            ISRS_input_adv = "",
            go_adv = 1
        )
        session$flushReact()
        expect_error(CP_allions_glob_adv())
    })
})
```

`CP_allions_glob_adv` is inside the server function and not returned. Do **not** call it by name unless it is accessible. Prefer:

```r
expect_error(session$flushReact())
```

only if that actually throws. Safer: extract validation into `resolve_advanced_tp_input(custom, selected, text, classes, C, Cl, Clmax, Br, Brmax)` in `CPions_utils.R` and unit-test that instead of testServer internals.

Add this unit test in `test-cpions-advanced.R` (this is the real RED for validation):

```r
test_that("resolve_advanced_tp_input blocks bad custom text and Br on PCA", {
    expect_error(
        CPxplorer:::resolve_advanced_tp_input(
            custom = TRUE, selected = "None", text = "-H+COOH",
            classes = "PCA", C = 10L, Cl = 6L, Clmax = 6L, Br = 0L, Brmax = 0L
        ),
        "COOH"
    )
    expect_error(
        CPxplorer:::resolve_advanced_tp_input(
            custom = FALSE, selected = "-Br+OH", text = "",
            classes = "PCA", C = 10L, Cl = 6L, Clmax = 6L, Br = 0L, Brmax = 0L
        ),
        "Br"
    )
    expect_identical(
        CPxplorer:::resolve_advanced_tp_input(
            custom = TRUE, selected = "None", text = " -H+OH ; -Cl+OH ",
            classes = "PCA", C = 10L, Cl = 6L, Clmax = 6L, Br = 0L, Brmax = 0L
        ),
        c("-H+OH", "-Cl+OH")
    )
    expect_identical(
        CPxplorer:::resolve_advanced_tp_input(
            custom = FALSE, selected = c("None", "-H+OH"), text = "",
            classes = "PCA", C = 10L, Cl = 6L, Clmax = 6L, Br = 0L, Brmax = 0L
        ),
        c("None", "-H+OH")
    )
})
```

Keep smoke `TP_custom_adv = FALSE` so default unchecked still works once the input exists.

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'cpions-advanced')"`

Expected: FAIL, `resolve_advanced_tp_input` not found.

- [ ] **Step 3: Implement resolver + UI + server wiring**

```r
resolve_advanced_tp_input <- function(custom, selected, text, classes, C, Cl, Clmax, Br, Brmax) {
    notations <- if (isTRUE(custom)) {
        parse_tp_list(text)
    } else {
        as.character(selected)
    }
    if (length(notations) == 0L) {
        notations <- "None"
    }
    err <- tp_feasibility_error(classes, notations, C, Cl, Clmax, Br, Brmax)
    if (!is.null(err)) {
        stop(err, call. = FALSE)
    }
    notations
}
```

UI in `R/CPions.R`:

- `choices = c("PCA", "PCO", "PCdiO", "PCtriO", "BCA")` for `Compclass_adv`
- `choices = tp_catalog_notations()` for `TP_adv` (cannot call package internals from UI if not in namespace; use `CPxplorer:::tp_catalog_notations()` inside `CPions()` or hardcode the same vector from `tp_catalog$notation` — calling the helper is required so new catalog rows appear automatically)
- After `TP_adv` selectInput:

```r
shiny::checkboxInput("TP_custom_adv", "Custom transformation product", value = FALSE),
shiny::conditionalPanel(
    condition = "input.TP_custom_adv == true",
    shiny::textInput(
        "TP_text_adv",
        "Transformation products",
        value = "",
        placeholder = "-H+OH; -Cl+OH; -2H+2OH"
    )
),
shiny::conditionalPanel(
    condition = "input.TP_custom_adv == false",
    shiny::selectInput(
        "TP_adv",
        "Transformation product",
        choices = CPxplorer:::tp_catalog_notations(),
        selected = "None",
        multiple = TRUE,
        selectize = TRUE
    )
)
```

Do not duplicate `TP_adv`: remove the old standalone `selectInput("TP_adv", ...)` and keep only the conditional one. Checkbox stays visible.

Server: replace `selectedTP_adv <- eventReactive(..., input$TP_adv)` with:

```r
selectedTP_adv <- shiny::eventReactive(input$go_adv, {
    CPxplorer:::resolve_advanced_tp_input(
        custom = isTRUE(input$TP_custom_adv),
        selected = input$TP_adv,
        text = input$TP_text_adv,
        classes = as.character(input$Compclass_adv),
        C = as.integer(input$Cmin_adv:input$Cmax_adv),
        Cl = as.integer(input$Clmin_adv:input$Clmax_adv),
        Clmax = as.integer(input$Clmax_adv),
        Br = as.integer(input$Brmin_adv:input$Brmax_adv),
        Brmax = as.integer(input$Brmax_adv)
    )
})
```

In `CP_allions_glob_adv`, `TP <- selectedTP_adv()` already. Wrap generation so `stop()` surfaces as a Shiny validation error:

```r
shiny::validate(shiny::need(TRUE, ""))
```

is wrong. Use:

```r
notations <- tryCatch(selectedTP_adv(), error = function(e) e)
shiny::validate(shiny::need(!inherits(notations, "error"), notations$message))
```

or let `eventReactive` throw and `validate` in the DT render. Preferred: `tryCatch` around `resolve_advanced_tp_input` inside `CP_allions_glob_adv` and `shiny::validate(shiny::need(is.null(err), err))` before the nested loops.

When `TP_custom_adv` is TRUE, `input$TP_adv` may be NULL; resolver must accept `selected = NULL`.

- [ ] **Step 4: Run tests**

Run:

```
Rscript -e "devtools::test(filter = 'cpions-advanced')"
Rscript -e "devtools::test(filter = 'cpions-app-smoke')"
```

Expected: both PASS. Smoke must set `TP_custom_adv = FALSE` and `TP_text_adv = ""`.

- [ ] **Step 5: Commit (only if the user asked)**

```bash
git add R/CPions.R R/CPions_utils.R tests/testthat/test-cpions-advanced.R tests/testthat/test-cpions-app-smoke.R
git commit -m "feat: custom TP input and PCdiO/PCtriO in Advanced settings"
```

---

### Task 9: Instructions

**Files:**
- Modify: `inst/instructions_CPions.md` Advanced settings section (lines 37–38)

**Interfaces:**
- Consumes: catalog notations and parent H formulas
- Produces: short user-facing docs rendered by `includeMarkdown`

- [ ] **Step 1: Replace the Advanced settings paragraph**

Replace:

```
## Advanced settings tab  
Mostly same initial parameters as Normal settings. In advanced settings, there is more flexibility to combine and mix the `Compound Class`, `Adduct`, `Charge`, and `Transformation product`.  
```

with:

```
## Advanced settings tab  
Mostly same initial parameters as Normal settings. Compound class, adduct, charge, and transformation product can be combined freely.  

__Compound Class__: `PCA` (alkane, H = 2C+2−Cl), `PCO` (mono-olefin, H = 2C−Cl), `PCdiO` (di-olefin, H = 2C−2−Cl), `PCtriO` (tri-olefin, H = 2C−4−Cl), `BCA` (bromo-chloro alkane, H = 2C+2−Cl−Br). PCdiO and PCtriO are only in Advanced settings. Br min/max apply to BCA.  

__Transformation product__: predefined notations (multiple allowed): `None`, `-Cl+OH`, `-H+OH`, `-2Cl+2OH`, `-2H+2OH`, `-2H+O`, `-H+SO4H`, `-H+C6H10O7`, `-2H+2O`, `-Br+OH`, `-2Br+2OH`, `-H+OCH3`, `-Cl+OCH3`, `-4H+2O`. `-Br+OH` / `-2Br+2OH` require BCA only.  

__Custom transformation product__: unchecked by default. When checked, the dropdown is replaced by a text field. Enter one or more notations separated by `;`, e.g. `-H+OH; -Cl+OH; -2H+2OH`. `-H` is loss of one H; `+OH` is gain of OH; `-2H` is loss of two H. Empty input is `None`. Invalid grammar or a TP that is impossible for the selected class and C/Cl/Br range blocks calculation (no table). `-2H+O` is valid (ketone).  
```

Keep the existing Skyline warning that `-Cl+OH` and `-H+OH` can give the same formula.

No HTML file exists; do not create one.

- [ ] **Step 2: Confirm the Instructions tab still points at this file**

`R/CPions.R` uses `system.file("instructions_CPions.md", package = "CPxplorer")`. No path change.

- [ ] **Step 3: Commit (only if the user asked)**

```bash
git add inst/instructions_CPions.md
git commit -m "docs: describe custom TPs and PCdiO/PCtriO in CPions instructions"
```

---

### Task 10: Final verification

**Files:** none new

- [ ] **Step 1: Run focused tests**

```
Rscript -e "devtools::test(filter = 'cpions-advanced')"
Rscript -e "devtools::test(filter = 'cpions-app-smoke')"
```

Expected: PASS.

- [ ] **Step 2: Run interference/skyline/normal only if shared helpers drifted**

If `build_skyline_transition_list` or `getAdduct_advanced` None-path changed:

```
Rscript -e "devtools::test(filter = 'cpions-skyline')"
Rscript -e "devtools::test(filter = 'cpions-normal')"
```

Expected: PASS. Do not claim CPquant coverage.

- [ ] **Step 3: Manual checklist**

- Catalog add-a-TP later = one `tribble` row + rerun `write_tp_formula_xlsx()`.
- Custom checkbox default unchecked.
- BCA `-Cl+OH` molecule Cl decremented.
- Glucuronidation Skyline name `PCA-C10_-H+C6H10O7`.
- Excel has all catalog rows including glucuronidation and Br TPs.

---

## Spec coverage

| Spec requirement | Task |
|---|---|
| Parser + group table | 1 |
| `tp_catalog` notation+name | 2 |
| Parent class catalog / PCdiO / PCtriO H counts | 2, 5, 8 |
| Apply deltas / feasibility / mixed PCA+BCA Br block / partial range | 3 |
| Refactor `getAdduct_advanced` / BCA Cl bug / O,S from deltas | 4 |
| Extra TPs through generation | 5 |
| Excel regenerated from catalog | 6 |
| Skyline parent carbon | 7 |
| Custom checkbox, text grammar, block submit | 8 |
| Instructions | 9 |
| Existing None snapshot unchanged | 4, 10 |
| Advanced-only olefins (no Normal adducts) | 8, 10 |
| TDD | each task Steps 1–4 |
