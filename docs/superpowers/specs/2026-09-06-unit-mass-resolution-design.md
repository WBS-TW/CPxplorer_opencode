# Unit Mass Resolution for CPions Interfering Ions

Date: 2026-09-06

## Problem

The Interfering ions tab always uses a numeric MS Resolution (default 20,000) and exact high-resolution m/z values. Users who work with low-resolution quadrupole instruments need a unit-mass mode: hide the resolution input, treat generated masses as integer (nominal) m/z, flag ions that share the same integer m/z as interfering, and export those integer masses to Skyline.

## Goals

- Add an unchecked-by-default checkbox named **Unit mass resolution** on the Interfering ions tab.
- When checked, hide the **MS Resolution** numeric input.
- When Calculate runs with the checkbox on, round `m/z` to nearest integer before interference calculation.
- Ions with the same integer m/z are flagged as interfering (`YES`).
- Skyline `mz` export uses the integer m/z already stored in the interference table. No extra Skyline control.
- Document the checkbox in `inst/instructions_CPions.md`.
- Cover the new behavior with focused interference tests. Existing high-res snapshots stay unchanged.

## Non-goals

- No Unit mass resolution checkbox on the Skyline tab.
- Do not change Normal settings or Advanced settings ion generation; those tables stay high-res.
- Do not add a separate nominal-mass column; `m/z` itself becomes integer in the interference table.
- No CPquant changes.
- Do not hand-edit generated `NAMESPACE` / `man/*.Rd`.

## Architecture

Keep the change local to the Interfering ions Calculate path in `R/CPions.R`, plus a small helper in `R/CPions_utils.R` if rounding needs a named, testable function. Reuse existing `compute_interference()`; do not add a `unit_mass` argument to it.

### UI

Interfering ions sidebar (`R/CPions.R`):

1. `shiny::checkboxInput("unit_mass_resolution", "Unit mass resolution", value = FALSE)`
2. Wrap the existing `numericInput("MSresolution", ...)` in `shiny::conditionalPanel(condition = "!input.unit_mass_resolution", ...)` so it disappears when the checkbox is on.
3. Leave `interfere_mode` and Calculate unchanged.

### Calculate path

On `input$go2`:

1. Load the Normal or Advanced ion table as today.
2. If `input$unit_mass_resolution` is TRUE, round the `m/z` column with `round(x, 0)` (R default, including half-to-even).
3. Call `compute_interference()`:
   - Unit-mass interference rule: after rounding, two ions interfere if and only if they share the same integer m/z (`difflag == 0` or `difflead == 0`). Neighbors at different integers do not interfere.
   - Implementation: after rounding, call `compute_interference()` with `ms_resolution = .Machine$integer.max` so only exact-zero deltas flag `YES`. Do not pass `1L`: after rounding, distinct integers have `difflag >= 1` and `reslag` equals at most the m/z value, so `ms_resolution = 1L` would mark every neighbor as `YES`.
   - When the checkbox is off, keep current behavior: `compute_interference(..., MSresolution())` with the numeric input (default 20000).
4. Store the result in `CP_allions_compl2` as today (plots, table, downloads, Skyline).

`MSresolution` eventReactive currently always reads `input$MSresolution`. When the checkbox is on, that input is hidden but still holds its last value. Skyline Interference-filtered mode already uses `MSresolution()` / `input$MSresolution`. After unit-mass Calculate, Skyline must use the integer `m/z` from `CP_allions_compl2` and treat same-integer masses as interfering. Passing the hidden 20000 value would re-evaluate high-res interference on already-rounded masses and mis-flag neighbors 1 Da apart as non-interfering (correct) but would also treat identical integer m/z as interfering (also correct, because delta == 0). For Interference-filtered Skyline, `has_ms_interference()` already treats `delta == 0` as interference, so the numeric resolution is irrelevant once masses are integers. Most intense Skyline currently reads Normal/Advanced tables (exact m/z), not `CP_allions_compl2`.

### Skyline export

User choice: auto from interference table.

- **Interference-filtered**: already reads `CP_allions_compl2()`. After unit-mass Calculate, `Precursor m/z` is the integer `m/z`. No Skyline UI change. Keep passing `MSresolution()`; same-integer interference is already handled by `delta == 0`.
- **Most intense**: currently reads `CP_allions_glob()` / `CP_allions_glob_adv()` (exact m/z). After any Interfering ions Calculate, Most intense must also use `CP_allions_compl2()` so unit-mass Calculate exports integer precursor m/z. If `CP_allions_compl2()` is still NULL, Most intense keeps using the Normal/Advanced table (current behavior).
- Rule: if `CP_allions_compl2()` exists, Skyline `mz` export (both Most intense and Interference-filtered) uses that table’s `m/z`, regardless of whether `skyline_mode` matches `interfere_mode`. Interference-filtered still requires Calculate first.
- Do not add a Skyline checkbox.

This means: after unit-mass Calculate, both Skyline Quant Ion modes export integer precursor m/z. After a later high-res Calculate (checkbox off), both modes export exact m/z from the new interference table. After generating ions but never clicking Calculate on Interfering ions, Most intense still exports exact m/z from Normal/Advanced.

### Helper

Prefer a tiny helper so tests do not need the Shiny app:

```r
apply_unit_mass_mz <- function(CP_allions) {
    dplyr::mutate(CP_allions, `m/z` = round(`m/z`, 0))
}
```

Unit-mass interference is then:

```r
compute_interference(apply_unit_mass_mz(CP_allions), ms_resolution = .Machine$integer.max)
```

Keep `compute_interference()` signature unchanged.

## Error handling

- Checkbox off: unchanged validation and Calculate flow.
- Checkbox on with no prior Normal/Advanced Submit: same as today (Calculate on empty/null tables).
- Hidden MS Resolution is not required and is not validated while the checkbox is on.
- Unchecking the checkbox shows MS Resolution again with its previous value; the user must click Calculate again to recompute high-res interference.

## Testing

Add focused tests in `tests/testthat/test-cpions-interference.R` (and Skyline if Most intense source table changes):

1. `apply_unit_mass_mz()` rounds `m/z` to nearest integer.
2. After unit-mass rounding + `compute_interference(.Machine$integer.max)`, ions that share an integer m/z are `YES`; neighbors at different integers are `NO`.
3. Existing `interference_fixture()` / high-res snapshot is unchanged.
4. Skyline built from a unit-mass interference table has integer `Precursor m/z`.
5. App smoke: Interfering ions tab still calculates with default (checkbox off). Optionally set the checkbox and Calculate; do not require a new snapshot for the full app.

Do not update existing high-res snapshots unless output actually changes (it should not).

## Instructions

Update `inst/instructions_CPions.md` Interfering ions section:

- Mention **Unit mass resolution** (unchecked by default).
- When checked, MS Resolution is hidden and m/z values are rounded to integers to reflect a low-resolution (quadrupole) mass spectrometer.
- Ions with the same nominal m/z are treated as interfering.
- Skyline export uses those integer masses after Calculate.

Keep the existing R=20,000 default text for the unchecked path.

## Verification

- `devtools::test(filter = "cpions-interference")`
- `devtools::test(filter = "cpions-skyline")` if Skyline source-table behavior changes
- `devtools::test(filter = "cpions-app-smoke")`
- `devtools::load_all(); CPions()` for UI: checkbox hides MS Resolution; Calculate produces integer m/z
