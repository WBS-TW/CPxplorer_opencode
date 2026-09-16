# Advanced Transformation Products Catalog and Custom Input

Date: 2026-09-03

## Problem

Advanced settings currently hardcodes transformation products (TPs) in three duplicated `case_when` branches (PCA, PCO, BCA) in `getAdduct_advanced()` and again for O/S counts in `generateInput_Envipat_advanced()`. Adding a TP requires editing several parallel tables. The packaged Excel file `inst/CPions_TP_formula.xlsx` already lists glucuronidation (`-H+C6H10O7`), which is not in the UI. BCA `-Cl+OH` / `-2Cl+2OH` currently do not decrement Cl.

Users also need a custom TP text input with the same notation as the dropdown (`-H+OH; -Cl+OH; -2H+2OH`), with grammar and chemical feasibility checks that block submit on failure.

This work also adds two Advanced-only parent classes: PCdiO (di-olefin) and PCtriO (tri-olefin).

## Goals

- One R catalog of predefined TPs (notation + display name). Parser computes element deltas.
- Same parser/applier for catalog TPs and custom text.
- Custom TP checkbox (unchecked by default) swaps the dropdown for a text input.
- Invalid grammar or infeasible chemistry blocks submit; no partial table.
- Extra predefined TPs: glucuronidation, carboxylic acid/ω-oxidation, Br-hydroxylation, methoxy, diketone.
- PCdiO and PCtriO in Advanced Compound Class only.
- `inst/CPions_TP_formula.xlsx` regenerated from the catalog.
- Document the new classes, TPs, and custom grammar in `inst/instructions_CPions.md`.
- Adding a TP later is one catalog row; Excel is rebuilt from that table.

## Non-goals

- No Normal-settings adducts for PCdiO/PCtriO (`[PCdiO+Cl]-` etc.).
- Excel is not the runtime source of TPs.
- No multi-charge ions, no new adduct types, no CPquant changes.
- Do not hand-edit generated `NAMESPACE` / `man/*.Rd`.

## Architecture

Keep generation inside `R/CPions_utils.R` and the Advanced UI in `R/CPions.R`. Extract TP and parent-class rules so `getAdduct_advanced()` no longer branches on TP strings.

### Catalog

Internal tibble `tp_catalog` with two columns:

- `notation`: exact dropdown / parser string (`None`, `-H+OH`, …)
- `name`: short display name (`H-hydroxylation`, …)

The parser derives net element deltas (`C, H, Cl, Br, O, S, F`) from `notation`. Developers never hand-edit atom counts.

To add a TP later:

1. Add one `tribble` row (`notation`, `name`).
2. Rebuild `inst/CPions_TP_formula.xlsx` from the catalog (same helper used in this change).
3. Add/adjust a focused test if the new notation introduces a new group token.

### Parent class catalog

Internal tibble `parent_class_catalog`:

| Class  | H formula              | Has Br | Scope              |
|--------|------------------------|--------|--------------------|
| PCA    | `2*C+2-Cl`             | no     | Normal + Advanced  |
| PCO    | `2*C-Cl`               | no     | Normal + Advanced  |
| PCdiO  | `2*C-2-Cl`             | no     | Advanced only      |
| PCtriO | `2*C-4-Cl`             | no     | Advanced only      |
| BCA    | `2*C+2-Cl-Br`          | yes    | Normal + Advanced  |

Advanced Compound Class choices become `c("PCA", "PCO", "PCdiO", "PCtriO", "BCA")`. Br min/max still apply only when BCA is selected (same as today for PCA/PCO).

### Pipeline in `getAdduct_advanced()`

1. Resolve parent atoms from `parent_class_catalog` (C/Cl crossing; BCA also Br, with `Br+Cl <= C`).
2. Parse TP to a delta vector (`None` is all zeros).
3. Add deltas to parent atoms; reject if any count is negative or the TP is impossible for the selected class/range (see Feasibility).
4. Set `Parent_Formula` from parent atoms and `Molecule_Formula` from transformed atoms via existing `create_formula()`.
5. Continue through existing `generateInput_Envipat_advanced()` using transformed C/H/Cl/Br/S/O/F, not TP-string `case_when`.

This also fixes BCA: `-Cl+OH` decrements Cl.

Skyline `Molecule List Name` currently special-cases PCA/PCO/BCA and extracts carbon count from `Molecule_Formula`. Extend it to any Advanced class: `{Class}-C{n}` or `{Class}-C{n}_{TP}`. Use the **parent** carbon count (pre-TP), not `Molecule_Formula`, so glucuronidation (`C+6`) still groups as `PCA-C10_...` rather than `PCA-C16_...`.

`Halo_perc` already uses `calculate_haloperc(Molecule_Formula)` in the advanced path; no class-specific change.

## UI (Advanced settings)

Keep the current TP dropdown. Add checkbox `Custom transformation product`, unchecked by default.

- Unchecked: existing `selectInput("TP_adv")` with `None` plus all catalog notations except that `None` remains first. Multiple selection stays.
- Checked: hide the dropdown; show `textInput` / `textAreaInput` with placeholder `-H+OH; -Cl+OH; -2H+2OH`. Split on `;`, trim whitespace. Empty input is `None`.

PCdiO and PCtriO appear in Compound Class. Submit still waits for Calculate (`go_adv`). Invalid grammar or infeasible TPs for the selected class and C/Cl/(Br) range block the run with `shiny::validate()` / `need()` (or equivalent) and a clear error; no partial table.

## Grammar and parser

Same notation as the catalog.

- Tokens: optional `-`/`+`, optional integer (default 1), then a group.
- Allowed groups: `H`, `Cl`, `Br`, `OH`, `O`, `SO4H`, `OCH3`, `C6H10O7`.
- Tokens concatenate: `-2H+2OH`.
- `;` separates TPs. `None` is allowed as a token or as empty custom input.
- Whitespace around `;` is ignored.

Parser output is a named integer delta vector. Worked examples:

| Notation        | C  | H  | Cl | Br | O | S | F |
|-----------------|----|----|----|----|---|---|---|
| `None`          | 0  | 0  | 0  | 0  | 0 | 0 | 0 |
| `-H+OH`         | 0  | 0  | 0  | 0  | 1 | 0 | 0 |
| `-Cl+OH`        | 0  | +1 | −1 | 0  | 1 | 0 | 0 |
| `-2Cl+2OH`      | 0  | +2 | −2 | 0  | 2 | 0 | 0 |
| `-2H+2OH`       | 0  | 0  | 0  | 0  | 2 | 0 | 0 |
| `-2H+O`         | 0  | −2 | 0  | 0  | 1 | 0 | 0 |
| `-H+SO4H`       | 0  | 0  | 0  | 0  | 4 | 1 | 0 |
| `-H+C6H10O7`    | +6 | +9 | 0  | 0  | 7 | 0 | 0 |
| `-2H+2O`        | 0  | −2 | 0  | 0  | 2 | 0 | 0 |
| `-Br+OH`        | 0  | +1 | 0  | −1 | 1 | 0 | 0 |
| `-2Br+2OH`      | 0  | +2 | 0  | −2 | 2 | 0 | 0 |
| `-H+OCH3`       | +1 | +2 | 0  | 0  | 1 | 0 | 0 |
| `-Cl+OCH3`      | +1 | +3 | −1 | 0  | 1 | 0 | 0 |
| `-4H+2O`        | 0  | −4 | 0  | 0  | 2 | 0 | 0 |

Group expansion used by the parser:

| Group     | C | H | Cl | Br | O | S | F |
|-----------|---|---|----|----|---|---|---|
| `H`       | 0 | 1 | 0  | 0  | 0 | 0 | 0 |
| `Cl`      | 0 | 0 | 1  | 0  | 0 | 0 | 0 |
| `Br`      | 0 | 0 | 0  | 1  | 0 | 0 | 0 |
| `OH`      | 0 | 1 | 0  | 0  | 1 | 0 | 0 |
| `O`       | 0 | 0 | 0  | 0  | 1 | 0 | 0 |
| `SO4H`    | 0 | 1 | 0  | 0  | 4 | 1 | 0 |
| `OCH3`    | 1 | 3 | 0  | 0  | 1 | 0 | 0 |
| `C6H10O7` | 6 | 10| 0  | 0  | 7 | 0 | 0 |

Sign and coefficient apply to the whole group (`-H+OH` = lose 1 H, gain 1 OH).

## Feasibility (block submit)

Fail the whole Advanced run if any of these is true:

1. Parse error: unknown group, missing sign, empty token after split, non-integer coefficient.
2. Class mismatch: any Br loss (`-Br`, `-nBr`) when **any** selected class has no Br (PCA, PCO, PCdiO, PCtriO). Mixed PCA+BCA with `-Br+OH` therefore blocks; the user must select BCA only.
3. After applying deltas, **every** parent in the selected C/Cl/(Br) range for a selected class has a negative atom count (the TP is impossible for that class+range). Example: Cl min=max=1 with `-2Cl+2OH`.
4. If some homologues in the range are valid and some are not (Cl min=1, max=15 with `-2Cl+2OH`), keep current filtering: drop impossible rows (`Cl > 0` after adduct, no negative atoms) and generate the rest. Block only when the TP leaves zero valid rows for a selected class.

Error messages must name the failing TP string and the reason (parse vs infeasible vs class).

Predefined dropdown TPs use the same rules. Example: `-Br+OH` with PCA selected blocks.

## Predefined catalog contents

| Notation       | Name                                      |
|----------------|-------------------------------------------|
| `None`         | None                                      |
| `-Cl+OH`       | Cl-hydroxylation                          |
| `-H+OH`        | H-hydroxylation                           |
| `-2Cl+2OH`     | Double Cl-hydroxylation                   |
| `-2H+2OH`      | Double H-hydroxylation                    |
| `-2H+O`        | Oxidation (ketone/aldehyde)               |
| `-H+SO4H`      | Sulfonation                               |
| `-H+C6H10O7`   | Glucuronidation                           |
| `-2H+2O`       | Carboxylic acid / ω-oxidation             |
| `-Br+OH`       | Br-hydroxylation                          |
| `-2Br+2OH`     | Double Br-hydroxylation                   |
| `-H+OCH3`      | Methoxylation                             |
| `-Cl+OCH3`     | Cl-methoxylation                          |
| `-4H+2O`       | Diketone                                  |

All except the Br TPs are valid on PCA, PCO, PCdiO, PCtriO, and BCA when Cl (and H) counts allow. Br TPs are BCA-only.

## Excel asset

`inst/CPions_TP_formula.xlsx` is regenerated from the catalog, not hand-maintained. Keep the current columns:

- Name of TP
- Transformation product
- General formula
- Example parent formula
- Example molecule formula
- Example adduct ion
- Example adduct ion mz (monoisotopic `[M-H]-`)
- Note

Example parent is PCA `C10H16Cl6` except Br TPs, which use BCA `C10H16Cl4Br2`. Example `[M-H]-` m/z is computed with the same formula path (enviPat / existing helpers). Preserve the existing note that `-Cl+OH` and `-H+OH` can yield the same formula at adjacent Cl counts.

A small internal helper may write this file so catalog and Excel cannot drift. The helper is for package maintenance, not a user-facing export.

## Instructions

Update `inst/instructions_CPions.md` Advanced settings section (keep entries short):

- Compound classes: PCA, PCO, PCdiO, PCtriO, BCA, with H formulas.
- Predefined TP list and what each notation means.
- Custom TP checkbox, grammar (`-H+OH; -Cl+OH`), and that invalid/infeasible input blocks submit.
- Keep the existing warning that some TPs can give the same molecular formula (`-Cl+OH` vs `-H+OH`).

If `inst/instructions_CPions.html` is still in use, update it to match the markdown.

## Testing

Follow TDD. Prefer `devtools::test(filter = "cpions-advanced")`. Do not update the existing advanced snapshot unless `TP = "None"` output actually changes (it should not).

Cover:

- Parser: valid tokens, `;` lists, `None`, whitespace, bad grammar.
- Apply-TP: PCA, PCO, BCA, PCdiO, PCtriO; Excel example formulas from `C10H16Cl6` (and BCA parent for Br TPs).
- Feasibility: `-Br+OH` on PCA fails; mixed PCA+BCA with `-Br+OH` fails; `-2Cl+2OH` with Cl min=max=1 fails; `-2H+O` on PCA `C10H16Cl6` yields `C10H14Cl6O`.
- Custom multi-TP string expands to the same formulas as selecting those catalog entries.
- App smoke: new checkbox default unchecked still submits with `TP_adv = "None"`.

Verify PCA/PCO/BCA (and PCdiO/PCtriO) molecule formulas against the Excel examples after regeneration.

## Error handling

- Custom empty string → `None`.
- Custom `" -H+OH ; -Cl+OH "` → two TPs.
- Unknown group e.g. `-H+COOH` → parse error, block.
- `-2H+O` is valid (ketone).
- Combining custom mode with the dropdown is impossible; checkbox is exclusive.

## Files expected to change

- `R/CPions_utils.R` — catalog, parser, apply-TP, `getAdduct_advanced()`, Skyline class naming, Excel rebuild helper if kept in-package
- `R/CPions.R` — Compound Class choices, custom TP checkbox/text input, submit validation
- `inst/CPions_TP_formula.xlsx`
- `inst/instructions_CPions.md` (and `.html` if still used)
- `tests/testthat/test-cpions-advanced.R` and possibly `test-cpions-app-smoke.R`

## Verification

- `devtools::test(filter = "cpions-advanced")`
- `devtools::test(filter = "cpions-app-smoke")` if Advanced inputs change
- Do not claim CPquant or Normal-settings coverage unless those tests are run
