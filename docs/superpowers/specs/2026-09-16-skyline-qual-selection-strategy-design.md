# Skyline Qual Ion Selection by Strategy

Date: 2026-09-16

## Problem

In the CPions Skyline tab, **Selection strategy** currently changes only Quan ion choice. Qual ions are always taken from leftover candidates in abundance order, and only if they do not interfere with already-selected ions. Interfering leftovers are skipped, so the requested Qual count can be short even when extra ions exist.

The three strategies should apply to Qual selection as well as Quan.

## Goals

- Apply the same three interference-filtered strategies to Qual ions:
  - Highest abundance
  - Balanced
  - Least interference
- After Quan is chosen, rank remaining ions for that molecule with the same strategy, then take up to `preferred_qual_n` Quals.
- If Balanced or Least interference cannot find enough non-interfering Quals, fill remaining slots with the next-best interfering leftovers and flag them `[INTERFERENCE]`.
- Highest abundance Quals may interfere; still keep them and flag `[INTERFERENCE]`.
- Quals must still not collide with already-selected ions from **other** molecules.
- Update Skyline UI help text and `inst/instructions_CPions.md` so strategy descriptions cover Qual, not only Quan.
- Cover the new Qual behavior with focused tests in `tests/testthat/test-cpions-skyline.R`.

## Non-goals

- Do not change **Most intense** Quant Ion mode. Strategy and Qual-count options still do not apply there.
- Do not change Skyline export column names.
- Do not change interference calculation, Normal/Advanced ion generation, or CPquant.
- Do not replace the sequential global interference check with a single per-molecule ranked slice (Quan first, next N Qual). That would ignore collisions with other molecules’ already-selected ions.
- Do not hand-edit generated `NAMESPACE` / `man/*.Rd`.

## Architecture

Keep the change local to Skyline interference-filtered selection in `R/CPions_utils.R`, plus UI copy in `R/CPions.R` and instructions in `inst/instructions_CPions.md`.

Current flow in `select_skyline_filtered_candidates()`:

1. `prepare_skyline_candidates()` annotates interference diagnostics.
2. `select_quan_by_strategy()` picks one Quan per `Molecule Name`.
3. Remaining candidates are walked in `Original Rel_ab Rank` order.
4. A leftover becomes Qual only if the molecule still needs Quals **and** `has_ms_interference()` is false against the growing selected m/z set.
5. Interfering leftovers are skipped; Qual count may stay below `preferred_qual_n`.

New flow:

1. Unchanged: prepare candidates and select Quan by strategy. All Quans are chosen before any Qual.
2. Process molecules in the same order Quan selection uses: `unique(candidates$Molecule Name)`.
3. For each molecule, take leftover candidates (not the chosen Quan).
4. Rank those leftovers with the same strategy used for Quan.
5. Walk that molecule’s ranked leftover list and keep up to `preferred_qual_n` Quals, using the Qual acceptance rules below. Later molecules see earlier molecules’ Quan+Qual in `selected_mz`.
6. Continue to recompute final transition interference and Note flags on the selected Quan+Qual set.

Do not pick Quan and Qual from one combined ranked list. Quan stays first via `select_quan_by_strategy()`. Qual ranking is applied only to leftovers. Do not walk all leftovers globally by `Original Rel_ab Rank`.

### Qual ranking

Reuse the Quan ranking rules on leftover candidates.

**Highest abundance (`abundance`)**

- Sort leftovers by descending `Rel_ab`, then `Precursor m/z`.
- Consider them in that order until `preferred_qual_n` Quals are kept or leftovers run out.

**Balanced (`balanced`)**

- Sort leftovers by descending `Rel_ab`, then `Precursor m/z`.
- Prefer candidates whose `.original_interference == "NO"`.
- If those are exhausted before reaching `preferred_qual_n`, continue with interfering leftovers in the same abundance order.

**Least interference (`interference`)**

- Sort leftovers by `.original_interference` (`NO` before `YES`), then `Interfering Candidate Count`, then `Resolution Needed` (missing treated as 0), then descending `Rel_ab`, then `Precursor m/z`.
- Walk that order until `preferred_qual_n` Quals are kept or leftovers run out.
- This already prefers non-interfering ions; remaining slots are filled by the next-best interfering leftovers.

### Qual acceptance vs already-selected ions

Two interference notions stay distinct:

- **Original interference** (`.original_interference`): from the full candidate pool / Interfering ions tab. Used for ranking and for fill-to-N fallback.
- **Collision with already-selected ions**: `has_ms_interference(candidate_mz, selected_mz, ms_resolution)` against Quan and Qual ions already accepted, including other molecules.

Acceptance rules while walking the ranked leftovers:

1. If the molecule already has `preferred_qual_n` Quals, stop considering that molecule.
2. If the candidate collides with an already-selected ion from **another** molecule, skip it. Do not keep it as a fill-to-N Qual.
3. Same-molecule leftovers that collide only with this molecule’s own already-selected ions (its Quan or an earlier Qual):
   - `abundance`: keep them (they are still the highest-abundance leftovers).
   - `balanced` / `interference`: keep them only as fill-to-N after leftovers that do not collide with this molecule’s selected ions are used, or when no such leftover remains.
4. Flag kept interfering Quals in `Note` with `[INTERFERENCE]`, same as Quan.

Implement the other-molecule skip by testing `has_ms_interference()` against selected m/z from other molecules only. Do not use the full `selected_mz` vector for that skip, or same-molecule Quan/Qual pairs would be dropped.

`preferred_qual_n = 0` still means no Quals.

If leftovers run out, the Qual count may still be below N. Fill-to-N uses interfering leftovers that were skipped only for original/same-molecule interference, not leftovers skipped because they collide with another molecule’s selected ions.

### Selected Reason

Keep existing Quan reasons unchanged.

Qual reasons:

- Non-interfering Qual: `Qual: non-interfering` (current text).
- Interfering Qual kept because the strategy allows or requires fill-to-N: `Qual: fallback interfering`.

`Ion Selection Rank` for Quals remains the candidate’s `Original Rel_ab Rank`.

### UI and instructions

`R/CPions.R` Skyline sidebar:

- Help text currently says the strategy selects the Quant ion. Change it so strategy applies to both Quan and Qual.
- Strategy labels can stay `Highest abundance Quan`, `Balanced`, `Least interference Quan` if the help/instructions make Qual coverage explicit. Prefer renaming labels to drop the Quan-only wording:
  - `Highest abundance`
  - `Balanced`
  - `Least interference`
- Keep input IDs (`skyline_strategy`, values `abundance` / `balanced` / `interference`).

`inst/instructions_CPions.md` Skyline section:

- Highest abundance: most abundant leftover Quals are kept even if they interfere; flagged `[INTERFERENCE]`.
- Balanced: highest-abundance non-interfering Quals first; if still short of N, fill with next-best interfering leftovers and flag them.
- Least interference: Quals ranked by interference then abundance; same fill-to-N fallback.
- Preferred number of Qual ions: target count per `Molecule Name`; may still be lower if too few leftovers remain after skipping collisions with other molecules’ selected ions.

If `inst/instructions_CPions.html` is still in use, update it to match the markdown.

## Testing

Keep existing Skyline snapshots/fixtures on `balanced` unless Qual picks actually change. Do not update snapshots just to make tests pass.

Add focused tests in `tests/testthat/test-cpions-skyline.R` using small constructed tibbles (same style as current strategy tests):

1. **Abundance + Qual N > 0**: Quals are the next-highest `Rel_ab` leftovers, including interfering ones; those Notes contain `[INTERFERENCE]`.
2. **Balanced fill-to-N**: non-interfering Quals first; if short of N, remaining Quals are interfering leftovers and are flagged.
3. **Least interference Qual order**: Qual order follows interference then `Rel_ab`; same fill-to-N fallback.
4. **Cross-molecule collision**: a leftover that collides with another molecule’s already-selected ion is not kept as Qual, even if that would help reach N.

Existing tests with `preferred_qual_n = 0L` stay Quan-only and should remain valid.

Verification:

- `devtools::test(filter = "cpions-skyline")`
- `devtools::test(filter = "cpions-app-smoke")` if Skyline UI copy/labels change
- Review any snapshot diffs before accepting them

## Error handling

- `ms_resolution` remains required for interference-filtered selection.
- Invalid/missing strategy still defaults to `balanced`.
- Invalid/missing `preferred_qual_n` still becomes `2L`, then `max(0L, as.integer(...))`.
- No new user-facing errors.

## Open decisions

None. Strategy applies to Qual as well as Quan; leftover ranking reuses Quan rules; Balanced/Least interference fill remaining Qual slots with interfering leftovers; cross-molecule collisions are never used as fill-to-N Quals.
