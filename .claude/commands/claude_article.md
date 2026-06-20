---
description: Standardize the paper-library filenames (English papers -> year_author_ShortTitle), renumber each folder 1.2.3..., then commit, push, and fetch origin.
argument-hint: "(no args) — optionally: dry-run"
allowed-tools: Bash, Read, Write, Edit, Glob, Grep, Workflow, Agent
---

# /claude_article — organize the LBM paper library

You are organizing a research-paper library (this repository). Apply the rules below to **every folder**, then commit + push + fetch. If the user passed `dry-run` in `$ARGUMENTS`, build and show the plan but DO NOT rename, commit, or push.

## Naming rules (agreed conventions)

1. **Chinese-named files** (filename contains any CJK Han character) → **keep the name unchanged**. They are already curated by the user.
2. **Already-curated English files** → **keep the name unchanged**. A file is "already curated" if it:
   - contains a `[tag]` such as `[CLBM]`, `[MRT-CM]`, `[MRT]`, **or**
   - already clearly contains the first author's surname (e.g. `Zhaoli Guo`, `Gehrke`, `Kramer-2017`, `Frapolli`, `Reyes Barraza2019`, `MBreuer`, `Heng Xiao`).
3. **Raw-titled / journal-code English papers** → **rename** to:
   ```
   <year>_<FirstAuthorSurname>_<ShortTitle>.pdf
   ```
   - `year`: 4-digit publication year.
   - `FirstAuthorSurname`: first author's surname, ASCII only, no spaces/accents (`Fröhlich`→`Frohlich`, `d'Humières`→`DHumieres`).
   - `ShortTitle`: 2–4 word CamelCase keyword from the title, ASCII letters/digits only (e.g. `WeaklyCompressibleLBM`, `RotatedFluxSolver`).
   - Only rename when you can **confidently** read all three from the PDF. If not confident, keep the name and report it.
4. **Code / data / image / note files** (`.py .dat .txt .f90 .png .ppt .md .tex .npz` …) → keep the body; they only participate in renumbering.
5. **Sequence numbering**: within each folder, prefix every entry with `1.`, `2.`, `3.`, … in the folder's **existing order** (preserve the user's curation; only fix/add numbers). Renamed papers keep their position.

## Excluded from renaming AND renumbering (do not touch internals)

These are interdependent code projects or an already-formatted review tree — renumbering would break file/path references:

- `Grid Gerneration/7.J Frohlich`
- `Grid Gerneration/8.NASA_Use_of_Poisson's_Equation_(GRAPE)`
- `Grid Gerneration/9.J_Frohlich_code`
- `PeriodicHill/WENO5-Interpolation/Total_literature_review`

Also never touch: `.git`, `.vscode`, `__pycache__`, `*.pyc`, `~$*`, `.DS_Store`, `.gitignore`.

## Procedure

1. **Inventory.** Walk the repo (skip the exclusions above). For each file record: folder, name, leading number, body, extension, whether it has CJK, whether it has a `[` tag. Classify deterministically: CJK→chinese, `[tag]`→curated, non-pdf→data, else→`english_pdf_candidate`.

2. **Read candidates.** For each `english_pdf_candidate`, gather evidence:
   - First try `pdfinfo "<abs>"` and `pdftotext -f 1 -l 3 -q "<abs>" -` (Poppler is installed).
   - If the text layer is empty/scanned, use the **Read tool** on the PDF with `pages: "1-3"` to view the rendered title page.
   - Decide `keep` vs `rename` per the rules; extract `year / author / shortTitle` when renaming.
   - **Parallelize** this: prefer a `Workflow` with one reader agent per candidate (structured output `{relpath, action, year, author, shortTitle, reason}`), or batched `Agent` calls. There are typically ~80 candidates.

3. **Build the plan.** Per folder, number files 1..N in existing order; set each body per the rules; compute `old → new`. Check for collisions. Write a `_plan.json` and print a human-readable table of only the **changed** rows, grouped by folder, plus any warnings (low confidence, numbering remnants).

4. **Show the plan to the user.** If `dry-run`, stop here.

5. **Execute.** Rename changed files with `os.rename` (works for tracked and untracked; git detects renames by content). Then:
   ```
   git add -A
   git commit -m "整理論文庫：英文論文改為 年_作者_標題、各資料夾重新編號"
   git push
   git fetch origin
   ```
6. **Report.** Summarize: # renamed, # renumbered, # kept, any files that need manual attention.

### Helper scripts (this repo)

Reusable Python helpers may exist at the repo root: `_inventory.py` (builds `_inventory.json` + extracts PDF text) and `_apply_rename.py` (`plan` / `apply` modes consuming `_inventory.json` + `_decisions.json`). Reuse/regenerate them as needed. These `_*.py` / `_*.json` scratch files are git-ignored and should not be committed.
