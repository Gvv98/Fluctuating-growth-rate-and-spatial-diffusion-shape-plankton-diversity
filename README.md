# MacroecologyMarine

Workspace for marine macroecology analyses.

## Structure

- `notebooks/figures.ipynb`: Python notebook for publication-style figures.
- `figures/`: generated figures to use in reports and manuscripts.
- `data/`: project data; large raw/intermediate data are ignored by default.
- `reports/`: Markdown reports for results and current project state.
- `scripts/`: Python wrappers and pipeline scripts.

## Figure workflow

Open `notebooks/figures.ipynb`, load or build the relevant data, and call
`save_figure(fig, "figure_name")`. The helper saves both PNG and PDF files in
`figures/`.

The notebook defaults follow `AGENTS.md`: seaborn, colorblind-friendly colors,
large line widths, no top/right spines, and Avenir Next when available.

If needed, install the plotting environment with:

```bash
python3 -m pip install -r requirements.txt
```

## GitHub

This directory is intended to be a local Git repository on branch `main`.
After creating the empty repository on GitHub, connect it with:

```bash
git remote add origin git@github.com:<user>/MacroecologyMarine.git
git push -u origin main
```
