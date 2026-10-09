# MacroecologyMarine

Workspace for marine macroecology analyses.

## Structure

- `notebooks/figures.ipynb`: Python notebook for publication-style figures.
- `notebooks/SimulationsCode.ipynb`: Python notebook for logistic model simulations.
- `notebooks/Marechiara_timeseries.ipynb`: Python notebook for MareChiara time series (LTER-MC) analysis.
- `notebooks/SAD_GRUMP.ipynb`: Python code for SAD from GRUMP dataset.
- `notebooks/Taylor_GRUMP.ipynb`: Python notebook for spatial patterns analysis from GRUMP dataset.
- `notebooks/patchiness_TARA.ipynb`: Python notebook for spatial patterns analysis from Tara dataset.
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

## Simulations workflow
Open 'notebooks/SimulationsCode.ipynb '. From there, you can either run a new set of simulations using the same parameters from the main text and Supplementary Material (SM), or load the pre-computed simulation outputs directly from the 'notebooks/Output' folder.

## Data analysis workflow
You can run the analyses for the MareChiara and GRUMP datasets from scratch by following their corresponding notebooks:

For the Tara dataset, since the analysis in 'patchiness_TARA.ipynb' can take a few minutes to complete, you have the option to bypass the computation and directly load the pre-processed dataset ('patch.feather').

Please note that the raw databases must be downloaded separately.

## GitHub

This directory is intended to be a local Git repository on branch `main`.
After creating the empty repository on GitHub, connect it with:

```bash
git remote add origin git@github.com:<user>/MacroecologyMarine.git
git push -u origin main
```
