## General attitude

- always prefer conciseness over grammar correctedness
- I am a scientist, I want you to behave in the same way
- no principle of authority, I can be wrong and I also make mistakes

## Comparing theory and simulations or data

- when comparing theoretical models with numerical simulations, always check if there is a mismatch
- if there is a mismatch report it and try to adress whether is a bug in the code, is something potentially wrong in the math, an expected mismatch given assumption, or something potentially interesting to be

## Coding and style

- use julia for simulating models and computationally expensive tasks
- use python as wrapper, to read/modify files, launch multiple codes in a pipeline
- use python for analyzing data
- use python for plotting
- use .md files for drafting reports and latex with pdf for finalized ones
- use revtex for latex

## Reports

- write reports about the current state of the project
- keep a track of the high-level task and results obtained
- write .md files in reports folder. Each file should reflect one or more results (something we understood) not a list of operations / tasks. Include math!
- insert to figures in these files, with clear captions, and modify the files when something new is done
- always use the $ notation for formulas, not the parenthesis.

## Paper and supplementary materials

- Never modify files in `my_notes`.

## Plotting style

- font text of axis names, label, legend should be large and of the same font size
- use Avenir Next font for all text in figures, including legends and axis labels, unless a figure-specific reason is documented. If Avenir Next is not available, use a similar sans-serif font.
- use large width of lines
- remove unnecessary elements (right and upper axis, squares around legends, grid, etc)
- use colorblind-friendly colors
- use seaborn for plotting in python

## Notes

- I am going to write personal notes in the folder called 'my_notes'. Never modify them. If you find a mistake inside, please report in the reports.
