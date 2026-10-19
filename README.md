# Hydraulics

An online textbook for the undergraduate Hydraulics course, Department of Civil Engineering, Jeonbuk National University.

**Read it at <https://daehakim-waterlab.github.io/hydraulics/>** · [PDF](https://daehakim-waterlab.github.io/hydraulics/hydraulics.pdf)

Built with [Jupyter Book](https://jupyterbook.org/) 1.x.

## Layout

- `_config.yml`: book settings (MathJax, code execution, PDF, footer "© Daeha Kim" and license).
- `_toc.yml`: chapter order. Chapter, section, figure and table numbers are written by hand in the chapter files.
- `intro.md`: the front page.
- `chapters/chNN_*.md`: one MyST Markdown notebook per chapter. Python code goes in `{code-cell} python` blocks and runs on every build, so figures and tables always come from the code. Tag a cell `hide-input` to fold its code under the figure, or `remove-input` to show only its output (e.g. a table).
- `_static/custom.css`: table styling for phones and dark mode.

## Build locally

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
jupyter-book build .
open _build/html/index.html
```

The PDF needs XeLaTeX (TeX Live); `jupyter-book build . --builder pdflatex` writes `_build/latex/hydraulics.pdf`. CI builds it for you.

## Publishing

`.github/workflows/deploy.yml` builds the HTML and the PDF on every push to `main` and deploys them to GitHub Pages. A code cell that raises an error, or any build warning, stops the deployment, and the previous version stays online. A failed PDF build does not stop it.

## License

[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). See [LICENSE](LICENSE).
