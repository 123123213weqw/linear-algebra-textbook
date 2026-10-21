# Linear Algebra Textbook (线性代数)

A Chinese-language introductory linear algebra textbook (《线性代数》) written in LaTeX,
covering the standard six-chapter outline of a university course:

| # | Chapter | Topics |
|---|---------|--------|
| 1 | 行列式 | definition, properties, cofactor expansion, Cramer's rule |
| 2 | 矩阵及其运算 | operations, inverse, adjugate, block matrices, elementary operations, rank |
| 3 | 向量组的线性相关性 | combinations, dependence, maximal independent sets, vector spaces |
| 4 | 线性方程组 | solvability criteria, fundamental solution systems, structure of solutions |
| 5 | 特征值与特征向量 | inner products, orthogonality, diagonalization, real symmetric matrices |
| 6 | 二次型 | standard forms, the law of inertia, positive definiteness |

Appendix A gives short answers or hints for selected exercises.

## Read the book (auto-built PDF)

Every push to `main` rebuilds the PDF and publishes it to GitHub Pages:

- **Book page:** <https://123123213weqw.github.io/linear-algebra-textbook/>
- **PDF (direct):** <https://123123213weqw.github.io/linear-algebra-textbook/linear-algebra.pdf>

## Build locally

The project compiles with [tectonic](https://tectonic-typesetting.github.io/) (XeTeX-based,
downloads packages on demand; no TeX Live install needed):

```bash
tectonic main.tex   # produces main.pdf
```

Any full XeLaTeX or pdfLaTeX toolchain also works, e.g. `latexmk -pdf main.tex`.
The document uses only standard packages, so it builds identically on macOS,
Linux, and CI.

## Repository layout

```
main.tex                     # document root: class, theorem environments, boxes
chapters/preface.tex         # 前言
chapters/ch1-determinants.tex
chapters/ch2-matrices.tex
chapters/ch3-vectors.tex
chapters/ch4-linear-systems.tex
chapters/ch5-eigen.tex            # in progress
chapters/ch6-quadratic.tex        # in progress
chapters/appendix-answers.tex     # in progress
tools/verify_examples.py      # numeric verification of the worked examples
site/index.html               # GitHub Pages landing page (embeds the PDF)
.github/workflows/build-pdf.yml
```

## Reviewing and giving feedback

This is a draft (v0.1) open for review. Please open a GitHub Issue per
comment — chapter, section, and item number (e.g. "Ch. 1, §1.4, Example 3")
plus a short description of the problem (typo, math error, pedagogical
suggestion). Numeric claims in the worked examples are checked by
`tools/verify_examples.py` (pure Python, no dependencies).

## Status

Draft v0.1 — Chapters 1–4 (determinants, matrices, vectors, linear systems) are
complete with definitions, theorems and proofs, worked examples, figures, Q&A
boxes, and end-of-chapter exercises. Chapters 5–6 and the answers appendix are
in progress and will be pushed as they are written. Every push re-runs the
numeric verification (39 checks so far) before the PDF is published.
