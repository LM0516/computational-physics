;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("preamble" "")))
   (TeX-run-style-hooks
    "latex2e"
    "chapters/01_ErrorAnalysis"
    "chapters/02_LinearSystems"
    "chapters/03_Interpolations"
    "chapters/04_RootsOfNonlinearEquations"
    "chapters/05_NumericalIntegration"
    "chapters/06_OrdinaryDifferentialEquations"
    "article"
    "art10"
    "preamble")
   (LaTeX-add-bibliographies
    "references"))
 :latex)

