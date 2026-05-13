.PHONY: all build clean view help

# Main target name (without extension)
MAIN := main

# Output PDF
PDF := $(MAIN).pdf

# LaTeX compiler options
LATEXMK := latexmk
LATEXMK_FLAGS := -pdf -interaction=nonstopmode

help:
	@echo "Makefile for LaTeX compilation"
	@echo "Available targets:"
	@echo "  make build     - Compile the LaTeX document"
	@echo "  make view      - Compile and open the PDF"
	@echo "  make clean     - Remove build artifacts"
	@echo "  make help      - Show this help message"

all: build

build: $(PDF)

$(PDF): $(MAIN).tex preamble.sty references.bib
	$(LATEXMK) $(LATEXMK_FLAGS) $(MAIN).tex

view: $(PDF)
	open $(PDF)

clean:
	$(LATEXMK) -c $(MAIN).tex
	rm -f $(MAIN).pdf
	rm -f *.loe

.SILENT: help
