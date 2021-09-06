# Run this script to convert LaTeX files into ODT and DOCX documents
# Requires pandoc and pandoc-crossref being installed

pandoc main.tex --filter=pandoc-crossref --bibliography=library.bib --csl=journal-of-evolutionary-biology.csl -o scherrer2021.odt
pandoc main.tex --filter=pandoc-crossref --bibliography=library.bib --csl=journal-of-evolutionary-biology.csl -o scherrer2021.docx
