# Copied from mlr-org/mlr3book
all: serve

.PHONY : help
help :
	@echo "serve   : Clear cache and serve book on http server."
	@echo "render     : Clear cache and render all formats."
	@echo "pdf     : Clear cache and render pdf."
	@echo "clean   : Remove auto-generated files."
	@echo "bibtex  : Reformats the bibtex file."

serve:
	quarto preview book/

clean:
	$(RM) -r book/_book book/.quarto book/site_libs;\
	find . -name "*.ps" -type f -delete;
	find . -name "*.dvi" -type f -delete;
	find . -type d -name "*_files" -exec rm -rf {} \;
	find . -type d -name "*_cache" -exec rm -rf {} \;

render:
	quarto render book/ --cache-refresh

pdf:
	quarto render book/ --to pdf --cache-refresh

bibtex:
	biber --tool --output-align --output-indent=2 --output-fieldcase=lower  book/library.bib -O book/library.bib
	rm book/library.bib.blg
