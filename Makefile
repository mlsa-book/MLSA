# Copied from mlr-org/mlr3book
all: install serve

.PHONY : help
help :
	@echo "install : Install book and dependencies."
	@echo "bookinstall : Install book without dependencies."
	@echo "serve   : Start a http server to serve the book."
	@echo "serverefresh   : Clear cache and start a http server to serve the book."
	@echo "render     : Render all formats."
	@echo "pdf     : Render book as pdf."
	@echo "html    : Render book as html."
	@echo "clean   : Remove auto-generated files."
	@echo "bibtex  : Reformats the bibtex file."

install:
	Rscript -e 'if (length(find.package("devtools", quiet = TRUE)) == 0) install.packages("devtools")' \
	        -e 'devtools::install_dev_deps(upgrade = "always")' \
			-e 'devtools::update_packages(upgrade = "always")' \
	        -e 'devtools::document()' \
			-e 'devtools::install()'

bookinstall:
	Rscript -e 'devtools::document()' \
			-e 'devtools::install()'

serve:
	Rscript -e 'renv::restore("book/", prompt = FALSE)'
	quarto preview book/

serveref:
	Rscript -e 'renv::restore("book/", prompt = FALSE)'
	quarto preview book/ --cache-refresh

clean:
	$(RM) -r book/_book book/.quarto book/site_libs;\
	find . -name "*.ps" -type f -delete;
	find . -name "*.dvi" -type f -delete;
	find . -type d -name "*_files" -exec rm -rf {} \;
	find . -type d -name "*_cache" -exec rm -rf {} \;

render:
	quarto render book/

html:
	quarto render book/ --to html

pdf:
	Rscript -e 'renv::restore("book/", prompt = FALSE)'
	quarto render book/ --to pdf

pdfref:
	Rscript -e 'renv::restore("book/", prompt = FALSE)'
	quarto render book/ --to pdf --cache-refresh


bibtex:
	biber --tool --output-align --output-indent=2 --output-fieldcase=lower book/book.bib -O book/book.bib
	rm book/book.bib.blg
