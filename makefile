all: README.html

README.html: README.Rmd
	- Rscript -e "rmarkdown::render('README.Rmd')"

clean:
	- rm -f README.html