all: preview

docs: red.mat
	quarto render --output-dir ../docs

preview: red.mat
	quarto preview --output-dir ../docs

red.mat:
	curl -L -o red.mat https://zenodo.org/records/21287414/files/red_1.mat

clean:
	rm -f red.mat .gitignore *.html
	rm -rf _site .quarto _freeze
