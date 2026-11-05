all: preview

docs: red.mat red_noise.mat
	quarto render --output-dir ../docs

preview: red.mat red_noise.mat
	quarto preview --output-dir ../docs

red.mat:
	curl -L -o red.mat https://zenodo.org/records/21287414/files/red_1.mat

red_noise.mat:
	curl -L -o red_noise.mat https://zenodo.org/records/21287414/files/red_noise.mat

clean:
	rm -f red.mat red_noise.mat .gitignore *.html
	rm -rf _site .quarto _freeze
