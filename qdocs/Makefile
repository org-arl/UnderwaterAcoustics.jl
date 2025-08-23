all: preview

docs: red.mat
	quarto render --output-dir ../docs

preview: red.mat
	quarto preview --output-dir ../docs

red.mat:
	cp ~/Data/ChannelModeling/UACR_Public/red/red.mat .

clean:
	rm -f red.mat .gitignore *.html
	rm -rf _site .quarto _freeze
