all: preview

docs: red.mat
	quarto render --output-dir ../docs

preview: red.mat
	quarto preview

red.mat:
	cp ~/Data/ChannelModeling/UACR_Public/red/red.mat .

clean:
	rm -f red.mat .gitignore
	rm -rf _site .quarto
