# GNP2TEX

*gnp2tex* is a tool to convert gnuplot plots to LaTeX style. It works using
**luaTeX** and **pdflatex** to transform a script written in gnuplot to a
pdf standalone file.

## Usage:

In order to run **gnp2tex**, first copy your file into /usr/bin/ to have the
command available,

```bash
sudo cp gnp2tex /usr/bin
```
and grant it executable permissions,

```bash
sudo chmod +x /usr/bin/gnp2tex
```
If you wanted to use it, just run it using,

```bash
gnp2tex nameFile.gn namePDF.pdf
```
The first argument is the name of the file to process, the second is the name
of the pdf to be created. The second argument is optional, if it is not set,
then **gnp2tex** will create a pdf file named **nameFile.pdf**. You can try
the script on the file **trialPlot.gn** and it should produced a plot containing
two harmonic functions.

