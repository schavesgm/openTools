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
The instructions on how to use it can be achieved using,

```bash
gnp2tex -h
```

