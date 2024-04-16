fn=Casimir_honeycomb_drive

rm $fn.aux
rm $fn.dvi
rm $fn.out
rm $fn.log
rm $fn.bib
rm $fn.bbl

xelatex $fn.tex  # generates aux,dvi etc
cp $fn.aux $fn.tex.aux

bibtex $fn.tex # generates a bbl
cp $fn.tex.bbl $fn.bbl

#xelatex $fn.tex
pdflatex $fn.tex

