# brew install --cask basictex
# sudo tlmgr update --self
# sudo tlmgr install collection-langcyrillic
# sudo tlmgr install tcolorbox
# sudo tlmgr install pdfcol
# sudo tlmgr install collection-latexextra
# sudo tlmgr install palatino
# sudo tlmgr install newpx pxfonts
# sudo tlmgr install mathpazo

# sudo tlmgr install cm-unicode
# sudo tlmgr install physics

# sudo updmap-sys




rm *.aux
rm *.log
rm *.out


# Список исходников
SOURCES=(
    chat-gravity-1.tex
    chat-gravity-2.tex
    sono.tex
    gravi.tex
    gravity.tex
    searl.tex
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    bibtex "$src"
    xelatex "$src"
    pdflatex "$src"
done



rm *.aux
rm *.log
rm *.out