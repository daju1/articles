# brew install --cask basictex
# sudo tlmgr update --self
# sudo tlmgr install collection-langcyrillic
# sudo tlmgr install tcolorbox
# sudo tlmgr install pdfcol
# sudo tlmgr install collection-latexextra
# sudo tlmgr install palatino
# sudo tlmgr install newpx pxfonts
# sudo tlmgr install mathpazo
# sudo updmap-sys




rm *.aux
rm *.log
rm *.out


# Список исходников
SOURCES=(
    chrono_casimir_en.tex
    chrono_casimir_ru.tex
    chrono_casimir_grebennikov_ru.tex
    chrono_casimir_grebennikov_en.tex
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    pdflatex "$src"
done



rm *.aux
rm *.log
rm *.out