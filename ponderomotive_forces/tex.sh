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
    Tamm_ponderomotive_forces.tex
    Tamm_ponderomotive_forces_update.tex
    Tamm_ponderomotive_forces_update-2.tex
    MenDrive_ponderomotive.tex
    MenDrive_ponderomotive-1.tex
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    pdflatex "$src"
done