rm *.aux
rm *.log
rm *.out
rm *.nav
rm *.toc

# Список исходников
SOURCES=(
    mendrive.tex
    mendrive_presentation_v5.tex
    # mendrive_presentation_v6_ru.tex
    mendrive_presentation_v6_en.tex
    mendrive_presentation_v7.tex
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    pdflatex "$src"
done