rm *.aux
rm *.log
rm *.out

# Список исходников
SOURCES=(
    some_colcars_description_and_study.tex
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    pdflatex "$src"
done