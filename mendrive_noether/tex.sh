rm *.aux
rm *.log
rm *.out

# Список исходников
SOURCES=(
    mendrive_noether.tex
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    pdflatex "$src"
done
