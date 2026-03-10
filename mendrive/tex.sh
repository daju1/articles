rm *.aux
rm *.log
rm *.out

# Список исходников
SOURCES=(
    mendrive.tex 
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    pdflatex "$src"
done