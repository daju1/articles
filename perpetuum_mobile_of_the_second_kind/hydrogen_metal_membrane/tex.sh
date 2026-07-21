rm *.aux
rm *.log
rm *.out

# Список исходников
SOURCES=(
    nickel_hydrogen_ring.tex
    nickel_hydrogen_ring_v02.tex
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    pdflatex "$src"
done
