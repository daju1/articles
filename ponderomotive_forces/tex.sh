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