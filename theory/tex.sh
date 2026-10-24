rm *.aux
rm *.log
rm *.out
rm *.toc


# Список исходников
SOURCES=(
    # Lazar_Braginsky_Drozdov_full_theory.tex
    # Drozdov_theory_full_exposition_v2.tex
    # Drozdov_theory_full_exposition_v3.tex
    # Drozdov_three_variants_calculation.tex
    # Drozdov_unified_model_article.tex
    Drozdov_unified_theory_detailed_ru.tex
    Drozdov_unified_theory_intro_ru.tex
    # insert_burgers_scale_v3.tex
    # insert_planck_sector_v3.tex
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
rm *.toc