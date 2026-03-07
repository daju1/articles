rm *.aux
rm *.log
rm *.out

# Список исходников
SOURCES=(
    # On_the_applicability_of_Lorentz_calibration_to_Nikolaev-s_electrodynamics.tex
    01.on_contradiction_arising_when_deriving_vector_potential_formula_in_Coulomb_calibration.tex
    02.on_classical_electrodynamics_contradiction_arising_when_deriving_wave_equation_Juggling_calibrations.tex
    03.do_lagging_potentials_equations_correspond_to_Lorentz_calibration.tex
    04.do_Lienard-Wichert_potentials_correspond_to_Lorentz_calibration.tex
    05.on_convective_derivative_of_electric_field_vector.tex
    06.fourth_Maxwell_equation.tex
    07.stationary_point_charge_Coulomb_field_divergence.tex
    08.is_everything_alright_with_Gauss_theorem.tex
    09.on_Lorentz_calibration_applicability_to_Nikolaev_electrodynamics.tex
    
    longitudinal_magnetic_permeability_and_chronometric_Fermi_field.tex
    longitudinal_magnetic_permeability_and_chronometric_Fermi_field_en.tex
    longitudinal_magnetic_permeability_and_chronometric_Fermi_field_uk.tex

    classical_electrodynamics_contradictions_Fermi_variational_principle_and_Nikolaev_electrodynamics.tex 
)

for src in "${SOURCES[@]}"; do
    xelatex "$src"
    pdflatex "$src"
done