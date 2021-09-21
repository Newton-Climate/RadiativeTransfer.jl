## 
## Using RadiativeTransfer.Absorption to perform absorption cross-sections
## 

using Revise
using BenchmarkTools
using RadiativeTransfer
using RadiativeTransfer.Absorption

## 
## STEP 1: Get the Hitran data from the transition par file
## 

hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)

## 
## STEP 2: Create a model from parameters
## 

modelCPU = make_hitran_model(hitran_data, Voigt(), wing_cutoff = 40, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=CPU())
modelGPU = make_hitran_model(hitran_data, Voigt(), wing_cutoff = 40, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=GPU())

grid = 6000:0.01:6400
# pressures = 250:250:1250
# temperatures = 100:75:400

# model = make_interpolation_model(hitran_data, Voigt(), ν_grid, pressures, temperatures, wing_cutoff = 40, CEF=ErfcHumliErrorFunctionVoigt())

## 
## STEP 3: Calculate the absorption cross section
## 

absorption_cross_section(modelCPU, (1e7/6400):0.01:(1e7/6000), 1000.1, 296.1, wavelength_flag=true)

value, derivs = absorption_cross_section(modelCPU, 6000:0.01:6400, 1000.1, 296.1, autodiff=true);
@benchmark absorption_cross_section(modelGPU, 6000:0.001:6400, 1000.1, 296.1)


hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)

# Create the model with parameters
modelCPU = make_hitran_model(hitran_data, Voigt(), wing_cutoff = 40, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=CPU())

# Compute the cross-section with autodifferentiation
value, derivs = absorption_cross_section(modelCPU, 6000:0.01:6400, 1000.1, 296.1, true, wavelength_flag=false);

absorption_cross_section(modelCPU, grid, 1000.1, 296.1) ≈ reverse(absorption_cross_section(modelCPU, reverse(1e7 ./ collect(grid)), 1000.1, 296.1, wavelength_flag=true))