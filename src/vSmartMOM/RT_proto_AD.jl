using Revise
using Plots
using Pkg
# Pkg.activate(".")
using RadiativeTransfer
using RadiativeTransfer.Absorption
using RadiativeTransfer.Scattering
using RadiativeTransfer.vSmartMOM

# Load parameters from file
parameters = vSmartMOM.parameters_from_yaml("src/vSmartMOM/ModelParameters/DefaultParameters.yaml")

# default_parameters

# Sets all the "specific" parameters
# parameters = vSmartMOM.default_parameters();
function runner(x, parameters=parameters)
    parameters.τAer_ref = [x[1]];
    @show parameters.p₀
    parameters.p₀ = [x[2]];
    @show parameters.p₀
    model = model_from_parameters(parameters);
    #model.params.architecture = RadiativeTransfer.Architectures.CPU()
    J = vSmartMOM.rt_run(model);
    @show J
    return J#; R_SFI[1,1,:]
end
# Generates all the derived attributes from above parameters
#model = model_from_parameters(parameters);
#model.params.
##model.params.architecture = RadiativeTransfer.Architectures.GPU()
#@time R_GPU, T_GPU, R_SFI, T_SFI = vSmartMOM.rt_run(model);
#@show R_GPU[1,1,:]
#@show R_SFI[1,1,:]
#plot(R_GPU[1,1,:])
#plot!(R_SFI[1,1,:])
#plot(R_GPU[1,1,:]./R_SFI[1,1,:])
# curr_parameters.

