using DataInterpolations
using NCDatasets

ils_Δ, ils_in, dispersion = read_ils_table("/net/fluo/data2/oco2_L1bScND_18688a_180105_B8100r_180206190633.h5", "src/InstrumentOperators/json/ils_oco2.json");
band = 1
footprint = 1
# Footprint and band index:
extended_dims  = [footprint,band]
dispPoly = Polynomial(view(dispersion, :, extended_dims...))
ν = dispPoly.(0:1015)
res = 0.001
grid_x = Float32(-0.45 / 1e3):Float32(res / 1e3):Float32(0.45 / 1e3)
ind_out = collect(1:1016);
oco2_kernel = VariableKernelInstrument5(ils_pixel, wl, ind_out)

ν_mod = collect(Float32, 758:res:771)
spectrum = zeros(Float32, length(ν_mod))

y_conv =  conv_spectra(oco2_kernel, ν_mod, spectrum)