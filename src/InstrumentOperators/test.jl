
oco = Dataset("/Users/cfranken/Downloads/oco2_L1bScND_18688a_180105_B8100r_180206190633.h5")
ils = oco.group["InstrumentHeader"]["ils_relative_response"][:]
ils_Δ = oco.group["InstrumentHeader"]["ils_delta_lambda"][:]

# Views:
band = 2
footprint = 4
ils2   = view(ils, :, :, footprint, band)
ils2_Δ = view(ils_Δ, :, :, footprint, band)
grid = -0.9 / 1e3:0.01 / 1e3:0.9 / 1e3;
ils_pixel = zeros(Float32, length(collect(grid)), 1016);

for i = 1:1016
    interp_linear = LinearInterpolation(ils2_Δ[:,i], ils2[:,i])
    ils_pixel[:,i] = interp_linear(grid);
end