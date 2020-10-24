"Convolve and resample"
function conv_spectra(m::FixedKernelInstrument, ν, spectrum)
    s = imfilter(spectrum, m.kernel)
    interp_cubic = CubicSplineInterpolation(ν, s)
    return interp_cubic(m.ν_out)
end;

function conv_spectra(m::VariableKernelInstrument5, ν, spectrum)
    FT = eltype(m.ν_out)
    # Define grid where to perform convolution:
    stride = 5 # Fixed now, can be changed later
    off = ceil(Int, size(m.kernel, 1) / 2)
    ind = off:stride:(length(ν) - off)

    knots = view(ν, ind)
    te = LinearInterpolation(m.ν_out, Float32.(m.ind_out))
    spec_out = zeros(FT, length(knots));
    for i in eachindex(knots)
        # Simple first, nearest neighbor ILS
        ind_fraction = round(Int, te(knots[i]));
        kernel = view(m.kernel, :, ind_fraction)
        for j in eachindex(kernel)
            spec_out[i] += kernel[j] * spectrum[ind[i] + j] 
        end
    end
    fin = Interpolations.LinearInterpolation(ν[ind], spec_out; extrapolation_bc=Interpolations.Flat())
    return fin(m.ν_out)
end;