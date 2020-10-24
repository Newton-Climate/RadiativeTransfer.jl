
"""
$(FUNCTIONNAME)(grid_x::AbstractRange, ils_in::Array{FT}, ils_Δ::Array{FT}, extended_dims::Array{Int}=[]) where {FT <: AbstractFloat}

Pre-compute the ILS table input as function of spectral distance from center converted to the modeling grid 
Input: grid_x, ils_in, ils_Δ, extended_dims
Output: Offset-Array with tabulated responses interpolated to grid_x
"""
function prepare_ils_table(grid_x::AbstractRange, ils_in::Array{FT}, ils_Δ::Array{FT}, extended_dims::Array{Int}=[]) where {FT <: AbstractFloat}
    @assert minimum(abs.(grid_x)) < eps(FT) "grid_x must include 0 (center pixel)"
    ind_0 = argmin(abs.(grid_x))
    axis_pixel = (-ind_0 + 1):(grid_x.len - ind_0)
    # Dimension of ILS table (at least 1D, first dimension needs to be across wavenumber/wavelength)
    dims  = size(ils_in);
    # number of spectral positions of ILS table
    n_x = dims[1]
    # Number of ILS per detector position (can be 1 if ILS is constant across detector grid)
    n_pos = dims[2]
    ils    = view(ils_in, :, :,  extended_dims...);
    ils_Δ_ = view(ils_Δ,  :, :,  extended_dims...);
    ils_pixel = zeros(FT, grid_x.len, n_pos);

    for i = 1:n_pos
        ind = findall(minimum(ils_Δ_[:,i]) .< grid_x .< maximum(ils_Δ_[:,i]));
        interp = Interpolations.LinearInterpolation(ils_Δ_[:,i], ils[:,i])
        ils_pixel[ind,i] = interp.(grid_x[ind]);
    end
    # normalize here
    return OffsetArray(ils_pixel ./ sum(ils_pixel, dims=1), axis_pixel, 1:n_pos)
end
# This can still derive min/max range automatically, need to double check
