
function truncate_phase(mod::δBGE, α, β, γ, δ, ϵ, ζ, C_sca, C_ext)
    @unpack l_max, Δ_angle =  mod
    # Obtain Gauss-Legendre quadrature points and weights for phase function:
    μ, w_μ = gausslegendre( length(β) );
    # Reconstruct phase matrix elements:
    f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄, P, P² = reconstruct_phase(α, β, γ, δ, ϵ, ζ, μ; returnLeg=true)
    # Find elements that exclude the peak (if wanted!)
    iμ = findall(x -> x<cosd(Δ_angle), μ)
    # Prefactor for P2:
    fac = zeros(l_max);
    for l=2:l_max-1
        fac[l+1] = sqrt(1 / ( ( l-1) * l * (l+1) * (l+2) ));
    end

    # Create subsets (for Ax=y weighted least-squares fits):
    y₁₁ = view(f₁₁, iμ)
    y₁₂ = view(f₁₂, iμ)
    y₃₄ = view(f₃₄, iμ)
    A   = view(P ,iμ,1:l_max)
    B   = fac' .* view(P²,iμ,1:l_max)

    # Weights (also avoid division by 0)
    minY = zeros(length(iμ)) .+ 1e-8;
    W₁₁ = Diagonal( w_μ[iμ] ./ max(abs.(y₁₁), minY) );
    W₁₂ = Diagonal( w_μ[iμ] ./ max(abs.(y₁₂), minY) );
    W₃₄ = Diagonal( w_μ[iμ] ./ max(abs.(y₃₄), minY) );
    
    # Julia backslash operator for least squares (just like Matlab);
    cl = ((W₁₁*A) \ (W₁₁*y₁₁))   # B in δ-BGR (β)
    γᵗ = ((W₁₂*B) \ (W₁₂*y₁₂))   # G in δ-BGE (γ)
    ϵᵗ = ((W₃₄*B) \ (W₃₄*y₃₄))   # E in δ-BGE (ϵ)
    # Integrate truncated function for later renormalization (here: fraction that IS still scattered):
    c₀ = ( w_μ' * (P[:,1:l_max] * cl) ) / 2

    # Compute truncated greek coefficients:
    βᵗ = cl / c₀                                    # Eq. 38a, B in δ-BGR (β)
    δᵗ = (δ[1:l_max] .- (β[1:l_max] .- cl)) / c₀    # Eq. 38b, derived from β
    αᵗ = (α[1:l_max] .- (β[1:l_max] .- cl)) / c₀    # Eq. 38c, derived from β
    ζᵗ = (ζ[1:l_max] .- (β[1:l_max] .- cl)) / c₀    # Eq. 38d, derived from β

    # To Do: Adjust ω̃ and extinction
    return αᵗ, βᵗ, γᵗ, δᵗ, ϵᵗ, ζᵗ, c₀
end

