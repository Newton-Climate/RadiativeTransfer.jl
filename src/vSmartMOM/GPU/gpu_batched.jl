using ForwardDiff
# batch Matrix inversion for GPU
# function batch_solve!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
#    temp = similar(A)
#    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
#    getri_strided_batched(A, temp, pivot) # inv stored in aux1
#    X = temp ⊠ B  
#    # synchronize()
# end

function batch_solve!(X::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    temp = similar(A)
    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
    synchronize()
    getri_strided_batched!(A, temp, pivot); # inv stored in aux1
    NNlib.batched_mul!(X, temp, B)
    synchronize()
end

# batch Matrix inversion for CPU
function batch_solve!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}, B::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views ldiv!(X[:,:,i], qr!(A[:,:,i]), B[:,:,i])
    end
end

function batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    pivot, info   = CUBLAS.getrf_strided_batched!(A, true);
    synchronize()
    getri_strided_batched!(A, X, pivot); # inv stored in aux1
    synchronize()
    return nothing
end

# Overload for Dual numbers
function batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT<:ForwardDiff.Dual}
    Atemp = ForwardDiff.value.(A)
    Xtemp = similar(Atemp);
    pivot, info   = CUBLAS.getrf_strided_batched!(Atemp, true);
    synchronize()
    getri_strided_batched!(Atemp, Xtemp, pivot); # inv stored in aux1
    synchronize()
    @show size(ForwardDiff.partials.(A))
    ForwardDiff.value.(X) .= Xtemp;
    for i=1:2
         -Xtemp ⊠ getindex.(ForwardDiff.partials.(A), i) ⊠ Xtemp;
         #getindex.(ForwardDiff.partials.(A), 1)
    end
    # need to pack into Dual number again somehow

    return nothing
end


# batch Matrix inversion for CPU
function batch_inv!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views X[:,:,i] = A[:,:,i]\I;
    end
end
