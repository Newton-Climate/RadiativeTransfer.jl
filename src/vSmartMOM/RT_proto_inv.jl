using Revise
using Plots
using Test
using RadiativeTransfer
using RadiativeTransfer.vSmartMOM
using ForwardDiff

n = 10;
nSpec = 5000;
RM = rand(n,n,nSpec)
RMgpu = CuArray(RM);
RMcpu = copy(RM);
x = [1.0,2.0,3.0,4.0,5.0]
function testInv(x, RM=RM)
    temp = x[1]* RM .+ x[2] .+ x[3]* RM .+ x[4]* RM .- x[5]* RM
    @show typeof(x[1])
    B = similar(temp);
    #@show ForwardDiff.npartials(ForwardDiff.partials.(temp))
    #@show ForwardDiff.partials.(temp,1)[1:10,1:10,1]
    @time vSmartMOM.batch_inv!(B,temp)
    #@show typeof(B)
    return B
end

RM = RMcpu;
dfdx_cpu = ForwardDiff.jacobian(testInv, x);
F_cpu = testInv(x);
RM = RMgpu;
#x = [1.0,2.0]
F_gpu = testInv(x);
dfdx_gpu = ForwardDiff.jacobian(testInv, x);

@test F_cpu ≈ Array(F_gpu)
@test dfdx_cpu ≈ Array(dfdx_gpu)
