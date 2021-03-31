using Revise
using Plots

using RadiativeTransfer
using RadiativeTransfer.vSmartMOM
using ForwardDiff

n = 20;
nSpec = 100;
RM = rand(n,n,nSpec)
x = [1.0,2.0]
function testInv(x, RM=RM)
    temp = x[1]* RM .+ x[2]
    B = similar(temp);

    vSmartMOM.batch_inv!(B,temp)
    @show typeof(B)
    return B
end


dfdx = ForwardDiff.jacobian(testInv, x);


