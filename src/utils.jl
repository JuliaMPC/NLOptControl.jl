


# generates Gauss points and weights on the interval [0,1]
# this function comes from NovcolN.jl which cites the RosettaCode page
function gauss(n)
    λ, Q = eig(SymTridiagonal(zeros(n), [ i / sqrt(4i^2 - 1) for i = 1:n-1 ]))
    x = (λ + 1) * 1/2
    w = [ 2*Q[1,i]^2 for i = 1:n ]/2
    return x, w
end


# change the interval

function Define_Tau(t0, tf)
