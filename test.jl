using ForwardDiff
using Plots
# example, interpolate f(x) = x^2 over 1<=x<=3 given
x_data = [1.0,2.0,3.0]
y_data = [1.0,4.0,9.0];
x0 = 1; xf = 3;
ns = 100;  # plotting points
x = Array(linspace(x0,xf,ns));
y = x.^2;

#MAX_CHUNK_SIZE = 10
#XLEN = MAX_CHUNK_SIZE
#YLEN = div(MAX_CHUNK_SIZE, 2) + 1
#X, Y = rand(XLEN), rand(YLEN)

function chebyquad!(y::Vector, x::Vector)
    tk = 1/length(x)
    for j = 1:length(x)
        temp1 = 1.0
        temp2 = 2x[j]-1
        temp = 2temp2
        for i = 1:length(y)
            y[i] += temp2
            ti = temp*temp2 - temp1
            temp1 = temp2
            temp2 = ti
        end
    end
    iev = -1.0
    for k = 1:length(y)
    #    y[k] *= tk
        y[k] = y[k]*tk

        if iev > 0
            y[k] += 1/(k^2-1)
        end
        iev = -iev
    end
    return y
end

#chebyquad(x) = (y = zeros(eltype(x), length(Y)); chebyquad!(y, x); return y)

#yc = chebyquad!(y_data,x_data)

d = x -> ForwardDiff.derivative(chebyquad!, x);

plot(x,y,label= "actual polynomial")
plot!(x,d(x))
