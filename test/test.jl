using JuMP, Ipopt
Ni=2; Nck=[10,5];
τ = [zeros(Float64,Nck[int],) for int in 1:Ni];
τ[1] = [-1.0,-0.927484,-0.763842,-0.525646,-0.236234,0.0760592,0.380665,0.647767,0.851225,0.971175];
τ[2] = [-1.0,-0.72048,-0.167181,0.446314,0.885792];

mdl = Model(solver = IpoptSolver());

@variable(mdl, tf_var)

@NLexpression(mdl, ts_JuMP[int=1:Ni,j=1:Nck[int]],  tf_var/2*τ[int][j] + tf_var/2)


x = [Array(Any,Nck[int],) for int in 1:Ni];
c_temp = [Array(Any,Nck[int],) for int in 1:Ni];

for int in 1:Ni
  for j in 1:Nck[int]
    x[int][j] = ts_JuMP[int,j];
  end
  N = length(x[int])
  mat = repmat(x[int],1,N);
  c_temp[int] = @NLexpression(mdl,[j=1:N],prod(mat[i,j] for i in 1:N))
end


#=
@variable(mdl, tf_var)
@unpack Ni, Nck, τ, ω = ps
const Ni_const = Ni  #TODO either way these can all be constaints ! --> simplify
const Nck_const = Nck
const τ_const = τ
const ω_const = ω
# calculate ts and ωₛ based off of tf
ts_JuMP, ωₛ_JuMP = create_intervals_JuMP(mdl,tf_var,Nck_const,Ni_const,τ_const,ω_const);

Nck=Nck_const; Ni=Ni_const;
D = [Matrix((Nck[int]+1),(Nck[int]+1)) for int in 1:Ni];
c = [Array(Any,Nck[int],) for int in 1:Ni];
DX = [Array(Any,Nck[int]+1,Nck[int]+1) for int in 1:Ni];

int=1;
x=ts_JuMP[int]

ss = [Symbol("Expr_D$(int)") for int in 1:Ni];
D_dict = Dict(); for int in 1:Ni; D_dict[ss[int]] = ss[int]; end


typeof(:($(ss1[1])) )


ss = [Int64(int) for int in 1:Ni];
D_dict = Dict();

for int in 1:Ni;
  D_dict[ss[int]] = [ss[int]];
end
@variable(m,D_expr[i=keys(D_dict)])

getvalue(y)

d = Dict();
d["foo"] = [:bar];
d["foo2"] = [:bar2];
m = Model();
x=@variable(m, [i = keys(d), d[i]]);
=#

A = [1 2 3;3 5 6; 9 8 7]
B = cumsum(A)
N = 3
C = zeros(N,N)
for i in 1:N
  for j in 1:N
    if i == 1
      C[i,j] = A[i,j]
    else
      C[i,j] = A[i,j]+C[i-1,j]
    end
  end
end
