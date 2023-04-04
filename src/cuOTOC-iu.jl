using Yao, BitBasis, Random, LinearAlgebra, Distributions, CUDA

include("OTOC.jl")

function cu_inf_otoc(op1, op2;t,hamiltonian, nbit::Int64,i::Int64,j::Int64)
	mat1 = CUDA.zeros(ComplexF64,(2^nbit,2^nbit,2^nbit))
	mat1 .+= cu(Matrix(Yao.mat(put(nbit, i=>op1)))) 
	
	mat2 = CUDA.zeros(ComplexF64,(2^nbit,2^nbit,2^nbit))
	mat2 .+= cu(UdagMatrix(t,hamiltonian)) * cu(Matrix(Yao.mat(put(nbit, j=>op2)))) * cu(UMatrix(t, hamiltonian))
	
	reglist = collect(basis(rand_state(nbit)))
	statelist = zeros(2^nbit,2^nbit)
	st = state.(arrayreg.(reglist))
	for i in 2^nbit
		statelist[:,i] += st[i]
	end
	statelist = cu(statelist)
	
	var = dot.(statelist,adjoint.(mat2) .* adjonin.(mat1) .* mat2 .* mat1 .* statelist)
	sum(var) / (2^nbit)

end
