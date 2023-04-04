

using Yao, BitBasis, Random, LinearAlgebra, Distributions


# export haimltonian, U Matrix, get_final_state, get_final_state_list, inf_otoc

include("hamiltonians.jl")

"""
    UMatrix(t)
Time evlution of haimltonian
"""

function UMatrix(t, hamiltonian)
    exp(Matrix(Yao.mat(hamiltonian)) * -im * t)
end

function UdagMatrix(t, hamiltonian)
    exp(Matrix(Yao.mat(hamiltonian)) * im * t)
end

function get_final_state(t, hamiltonian;initstate)
    finalstate = UMatrix(t, hamiltonian) * initstate
    return finalstate
end


function  get_final_state_list(t, hamiltonian,nbit::Int64)
    reglist = collect(basis(rand_state(nbit)))
    statelist = [state(arrayreg(r)) for r in reglist]
    [U(t, hamiltonian) * statelist[i] for i in 1:2^nbit]  
end

# get_final_reg
get_final_reg = finalstate -> arrayreg(finalstate)



"""
	inf_otoc(op1, op2;t,hamiltonian, nbit::Int64,i::Int64,j::Int64)
Calculate infinite temperature for Location Pauli Matrix X, Y, Z
	
		F(t) = < W(t)dagger V(0)dagger W(t) V(0)>_beta=0
	
	for one configured hamiltonina
"""
function inf_otoc(op1, op2;t,hamiltonian, nbit::Int64,i::Int64,j::Int64)
	mat1 = Matrix(mat(put(nbit, i=>op1)))
	mat2 = UdagMatrix(t,hamiltonian) * Matrix(mat(put(nbit, j=>op2))) * UMatrix(t, hamiltonian)
	
	reglist = collect(basis(rand_state(nbit)))
	statelist                                                                                                  = [state(arrayreg(r)) for r in reglist]

	var = [dot(state,mat2'*mat1'*mat2*mat1*state) for state in statelist]
	sum(var) / (2^nbit)

end


