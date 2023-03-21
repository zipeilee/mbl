

using Yao, BitBasis, Random, LinearAlgebra, Distributions


# export haimltonian, U Matrix, get_final_state, get_final_state_list, inf_otoc



"""
	hamiltonian(nbit::Int, periodic::Bool=true)

1D Heisenberg hamiltonian with Dzyaloshinkii-Moriya and random-field
h_c = 9.45

"""
function the_hamiltonian(nbit::Int, W, periodic::Bool=false)
	sx = i -> 0.5 * put(nbit,i=>X)
	sy = i -> 0.5 * put(nbit,i=>Y)
	sz = i -> 0.5 * put(nbit,i=>Z)
	
	hopping_term = map(1:(periodic ? nbit : nbit-1)) do i
		j = i%nbit + 1
		sx(i)*sx(j)+sy(i)*sy(j)+sz(i)*sz(j) + sy(i)*sz(j) - sz(i)*sy(j) -sx(i)*sz(j) + sz(i)*sx(j) + sx(i)*sy(j) - sy(i)*sx(j)
	end |> sum
	
	Wi = rand(-W:0.0001:W, nbit)

	disorder = map(1:nbit) do i
		Wi[i] * sx(i)
	end |> sum

	hopping_term + disorder
end


function phenomenological(nbit::Int, h, J, xi)
end



"""
	xxz(nbit::Int, J1, J2, W, periodic::Bool=false)

Random Field XXZ model 
Wc = 5
J2 = 0.2
Ref: Science Bulletin 62 (2017) 707-711
"""


function xxz(nbit::Int, J1, J2, W, periodic::Bool=false)
	sx = i -> 0.5 * put(nbit,i=>X)
	sy = i -> 0.5 * put(nbit,i=>Y)
	sz = i -> 0.5 * put(nbit,i=>Z)

	hopping_term = map(1:(periodic ? nbit : nbit-1)) do i
		j = i%nbit + 1
		J1 * (sx(i)*sx(j)+sy(i)*sy(j))+J2 * (sz(i)*sz(j)) 
	end |> sum

	Wi = rand(-W:0.1:W, nbit)

	disorder = map(1:nbit) do i
		Wi[i] * sx(i)
	end |> sum

	hopping_term + disorder
	
end




"""
	Ising Model in NMR system to test TOTC
	Ref:PhysRevX.7,031011(2017)
"""
function nmr_ising(nbit::Int, g, h, periodic::Bool=false)

	hopping = map(1:(periodic ? nbit : nbit-1)) do i
		repeat(nbit, Z, (i, i%nbit+1))
	end |> sum
	-hopping + g*sum(map(i->put(nbit, i=>X), 1:nbit)) + h*sum(map(i->put(nbit, i=>Z), 1:nbit))
end


"""
	1D-Heisenberg model with random field
	hc = 3.75
	Ref: 1411.0660

"""
function rand_heisenberg(nbit::Int64, h;periodic::Bool=false)
	sx = i -> 0.5 * put(nbit,i=>X)
	sy = i -> 0.5 * put(nbit,i=>Y)
	sz = i -> 0.5 * put(nbit,i=>Z)
	
	hopping = map(1:(periodic ? nbit : nbit-1)) do i
		j = i%nbit + 1
		(sx(i)*sx(j)+sy(i)*sy(j)+sz(i)*sz(j))
	end |> sum

	hi = rand(Uniform(-h,h),nbit)

	disoder = map(1:nbit) do i
		hi[i] * sx(i)
	end |> sum

	hopping + disoder
end


"""
    UMatrix(t)
Time evlution of haimltonian
"""

function UMatrix(t, hamiltonian)
    exp(Matrix(mat(hamiltonian)) * -im * t)
end

function UdagMatrix(t, hamiltonian)
    exp(Matrix(mat(hamiltonian)) * im * t)
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
	inv_otoc(op1, op2;t,hamiltonian, nbit::Int64,i::Int64,j::Int64)
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


