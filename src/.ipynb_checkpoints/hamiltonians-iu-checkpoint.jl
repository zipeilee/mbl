using Yao, BitBasis, Random, LinearAlgebra, Distributions


# hamiltonians

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

	Wi = rand(-W:0.00001:W, nbit)

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

