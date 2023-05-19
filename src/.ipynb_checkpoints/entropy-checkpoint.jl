using LinearAlgebra
using SparseArrays
using OMEinsum

export my_ptrace

"""
 partical trace using einsum
 inspired by zipei
 create by ylxdzsw
"""
function my_ptrace(x, sys, dims)
    simplified_dims = prod(dims[sys+1:end]), dims[sys], prod(dims[1:sys-1])
    r = ein"abcdbf->acdf"(reshape(x, simplified_dims..., simplified_dims...))
    reshape(r, prod(dims) ÷ dims[sys], prod(dims) ÷ dims[sys])
end