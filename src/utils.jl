"""
    triangle(n::Int)

Get triangular numbers T(n) = n(n+1)/2
Uses division so ensure return Int.
"""
triangle(n::Int)::Int = n * (n + 1) ÷ 2


"""
    get_pt_idx(gt::Int, locus::Int)

Return the phenotype index (1-3) of the genotype given the locus. Used in several places.

"""
get_pt_idx(gt::Int, locus::Int) = mod1(ceil(Int, gt / 3^(locus - 1)), 3)



"""
    overshoot(W::Array, dW::Array)

Test for overshooting (i.e. `W + dW ≮ 0`) without allocating temporaries.
"""
function overshoot(W::Array, dW::Array)::Bool
	for i in eachindex(W)
		W[i] + dW[i] < 0 && return true
	end #for

	false
end #fn


"""
    ensure_positive!(W::Array, name::String)

Check if any entry in `W` in negative. If so, print an alert and set that entry to 0.
"""
function ensure_positive!(W::Array, name::String)
	for i in eachindex(W)
		if W[i] < 0
			println("$name[$i] = ", W[i])
			W[i] = 0
		end #if
	end #for

	nothing
end #fn

#=
"Take the array of alleles and return the number of genotypes (single locus version)."
function get_num_GTs(X::Array{Float64,1})::Int
	prod(triangle(length(X)))
end #fn


"Take the array of alleles and return the number of genotypes (multiple locus version)."
function get_num_GTs(X::Array{Array{Float64,1}, 1})::Int
	prod((triangle(length(xs)) for xs in X))
end #fn
=#
