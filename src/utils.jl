# triangular numbers T(n) = n(n+1)/2
# uses division so ensure return Int
triangle(n::Int)::Int = div(n * (n + 1), 2)


# test for overshooting
# don't allocate temporaries
function overshoot(W::Array, dW::Array)::Bool
	for i in eachindex(W)
		if W[i] + dW[i] < 0
			return true
		end #if
	end #for

	false
end #fn


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
# take the array of alleles and return the number of genotypes
# single locus version
function get_num_GTs(X::Array{Float64, 1})::Int
	prod(triangle(length(X)))
end #fn


# multiple locus version
function get_num_GTs(X::Array{Array{Float64, 1}, 1})::Int
	prod((triangle(length(xs)) for xs in X))
end #fn
=#
