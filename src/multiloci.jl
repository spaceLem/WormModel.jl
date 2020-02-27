"""
    get_alleles(W::Array{Int,2})

Get frequency of just the resistant alleles from the population `W`.
"""
get_alleles(W::AbstractArray{Int,2}) = g2a(vec(sum(W, dims=2)))[2:2:end]


"""
    prop_res(W::Array{Int,2}, gt_fitness_drug::Array{Float64,2})

Calculate the proportion of resistance in the Worm population `W`.
"""
function prop_res(W::AbstractArray{Int,2},
				  gt_fitness_drug::Array{Float64,2})

	n = size(gt_fitness_drug, 2)
	s = sum(W)
	if s > 0
		pr = [sum(W .* gt_fitness_drug[:, i]) / s for i in 1:n]
	else
		pr = zeros(n)
	end
	pr
end #fn


