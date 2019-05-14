using LinearAlgebra
using StaticArrays

include("utils.jl")

#= old code for n >= 2 alleles
function a2g1(Ad::Array{Float64,1})
	n = length(Ad)

	Gf = zeros(triangle(n))
	idx = 0
	for i = 1:n
		Gf[idx += 1] = Ad[i]^2 # i == j
		for j = i+1:n
			Gf[idx += 1] = 2 * Ad[i] * Ad[j] # i != j
		end #for
	end #for

	Gf
end #fn
=#

function unmixed_a2g(Ad::Array{Float64,1})
	num_loci = length(Ad) ÷ 2
	num_gts = 3^num_loci
	Gd = ones(num_gts)

	# for each locus in turn, find gt distribution
	# then multiply all gt indices that have that gt
	for i in 1:num_loci
		freqs = @MVector [Ad[2i-1], 0, Ad[2i]]

		for j in 1:num_gts
			idx = mod1(ceil(Int, j / 3^(i-1)), 3)
			Gd[j] *= freqs[idx]
		end #for
	end #for

	normalize!(Gd, 1)
end #fn

# turn allele distribution into genotypes distribution (1-dim)
function a2g1(Ad::Array{Float64, 1})
	@MVector [Ad[1]^2, 2 * Ad[1] * Ad[2], Ad[2]^2]
end #fn


# turn allele distribution into genotype distribution
function a2g(Ad::Array{Float64,1})
	num_loci = length(Ad) ÷ 2
	num_gts = 3^num_loci
	Gd = ones(num_gts)

	# for each locus in turn, find gt distribution
	# then multiply all gt indices that have that gt
	for i in 1:num_loci
		#freqs = a2g1(Ad[2i - 1:2i])
		freqs = @MVector [Ad[2i-1]^2, 2 * Ad[2i-1] * Ad[2i], Ad[2i]^2]

		for j in 1:num_gts
			idx = mod1(ceil(Int, j / 3^(i-1)), 3)
			Gd[j] *= freqs[idx]
		end #for
	end #for

	normalize!(Gd, 1)
end #fn


# get matrix for g2a
#= Something like A = [
	2 1 1 0 0 0
	0 1 0 2 1 0
	0 0 1 0 1 2]'
=#
function g2a_mat(x::Int)
	A = zeros(Int, triangle(x), x)
	row = 0
	for offset = 1:x
		A[row += 1, offset] = 2
		for col = offset + 1:x
			A[row += 1, offset] = 1
			A[row, col] = 1
		end #for
	end #for

	A
end #fn


# turn genotype frequency into alleles distribution
function g2a(Gf::Array{Int64,1})
	num_loci = round(Int, log(3, length(Gf)))
	num_gts = 3^num_loci

	Ad = zeros(2num_loci)
	if sum(Gf) == 0
		return Ad
	end #if

	A = [2 1 0; 0 1 2]

	for i in 1:num_loci
		rng = (-1:0) .+ 2i
		for j in 1:num_gts
			idx = mod1(ceil(Int, j / 3^(i - 1)), 3)
			Ad[rng] += A[:, idx] * Gf[j]
		end #for
		Ad[rng] /= sum(Ad[rng])
	end #for

	Ad
end #fn


#= something like
	A = [1,1,2]
	A = [1,1,1,2,2,3]
	A = [1,1,1,1,2,2,2,3,3,4]
=#
function gt_fitness1(num_alleles::Int)
	A = zeros(Int, triangle(num_alleles))
	idx = 0
	for i in 1:num_alleles, j in 1:(num_alleles + 1 - i)
		A[idx += 1] = i
	end #for
	A
end #fn


function calc_gt_fitness(allele_fitness::Array{Float64, 1}, dominance::Array{Float64, 1})
	num_loci = length(allele_fitness) ÷ 2
	num_gts = 3^num_loci

	gt_fitness = zeros(num_gts)

	for i in 1:num_loci, j in 1:num_gts
		idx = mod1(ceil(Int, j / 3^(i - 1)), 3)
		if idx == 1
			gt_fitness[j] += allele_fitness[2i - 1]
		elseif idx == 3
			gt_fitness[j] += allele_fitness[2i]
		else
			gt_fitness[j] += (1.0 - dominance[i]) * allele_fitness[2i - 1]
			gt_fitness[j] += dominance[i] * allele_fitness[2i]
		end #if
	end #for
	gt_fitness ./= num_loci

	gt_fitness
end #fn


function prop_res(W, gt_fitness_drug)
	s = sum(W)
	s > 0 ? sum(W .* gt_fitness_drug) / s : 0.0
end #fn

function prop_sus(W, gt_fitness_drug)
	s = sum(W)
	s > 0 ? 1 - sum(W .* gt_fitness_drug) / s : 0.0
end #fn



function ingest!(K_establish, K_ingest, num_treated)
	num_gts, num_sheep = size(K_establish)
	num_refugia = num_sheep - num_treated

	treated = 1:num_treated
	refugia = num_treated + 1:num_sheep

	for gt in 1:num_gts
		# treated = 1, refugia = 2
		if num_treated > 0
			K_establish[gt, treated] .= rand(Multinomial(K_ingest[gt, 1], num_treated))
		end #if

		if num_treated < num_sheep
			K_establish[gt, refugia] .= rand(Multinomial(K_ingest[gt, 2], num_refugia))
		end #if
	end #for
end #fn



#=
function print_gt_dist(gt_dist::Array{Float64,1}, alleles::Array{Int64,1})
	nA = length(alleles)
	alphabet = 'A':'Z'
	cpa = 1 ./ cumprod(alleles)
	for i in 1:length(gt_dist)
		for j in 1:nA
			x = mod1(ceil(Int, , alleles[i])
			print(alphabet[j], x)
		end #for
		println("\t= ", gt_dist[i])
	end #for
end #fn
=#

