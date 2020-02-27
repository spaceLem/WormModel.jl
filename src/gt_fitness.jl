"""
    calc_gt_fitness(P; fit::Symbol = :drug)

Calculate genotype fitness given parameters `P`, and choose `fit` is `:drug` or `:establish`. This only needs to be done once, at the start of the simulation.
"""
function calc_gt_fitness(P; fit::Symbol = :drug)

	allele_fitness = P[fit == :drug ? :fitness_drug : :fitness_establish]
	dominance = P[:dominance]
	num_loci = P[:num_loci]
	num_gts = 3^num_loci

	# array for phenotype [SS, SR, RR]
	pt = @MVector zeros(3)

	if P[:multidrug] && fit == :drug
		# need every drug combination
		gt_fitness = ones(num_gts, 2^num_loci)
	elseif fit == :drug
		# just one drug to worry about, so average across
		gt_fitness = zeros(num_gts, 1)
	else
		# establishment cost
		gt_fitness = ones(num_gts, 1)
	end #if

	for locus in 1:num_loci
		dom = dominance[locus]
		Sx, Rx = allele_fitness[2locus .+ (-1:0)]

		pt .= @MVector [Sx, (1.0 - dom) * Sx + dom * Rx, Rx]

		for gt in 1:num_gts
			idx = get_pt_idx(gt, locus)

			if fit == :drug
				if P[:multidrug]
					# binary select drugs (last result = no drugs)
					for L in 1:size(gt_fitness, 2)
						if mod(L ÷ 2^(locus -1), 2) == 1
							gt_fitness[gt, L] *= pt[idx]
						end #if
					end #for
				else
					# just one drug to worry about, so average across loci
					gt_fitness[gt, 1] += pt[idx] / num_loci
				end #if
			else
				# establishment cost
				gt_fitness[gt, 1] *= pt[idx]
			end #if
		end #for
	end #for

	gt_fitness
end #fn


#=
"""
    gt_fitness1(num_alleles::Int)

Generate the genotype fitness array for `num_alleles` assuming left dominance:

    1: [1,1,2]
    2: [1,1,1,2,2,3]
    3: [1,1,1,1,2,2,2,3,3,4]
"""
function gt_fitness1(num_alleles::Int)
	A = zeros(Int, triangle(num_alleles))
	idx = 0
	for i in 1:num_alleles, j in 1:(num_alleles + 1 - i)
		A[idx += 1] = i
	end #for
	A
end #fn
=#



