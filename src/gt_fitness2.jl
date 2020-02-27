"""
    calc_gt_fitness(P; fit::Symbol = :drug)

Calculate genotype fitness given parameters `P`, and choose `fit` is `:drug` or `:establish`. This only needs to be done once, at the start of the simulation.
"""
function calc_gt_fitness2(P; fit::Symbol = :drug)
	loci = P[:loci]          # [2, 1]
	num_drugs = length(loci) # 2
	num_loci  = sum(loci)    # 3
	num_gts   = 3^num_loci   # 27
	dominance = P[:dominance]

	# get the allele fitnesses
	allele_fitness = P[fit == :drug ? :fitness_R : :fitness_establish]
	base_fitness = fit == :drug ? P[:fitness_S] : ones(num_loci)

	# less work if establishing
	# some duplication here, may be able to remove it later
	if fit == :establish
		gt_fitness = ones(num_gts, 1)

		for locus in 1:num_loci
			dom = dominance[locus]
			Sx = 1.0
			Rx = allele_fitness[locus]
			pt = [Sx, (1 - dom) * Sx + dom * Rx, Rx]
			gt = repeat(pt, inner=3^(locus - 1), outer=3^(num_loci - locus))
			gt_fitness .*=  gt
		end #for

		return gt_fitness
	end #if

	# store the fitness for each genotype with each drug combination
	gt_fitness1 = zeros(num_gts, fit == :drug ? num_drugs : 1)
	gt_fitness = ones(num_gts, fit == :drug ? 2^num_drugs : 1)

	offset = 0

	for drug in 1:num_drugs, locus in 1:loci[drug]
		lidx = locus + offset

		dom = dominance[lidx]
		Sx = base_fitness[drug]
		Rx = allele_fitness[lidx]
		pt = [Sx, (1 - dom) * Sx + dom * Rx, Rx]
		gt = repeat(pt, inner=3^(lidx - 1), outer=3^(num_loci - lidx))
		gt_fitness1[:, 2^(drug-1)] .+=  gt / loci[drug]

		if locus == loci[drug]
			offset += loci[drug]
		end #if
	end #for

	# calculate drug combinations
	for i in 1:2^num_drugs, drug in 1:num_drugs
		# basically look at drug combination in binary and x those drugs together
		if mod(i ÷ 2^(drug -1), 2) == 1
			gt_fitness[:, i] .*= gt_fitness1[:, drug]
		end #if
	end #for

	gt_fitness
end #fn

