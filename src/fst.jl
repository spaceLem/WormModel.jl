include("multiloci.jl")

function calc_fst(W, num_loci, num_treated)
	FST = zeros(num_loci)

	num_gts, num_sheep = size(W)

	treated = 1:num_treated
	refugia = (num_treated + 1):num_sheep

	# divide the groups up
	gts1 = vec(sum(W[:, treated], dims=2))
	gts2 = vec(sum(W[:, refugia], dims=2))
	gts  = gts1 + gts2

	# allele frequencies in each group
	TAD = g2a(gts1)
	RAD = g2a(gts2)
	AD  = g2a(gts)

	# total number of worms in each group
	sum_gts1 = sum(gts1)
	sum_gts2 = sum(gts2)
	sum_gts  = sum_gts1 + sum_gts2

	for i in 1:num_loci
		# subgroup allele frequencies
		p1, q1 = TAD[(2i - 1) : 2i]
		p2, q2 = RAD[(2i - 1) : 2i]
		p, q   = AD[(2i - 1) : 2i]

		Hexp  = 2p * q
		Hexp1 = 2p1 * q1
		Hexp2 = 2p2 * q2

		H_S = (Hexp1 * sum_gts1 + Hexp2 * sum_gts2) / sum_gts
		H_T = Hexp

		FST[i] = H_T != 0.0 ? 1.0 - H_S / H_T : 0.0
	end #for

	FST
end #fn



function FST_strings(num_worms::Int = 5,
					 dna_length::Int = 10,
					 num_generations::Int = 10)

	# create a gene pool
	X = zeros(Int, num_worms, dna_length)

	# add some mutations
	for i in eachindex(X)
		if rand() < 0.1
			X[i] = true
		end #if
	end #for

	display(X)
	X2 = zeros(Int, size(X))

	# for each generation
	for i in 1:num_generations
		X2 .= X
		# choose parents and offspring
		m, f, c = rand(1:num_worms, 3)

		# start with one of the parents
		p = rand(Bool)

		# zip along one of the parents' DNA strands and copy it
		for j in 1:dna_length
			# randomly test for crossover (small probability)
			# in which case swap parents
			if rand() < 1e2
				p = !p
			end
			X2[c,j] = p ? X[m,j] : X[c,j]

			# sometimes add a random mutation
			if rand() < 1e3
				X[m,j] = rand(0:1)
			end
		end
		X .= X2
	end

	X
end #fn

