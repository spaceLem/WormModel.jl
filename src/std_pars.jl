function std_pars(num_loci::Int=1)

	# all ∈ [0,1]
	distn0            = repeat([0.9, 0.1], outer=num_loci)
	fitness_drug      = repeat([0.05, 0.95], outer=num_loci)
	fitness_establish = repeat([1.0, 0.8], outer=num_loci)
	dominance         = zeros(num_loci)

	P = Dict(# genetic parameters
			 "num_loci"          => num_loci,          # no. of loci
			 "distn0"            => distn0,            # initial allele distribution
			 "fitness_drug"      => fitness_drug,      # prob surviving drug
			 "fitness_establish" => fitness_establish, # prob establishing
			 "dominance"         => dominance, # 0.0 to 1.0 (rec to dom)

			 # treatment parameters
			 "p" => 1.0, # proportion of sheep treated ∈ [0,1]
			 "ϕ" => 0.2, # fraction of between-group transmission ∈ [0, 0.5]
			 "freq_treatments" => 6, # e.g. treat 6x per year

			 # model setup parameters
			 "T"         => 5.0,       # timescale in years
			 "num_itns"  => 10,        # no. of iterations
			 "num_recs"  => 20 * 12,   # record resolution (num per itn)
			 "num_sheep" => 20,        # no. of sheep in the flock
			 "worms0"    => 10^3,      # initial no. of worms per sheep

			 # model rate parameters
			 "ρ"  => 1e2,  # larvae output for females
			 "μL" => 2.0,  # larvae mortality rate
			 "μW" => 1e-2, # density independent mortality rate
			 "νW" => 1e-4, # density dependent mortality rate
			 "β"  => 1e-3, # ingestion (infection) rate
			 # clumping = 1.0  # see Cornell (2003)

			 "debug" => false,
			 "resample" => true)

end #fn

