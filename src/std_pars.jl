"""
    std_pars(num_loci::Int=1)

Create a dictionary with all the parameters based on the number of loci `num_loci`.
"""
function std_pars(num_loci::Int=1)
	Dict(# genetic parameters
		 :num_loci          => num_loci, # no. of loci
		 :distn0            => repeat([0.9, 0.1], outer=num_loci),   # initial allele distribution
		 :fitness_drug      => repeat([0.05, 0.95], outer=num_loci), # prob surviving drug
		 :fitness_establish => repeat([1.0, 0.8], outer=num_loci),   # prob establishing
		 :dominance         => zeros(num_loci),                      # 0.0 to 1.0 (rec to dom)

		 # treatment parameters
		 :multidrug => false, # whether 1 locus/1 drug, or N loci/1 drug
		 :drug_set => collect(1:num_loci), # by default use all drugs

		 :p => 1.0, # proportion of sheep treated ∈ [0,1]
		 :ϕ => 0.2, # mixing proportion ∈ [0,1]
		 :freq_treatments => 4, # e.g. treat 4x per year
		 # how to treat multiple drugs
		 # 1 = ?
		 # 2 = ?
		 :regimen => 1,


		 # model setup parameters
		 :T         => 10.0,    # timescale in years
		 :num_itns  => 10,      # no. of iterations
		 :num_recs  => 20 * 12, # record resolution (num per itn)
		 :num_sheep => 20,      # no. of sheep in the flock
		 :worms0    => 10^3,    # initial no. of worms per sheep

		 # model rate parameters
		 :ρ  => 1e2,  # larvae output for females
		 :μL => 2.0,  # larvae mortality rate
		 :μW => 1e-2, # density independent mortality rate
		 :νW => 1e-4, # density dependent mortality rate
		 :β  => 1e-3, # ingestion (infection) rate

		 # not yet implemented
		 # clumping = 1.0  # see Cornell (2003)

		 # switches
		 :cross => true,    # allow heterozygotes
		 :debug => false,   # turn on extra information
		 :progress => true # show progress meter
		 )
	
end #fn

