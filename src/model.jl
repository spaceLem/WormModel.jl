using Distributions
using LinearAlgebra
using ProgressMeter

include("utils.jl")
include("multiloci.jl")
include("fst.jl")

# main model
function model(P::Dict{String, Any})

	#for (key, value) in P println("$key = $value") end

	# meta parameters
	debug = P["debug"];
	resample = P["resample"]

	T         = P["T"] * 365.0 # timescale years -> days?
	num_itns  = P["num_itns"]  # no. of iterations
	num_recs  = P["num_recs"]  # record resolution (num per itn)
	num_sheep = P["num_sheep"] # no. of sheep in the flock
	num_loci  = P["num_loci"]  # no. of loci
	worms0    = P["worms0"]    # initial no. of worms per sheep


	# (daily) model parameters
	ρ  = P["ρ"]  # larvae output for females
	μL = P["μL"] # larvae mortality rate
	μW = P["μW"] # density independent mortality rate
	νW = P["νW"] # density dependent mortality rate
	β  = P["β"]  # ingestion (infection) rate
	# clumping = P["clumping"]  # see Cornell (2003)

	# treatment
	p = P["p"] # proportion treated ∈ [0, 1]
	ϕ = P["ϕ"] # mixing ∈ [0, 0.5]
	freq_treatments = P["freq_treatments"] # e.g. treat 6x per year

	mix1 = [1 - (1 - p) * ϕ, (1 - p) * ϕ]'
	mix2 = [p * ϕ, 1 - p * ϕ]'

	num_treated = floor(Int, p * num_sheep)              # ∈ [0, num_sheep]
	num_refugia = num_sheep - num_treated                # ∈ [0, num_sheep]
	treated     = 1:num_treated                          # 1 .. pN
	refugia     = (num_treated + 1):num_sheep            # pN+1 .. N
	treat_times = 365normalize(ones(freq_treatments), 1) # how long until next treatment

	# each locus gives a triangular number of combinations
	num_gts = 3^num_loci

	# these need to match the no. of genotypes
	gt_fitness_drug      = calc_gt_fitness(P["fitness_drug"], P["dominance"])
	gt_fitness_establish = calc_gt_fitness(P["fitness_establish"], P["dominance"])

	μ_drug = (1.0 .- gt_fitness_drug) .* 365 ./ freq_treatments

	# Storage for results, need:
	# - worm & larvae burden (2)
	# - proportion {susceptible, resistant} worms on {treated, refugia, global} (6)
	# - allele distribution in treated / refugia / global zones (6x loci)
	# - FST and before/after FST (2x loci)
	records   = zeros(8num_loci + 8, num_recs, num_itns)
	rec_width = T / num_recs

	# calculate initial genotype distribution
	gt_distn0 = resample ? a2g(P["distn0"]) : unmixed_a2g(P["distn0"])

	# Variables
	W = zeros(Int, num_gts, num_sheep) # Worms in each sheep
	L = zeros(Int, num_gts, 2)         # Larvae (treated / untreated)
	E = zeros(num_gts, 2)              # Egg genotype distn

	# Temporary variables for calculating FST between treatments
	Wtmp   = zeros(Int, size(W))
	FSTtmp = zeros(num_loci)

	# convenient to store the number of worms like this
	num_worms = zeros(Int, 1, num_sheep)

	# containers for event rates ()
	birth_L = zeros(size(L))
	death_L = zeros(size(L))
	death_W = zeros(size(W))
	# ingestion removed from L then randomly assigned to W
	ingest = zeros(Float64, size(L))

	# temporary arrays for changes
	dW = zeros(Int, size(W))
	dL = zeros(Int, size(L))

	K_birth_L   = zeros(Int, size(L))
	K_death_L   = zeros(Int, size(L))
	K_death_W   = zeros(Int, size(W))
	K_ingest    = zeros(Int, size(L))
	K_establish = zeros(Int, size(W))

	# keep an eye on dt
	# watch_dt = Array{Float64}(0)

	# Set up a handy progress bar
	prog = Progress(num_itns * num_recs, 1)

	tmp_gts = zeros(num_gts)

	for itn in 1:num_itns

		# populate worms with initial genotype frequencies
		W .= rand(Multinomial(worms0, gt_distn0), num_sheep)
		L .= 0
		t  = 0.0
		Wtmp .= 0
		FSTtmp .= 0.0

		dt  = 0.5
		rec = 0
		t_next_rec = 0.0
		treat_pos = 1
		t_next_treat = treat_times[treat_pos]

		while t < T
			# keep records
			while t > t_next_rec && rec < num_recs
				next!(prog)
				t_next_rec += rec_width

					#@show g2a(vec(sum(W[:, treated], dims=2)))
					#@show g2a(vec(sum(W[:, refugia], dims=2)))
					#@show calc_fst(W, num_loci, num_treated)
					#@show FSTtmp

				records[:, rec += 1, itn] = [
					# Worm burden
					sum(W) / num_sheep
					sum(L) / num_sheep
					# proportion susceptible / resistant worms on treated / refugia / global
					prop_sus(W[:, treated], gt_fitness_drug) # Treated S
					prop_res(W[:, treated], gt_fitness_drug) # Treated R
					prop_sus(W[:, refugia], gt_fitness_drug) # Refugia S
					prop_res(W[:, refugia], gt_fitness_drug) # Refugia R
					prop_sus(W, gt_fitness_drug) # Global S
					prop_res(W, gt_fitness_drug) # Global R
					# Allele distribution in treated / refugia / global zones
					g2a(vec(sum(W[:, treated], dims=2))) # Treated (x2 loci)
					g2a(vec(sum(W[:, refugia], dims=2))) # Refugia (x2 loci)
					g2a(vec(sum(W, dims=2)))             # Global (x2 loci)
					# FST
					calc_fst(W, num_loci, num_treated) # (x loci)
					FSTtmp                             # (x loci)
				]
			end #while

			# Event: treat with drug
			while t > t_next_treat
				Wtmp .= W
				treat_pos = mod1(treat_pos + 1, length(treat_times))
				t_next_treat += treat_times[treat_pos]
				W[:, treated] .= rand.(Binomial.(W[:, treated], gt_fitness_drug))
				FSTtmp .= calc_fst([W Wtmp], num_loci, num_sheep)
			end #while


			# --- Event Rates ---

			# calc egg genotype distribution
			E .= 0
			num_worms .= sum(W, dims=1)
			for sheep in 1:num_sheep
				# need at least 2 worms to reproduce
				if num_worms[sheep] >= 2
					tmp_gts .= if resample
						num_worms[sheep] .* a2g(g2a(W[:, sheep]))
					else
						W[:, sheep]
					end #if
					E .+= (sheep <= num_treated ? mix1 : mix2) .* tmp_gts
				end #if
			end #for

			# This is where the eggs are mixed up on pasture
			# note: vec() necessary to get the dimensions right if L is 1-dim
			birth_L .= ρ * E
			death_L .= μL * L
			#death_W .= (μW .+ μ_drug .+ νW * num_worms) .* W
			death_W .= (μW .+ νW * num_worms) .* W
			ingest  .= β .* L

			# --- Apply Timestep ---

			# quick hack, fix this later
			dt *= 1.2

			dW .= 0
			dL .= 0

			while true
				K_birth_L .= rand.(Poisson.(birth_L * dt))
				K_death_L .= rand.(Poisson.(death_L * dt))
				K_death_W .= rand.(Poisson.(death_W * dt))
				K_ingest  .= rand.(Poisson.(ingest * dt))

				# calculation of K_establish is more involved
				ingest!(K_establish, K_ingest, num_treated)
				K_establish .= rand.(Binomial.(K_establish, gt_fitness_establish))

				dL .= + K_birth_L .- K_death_L .- K_ingest
				dW .= + K_establish .- K_death_W

				# overshoot tests if change leaves variable < 0
				if overshoot(W, dW) || overshoot(L, dL)
					dt /= 2
					continue
				end #if

				break
			end #while

			W .+= dW
			L .+= dL

			# ensure_positive!(W, "W")
			# ensure_positive!(L, "L")

			# push!(watch_dt, dt)

			t += dt
		end #while

		# finish recording if necessary
		while rec < num_recs
			next!(prog)

			# treated / refugia allele dists

			records[:, rec += 1, itn] = [
				# Worm burden
				sum(W) / num_sheep
				sum(L) / num_sheep
				# proportion susceptible / resistant worms on treated / refugia / global
				prop_sus(W[:, treated], gt_fitness_drug) # Treated S
				prop_res(W[:, treated], gt_fitness_drug) # Treated R
				prop_sus(W[:, refugia], gt_fitness_drug) # Refugia S
				prop_res(W[:, refugia], gt_fitness_drug) # Refugia R
				prop_sus(W, gt_fitness_drug) # Global S
				prop_res(W, gt_fitness_drug) # Global R
				# Allele distribution in treated / refugia / global zones
				g2a(vec(sum(W[:, treated], dims=2))) # Treated (x2 loci)
				g2a(vec(sum(W[:, refugia], dims=2))) # Refugia (x2 loci)
				g2a(vec(sum(W, dims=2)))             # Global (x2 loci)
				# FST
				calc_fst(W, num_loci, num_treated) # (x loci)
				FSTtmp                             # (x loci)
			]
		end #while
	end #for

	t = T * (0 : num_recs - 1) / (365 * num_recs)

	# println("dt stats: ", signif.(quantile(watch_dt), 3))


	#=
	records_dict = Dict("S"  => dropdims(sum(records[1:2:4num_loci, :, :], dims=1), dims=1) / num_loci,
						"R"  => dropdims(sum(records[2:2:4num_loci, :, :], dims=1), dims=1) / num_loci,
						"W"  => records[4num_loci + 1, :, :],
						"L"  => records[4num_loci + 2, :, :],
						"pS" => records[4num_loci + 3, :, :],
						"pR" => records[4num_loci + 4, :, :],
						"F"  => records[4num_loci + 7:5num_loci + 6, :, :],
						"Fb" => records[4num_loci + 7:end, :, :])
						=#

	return records, t
end #fn
