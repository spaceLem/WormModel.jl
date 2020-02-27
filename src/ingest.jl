"""
    function ingest!(K_establish::Array{Int,2},
                     K_ingest::Array{Int,2},
                     num_treated::Int)

Take the number of worms ingested `K_ingest` and the treatment number `num_treated` and calculate how many worms will establish in `K_establish`.
"""
function ingest!(K_establish::Array{Int,2},
				 K_ingest::Array{Int,2},
				 num_treated::Int)

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


