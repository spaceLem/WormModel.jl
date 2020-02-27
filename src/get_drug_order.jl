"""
    get_drug_order(P)

Provide the treatment schedule for all the drugs described by `P`.
"""
function get_drug_order(P)
	num_loci = P[:num_loci]
	drug_set = P[:drug_set]
	num_drugs = length(drug_set)

	num_treatments = round(Int, P[:T] * P[:freq_treatments])

	if P[:multidrug]
		regimen = P[:regimen]

		drug_order = if regimen == 1 # 1 1 1 ... 2 2 2 ...
			repeat(drug_set, inner = num_treatments ÷ num_drugs)
		elseif regimen == 2 # 1 2 1 2 1 2 ...
			repeat(drug_set, num_treatments ÷ num_drugs)
		elseif regimen == 3
			repeat([1,1,1,1,2], num_treatments)[1:num_treatments]
		elseif regimen == 0 # 0 0 0 ... 0
			# 2^num_loci is use all drug_set together
			repeat([sum(2 .^(drug_set .- 1))], num_treatments)
		end #if

	else
		# just use drug no. 1 repeatedly
		drug_order = repeat([1], num_treatments)
	end #if

	drug_order
end #fn

