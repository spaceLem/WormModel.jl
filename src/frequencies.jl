"""
    unmixed_a2g(Ad::Array{Float64,1})

Take an array of allele distributions `Ad` and return the genotype frequencies.
Assumes no heterozygotes are possible.
"""
function unmixed_a2g(Ad::Array{Float64,1})
	num_loci = length(Ad) ÷ 2
	num_gts = 3^num_loci
	Gd = ones(num_gts)

	freqs = @MVector zeros(3)

	# for each locus in turn, find gt distribution
	# then multiply all gt indices that have that gt
	for locus in 1:num_loci
		freqs .= @MVector [Ad[2locus-1], 0, Ad[2locus]]

		for gt in 1:num_gts
			idx = get_pt_idx(gt, locus)
			Gd[gt] *= freqs[idx]
		end #for
	end #for

	normalize!(Gd, 1)
end #fn


"""
    a2g1(Ad::Array{Float64,1})

Take an array of allele distributions `Ad` and return the genotype frequencies (1-dim).
"""
function a2g1(Ad::Array{Float64,1})
	@SVector [Ad[1]^2,
			  2 * Ad[1] * Ad[2],
			  Ad[2]^2]
end #fn


"""
    a2g(Ad::Array{Float64,1})

Take an array of allele distributions `Ad` and return the genotype frequencies (N-dim).
"""
function a2g(Ad::Array{Float64,1})
	num_loci = length(Ad) ÷ 2
	num_gts = 3^num_loci
	Gd = ones(num_gts)

	# for each locus in turn, find gt distribution
	# then multiply all gt indices that have that gt
	for locus in 1:num_loci
		#freqs = a2g1(Ad[2i - 1:2i])
		p, q = Ad[2locus-1], Ad[2locus]
		freqs = @SVector [p^2, 2p*q, q^2]

		for gt in 1:num_gts
			idx = get_pt_idx(gt, locus)
			Gd[gt] *= freqs[idx]
		end #for
	end #for

	normalize!(Gd, 1)
end #fn


"""
    g2a_mat(num_loci::Int)

Get matrix for g2a with `num_loci` loci. Looks like:
```
A =
3×2 Array{Int64,2}:
 2 0
 1 1
 0 2

A =
6×3 Array{Int64,2}:
 2  0  0
 1  1  0
 1  0  1
 0  2  0
 0  1  1
 0  0  2
```
"""
function g2a_mat(num_loci::Int)
	A = zeros(Int, triangle(num_loci), num_loci)
	row = 0
	for offset = 1:num_loci
		A[row += 1, offset] = 2
		for col = offset + 1:num_loci
			A[row += 1, offset] = 1
			A[row, col] = 1
		end #for
	end #for

	A
end #fn


"""
    g2a(Gf::Array{Int,1})

Turn genotype frequency `Gf` into allele distribution `Ad`.
"""
function g2a(Gf::AbstractArray{Int,1})
	num_loci = round(Int, log(3, length(Gf)))
	num_gts = length(Gf)

	Ad = zeros(2num_loci)

	# effectively prevents divide-by-zero
	if sum(Gf) == 0
		return Ad
	end #if

	#A = [2 1 0; 0 1 2]
	A = @SArray [2 1 0; 0 1 2]

	for locus in 1:num_loci
		Adr = view(Ad, 2locus-1:2locus)
		for gt in 1:num_gts
			idx = get_pt_idx(gt, locus)
			Adr .+= view(A, :, idx) .* Gf[gt]
			#Adr .+= A[:, idx] .* Gf[j]
		end #for
		Adr ./= sum(Adr)
	end #for

	Ad
end #fn



