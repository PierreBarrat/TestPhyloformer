### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 16e70742-9359-11ed-2627-b15a0591189a
begin
	using Pkg; Pkg.activate(".")
	using Chain
	using CSV
	using DataFrames
	using Logging
	using LoggingExtras
	using Plots
	using TreeKnit
	using TreeTools
end

# ╔═╡ ec7b6daa-5172-4e63-8d1c-cf145cf25248
md"""
The inferred tree is compared to a slightly modified "real" tree where: 
- the length of each branch is set to the number of mutations that occurred on it; 
- branches with no mutations are removed from the tree, creating polytomies. 

I believe this makes for a more fair comparison: the inference technique does not have to find splits that are impossible to find from the data (0 mutations), and should find branch length that are those observed (number of mutations. )
"""

# ╔═╡ e37fc29d-87a9-42dc-a79d-9d1bebacbe44
n = 25

# ╔═╡ df1f46ed-d6fb-4522-864a-7f90af454b23
datdir = "data_n$n"

# ╔═╡ 0e9641e0-2e03-46cc-b5d5-aef0a2a4cf31
md"## Functions"

# ╔═╡ 76630149-b555-4c70-8020-9e654348c02d
function is_valid_folder(f)
	contents = readdir(f)
	files = ["alignment.fasta", "inferred_iqtree.nwk", "real_tree_mut_length.nwk", "real_tree.nwk", "inferred_phyloformer.nwk"]
	return mapreduce(in(contents), &, files)
end

# ╔═╡ 944ca12a-3881-4d6e-95d4-552a5a9e4d35
function parse_folder(f)
	m = match(r"data_n[0-9]+/n([0-9]+)_L([0-9]+)_alpha([0-9]+\.[0-9]+)_μ([0-9]+\.[0-9]+)_([0-9]+)", f)
	return (
		n = parse(Int, m.captures[1]), 
		L = parse(Int, m.captures[2]),
		α = parse(Float64, m.captures[3]),
		μ = parse(Float64, m.captures[4]),
		rep = parse(Int, m.captures[5]),
	)
end

# ╔═╡ d8a78390-7c99-48ed-9fde-52fb71817234
function spr_distance(t1, t2)
	oa = OptArgs(γ=1.1, nMCMC=250)
	mccs = with_logger(MinLevelLogger(current_logger(), Logging.Info)) do 
		run_treeknit(t1, t2, oa)
	end
	return @chain mccs getproperty(_, :mccs) values first length _-1
end

# ╔═╡ ebde6999-89df-4341-9d2e-f13d1b116461
md"""
### Custom RF distance

We have a potentially unresolved real tree $T_r$ and an inferred tree $T_i$. The modified RF distance will penalize: 
1. the splits $Q_i$ that are in $T_i$ and are not *compatible* with $T_r$; 
2. the splits $Q_r$ that are in $T_r$ but not in $T_i$. 

The point is that splits in $T_i$ not in $T_r$ may be simply unresolved in $T_r$, in which case they should not be considered "mistakes". 

In brief, 

$$d = \frac{\vert Q_i \vert + \vert Q_r \vert}{Z},$$

where $Z$ is a normalization (number of splits in both trees). 
"""

# ╔═╡ 0280557e-13f9-4afa-8c58-499e16ee67da
function custom_rf(tr, ti)
	Sr = SplitList(tr)
	Si = SplitList(ti)

	@assert Sr.leaves == Si.leaves

	Qi = count(s -> !iscompatible(s, Sr), Si)
	Qr = count(s -> !in(s, Si), Sr)

	# -2 because of the irrelevant root split
	return (Qi + Qr)/(length(Sr) + length(Si) - 2) 
end

# ╔═╡ f1ac56a3-f852-42e8-b0e3-6e14a8c2e60b
"""
	measures!(real_tree, inf_tree)

Return named tuple. Modifies `inf_tree` (null branches).
"""
function measures!(real_tree, inf_tree)
	return (
		spr_distance = spr_distance(real_tree, inf_tree),
		rf_distance = TreeTools.RF_distance(real_tree, inf_tree, normalize=true),
		modified_rf_distance = custom_rf(real_tree, inf_tree),
	)
end

# ╔═╡ 3f7dc878-a2aa-4c07-83d9-75f91546dbd6
dat_iqtree = let
	dat = []
	
	for f in Iterators.filter(is_valid_folder, readdir(datdir, join=true))
		real_tree = read_tree(f * "/real_tree_mut_length.nwk")
		iqtree_tree = read_tree(f * "/inferred_iqtree.nwk")
		
		measures = measures!(real_tree, iqtree_tree)
		
		P = parse_folder(f)
		row = Dict(
			:n => P.n,
			:μ => P.μ,
			:α => P.α,
			:L => P.L,
			:rep => P.rep,
		)

		row[:α_observed] = TreeTools.resolution_index(real_tree)
		for (k,v) in pairs(measures)
			row[k] = v
		end
		row
		push!(dat, row)
	end
	DataFrame(dat)
end

# ╔═╡ dc3b51b1-2405-4cb1-8923-e2d390f3ff3a
dat_phylo = let
	dat = []
	
	for f in Iterators.filter(is_valid_folder, readdir(datdir, join=true))
		real_tree = read_tree(f * "/real_tree_mut_length.nwk")
		phyloformer_tree = read_tree(f * "/inferred_phyloformer.nwk")
		
		measures = measures!(real_tree, phyloformer_tree)
		
		P = parse_folder(f)
		row = Dict(
			:n => P.n,
			:μ => P.μ,
			:α => P.α,
			:L => P.L,
			:rep => P.rep,
		)

		row[:α_observed] = TreeTools.resolution_index(real_tree)
		for (k,v) in pairs(measures)
			row[k] = v
		end
		row
		push!(dat, row)
	end
	DataFrame(dat)
end

# ╔═╡ 2a0bade7-7dc4-4add-b069-ee41664eb3a1
let
	CSV.write("measures_iqtree_n$n.csv", dat_iqtree)
	CSV.write("measures_phylo_n$n.csv", dat_phylo)
end

# ╔═╡ ab8935a5-a0e9-474b-85bb-5c037443c29f
md"# Tests"

# ╔═╡ Cell order:
# ╠═16e70742-9359-11ed-2627-b15a0591189a
# ╟─ec7b6daa-5172-4e63-8d1c-cf145cf25248
# ╠═e37fc29d-87a9-42dc-a79d-9d1bebacbe44
# ╠═df1f46ed-d6fb-4522-864a-7f90af454b23
# ╠═3f7dc878-a2aa-4c07-83d9-75f91546dbd6
# ╠═dc3b51b1-2405-4cb1-8923-e2d390f3ff3a
# ╠═2a0bade7-7dc4-4add-b069-ee41664eb3a1
# ╟─0e9641e0-2e03-46cc-b5d5-aef0a2a4cf31
# ╠═76630149-b555-4c70-8020-9e654348c02d
# ╠═944ca12a-3881-4d6e-95d4-552a5a9e4d35
# ╠═f1ac56a3-f852-42e8-b0e3-6e14a8c2e60b
# ╠═d8a78390-7c99-48ed-9fde-52fb71817234
# ╟─ebde6999-89df-4341-9d2e-f13d1b116461
# ╠═0280557e-13f9-4afa-8c58-499e16ee67da
# ╠═ab8935a5-a0e9-474b-85bb-5c037443c29f
