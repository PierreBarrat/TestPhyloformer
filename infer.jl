### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 547849af-bf8b-42cc-99da-8e1f61406f64
begin
	using Pkg; Pkg.activate(".")
	using TreeTools
	using Chain
end

# ╔═╡ f6473643-06d2-4c92-bf78-9a6f85dcf5fb
md"""
## Note

This is meant to be ran as a script in the right conda environment for phyloformer.
"""

# ╔═╡ 6d736090-9356-11ed-122a-ebc9609d61e8
dir = ARGS[1]

# ╔═╡ e9d93ec7-f618-4dbe-b8aa-813ef136aeaa
for f in readdir(dir, join=true)
	try
		run(`iqtree.jl --aln $(f*"/alignment.fasta") -o $(f*"/inferred_iqtree.nwk") --cleanup`)
		tree = read_tree(f*"/inferred_iqtree.nwk")
		# TreeTools.root!(tree, @chain tree["outgroup"] ancestor label)
		TreeTools.root!(tree; method = :midpoint)
		write(f*"/inferred_iqtree.nwk", tree)
	catch err
	end
end

# ╔═╡ 3d0fd6cb-2ae3-4c4c-a6f4-426d36b54422
for f in readdir(dir, join=true)
	cmd = Cmd([
		"python",
		homedir() * "/Documents/BaleLabo/Code/External/Phyloformer/predict.py",
		f,
	])
	run(cmd)
	tree = read_tree(f*"/predicted_alignment.nwk")
	# TreeTools.root!(tree, @chain tree["outgroup"] ancestor label)
	TreeTools.root!(tree; method = :midpoint)
	write(f*"/inferred_phyloformer.nwk", tree)
end

# ╔═╡ Cell order:
# ╠═547849af-bf8b-42cc-99da-8e1f61406f64
# ╠═f6473643-06d2-4c92-bf78-9a6f85dcf5fb
# ╠═6d736090-9356-11ed-122a-ebc9609d61e8
# ╠═e9d93ec7-f618-4dbe-b8aa-813ef136aeaa
# ╠═3d0fd6cb-2ae3-4c4c-a6f4-426d36b54422
