### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 06591e16-928b-11ed-3928-6fa568fc6275
begin
	using Pkg; Pkg.activate(".")
	using Chain
	using Roots
	using TreeAlgs
	using TreeTools
	using BackwardCoalescent
end

# ╔═╡ fdffe5d3-1a93-4796-8e10-14821f71812b
function outfolder(n, L, μ, α, rep)
	"data_n$(n)/n$(n)_L$(L)_alpha$(round(α, sigdigits=2))_μ$(round(μ, sigdigits=2))_$(rep)"
end

# ╔═╡ 5e19c415-6f0d-4524-975e-97c93d6f3cb1
begin
	n = 25
	N = 1
	L = 500
	αvals = [0.8] # ratio of internal nodes to leaves
	n_trees = 50
end

# ╔═╡ 26765721-552b-45cb-83ab-d022ea5f30c3
outfolder(n, L, 0.03, αvals[1], 1)

# ╔═╡ 7381df1e-45c5-4e21-b9a0-9a11163ac993
genealogy(KingmanCoalescent(5, N))

# ╔═╡ 7391b3ac-10b7-4ef4-a699-28ef315ddcd5
coalescent = KingmanCoalescent(n, N)

# ╔═╡ a0cf0743-864d-42d3-870c-26eb250aaec7
md"## Functions"

# ╔═╡ a90579d5-08f2-49bf-8c4f-b18b1b6fb1c9
md"""
## Mutation rate

I aim for a ratio of internal nodes to leaves similar to what we have in flu: $\alpha \sim 0.7$. 
For a given $\mu$, the probability that a branch of length $t$ is removed $e^{-\mu Lt}$ (no mutations).

Let us assume that the branches above internal nodes have times $t_1 \ldots t_n$ (in reality $n-1$ but $n\gg 1$). The number of branches that will be removed is 

$$\sum_{i=1}^n e^{-\mu L t_i} = (1-\alpha) n$$

by definition of $\alpha$ (we removed a fraction $1-\alpha$ of the $n$ internals). So we have to pick a corresponding $\mu$. 

We can trivially constraint the search by noting that 

$$\begin{align}
&n e^{-\mu L t_{max}} \leq (1-\alpha) n = \sum_{i=1}^n e^{-\mu L t_i} \leq n e^{-\mu L  t_{min}} \\
&\rightarrow -\frac{\log(1-\alpha)}{t_{max}} \leq \mu L \leq -\frac{\log(1-\alpha)}{t_{min}}
\end{align}$$
"""

# ╔═╡ 75a62135-3215-4ba1-88e1-806065831a06
function pick_mu(tree, α)
	times = @chain tree internals map(branch_length, _) skipmissing collect
	n = length(times)
	f(μ) = sum(t -> exp(-μ*t), times) - (1-α)*n

	tmin, tmax = extrema(times)
	μ_range =  (-log(1-α)/tmax, -log(1-α)/tmin)

	μ = try
		find_zero(f, μ_range, Bisection())
	catch err
		@info f.(μ_range)
		error(err)
	end
	return μ 
end

# ╔═╡ a37cd1cd-a7be-4474-a143-17c71729a157
function test_remove(tree, μ)
	n = length(leaves(tree))
	z = 0
	for x in internals(tree)
		if !isroot(x) && rand() < exp(-μ * branch_length(x))
			z += 1
		end
	end
	return z/n
end

# ╔═╡ 5693e3e5-1969-4c51-b647-093a3b59325c
md"""
Below: function to set tree branch length to the number of mutations that occured on the branch. This is the best a tree builder can do (without a prior)
"""

# ╔═╡ ba68356d-ab8d-48c3-b596-5ad7e19955e4
function mutation_branch_length!(tree)
	@chain tree nodes(_; skiproot=true) map(_) do n
		seq_c = n.data[:seq]
		seq_a = ancestor(n).data[:seq]
		h = count(x -> x[1] != x[2], zip(seq_c, seq_a)) # hamming
		branch_length!(n, h)
	end

	TreeTools.delete_null_branches!(tree; threshold = 1/2)
	return tree
end

# ╔═╡ 7da5844a-1952-4536-9054-4f51454ebe02
for α in αvals, rep in  1:n_trees
	tree = @chain coalescent genealogy convert(Tree{MiscData}, _)
	label!(tree, tree.root, "root")

	# add_outgroup!(tree)
	
	μ = pick_mu(tree, α) / L
	TreeAlgs.Evolve.evolve!(tree, L, μ)

	mkpath(outfolder(n, L, μ, α, rep))
	
	fasta_file = outfolder(n, L, μ, α, rep) * "/alignment.fasta"
	TreeAlgs.Sequences.write_fasta(fasta_file, tree, :seq, root=false)

	tree_file = outfolder(n, L, μ, α, rep) * "/real_tree.nwk"
	write(tree_file, tree)

	tree_file_mut_length = outfolder(n, L, μ, α, rep) * "/real_tree_mut_length.nwk"
	_t = mutation_branch_length!(copy(tree))
	write(tree_file_mut_length, _t)
end

# ╔═╡ d67291b4-8bc4-433c-9603-bcbc8526042e
function add_outgroup!(tree::Tree)
	height = distance(first(leaves(tree)), tree.root) # distance from leaves to root
	outgroup = TreeNode(
		tau = height/2, 
		label = "outgroup",
		data = MiscData(),
	)
	new_root = TreeNode(label = "root", data = MiscData())
	TreeTools.graftnode!(new_root, tree.root, tau = height/2)
	label!(tree, tree.root, "NODE_0")
	node2tree!(tree, new_root)
	graft!(tree, outgroup, tree.root)

	return nothing
end

# ╔═╡ d6287502-3a5b-4283-bb57-d4d862923ca7
md"# Tests"

# ╔═╡ Cell order:
# ╠═06591e16-928b-11ed-3928-6fa568fc6275
# ╠═fdffe5d3-1a93-4796-8e10-14821f71812b
# ╠═26765721-552b-45cb-83ab-d022ea5f30c3
# ╠═7381df1e-45c5-4e21-b9a0-9a11163ac993
# ╠═5e19c415-6f0d-4524-975e-97c93d6f3cb1
# ╠═7391b3ac-10b7-4ef4-a699-28ef315ddcd5
# ╠═7da5844a-1952-4536-9054-4f51454ebe02
# ╟─a0cf0743-864d-42d3-870c-26eb250aaec7
# ╟─a90579d5-08f2-49bf-8c4f-b18b1b6fb1c9
# ╠═75a62135-3215-4ba1-88e1-806065831a06
# ╠═a37cd1cd-a7be-4474-a143-17c71729a157
# ╟─5693e3e5-1969-4c51-b647-093a3b59325c
# ╠═ba68356d-ab8d-48c3-b596-5ad7e19955e4
# ╠═d67291b4-8bc4-433c-9603-bcbc8526042e
# ╟─d6287502-3a5b-4283-bb57-d4d862923ca7
