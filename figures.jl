### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ e4fced54-c7cf-11ed-2b5d-19dde8519a7b
begin
	using Pkg; Pkg.activate(".")
	using CSV
	using DataFrames
	using Plots
	using StatsBase
end

# ╔═╡ 2ecdfc10-9138-4918-8194-d6f28ec96f20
n = 25

# ╔═╡ ff709c3b-a823-40c6-ad90-777bb4e863f0
dat_iqtree = CSV.read("measures_iqtree_n$(n).csv", DataFrame)

# ╔═╡ 624d0256-6209-4ca9-a8b0-be52c2a2f759
dat_phylo = CSV.read("measures_phylo_n$(n).csv", DataFrame)

# ╔═╡ 34df710f-6bce-49fd-9404-f4e1bf030b45
sort(dat_phylo, [:spr_distance])

# ╔═╡ e78a8084-b751-45fe-9fa8-9988f7d8b2b9
md"# RF distance"

# ╔═╡ 61676039-3ed8-49fb-9e9c-4b801c3d317d
plts_rf = let
	plts = Dict()
	bins = 0:0.1:1
	bin_centers = (bins[2:end] + bins[1:end-1])/2
	
	for (k, df) in pairs(groupby(dat_phylo, :α))
		α = k.α
		h = fit(Histogram, df.modified_rf_distance, bins)
		plts[α] = plot(bin_centers, h.weights, label = "phyloformer",)
	end
	for (k, df) in pairs(groupby(dat_iqtree, :α))
		α = k.α
		h = fit(Histogram, df.modified_rf_distance, bins)
		plot!(plts[α], bin_centers, h.weights, label = "iqtree",)
	end
		
	for (α, p) in plts 
		plot!(p, xlabel = "RF (normalized)")
	end
	plts
end

# ╔═╡ 2db5bb89-dcc6-4dd5-971a-6cffacc4124f
plts_rf[0.8]

# ╔═╡ 5ed7fd51-82bf-4d06-8a47-afb616ebaf31
savefig(plts_rf[0.8], "RF_distance_histogram.png")

# ╔═╡ 9fefe9e2-26fe-47b1-a264-242c1593b697
md"# SPR distance"

# ╔═╡ 0cc07eb1-1f22-466b-95d6-1f8990f2bdca
plts_spr = let
	plts = Dict()
	bins = -0.5:1:12.5
	bin_centers = (bins[2:end] + bins[1:end-1])/2
	
	for (k, df) in pairs(groupby(dat_phylo, :α))
		α = k.α
		h = fit(Histogram, df.spr_distance, bins)
		plts[α] = plot(bin_centers, h.weights, label = "phyloformer",)
	end
	for (k, df) in pairs(groupby(dat_iqtree, :α))
		α = k.α
		h = fit(Histogram, df.spr_distance, bins)
		plot!(plts[α], bin_centers, h.weights, label = "iqtree",)
	end
		
	for (α, p) in plts 
		plot!(p, xlabel = "SPR distance")
	end
	plts
end

# ╔═╡ 4e73e972-f10e-46df-a39f-91aff0367cea
plts_spr[0.8]

# ╔═╡ bd4b0614-d152-4d2d-9f61-6816c6592dd0
savefig(plts_spr[0.8], "SPR_distance_histogram.png")

# ╔═╡ Cell order:
# ╠═e4fced54-c7cf-11ed-2b5d-19dde8519a7b
# ╠═2ecdfc10-9138-4918-8194-d6f28ec96f20
# ╠═ff709c3b-a823-40c6-ad90-777bb4e863f0
# ╠═624d0256-6209-4ca9-a8b0-be52c2a2f759
# ╠═34df710f-6bce-49fd-9404-f4e1bf030b45
# ╠═e78a8084-b751-45fe-9fa8-9988f7d8b2b9
# ╠═61676039-3ed8-49fb-9e9c-4b801c3d317d
# ╠═2db5bb89-dcc6-4dd5-971a-6cffacc4124f
# ╠═5ed7fd51-82bf-4d06-8a47-afb616ebaf31
# ╠═9fefe9e2-26fe-47b1-a264-242c1593b697
# ╠═0cc07eb1-1f22-466b-95d6-1f8990f2bdca
# ╠═4e73e972-f10e-46df-a39f-91aff0367cea
# ╠═bd4b0614-d152-4d2d-9f61-6816c6592dd0
