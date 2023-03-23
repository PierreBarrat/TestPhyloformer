# Installing

- Install [julia](https://julialang.org/)
- Clone this repo and navigate to it
- start julia and enter the following
  ```
  using Pkg
  Pkg.activate(".") # activates the local environment
  Pkg.instantiate() # download and precompile the needed packages using Manifest.toml
  ```

*Note*: the scripts `make_trees_kingman.jl`, `figures.jl` and `measures.jl` are actually [Pluto](https://github.com/fonsp/Pluto.jl) notebooks. You can open them in interactive mode in this way
```
using Pkg
Pkg.activate(".") # the environment should have Pluto installed
using Pluto; Pluto.run() # opens in a browser --> open the notebooks
```
This may be better if you want to edit them. 

# Generating trees and alignments

Run `julia make_trees_kingman.jl`. This will add new trees to the `data_n25` folder. You can change `n` (number of leaves) by editing the script. 

Details: 
- the trees are sampled from the Kingman coalescent, using by default `n=25` leaves. The population size `N` just gives an overall scale for branch length, so it is not very relevant. 
- For each tree, a mutation rate `m` is chosen such that the resolution index is about 0.8. The resolution index is defined by `R = (I-1) / (L-2)` where `L` is the number of leaves and `I` the number of internal nodes. 
- Sequences are sampled using the JC69 model, a sequence length of 500 and the mutation rate chosen above. The code for this is [here](https://github.com/PierreBarrat/TreeAlgs/tree/main/src/Evolve)

# Inferring trees

Run `julia infer.jl data_n25` (or any other data folder you generated).

Trees are inferred using 
- `iqtree`, with pretty much default settings. The script `iqtree.jl` is a wrapper (which is completely overkill in this case). 
- `phyloformer`

# Measures

I looked at the RF distance (corrected to account for not fully resolved trees, see `measures.jl`) and the SPR distance (using [TreeKnit](https://github.com/PierreBarrat/TreeKnit.jl)). Run `julia measures.jl` and `julia figures.jl` to generate the data tables and figures. 


