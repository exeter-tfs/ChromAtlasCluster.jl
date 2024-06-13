# ChromAtlasCluster.jl


Julia project for clustering single-cell chromatin accessibility atlas on `slade.ex.ac.uk`. Data taken on chromatin accessibility of 222 cell types/tissues.

Zhang, K., Hocker, J. D., Miller, M., Hou, X., Chiou, J., Poirion, O. B., Qiu, Y., Li, Y. E., Gaulton, K. J., Wang, A., Preissl, S., & Ren, B. (2021). A single-cell atlas of chromatin accessibility in the human genome. *Cell*, 184(24), 5985–6001.e19. https://doi.org/10.1016/j.cell.2021.10.024


#### All data is against hg38

### Prerequistes

Julia ≥ 1.10

## Installation
```bash
git clone https://github.com/exeter-tfs/ChromAtlasCluster.jl
cd ChromAtlasCluster.jl
julia
```
Within julia activiate the local 
```julia
] # to enter into Pkg mode
activate .
instantiate ## for first time installation