function plotregion(chrom, loc, cell, meta::DataFrame; norm=false, kwargs...)
   
    if norm
        P = normbigwigdata(chrom, loc, bwfile(cell, meta), normloc=loc)
    else
        P = bigwigdata(chrom, loc, bwfile(cell, meta))
    end
    plotregion(chrom, loc, cell, P; norm=norm, kwargs...)
end

function plotregion(chrom, loc, cell, P; dx=ifelse(length(loc) < 2000, 10, 1), coordscale=1e+6, norm=false, kwargs...)
    lxp = (first(loc):dx:last(loc))/coordscale
    plot(lxp, avh(P, dx), lab="", fill=0, grid=false; kwargs...)
    xlims!(first(lxp), last(lxp))
    ylm = last(ylims())
    if norm
        ylims!(0, 1)
        yticks!([0, 1])
    else
        yticks!([0, 10*div(ylm, 10)])
    end


    ap = first(lxp)*0.95 + last(lxp)*0.05
    annotate!((ap, ylm/2, text(cell, font("helvetica", 10, :left))))
    plot!(size=(900, 100))

end

function clusterlines!(KMK)
    cs = KMK.KM.counts[KMK.KSI] 
    hline!(cumsum(cs) .+ .5, c=:white, lab="")
    yticks!(cumsum(cs) .- cs./2, string.(1:KMK.k))
end

function plot_cluster_centres_group(KMK; markregion=nothing, kwargs...)
    pc = plot(; kwargs...)
    plot!(KMK.loc, KMK.KM.centers[:, KMK.KSI], xticks=false, lab=string.("C", (1:KMK.k)'), ylabel="Accessibility", title="Clustering of Human scATAC-seq Atlas")
    
    markregion!(markregion, c=:black, ls=:dash)
    
    ccs = KMK.KM.counts[KMK.KSI]
    ps = bar(1:KMK.k, ccs, c=1:KMK.k,  xlabel="Cluster Size", ylabel="Cluster", title="Cluster Sizes", leg=false, orientation=:horizontal, ylims=(0, KMK.k +1), yticks=(1:KMK.k, string.(1:KMK.k)), xlims=(0, 1.2*maximum(KMK.KM.counts)))
    annotate!([(ccs[i] + 2, i, text(string(ccs[i]), font("helvetica", 12, :left))) for i = 1:KMK.k], grid=false)
    plot(pc, ps, size=(700, 300), layout=grid(1, 2, widths=[0.7, 0.3]), bottom_margin=5mm, left_margin=3mm, titlefont=font("helvetica", 12), fontfamily="helvetica", grid=false)
end


function plot_cluster_centres_grid(KMK, meta; markregion=nothing, showcells=4, kwargs...)
    
    
    phs = Plots.Plot[]
    m = 0
    for i = 1:KMK.k
        
        t = string("C", i, ": ", KMK.KM.counts[KMK.KSI[i]], " cells")
        m = max(m, maximum(KMK.KM.centers[:, KMK.KSI[i]]))
        p = plot(KMK.loc, KMK.KM.centers[:, KMK.KSI[i]], xticks=false, lab="", grid=false, fill=0, fillalpha=0.05, title=t)
        markregion!(markregion, c=:black, ls=:dash)
        push!(phs, p)
    end
    
    for (i, p) in zip(1:KMK.k, phs)
        plot!(p, ylims=(0, 1.2*m))
        if showcells > 0
            xl = xlims()
            cells = meta.Cell[KMK.A .== i]
            cells = cells[1:min(length(cells), showcells)]
            annotate!(((xl[1] + xl[2])/2, m*1.1, text(join(cells, ", "), font("helvetica", 6))))
        end
    end
    
    plot(phs..., size=(1200, 500), fontfamily="helvetica", link=:y, titlefont=font(10, "helvetica"); kwargs...)
end

function cluster_heatmap(KMK; markregion=nothing, coordscale=1e+6, kwargs...)
    p = plot(fontfamily="helvetica")
    
    heatmap!(KMK.loc/coordscale, 1:size(KMK.X, 2), KMK.X[:, KMK.SI]', ylabel="Cluster", yflip=true; kwargs...)
    !isnothing(markregion) && markregion!(markregion, c=:white, ls=:dash, coordscale=coordscale)
    clusterlines!(KMK)
    p
end



# markregion!(::Nothing) = nothing
function markregion!(cl;  coordscale=1, kwargs...)
    chrom, loc = cl
    vline!([loc[1], loc[end]]./coordscale, c=:white, lab=""; kwargs...)
end
