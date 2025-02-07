
"""
    gethostdir()

    Chrom atlas data directory on slade and penrose.
"""
function gethostdir()

    host = gethostname()
    if host == "slade"
        return "/slade/projects/Research_Project-MRC158833/tfs_chromatin/chromatlas"
    elseif host == "penrose"
        return "/penrose/projects/data/chromatlas"
    else
        error("$host not recognosied")
    end
end


"""
    loadchromatlasmeta(; folder=gethostdir())

    Load data frame of bed files and bigwig files of chrom atlas
"""
function loadchromatlasmeta(; folder==gethostdir())
    bedfiles = glob("*.bed", joinpath(folder, "peaks"))
    samples = replace.(basename.(bedfiles), ".bed" => "")
    totalpeaks = @showprogress map(countlines, bedfiles);
    bigwigs = joinpath.(folder, "bigwig", string.(samples, ".bw"))

    @assert all(isfile, bigwigs)
    bwdict = Dict(s => bw for (s, bw) in zip(samples, bigwigs));
    meta = DataFrame(Study="ChromAtlas", Cell=samples, TotalPeaks=totalpeaks, BedFile=bedfiles, BigWigFile=bigwigs)
    meta, bwdict
end


#### helper functions
function cellind(cell, cells)
    i = findfirst(c -> cell .== c, cells)
    if isnothing(i)
        error("Cell: $cell not found")
    end
    i
end
bwfile(cell, meta) = meta.BigWigFile[cellind(cell, meta.Cell)]


"""
    bigwigdata(chrom, loc, bwf) 

    extract the bigwig data for `chrom` and `loc`
    
"""
function bigwigdata(chrom, loc, bwf)
    reader = open(BigWig.Reader, bwf)
    data = zeros(length(loc))
    
    for record in eachoverlap(reader, Interval(chrom, loc))
        
        s = BigWig.chromstart(record)
        e = BigWig.chromend(record)
        v = BigWig.value(record)
        rind = (s:e) .- first(loc) .+ 1
        ind = intersect(rind, 1:length(data))
        data[ind] .= v
        
    end
    close(reader)
    data
end


"""
    normbigwigdata(chrom, loc, bwf; ex=100_000 normloc=(first(loc) - ex):(last(loc) + ex))

    extract the bigwig data for `chrom` and `loc`, normalising to max data +/- `ex` bp around the interval
    
"""
function normbigwigdata(chrom, loc, bwf; ex=100_000, normloc=(first(loc) - ex):(last(loc) + ex), applynorm=true)
    reader = open(BigWig.Reader, bwf)
    data = zeros(length(loc))

    m = 0
    for record in eachoverlap(reader, Interval(chrom, normloc))
        
        s = BigWig.chromstart(record)
        e = BigWig.chromend(record)
        v = BigWig.value(record)
        rind = (s:e) .- first(loc) .+ 1
        ind = intersect(rind, 1:length(data))
        if !isempty(ind)
            data[ind] .= v
        end
        m = max(m, v)
        
    end
    close(reader)
    applynorm && (data ./= m)
    data
end

"""
    avh(H::Vector{T}, δ) where {T}
    average a vector into bins of δ in width
"""
function avh(H::Vector{T}, δ) where {T}
    n = length(H)

    w = cld(n, δ)
    AH = zeros(Float64, w)
    te = zeros(Int, w)
    for i = 1:n
        wi = cld(i, δ)
        te[wi] += 1
        AH[wi] += H[i]
    end
    if te[end] < 0.75δ
        AH[end] += AH[end-1]
        te[end] += te[end-1]
    end
    AH./te
end

bigwigdatasets(chrom, loc, bwfs) = @showprogress [bigwigdata(chrom, loc, bwf) for bwf in bwfs]
normbigwigdatasets(chrom, loc, bwfs; kwargs...) = @showprogress [normbigwigdata(chrom, loc, bwf; kwargs...) for bwf in bwfs]



maxnorm(x) = x/maximum(x)

function clusterall(chrom, loc, meta; k =5, dx=1, normfun=identity, ex=100_000, maxnorm=true)
    data = normbigwigdatasets(chrom, loc, meta.BigWigFile, ex=ex, applynorm=maxnorm)
    KMK = clusterdata(data, k=k, dx=dx, normfun=normfun)
    (chrom=chrom, loc=loc, KMK...)
end

function clusterdata(data; k = 5, dx = 1, normfun=identity, seed=1618)
    SP = replace(mapreduce(p -> avh(normfun(p), dx), hcat, data), NaN => 0)
    @time KMK = kmeansorder(SP, k, seed);
    KMK
end



clustertable(KMK, meta) = DataFrame(Cell=meta.Cell, Cluster=KMK.A)
clustersummary(KMK, meta) = combine(groupby(clustertable(KMK, meta), :Cluster), nrow => :count, :Cell => (x -> join(x, ", ")) => :Cell)