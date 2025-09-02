using Nucleus

using ImmuneReceptor

# ---- #

const S1_ = String[]

const S2_ = String[]

const S3_ = String[]

for fa in (
    joinpath(ImmuneReceptor.IN, "gliph", "db", fa) for fa in (
        "rubelt-naive-CD4.fa",
        "rubelt-naive-CD8.fa",
        "tcrab-naive-refdb.fa",
        "warren-naive.fa",
    )
)

    for s1 in eachsplit(read(fa, String), '>'; keepempty = false)

        s2, s3 = split(s1, '\n')

        @assert s3[1] === 'C'

        @assert s3[end] === 'F'

        push!(S3_, s3)

        s4, s5, _ = split(s2, ','; limit = 3)

        push!(S1_, s4)

        push!(S2_, s5)

    end

end

Nucleus.Plotly.writ(
    joinpath(ImmuneReceptor.OU, "1.html"),
    (Dict("type" => "histogram", "x" => map(lastindex, S3_)),),
    Dict(
        "title" => Dict("text" => "$(lastindex(unique(S3_))) unique CDR3s"),
        "yaxis" => Dict("title" => Dict("text" => "Count")),
        "xaxis" => Dict("title" => Dict("text" => "Number of amino acids in CDR3")),
    ),
)

# ---- #

const U1_ = unique(S1_)

const U2_ = unique(S2_)

const D1 = Dict(st => nd for (nd, st) in enumerate(U1_))

const D2 = Dict(st => nd for (nd, st) in enumerate(U2_))

const I = zeros(Int, lastindex(U1_), lastindex(U2_))

for nd in eachindex(S1_)

    I[D1[S1_[nd]], D2[S2_[nd]]] += 1

end

Nucleus.HeatPlot.writ(
    joinpath(ImmuneReceptor.OU, "2.html"),
    U1_,
    U2_,
    I,
    Dict(
        "title" => Dict("text" => "$(lastindex(S1_)) TCRs"),
        "yaxis" => Dict("title" => Dict("text" => "7 V gene")),
        "xaxis" => Dict("title" => Dict("text" => "7 J gene")),
    ),
)
