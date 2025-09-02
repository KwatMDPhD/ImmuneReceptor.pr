using Nucleus

using ImmuneReceptor

# ---- #

const V1_ = String[]

const J1_ = String[]

const C1_ = String[]

for fa in (joinpath(ImmuneReceptor.IN, "gliph", "db", fa) for fa in (
    #"rubelt-naive-CD4.fa",
    #"rubelt-naive-CD8.fa",
    #"tcrab-naive-refdb.fa",
    "warren-naive.fa",
))

    for s1 in eachsplit(read(fa, String), '>'; keepempty = false)

        s2, s3 = split(s1, '\n')

        @assert s3[1] === 'C'

        @assert s3[end] === 'F'

        push!(C1_, s3)

        s4, s5, _ = split(s2, ','; limit = 3)

        push!(V1_, s4)

        push!(J1_, s5)

    end

end

# ---- #

const V2_::Vector{String}, J2_::Vector{String}, C2_ = eachcol(
    filter!(
        an_ ->
            ImmuneReceptor.is_t_gene(an_[1]) &&
                ImmuneReceptor.is_t_gene(an_[2]) &&
                ImmuneReceptor.is_cdr3(an_[3]),
        Nucleus.Table.rea(
            joinpath(ImmuneReceptor.IN, "su", "41586_2017_BFnature22976_MOESM2_ESM.xlsx"),
            "Raw",
        )[
            !,
            ["TCRBV", "TCRBJ", "CDR-H3"],
        ],
    ),
)

# ---- #

for (nd, (vg_, jg_, U)) in
    ((1, ImmuneReceptor.make_vj(V1_, J1_)), (2, ImmuneReceptor.make_vj(V2_, J2_)))

    ImmuneReceptor.write_vj(joinpath(ImmuneReceptor.OU, "$nd.vj.html"), vg_, jg_, U)

end

# ---- #

Nucleus.Plotly.writ(
    joinpath(ImmuneReceptor.OU, "length.html"),
    (
        Dict(
            "name" => 1,
            "type" => "histogram",
            "x" => map(lastindex, rand(C1_, lastindex(C2_))),
        ),
        Dict("name" => 2, "type" => "histogram", "x" => map(lastindex, C2_)),
    ),
    Dict(
        "yaxis" => Dict("title" => Dict("text" => "Count")),
        "xaxis" => Dict("title" => Dict("text" => "Length")),
    ),
)

# ---- #

const U1_ = unique(C1_)

const U2_ = unique(C2_)

# ---- #

const D1_ = ImmuneReceptor.make_distance(U1_)

const D2_ = ImmuneReceptor.make_distance(U2_)

# ---- #

Nucleus.Plotly.writ(
    joinpath(ImmuneReceptor.OU, "distance.html"),
    (
        Dict("name" => 1, "type" => "histogram", "x" => rand(D1_, lastindex(D2_))),
        Dict("name" => 2, "type" => "histogram", "x" => D2_),
    ),
    Dict(
        "yaxis" => Dict("title" => Dict("text" => "Count")),
        "xaxis" => Dict("title" => Dict("text" => "Distance")),
    ),
)

# ---- #

const M1_ = ImmuneReceptor.get_motif(C1_, 2, 3)

const M2_ = ImmuneReceptor.get_motif(C2_, 2, 3)

# ---- #


