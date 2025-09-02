using Nucleus

using ImmuneReceptor

# ---- #

function is_t(an::AbstractString)

    startswith(an, 'T')

end

function is_t(::Any)

    false

end

function is_c(st)

    st[1] === 'C' && st[end] === 'F'

end

# ---- #

const A = filter!(
    an_ -> is_t(an_[1]) && is_t(an_[2]) && is_c(an_[3]),
    Nucleus.Table.rea(
        joinpath(ImmuneReceptor.IN, "su", "41586_2017_BFnature22976_MOESM2_ESM.xlsx"),
        "Raw",
    )[
        !,
        ["TCRBV", "TCRBJ", "CDR-H3"],
    ],
)

# ---- #

const S1_::Vector{String}, S2_::Vector{String}, S3_ = eachcol(A)

# ---- #

const UN_ = unique(S3_)

# ---- #

const HA_ = ImmuneReceptor.make(UN_)

Nucleus.Plotly.writ(
    joinpath(ImmuneReceptor.OU, "4.html"),
    (Dict("type" => "histogram", "x" => HA_),),
    Dict(
        "yaxis" => Dict("title" => Dict("text" => "Count")),
        "xaxis" => Dict("title" => Dict("text" => "Hamming distance")),
    ),
)

# ---- #
# Do local

# ---- #
# Group

# ---- #
# Score

unique(S3_)
