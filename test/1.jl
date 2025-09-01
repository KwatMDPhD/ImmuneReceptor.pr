using Nucleus

using ImmuneReceptor

# ---- #

const FA = joinpath(ImmuneReceptor.IN, "gliph", "db", "rubelt-naive-CD4.fa")

const DI = ImmuneReceptor.rea(FA)

# ---- #

const ST_ = collect(keys(DI))

const UM_ = map(lastindex, values(DI))

const IN_ = sortperm(UM_)

Nucleus.Plotly.writ(
    joinpath(ImmuneReceptor.OU, "1.html"),
    (Dict("y" => UM_[IN_], "x" => ST_[IN_]),),
    Dict(
        "title" => Dict("text" => "Title"),
        "yaxis" => Dict("title" => Dict("text" => "Number of TODOs")),
        "xaxis" => Dict("title" => Dict("text" => "Regions")),
    ),
)
