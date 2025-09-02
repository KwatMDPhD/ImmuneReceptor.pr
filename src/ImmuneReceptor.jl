module ImmuneReceptor

const PK = pkgdir(ImmuneReceptor)

const IN = joinpath(PK, "in")

const OU = joinpath(PK, "ou")

# ----------------------------------------------------------------------------------------------- #

using ProgressMeter: @showprogress

using Nucleus

# =============================================================================================== #
# Read
# =============================================================================================== #

function is_t_gene(st::AbstractString)

    startswith(st, 'T')

end

function is_t_gene(::Any)

    false

end

function is_cdr3(st)

    st[1] === 'C' && st[end] === 'F'

end

# =============================================================================================== #
# VJ
# =============================================================================================== #

function make_vj(vg_, jg_)

    u1_ = unique(vg_)

    u2_ = unique(jg_)

    d1 = Dict(un => nd for (nd, un) in enumerate(u1_))

    d2 = Dict(un => nd for (nd, un) in enumerate(u2_))

    U = zeros(Int, lastindex(u1_), lastindex(u2_))

    for nd in eachindex(vg_)

        U[d1[vg_[nd]], d2[jg_[nd]]] += 1

    end

    u1_, u2_, U

end

function write_vj(ht, vg_, jg_, U)

    Nucleus.HeatPlot.writ(
        ht,
        vg_,
        jg_,
        U,
        Dict(
            # TODO: Size
            "yaxis" => Dict("title" => Dict("text" => "V gene")),
            "xaxis" => Dict("title" => Dict("text" => "J gene")),
        ),
    )

end

# =============================================================================================== #
# CDR3
# =============================================================================================== #

function make_hamming_distance(s1, s2)

    sum(c1 != c2 for (c1, c2) in zip(s1, s2))

end

function make_todo_distance(s1, s2)

end

function make_distance(st_)

    um = lastindex(st_)

    di_ = Vector{Int}(undef, div(um * (um - 1), 2))

    i1 = 0

    @showprogress for i2 in 1:um, i3 in (i2 + 1):um

        s1 = st_[i2]

        s2 = st_[i3]

        if lastindex(s1) == lastindex(s2)

            di_[i1 += 1] = make_hamming_distance(s1, s2)

        end

    end

    resize!(di_, i1)

end

# =============================================================================================== #
# Motif
# =============================================================================================== #

function get_motif(c1::AbstractString, u1)

    mo_ = String[]

    c2 = c1[4:(end - 3)]

    u2 = lastindex(c2)

    i1 = 0

    i2 = i1 + u1 - 1

    while i2 < u2

        push!(mo_, c2[(i1 += 1):(i2 += 1)])

    end

    mo_

end

function get_motif(mo__, um)

    di = Dict{String, Vector{Int}}()

    for nd in eachindex(mo__)

        for mo in mo__[nd]

            if !haskey(di, mo)

                di[mo] = Int[]

            end

            push!(di[mo], nd)

        end

    end

    [st for (st, in_) in di if um <= lastindex(in_)]

end

end
