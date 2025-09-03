using ProgressMeter: @showprogress

using Random: seed!

using Nucleus

using ImmuneReceptor

# ---- #

const V1_ = String[]

const J1_ = String[]

const C1_ = String[]

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

        s4, s5, _ = split(s2, ','; limit = 3)

        @assert ImmuneReceptor.is_cdr3(s3)

        push!(V1_, s4)

        push!(J1_, s5)

        push!(C1_, s3)

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

for (nd, (s1_, s2_, U)) in
    ((1, ImmuneReceptor.make_vj(V1_, J1_)), (2, ImmuneReceptor.make_vj(V2_, J2_)))

    Nucleus.HeatPlot.writ(
        joinpath(ImmuneReceptor.OU, "$nd.vj.html"),
        s1_,
        s2_,
        U,
        Dict(
            "yaxis" => Dict("title" => Dict("text" => "V gene")),
            "xaxis" => Dict("title" => Dict("text" => "J gene")),
        ),
    )

end

# ---- #

ImmuneReceptor.writ(
    joinpath(ImmuneReceptor.OU, "length.html"),
    "Length",
    map(lastindex, C1_),
    map(lastindex, C2_),
)

# ---- #

const U1_ = unique(C1_)

const U2_ = unique(C2_)

# ---- #

# TODO: Precompute
const D1_ = ImmuneReceptor.make_distance(U1_)

const D2_ = ImmuneReceptor.make_distance(U2_)

# ---- #

ImmuneReceptor.writ(joinpath(ImmuneReceptor.OU, "distance.html"), "Distance", D1_, D2_)

# ---- #

seed!(20250902)

const UM = 1000

# TODO: Match length
const CD__ = map(_ -> rand(U1_, lastindex(U2_)), 1:UM)

# ---- #

const MO__ = map(cd -> ImmuneReceptor.get_motif(cd, 2), U2_)

const M1_ = reduce(vcat, MO__)

const M2_ = ImmuneReceptor.get_motif(MO__, 3)

const M3_ = String[]

# ---- #

@showprogress for st in M2_

    um = count(==(st), M1_)

    @assert 3 <= um

    um_ = map(CD__) do cd_

        count(==(st), reduce(vcat, map(cd -> ImmuneReceptor.get_motif(cd, 2), cd_)))

    end

    pv = count(>=(um), um_) / UM

    if (iszero(pv) ? 1 / UM : pv) <= 0.01

        push!(M3_, st)

    end

end

@info "" M3_

# ---- #

@showprogress for cd_ in CD__

    di = Dict{String, Int}()

    for st_ in map(cd -> ImmuneReceptor.get_motif(cd, 2), cd_)

        for st in st_

            if !haskey(di, st)

                di[st] = 0

            end

            di[st] += 1

        end

    end

  # TODO: Track di

end

  for st in M2_

      pv = if haskey(di, st)

          1 / UM

      else

          di[st] / UM

      end

      if (iszero(pv) ? 1 / UM : pv) <= 0.01

          push!(M3_, st)

      end

  end

@info "" M3_
