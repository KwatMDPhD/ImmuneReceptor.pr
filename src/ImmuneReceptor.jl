module ImmuneReceptor

const PK = pkgdir(ImmuneReceptor)

const IN = joinpath(PK, "in")

const OU = joinpath(PK, "ou")

# ----------------------------------------------------------------------------------------------- #

using ProgressMeter: @showprogress

function make(s1, s2)

    sum(c1 != c2 for (c1, c2) in zip(s1, s2))

end

for (s1, s2, re) in (
    ("AA", "AA", 0),
    ("AA", "BA", 1),
    ("AA", "AB", 1),
    ("AAA", "AAA", 0),
    ("AAA", "BAA", 1),
    ("AAA", "BBA", 2),
    ("AAA", "BBB", 3),
    ("AAA", "ABB", 2),
    ("AAA", "AAB", 1),
)

    @assert make(s1, s2) === re

end

function make(un_)

    # TODO: Preallocate
    ha_ = Int[]

    @showprogress for s1 in un_, s2 in un_

        if s1 == s2 || lastindex(s1) != lastindex(s2)

            continue

        end

        push!(ha_, make(s1, s2))

    end

    ha_

end

end
