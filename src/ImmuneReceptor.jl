module ImmuneReceptor

const PK = pkgdir(ImmuneReceptor)

const IN = joinpath(PK, "in")

const OU = joinpath(PK, "ou")

# ----------------------------------------------------------------------------------------------- #

function ge(fa)

    s1, s2 = split(fa, '\n')

    s3 = split(s1, ';'; limit = 2)[1]

    s4, s5 = rsplit(s3, ','; limit = 2)

    @assert s5 == s2

    s4, s5

end

function rea(fa)

    di = Dict{String, Vector{String}}()

    for s1 in eachsplit(read(fa, String), '>'; keepempty = false)

        s2, s3 = ge(s1)

        if !haskey(di, s2)

            di[s2] = String[]

        end

        push!(di[s2], s3)

    end

    di

end

end
