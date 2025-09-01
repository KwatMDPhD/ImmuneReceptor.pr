using Nucleus

using ImmuneReceptor

# ---- #

const A = Nucleus.Table.rea(joinpath(ImmuneReceptor.IN, "vdjdb.human_tcrb.tsv.gz"))

# ---- #

unique!(A[!, "MHC B"])
