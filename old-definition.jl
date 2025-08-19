    subspaces::Vector{SubSpace}
    subspaceinfo::SubSpaceInfo    # Info object containing references to all the indexing of outer and inner subspaces
    I_op::Vector{Vector{Int}}     
    I_ensemble_op::Vector{Vector{Is}}

    # Abstract operators
    operatortypes::Vector{OperatorType}
    operator_names::Vector{String}
    operatortypes_commutator_mat::Matrix{Bool}

    # Parameter fields:
    vars::Vector{Parameter}
    how_many_by_ensemble::Dict{Int,Int}   # for ensemble subspaces, how many variables do we have 
    where_by_ensemble::Dict{Int,Vector{Vector{Int}}}   # for ensemble subspaces, where are our variables  
    where_by_ensemble_var::Vector{Vector{Vector{Int}}}
    where_by_time::Vector{Int}
    where_const::Vector{Int}
    
    subspace_by_ind::Vector{Int}
    c_one::CAtom