""" 
    Parameter(var_name::String, var_of_t::Bool, var_of_ensemble::Bool, var_ensemble_index::Int=0; var_val::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="")

Parameter is a struct that represents a variable in the state space, and information of how to access and print it.
"""
mutable struct Parameter
    var_name::String
    var_str_fun::Function
    var_latex_fun::Function
    var_val::Union{Nothing,Number,Vector{Number},Function}
    var_of_t::Bool
    var_of_ensemble::Bool
    var_indexes::
    var_ensemble_index::Int
    var_subspace_indexes::Vector{Int}

    function Parameter(var_name::String, var_of_t::Bool, var_of_ensemble::Bool, var_ensemble_index::Int=0; var_val::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="", var_subspace_indexes::Vector{Int}=Int[])
        var_suffix_str = str2sub(var_suffix)
        var_str = ""
        if haskey(var_substitution, var_name)
            var_str *= var_substitution[var_name] * var_suffix_str
        else
            var_str *= var_name * var_suffix_str
        end
        var_suffix_latex = ""
        if length(var_suffix) > 0
            var_suffix_latex *= "_{" * var_suffix * "}"
        end
        var_latex = ""
        if haskey(var_substitution_latex, var_name)
            var_latex = var_substitution_latex[var_name] * var_suffix_latex
        else
            var_latex = var_name * var_suffix_latex
        end
        if var_of_t
            add_index(i) = i > 0 ? "_"*string(i) : ""
            var_str_fun(t_ind::Int=0) = var_str * "(t"*str2sub(t_ind)*")"
            var_latex_fun(t_ind::Int=0) = var_latex * "(t"*add_index(t_ind)*")"
        else
            var_str_fun(t_ind::Int=0) = var_str
            var_latex_fun(t_ind::Int=0) = var_latex 
        end
        var_name_mod = var_name
        if length(var_suffix) > 0
            var_name_mod *= "_" * var_suffix
        end
        new(var_name_mod, var_str_fun, var_latex_fun, var_val, var_of_t, var_of_ensemble, var_ensemble_index, var_subspace_indexes)
    end
end
function Base.show(io::IO, p::Parameter)
    print(io, "par: ", p.var_str_fun())
end

struct ParameterInfo 
    vars::Vector{Parameter} 
    time_var_inds::Vector{Int}
    c_one::CAtom


    function ParameterInfo(vars::Vector{Parameter}) 
        time_var_inds = Int[]
        for i in 1:length(vars)
            if vars[i].var_of_t
                push!(time_var_inds, i)
            end
        end
        return new(vars, time_var_inds) 
    end
end

