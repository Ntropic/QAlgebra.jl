@inline function depends_on_time(cfun::CFunction, param_info::ParameterInfo)::Bool
    return depends_on_inds(cfun, param_info.time_var_inds)  
end
@inline function depends_on_time(qfun::QAtomProduct, param_info::ParameterInfo)::Bool 
    return depends_on_time(qfun.coeff, param_info)
end
@inline function depends_on_time(qfun::QFunction, param_info::ParameterInfo)::Bool 
    return depends_on_time.(each_coeff(qfun), Ref(param_info)) 
end
@inline function depends_on_time(qfun::QAtom, param_info::ParameterInfo)::Bool 
    return error("Depends_on_time is by design not implemented for QAtoms. Shouldn#t be called on QAtom.")
end

abstract type TFunction end

struct TAtom 
    coeff_fun::CFunction 
    time_ind::Int 
    time_dependent::Bool 
    param_info::ParameterInfo
    function TAtom(coeff_fun::CFunction, time_ind::Int=0, depends_on_time::Bool=true, param_info::ParameterInfo) 
        new(coeff_fun, time_ind, depends_on_time, param_info)
    end
    function TAtom(coeff_fun::CFunction, time_ind::Int, param_info::ParameterInfo) 
        time_dependent = depends_on_time(coeff_fun, param_info) 
        if !time_dependent
            time_ind = 0
        end
        new(coeff_fun, time_ind, time_dependent, param_info)
    end
end
function TAtom(coeff_fun::CFunction, time_ind::Int, statespace::StateSpace)
    return TAtom(coeff_fun, time_ind, statespace.param_info)
end
function TAtom(coeff_fun::CFunction, time_ind::Int, depends_on_time::Bool=true, statespace::StateSpace)
    return TAtom(coeff_fun, time_ind, depends_on_time, statespace.param_info)
end
depends_on_time(coeff_fun::TAtom, param_info::ParameterInfo)::Bool =  coeff_fun.time_dependent
depends_on_time(coeff_fun::TAtom)::Bool = coeff_fun.time_dependent
modify_coeff(coeff_fun::TAtom, new_coeff_fun::CFunction)::TAtom = TAtom(new_coeff_fun, coeff_fun.new_time_ind, coeff_fun.param_info)
Base.isless(a::TAtom, b::TAtom) = a.time_ind < b.time_ind
function unifiable(t1::TAtom, t2::TAtom)::Tuple{Bool, Int}
    if t1.time_ind == 0 
        return true, t2.time_ind
    elseif t2.time_ind == 0
        return true, t1.time_ind
    else
        return t1.time_ind == t2.time_ind , t1.time_ind 
    end
end

struct TProd <: TFunction
    atoms::Vector{TAtom}
    function TProd(atoms::Vector{TAtom})
        isempty(atoms) && throw(ArgumentError("TProd must have at least one atom"))
        if length(atoms) == 1
            return atoms[1]  # return TAtom directly
        end
        return new(sort(atoms))
    end
    function TProd(atoms::Vector{TAtom}, ::Val{:sorted})
        isempty(atoms) && throw(ArgumentError("TProd must have at least one atom"))
        if length(atoms) == 1
            return atoms[1]  # return TAtom directly
        end
        return new(atoms)
    end
end
depends_on_time(t::TProd)::Bool = true 
depends_on_time(t::TProd, param_info::ParameterInfo)  
time_indices(prod::TProd) = getfield.(prod.atoms, :time_ind)

#### Algebraic operations +, -, *, /, ^, etc. on TFunction 
# Are already defined on CFunction, here we need to check that the elements are unifiable. otherwise we return a 
import Base: +, -, *, /, ^, ==

function *(t1::TAtom, t2::TAtom)::TFunction
    do_unify, new_time_index = unifiable(t1, t2) 
    if do_unify
        coeff = t1.coeff_fun * t2.coeff_fun
        return TAtom(coeff, new_time_index, t1.param_info) 
    else
        return TProd([t1, t2])  # fallback to a product
    end
end
function merge_sorted(a::Vector{TAtom}, b::Vector{TAtom})
    out = Vector{TAtom}()
    sizehint!(out, length(a) + length(b))

    i = j = 1
    while i <= length(a) && j <= length(b)
        if a[i].time_ind == b[j].time_ind
            atom = a[i] * b[j]
            if !iszero(atom.coeff_fun)
                push!(out, atom)
            end
            i += 1
            j += 1
        elseif a[i] < b[j]
            push!(out, a[i]); i += 1
        else
            push!(out, b[j]); j += 1
        end
    end

    # append the rest
    while i <= length(a)
        push!(out, a[i]); i += 1
    end
    while j <= length(b)
        push!(out, b[j]); j += 1
    end

    return out
end
function *(p1::TProd, p2::TProd)::TFunction
    return TProd(merge_sorted(p1.atoms, p2.atoms))
end

function *(p1::TProd, p2::TAtom)::TFunction
    return TProd(merge_sorted(p1.atoms, [p2]))
end
*(p1::TAtom, p2::TProd) = p2 * p1

-(a::TAtom) = modify_coeff(a, -q.coeff_fun) 
-(a::TProd) = modify_coeff(a, -q.coeff_fun) 


