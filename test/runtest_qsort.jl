using Test
using QAlgebra: sorted_push, sorted_push!, sorted_push_unique, sorted_push_unique!, sorted_append, sorted_append!, sorted_append_unique, sorted_append_unique!, sorted_append_unify!, sorted_append_unify, sorted_push_unify_branches, sorted_append_unify_branches
using ComplexRationals

@testset "sorted push / append functions" begin
    @test sorted_push([1,3,5], 4) == [1,3,4,5]
    @test sorted_push([1,3,5], 7) == [1,3,5,7]
    @test sorted_push_unique([1,3,5], 3) == [1,3,5]
    @test sorted_push_unique([1,3,5], 2) == [1,2,3,5]

    v = [1,3,5]; sorted_push!(v,4); @test v == [1,3,4,5]

    @test sorted_append([1,3,5],[2,4,6]) == [1,2,3,4,5,6]

    u = [1,3,5]
    sorted_append_unique!(u, [3,4,5,6])
    @test u == [1,3,4,5,6]

    @test sorted_append_unique([1,3,5],[3,4,5,6]) == [1,3,4,5,6]
end

# Minimal Atom type for testing unify behavior
struct Atom
    idx::Int
    exp::Int
end

Base.isless(a::Atom, b::Atom) = a.idx < b.idx
Base.:(==)(a::Atom, b::Atom) = a.idx == b.idx && a.exp == b.exp

# Unifier: combine exponents for equal indices
unify_addexp(a::Atom, b::Atom) = Atom(a.idx, a.exp + b.exp)

@testset "sorted_append_unify! basic merge" begin
    dest = Atom[Atom(1, 1), Atom(3, 2), Atom(6, 1)]
    src  = Atom[Atom(2, 5), Atom(3, 3), Atom(5, 4)]

    sorted_append_unify!(dest, src, unify_addexp; lt=isless)

    @test dest == Atom[
        Atom(1, 1),
        Atom(2, 5),
        Atom(3, 5),
        Atom(5, 4),
        Atom(6, 1),
    ]
end

@testset "sorted_append_unify (pure version)" begin
    dest = Atom[Atom(1, 1), Atom(3, 2)]
    src  = Atom[Atom(2, 5), Atom(3, 3)]

    dest_orig = copy(dest)
    src_orig  = copy(src)

    out = sorted_append_unify(dest, src, unify_addexp; lt=isless)

    @test dest == dest_orig
    @test src == src_orig

    @test out == Atom[
        Atom(1, 1),
        Atom(2, 5),
        Atom(3, 5),
    ]

    @test out !== dest
    @test out !== src
end

branch_unify(a::Atom, b::Atom) = [
    (ComplexRational(1, 0, 2), Atom(a.idx, a.exp + b.exp)),
    (ComplexRational(1, 0, 2), Atom(a.idx, a.exp - b.exp))
]

@testset "sorted_push_unify_branches" begin
    base = Atom[Atom(2, 3)]
    branches = sorted_push_unify_branches(base, Atom(2, 1); lt=isless, unify=branch_unify)
    @test length(branches) == 2
    weights = first.(branches)
    vectors = last.(branches)
    @test all(w -> w == ComplexRational(1, 0, 2), weights)
    @test Atom(2, 4) in (v[1] for v in vectors)
    @test Atom(2, 2) in (v[1] for v in vectors)
end

@testset "sorted_append_unify_branches" begin
    dest = Atom[Atom(1, 1), Atom(3, 2)]
    src  = Atom[Atom(2, 5), Atom(3, 3)]
    branches = sorted_append_unify_branches(dest, src; lt=isless, unify=branch_unify)
    @test length(branches) == 2
    weights = first.(branches)
    vectors = last.(branches)
    @test all(w -> w == ComplexRational(1, 0, 2), weights)
    expected = [
        Atom[Atom(1,1), Atom(2,5), Atom(3,5)],
        Atom[Atom(1,1), Atom(2,5), Atom(3,-1)]
    ]
    @test all(vec -> vec in expected, vectors)
end
