using QAlgebra: sorted_push, sorted_push!, sorted_push_unique, sorted_push_unique!, sorted_append, sorted_append!, sorted_append_unique, sorted_append_unique!

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

