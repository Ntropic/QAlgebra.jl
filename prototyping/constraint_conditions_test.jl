using QAlgebra
import QAlgebra.QExpressions: ConstrainedIndexBlock, neq_expand
using QAlgebra.QSpaces: SubSpaceIndex

function demo_block()
    indexes = SubSpaceIndex[
        SubSpaceIndex(1, 2, 102),
        SubSpaceIndex(1, 4, 104),
    ]
    constraints = BitVector[
        BitVector([false, true, true, true]),
        BitVector([false, true, false, true]),
    ]
    return ConstrainedIndexBlock(1, 4, 1, indexes, constraints)
end

function main()
    block = demo_block()
    where_defined = BitVector([true, false, true, false])

    println("Original block:")
    for (idx, row) in zip(block.indexes, block.constraints)
        println("  index expanded=$(idx.expanded) -> $(collect(row))")
    end
    println("where_defined = $(collect(where_defined))\n")

    branches = neq_expand(block, where_defined)
    println("Generated $(length(branches)) inequality branches:\n")
    for (i, (blk, defined, actions)) in enumerate(branches)
        println("Branch $i:")
        println("  where_defined = $(collect(defined))")
        for (removed, col) in actions
            println("  equality: $(removed.expanded) -> column $(col)")
        end
        for (idx, row) in zip(blk.indexes, blk.constraints)
            println("  index expanded=$(idx.expanded) -> $(collect(row))")
        end
        println()
    end
end

main()
