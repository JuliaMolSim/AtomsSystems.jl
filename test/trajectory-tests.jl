using AtomsBase
using AtomsBaseTesting
using AtomsSystems
using AtomsSystems.AtomsTrajectories
using Unitful
using Test



@testset "Trajectory" begin
    ref = make_test_system()
    sys = generic_system(ref.system)
    sys2 = deepcopy( sys )
    translate_system!(sys2, [1., 2., 3.]u"Å")
    trj = ConstantVolumeTrajectory([sys, sys2])
    @test length(trj) == 2
    @test all( species(trj[1], :) .=== species(sys, :) )
    @test all( position(trj[1], :) .≈ position(sys, :) )
    @test all( position(trj[2], :) .≈ position(sys2, :) )
    @test all( velocity(trj[1], :) .≈ velocity(sys, :) )
    @test all( velocity(trj[2], :) .≈ velocity(sys2, :) )
    @test cell(trj[1]) == cell(sys)
    @test cell(trj[1]) == cell(sys2) 
    @test all( AtomsSystems.distance(trj, 1, 2, 1:2) .≈ AtomsSystems.distance(trj, 1, 2, :) )
    @test all( bond_angle(trj, 1, 2, 3, 1:2) .≈ bond_angle(trj, 1, 2, 3, :) )
    @test all( dihedral_angle(trj, 1, 2, 3, 4, 1:2) .≈ dihedral_angle(trj, 1, 2, 3, 4, :) )
end