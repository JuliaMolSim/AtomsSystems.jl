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
    @test AtomsSystems.distance(trj[1], 1, 2) ≈ AtomsSystems.distance(sys, 1, 2)
    @test AtomsSystems.distance(trj[2], 1, 2) ≈ AtomsSystems.distance(sys2, 1, 2)
    @test AtomsSystems.bond_angle(trj[1], 1, 2, 3) ≈ AtomsSystems.bond_angle(sys, 1, 2, 3)
    @test AtomsSystems.bond_angle(trj[2], 1, 2, 3) ≈ AtomsSystems.bond_angle(sys2, 1, 2, 3)
    @test AtomsSystems.dihedral_angle(trj[1], 1, 2, 3, 4) ≈ AtomsSystems.dihedral_angle(sys, 1, 2, 3, 4)
    @test AtomsSystems.dihedral_angle(trj[2], 1, 2, 3, 4) ≈ AtomsSystems.dihedral_angle(sys2, 1, 2, 3, 4)

    traj = VariableVolumeTrajectory(sys)
    push!(traj, sys2)
    @test length(traj) == 2
    @test all( species(traj[1], :) .=== species(sys, :) )
    @test all( position(traj[1], :) .≈ position(sys, :) )
    @test all( position(traj[2], :) .≈ position(sys2, :))
    @test all( velocity(traj[1], :) .≈ velocity(sys, :) )
    @test all( velocity(traj[2], :) .≈ velocity(sys2, :))
    @test cell(traj[1]) == cell(sys)
    @test cell(traj[2]) == cell(sys2)
    @test AtomsSystems.distance(traj[1], 1, 2) ≈ AtomsSystems.distance(sys, 1, 2)
    @test AtomsSystems.distance(traj[2], 1, 2) ≈ AtomsSystems.distance(sys2, 1, 2)
    @test AtomsSystems.bond_angle(traj[1], 1, 2, 3) ≈ AtomsSystems.bond_angle(sys, 1, 2, 3)
    @test AtomsSystems.bond_angle(traj[2], 1, 2, 3) ≈ AtomsSystems.bond_angle(sys2, 1, 2, 3)
    @test AtomsSystems.dihedral_angle(traj[1], 1, 2, 3, 4) ≈ AtomsSystems.dihedral_angle(sys, 1, 2, 3, 4)
    @test AtomsSystems.dihedral_angle(traj[2], 1, 2, 3, 4) ≈ AtomsSystems.dihedral_angle(sys2, 1, 2, 3, 4)  
end