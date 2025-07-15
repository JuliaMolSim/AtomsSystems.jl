using Chemfiles
using AtomsSystems
using AtomsSystems.AtomsTrajectories
using AtomsBase
using Unitful
using Test


@testset "Chemfiles extension tests" begin
    @testset "Systems" begin
        traj = Chemfiles.Trajectory("test_data/methanol.xyz")
        frame = read(traj)  
        close(traj) 
        
        
        sys = SimpleSystem(frame)
        
        # Check the number of atoms
        @test length(sys) == length(frame)
        
        # Check species and positions
        @test all( zip(sys, frame)  ) do (atom, atom_frame)
            AtomsBase.atomic_number(atom) == Chemfiles.atomic_number(atom_frame)
        end
        
        new_pos = ustrip.( u"Å", position_as_matrix(sys, :))
        @test new_pos ≈ Chemfiles.positions(frame)
        
        # Test for cell
        set_cell!(frame, UnitCell([10.0, 11.0, 12.0], [75., 80., 85.]))
        sys = CellSystem(frame)
        
        new_cell = ustrip.( u"Å", cell_matrix(sys) )
        @test new_cell ≈ Chemfiles.matrix(Chemfiles.UnitCell(frame))' # Transpose to match the Chemfiles matrix shape

        # Check the number of atoms
        @test length(sys) == length(frame)
        
        # Check species and positions
        @test all( zip(sys, frame)  ) do (atom, atom_frame)
            AtomsBase.atomic_number(atom) == Chemfiles.atomic_number(atom_frame)
        end
        
        new_pos = ustrip.( u"Å", position_as_matrix(sys, :))
        @test new_pos ≈ Chemfiles.positions(frame)
    end
    @testset "Trajectories" begin
        traj_name = "test_data/methanol-traj.xyz"
        ctraj = ConstantVolumeTrajectory(traj_name)
        @test length(ctraj) == 2
        @test length(ctraj[1]) == 6 

        vtraj = VariableVolumeTrajectory(traj_name)
        @test length(vtraj) == 2
        @test length(vtraj[1]) == 6 
    end
end

