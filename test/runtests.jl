using AtomsBase
using AtomsBaseTesting
using AtomsSystems
using Rotations
using Unitful
using Test


include("Aqua.jl")

@testset "AtomsSystems.jl" begin
    # Write your tests here.
    ref = make_test_system()
    @testset "SimpleSystem" begin
        sys = SimpleSystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( species(sys, :) .=== species( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy] 
    end
    @testset "SimpleVelocitySystem" begin
        sys = SimpleVelocitySystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .=== species( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy]
    end
    @testset "AtomicPropertySystem" begin
        sys = AtomicPropertySystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .=== species( ref.system, :) )
        @test all( mass(sys, :) .== mass( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy]
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
    end
    @testset "CellSystemSystem" begin
        sys = CellSystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .=== species( ref.system, :) )
        @test all( mass(sys, :) .== mass( ref.system, :) )
        @test isa(cell(sys), PeriodicCell)
        @test all( cell_vectors(sys) .≈ ref.cell_vectors )
        @test all( sys[:periodicity] .== ref.periodicity )
        @test all( sys[:cell_vectors] .≈ ref.cell_vectors )
        @test_throws KeyError sys[:dummy]
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
    end
    @testset "GeneralSystem" begin
        sys = generic_system(ref.system)
        @test isa(sys, AtomsSystems.GeneralSystem)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .=== species( ref.system, :) )
        @test all( mass(sys, :) .== mass( ref.system, :) )
        @test isa(cell(sys), PeriodicCell)
        @test all( cell_vectors(sys) .≈ ref.cell_vectors )
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
        @test all( x-> haskey(sys, x), keys(ref.sysprop) )
        @test all( x-> haskey(ref.sysprop, x), keys(sys) )
        @test all( [all(sys[k] .== v) for (k,v) in pairs(ref.sysprop)] )
        @test_throws KeyError sys[:dummy]
    end
    @testset "Utils" begin
        sys = SimpleSystem(ref.system)
        # rotation tests
        q = rand(QuatRotation)
        sys2 = rotate_system(sys, q)
        @test all( i-> position(sys2, i) ≈ q * position(sys, i), 1:length(sys) )
        @test bond_angle(sys, 1, 2, 3) ≈ bond_angle(sys2, 1, 2, 3)
        @test dihedral_angle(sys, 1,2,3,4) ≈ dihedral_angle(sys2, 1,2,3,4)
        @test distance_vector(sys2, 1, 2) ≈ q * distance_vector(sys, 1, 2)

        # translation tests
        cms = center_of_mass(sys)
        sys2 = translate_system(sys, -cms)
        @test all( i-> position(sys2, i) ≈ position(sys, i) - cms, 1:length(sys) )
        @test bond_angle(sys, 1, 2, 3) ≈ bond_angle(sys2, 1, 2, 3)
        @test dihedral_angle(sys, 1,2,3,4) ≈ dihedral_angle(sys2, 1,2,3,4)
        @test distance_vector(sys2, 1, 2) ≈ distance_vector(sys, 1, 2)

        
        sys3 = add_systems(sys, sys2)
        @test isa(sys3, typeof(sys))
        @test length(sys3) == length(sys) + length(sys2)
        
        # now with cell
        sys = generic_system(ref.system)

        fp = fractional_coordinates(sys, :)
        @test length(fp) == length(sys)
        fpm = fractional_coordinates_as_matrix(sys, :)
        @test size(fpm) == (3, length(sys))

        clm = cell_matrix(sys)
        clv = cell_vectors(sys)
        @test size(clm) == (3,3)
        @test all( x -> all(x[1] .≈ x[2]), zip(clv, eachcol(clm)) )
        icell = inv_cell(sys)
        tmp = icell * clm

        # Repeat system
        csys = CellSystem(ref.system)
        sys123 = repeat(csys, (1, 2, 3))
        @test length(sys123) == 6 * length(csys)
        @test all( species( sys123, :) .=== repeat( species(csys,:), 6) )
        c1 = cell_vectors(csys)
        c2 = cell_vectors(sys123)
        @test c1[1] ≈   c1[1]
        @test c2[2] ≈ 2*c1[2]
        @test c2[3] ≈ 3*c1[3]
        
        # wrap coordinates
        pbc = [true, true, true]
        cvec = cell_vectors(ref.system)
        c = PeriodicCell(cell_vectors=cvec, periodicity=pbc)
        csys = CellSystem(SimpleSystem(ref.system), c)
        sys4 = wrap_coordinates!(csys)
    end
    @testset "SimpleAtom" begin
        sys = generic_system(ref.system)
        va = sys[:]
        @test all( k -> k in atomkeys(sys), atomkeys(va) )
        @test all( k -> k in atomkeys(va), atomkeys(sys) )
        @test all( mass(va, :) .≈ mass(sys, :) )
        @test all( species(va, :) .== species(sys, :) )
        @test all( position(va, :) .≈ position(sys, :) )
        @test all( velocity(va, :) .≈ velocity(sys, :) )
        @test hasatomkey(va, :charge)
        @test cell(va) isa IsolatedCell
        @test all( periodicity(va) .== (false, false, false) )
        @test all( cell_vectors(va) .≈ cell_vectors(cell(va)) )
        @test mass(va, 2) ≈ mass(sys, 2)
        @test species(va, 2) === species(sys, 2)
        @test position(va, 2) ≈ position(sys, 2)
        @test velocity(va, 2) ≈ velocity(sys, 2)

        sa = SimpleAtom( va[1] ; mark=4 )
        @test species(sa) === Species(va[1])
        @test position(sa) ≈ position(va[1])
        @test velocity(sa) ≈ velocity(va[1])
        @test mass(sa) == mass(va[1])
        @test haskey(sa, :mark)
        @test sa[:mark] == 4

        at = AtomsBase.Atom( :H, [0.0, 0.0, 0.0]u"Å", charge=-1.0u"q" )
        sa = SimpleAtom(at)
        @test species(sa) === species(at)
        @test position(sa) ≈ position(at)
        @test haskey(sa, :charge)
        @test sa[:charge] == -1.0u"q"

        @test AtomsBase.n_dimensions(sa) == 3

        sp = ChemicalSpecies(:C13; atom_name=:myC)
        sa = SimpleAtom(sp, [1.0, 2.0, 3.0]u"Å")
        @test species(sa) === sp
        @test AtomsBase.atom_name(sa) == :myC
        @test AtomsBase.atomic_symbol(sa) == :C13
        @test AtomsBase.element_symbol(sa) == :C
        @test AtomsBase.atomic_number(sa) == 6
        @test AtomsBase.atomic_mass(sa) ≈ AtomsBase.mass(sp) 

        sa = SimpleAtom( :H, [0.0, 0.0, 0.0]u"Å" )
        @test species(sa) === ChemicalSpecies(:H)
        @test position(sa) ≈ [0.0, 0.0, 0.0]u"Å"
        sa = SimpleAtom( :O => [1.0, 0.0, 0.0]u"Å" )
        @test species(sa) === ChemicalSpecies(:O)
        @test position(sa) ≈ [1.0, 0.0, 0.0]u"Å"
        sa = SimpleAtom( ChemicalSpecies(:C), [1.0, 2.0, 3.0]u"Å", [0.1, 0.2, 0.3]u"Å/s"; mass = 12.0u"u", charge = -1.0u"q" )
        @test species(sa) === ChemicalSpecies(:C)
        @test position(sa) ≈ [1.0, 2.0, 3.0]u"Å"
        @test velocity(sa) ≈ [0.1, 0.2, 0.3]u"Å/s"
        @test mass(sa) == 12.0u"u"
        @test sa[:charge] == -1.0u"q"
        ab = AtomsBase.Atom( :O, [1.0, 0.0, 0.0]u"Å")
        sa2 = SimpleAtom(ab)
        @test species(sa2) === ChemicalSpecies(:O)
        @test position(sa2) ≈ [1.0, 0.0, 0.0]u"Å"
        ab = AtomsBase.Atom(ChemicalSpecies(:C), [1.0, 2.0, 3.0]u"Å", [0.1, 0.2, 0.3]u"Å/s"; mass = 12.0u"u", charge = -1.0u"q")
        sa3 = SimpleAtom(ab)
        @test species(sa3) === ChemicalSpecies(:C)
        @test position(sa3) ≈ [1.0, 2.0, 3.0]u"Å"
        @test velocity(sa3) ≈ [0.1, 0.2, 0.3]u"Å/s"
        @test mass(sa3) == 12.0u"u"
        @test sa3[:charge] == -1.0u"q"
    end

    @testset "Views" begin
        sys = generic_system(ref.system)
        @testset "SimpleSystemView" begin
            sys1 = SimpleSystem(sys)
            sv = SimpleSystemView(sys1, 1:2)
            AtomsBase.set_species!(sv, 2, ChemicalSpecies(:Al))
            translate_system!(sv, [1., 2., 3.]u"Å")
            @test all( species(sv, :) .=== species(sys1, 1:2) )
            @test all( position(sv, :) .≈ position(sys1, 1:2) )
            @test isa(cell(sv), IsolatedCell)
            @test all( sv[:] .== sys1[1:2] )
            @test_throws KeyError sv[:dummy]
        end
        @testset "SimpleVelocitySystemView" begin
            sys1 = SimpleVelocitySystem(sys)
            sv = SimpleVelocitySystemView(sys1, 1:2)
            AtomsBase.set_species!(sv, 2, ChemicalSpecies(:U))
            AtomsBase.set_position!(sv, 1, [1., 2., 3.]u"Å")
            AtomsBase.set_velocity!(sv, 2, [1., 2., 3.]u"Å/s")
            @test all( species(sv, :) .=== species(sys1, 1:2) )
            @test all( position(sv, :) .≈ position(sys1, 1:2) )
            @test all( velocity(sv, :) .≈ velocity(sys1, 1:2) )
            @test isa(cell(sv), IsolatedCell) 
            @test all( sv[:] .== sys1[1:2] )
            @test_throws KeyError sv[:dummy]
        end
        @testset "AtomicPropertySystemView" begin
            ap = AtomicPropertySystem(sys)
            av = AtomicPropertySystemView(ap, 1:2)
            @test all( species(av, :) .=== species(ap, 1:2) )
            @test all( position(av, :) .≈ position(ap, 1:2) )
            @test all( velocity(av, :) .≈ velocity(ap, 1:2) )
            @test isa(cell(av), IsolatedCell)
            @test all( av[:] .== ap[1:2] ) 
            @test_throws KeyError av[:dummy]
        end
        @testset "CellSystemView" begin
            cs = CellSystem(sys)
            cv = CellSystemView(cs, 1:2)
            @test all( species(cv, :) .=== species(sys, 1:2) )
            @test all( position(cv, :) .≈ position(sys, 1:2) )
            @test all( velocity(cv, :) .≈ velocity(sys, 1:2) )
            @test isa(cell(cv), PeriodicCell)
            @test all( cell_vectors(cv) .≈ cell_vectors(cs) )
            @test all( periodicity(cv) .== periodicity(cs) )
            @test all( cv[:] .== cs[1:2] )
            @test_throws KeyError cv[:dummy]
        end
    end

end
