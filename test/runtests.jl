using AtomsBase
using AtomsBaseTesting
using AtomsSystems
using LinearAlgebra: diagm, norm
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
        @test haskey(sys, :cell_vectors)
        @test all( x -> in(x, (:position, :species)), atomkeys(sys) )
        @test hasatomkey(sys, :velocity) == false

        ss = SimpleSystem(:H, [0.0, 0.0, 0.0]u"Å")
        @test species(ss,1) === ChemicalSpecies(:H)
        @test position(ss, 1) ≈ [0.0, 0.0, 0.0]u"Å"
        @test length(ss) == 1
        @test element_symbol(ss, 1) == :H

        ss = SimpleSystem(sys, ChemicalSpecies(:H))
        @test length(ss) == count( species(sys, :) .== ChemicalSpecies(:H) )
        @test all( species(ss, :) .== ChemicalSpecies(:H) )

        ss = SimpleSystem(sys, ChemicalSpecies(:H), ChemicalSpecies(:He))
        @test length(ss) == count( species(sys, :) .== ChemicalSpecies(:H) ) +
                            count( species(sys, :) .== ChemicalSpecies(:He) )
        
        ss = SimpleSystem(sys)
        deleteat!(ss, 3)
        @test length(ss) == length(sys) - 1
        @test all( position(ss, 3) .≈ position(sys, 4) )
        @test species(ss, 3) === species(sys, 4)
        
    end
    @testset "SimpleVelocitySystem" begin
        sys = SimpleVelocitySystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .=== species( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy]
        @test all( x -> in(x, (:position, :species, :velocity)), atomkeys(sys) )
        @test hasatomkey(sys, :velocity) == true
        AtomsBase.set_velocity!(sys, 1, [1.0, 2.0, 3.0]u"Å/s")
        @test velocity(sys, 1) ≈ [1.0, 2.0, 3.0]u"Å/s"

        sv = SimpleVelocitySystem(:H, [0.0, 0.0, 0.0]u"Å", [1.0, 2.0, 3.0]u"Å/s")
        @test species(sv,1) === ChemicalSpecies(:H)
        @test position(sv, 1) ≈ [0.0, 0.0, 0.0]u"Å"
        @test velocity(sv, 1) ≈ [1.0, 2.0, 3.0]u"Å/s"
        @test length(sv) == 1
        @test element_symbol(sv, 1) == :H

        sv = SimpleVelocitySystem(sys, ChemicalSpecies(:H))
        @test length(sv) == count( species(sys, :) .== ChemicalSpecies(:H) )
        @test all( species(sv, :) .== ChemicalSpecies(:H) )

        sv = SimpleVelocitySystem(sys, ChemicalSpecies(:H), ChemicalSpecies(:He))
        @test length(sv) == count( species(sys, :) .== ChemicalSpecies(:H) ) +
                            count( species(sys, :) .== ChemicalSpecies(:He) )

        sv = SimpleVelocitySystem(sys)
        deleteat!(sv, 3)
        @test length(sv) == length(sys) - 1
        @test all( position(sv, 3) .≈ position(sys, 4) )
        @test all( velocity(sv, 3) .≈ velocity(sys, 4) )
        @test species(sv, 3) === species(sys, 4)

        ss = SimpleSystem(sys)
        sv = generic_system(ss, velocity(sys, :))
        @test all( position(sv, :) .≈ position(sys, :) )
        @test all( velocity(sv, :) .≈ velocity(sys, :) )
        @test all( species(sv, :) .=== species(sys, :) )   
    end
    @testset "AtomicPropertySystem" begin
        sys = AtomicPropertySystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .=== species( ref.system, :) )
        @test all( mass(sys, :) .≈ mass( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy]
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        @test keys(sys) == (:cell_vectors, :periodicity)
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end

        ss = AtomicPropertySystem(sys[:])
        @test all( position(ss, :) .≈ position(sys, :) )
        @test all( velocity(ss, :) .≈ velocity(sys, :) )
        @test all( species(ss, :) .=== species(sys, :) )
        @test all( mass(ss, :) .≈ mass(sys, :) )

        sv = SimpleVelocitySystem(sys)
        ss = AtomicPropertySystem(sv)
        @test all( position(ss, :) .≈ position(sv, :) ) 
        @test all( velocity(ss, :) .≈ velocity(sv, :) )
        @test all( species(ss, :) .=== species(sv, :) )

        ss = AtomicPropertySystem(sv, 1:2)
        @test all( position(ss, :) .≈ position(sv, 1:2) )
        @test all( velocity(ss, :) .≈ velocity(sv, 1:2) )
        @test all( species(ss, :) .=== species(sv, 1:2) )

        ss = AtomicPropertySystem(sv, ChemicalSpecies(:H))
        @test length(ss) == count( species(sv, :) .== ChemicalSpecies(:H) ) 
        @test all( species(ss, :) .=== ChemicalSpecies(:H) )
        
        gsys = generic_system(ref.system)
        ss = AtomicPropertySystem(gsys, ChemicalSpecies(:H))
        @test length(ss) == count( species(gsys, :) .== ChemicalSpecies(:H) )
        @test all( species(ss, :) .=== ChemicalSpecies(:H) )

        ss = AtomicPropertySystem(gsys, ChemicalSpecies(:H), ChemicalSpecies(:He))  
        @test length(ss) == count( species(gsys, :) .== ChemicalSpecies(:H) ) +
                            count( species(gsys, :) .== ChemicalSpecies(:He) )
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


        sc = CellSystem(sys, IsolatedCell(3))
        @test isa(sc, AtomsSystems.AbstractIsolatedSystem)

        sc = CellSystem(sys[1:2])
        @test isa(sc, AtomsSystems.AbstractIsolatedSystem)
        @test all( position(sc, :) .≈ position(sys, 1:2) )
        @test all( species(sc, :) .=== species(sys, 1:2) )
        @test all( mass(sc, :) .≈ mass(sys, 1:2) )

        sc = CellSystem(sys, ChemicalSpecies(:H))
        @test length(sc) == count( species(sys, :) .== ChemicalSpecies(:H) )
        @test all( species(sc, :) .=== ChemicalSpecies(:H) )

        sc = CellSystem(sys, ChemicalSpecies(:H), ChemicalSpecies(:He))
        @test length(sc) == count( species(sys, :) .== ChemicalSpecies(:H) ) +
                            count( species(sys, :) .== ChemicalSpecies(:He) )

        sc = CellSystem(sys)
        translate_system!(sc, [1.0, 2.0, 3.0]u"Å")
        ts = add_systems(sys, sc)
        length(ts) == 2*length(sys)
        @test all( position(ts, length(sys)+1:2*length(sys)) .≈ position(sc, :) )

        sc = CellSystem(sys)
        AtomsBase.set_periodicity!(sc, (true, true, true))
        @test all( periodicity(sc) .== (true, true, true) )
        AtomsBase.set_cell_vectors!(sc, [1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å")
        cv = cell_vectors(sc)
        @test cv[1] ≈ [1.0, 0.0, 0.0]u"Å"
        @test cv[2] ≈ [0.0, 1.0, 0.0]u"Å"
        @test cv[3] ≈ [0.0, 0.0, 1.0]u"Å"


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

        sys = generic_system"""
            H 0.0 0.0 0.0
            He 1.0 2.0 3.0
            O 4.0 5.0 6.0
        """
        @test length(sys) == 3
        @test species(sys, 1) === ChemicalSpecies(:H)
        @test position(sys, 1) ≈ [0.0, 0.0, 0.0]u"Å"
        @test species(sys, 2) === ChemicalSpecies(:He)
        @test position(sys, 2) ≈ [1.0, 2.0, 3.0]u"Å"
        @test species(sys, 3) === ChemicalSpecies(:O)
        @test position(sys, 3) ≈ [4.0, 5.0, 6.0]u"Å"

        sys = generic_system([:H => [0.0, 0.0, 0.0]u"Å", :O => [1.0, 0.0, 0.0]u"Å"]; energy = 10.0u"eV")
        @test length(sys) == 2
        @test species(sys, 1) === ChemicalSpecies(:H)
        @test position(sys, 1) ≈ [0.0, 0.0, 0.0]u"Å"
        @test species(sys, 2) === ChemicalSpecies(:O)
        @test position(sys, 2) ≈ [1.0, 0.0, 0.0]u"Å"
        @test sys[:energy] ≈ 10.0u"eV"

        rsys = generic_system(ref.system) 
        sys = generic_system(rsys, 1:3; energy = 10.0u"eV")
        @test length(sys) == 3
        @test all( position(sys, :) .≈ position(rsys, 1:3) )
        @test all( species(sys, :) .=== species(rsys, 1:3) )
        @test all( mass(sys, :) .≈ mass(rsys, 1:3) )
        @test sys[:energy] ≈ 10.0u"eV"

        sys = generic_system(rsys, ChemicalSpecies(:H), ChemicalSpecies(:He))
        @test length(sys) == count( species(rsys, :) .== ChemicalSpecies(:H) ) +
                            count( species(rsys, :) .== ChemicalSpecies(:He) )  
        
        sys = generic_system(species(rsys, :), position(rsys, :), velocity(rsys, :))
        @test length(sys) == length(rsys)
        @test all( position(sys, :) .≈ position(rsys, :) )
        @test all( velocity(sys, :) .≈ velocity(rsys, :) )
        @test all( species(sys, :) .=== species(rsys, :) )

        sys = generic_system(species(rsys, :), position(rsys, :))
        @test length(sys) == length(rsys)
        @test all( position(sys, :) .≈ position(rsys, :) )
        @test all( species(sys, :) .=== species(rsys, :) )

        rsys = generic_system(ref.system)
        ss = generic_system(ref.system, velocity(rsys, :))
        @test all( position(ss, :) .≈ position(rsys, :) )
        @test all( velocity(ss, :) .≈ velocity(rsys, :) )
        @test all( species(ss, :) .=== species(rsys, :) )
        @test Set( atomkeys(ss) ) == Set( atomkeys(rsys) )
    end
    @testset "Get Started test" begin
        sys = generic_system"""
            O     -2.1   0.6    0.0
            H     -1.4   0.4    0.6
            H     -1.8   1.3   -0.6
        """
        sys = generic_system(sys; cell_vectors = [
            [5., 0., 0.]u"Å",
            [0., 5., 0.]u"Å",
            [0., 0., 5.]u"Å"],
        )
        sys = generic_system(sys; periodicity=(false, true, false))
        sys = generic_system(sys; energy=1.2u"eV", label="my water") 
        
        @test length(sys) == 3
        @test species(sys, 1) === ChemicalSpecies(:O)
        @test species(sys, 2) === ChemicalSpecies(:H)
        @test species(sys, 3) === ChemicalSpecies(:H)

        @test periodicity(sys) == (false, true, false)
        @test all( cell_vectors(sys) .≈ [
            [5.0, 0.0, 0.0]u"Å",
            [0.0, 5.0, 0.0]u"Å",
            [0.0, 0.0, 5.0]u"Å",
        ] )
        @test sys[:energy] ≈ 1.2u"eV"
        @test sys[:label] == "my water"


        atoms = [
            SimpleAtom(:O, [-2.1, 0.6, 0.0]u"Å"; mass=16.5u"u" )
            SimpleAtom(:H, [-1.4, 0.4, 0.6]u"Å"; mass=2.3u"u"  )
            SimpleAtom(:H, [-1.8, 1.3, -0.6]u"Å"; mass=3.3u"u" )
        ]
        atoms = [
            SimpleAtom( atoms[1]; velocity=[0.1, 0.2, 0.0]u"Å/fs" , charge=-1.0u"q" )
            SimpleAtom( atoms[2]; velocity=[-0.2, 0.0, 0.1]u"Å/fs", charge=1.0u"q"  )
            SimpleAtom( atoms[3]; velocity=[0.0, -0.1, 0.2]u"Å/fs", charge=0.0u"q"  )
        ]
        sys = generic_system(atoms)
        @test length(sys) == 3
        @test species(sys, 1) === ChemicalSpecies(:O)
        @test species(sys, 2) === ChemicalSpecies(:H)
        @test species(sys, 3) === ChemicalSpecies(:H)
        @test Set( atomkeys(sys) ) == Set( (:position, :species, :velocity, :mass, :charge) )
        @test Set( atomkeys(sys) ) == Set( atomkeys(atoms) )

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
        gsys = generic_system(sys; label="test rotation")
        rotate_system!(gsys, q)
        @test all( position(gsys, :) .≈ position(sys2, :) )


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
        @test fractional_coordinates_as_matrix(sys, 2) ≈ fractional_coordinates(sys, 2)

        clm = cell_matrix(sys)
        clv = cell_vectors(sys)
        @test size(clm) == (3,3)
        @test all( x -> all(x[1] .≈ x[2]), zip(clv, eachcol(clm)) )
        icell = inv_cell(sys)
        tmp = icell * clm

        # isolated cell
        sys = SimpleSystem(ref.system)
        icell = inv_cell(sys)
        @test diagm(ones(3))*unit(icell[1,1]) == icell
        @test fractional_coordinates( cell(sys), position(sys, 1) ) ≈ fractional_coordinates(sys, 1)  
        @test fractional_coordinates_as_matrix( cell(sys), position(sys, :) ) ≈ fractional_coordinates_as_matrix(sys, :)
        

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
        sys333 = repeat(csys, 3)
        c2 = cell_vectors(sys333)
        @test c2[1] ≈ 3*c1[1]
        @test c2[2] ≈ 3*c1[2]
        @test c2[3] ≈ 3*c1[3]
        
        # wrap coordinates
        pbc = [true, true, true]
        cvec = cell_vectors(ref.system)
        c = PeriodicCell(cell_vectors=cvec, periodicity=pbc)
        csys = CellSystem(SimpleSystem(ref.system), c)
        sys4 = wrap_coordinates!(csys)

        # isolated cell wraps should not change coordinates
        ssys = SimpleSystem(ref.system)
        sys5 = wrap_coordinates!(ssys)
        @test all( position(ssys, :) .≈ position(sys5, :) )
        @test wrap_coordinates!( cell(ssys), position(ssys, 1) ) ≈ position(ssys, 1)

        # distances
        sys = SimpleSystem(ref.system)
        @test distance_vector(sys, 1, 2) ≈ position(sys, 2) - position(sys, 1)
        distance(sys, 1, 2) ≈ norm( position(sys, 2) - position(sys, 1) )
        dis = distance(sys, sys)
        @test dis[1, 2] ≈ distance(sys, 1, 2)
        @test dis[2, 3] ≈ distance(sys[1:3], 2, 3)
        @test dis[3, 1] ≈ distance(sys[1:3], 3, 1)

        sys = generic_system(ref.system)
        dis = distance(sys, 1, :)
        @test dis[2] ≈ distance(sys, 1, 2)
        @test dis[3] ≈ distance(sys, 1, 3)
        dis2 = distance(sys, sys)
        @test dis2[1, 3] ≈ distance(sys, 1, 3)
        @test dis2[2, 4] ≈ distance(sys, 2, 4)
         
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
        @test species(sa) === species(va[1])
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
        @test AtomsBase.mass(sa) ≈ AtomsBase.mass(sp) 

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
            @test Set( atomkeys(sv) ) == Set( (:species, :position) )
            
            sv = system_view(sys1, 1:4)
            sv1 = system_view(sv, 1:2)
            @test all( species(sv1, :) .=== species(sv, 1:2) )
            @test all( position(sv1, :) .≈ position(sv, 1:2) )   

            sv = system_view(sys1, 2)
            @test species(sv, 1) === species(sys1, 2)
            @test position(sv, 1) ≈ position(sys1, 2)
            @test length(sv) == 1
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
            @test Set( atomkeys(sv) ) == Set( (:species, :position, :velocity) )

            sv = system_view(sys1, 1:4)
            sv1 = system_view(sv, 1:2)
            @test all( species(sv1, :) .=== species(sv, 1:2) )
            @test all( position(sv1, :) .≈ position(sv, 1:2) )
            @test all( velocity(sv1, :) .≈ velocity(sv, 1:2) )

            sv = system_view(sys1, 2)
            @test species(sv, 1) === species(sys1, 2)
            @test position(sv, 1) ≈ position(sys1, 2)
            @test velocity(sv, 1) ≈ velocity(sys1, 2)
            @test length(sv) == 1
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
            @test Set( atomkeys(av) ) == Set( atomkeys(ap) )

            av = system_view(ap, 1:4)
            av1 = system_view(av, 1:2)
            @test all( species(av1, :) .=== species(av, 1:2) )
            @test all( position(av1, :) .≈ position(av, 1:2) )
            @test all( velocity(av1, :) .≈ velocity(av, 1:2) )
            @test all( mass(av1, :) .≈ mass(av, 1:2) )

            av = system_view(ap, 2)
            @test species(av, 1) === species(ap, 2)
            @test position(av, 1) ≈ position(ap, 2)
            @test velocity(av, 1) ≈ velocity(ap, 2)
            @test mass(av, 1) ≈ mass(ap, 2)
            @test length(av) == 1
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

            cv = system_view(cs, 1:4)
            cv1 = system_view(cv, 1:2)
            @test all( species(cv1, :) .=== species(cv, 1:2) )
            @test all( position(cv1, :) .≈ position(cv, 1:2) )
            @test all( velocity(cv1, :) .≈ velocity(cv, 1:2) )
            @test all( mass(cv1, :) .≈ mass(cv, 1:2) )
            @test all( cell_vectors(cv1) .≈ cell_vectors(cv) )
            @test all( periodicity(cv1) .== periodicity(cv) )

            cv = system_view(cs, ChemicalSpecies(:H))
            @test length(cv) == count( species(cs, :) .== ChemicalSpecies(:H) )
            @test all( species(cv, :) .=== ChemicalSpecies(:H) )
            @test cell(cv) == cell(cs)

            cv = system_view(cs, 2)
            @test species(cv, 1) === species(cs, 2)
            @test position(cv, 1) ≈ position(cs, 2)
            @test velocity(cv, 1) ≈ velocity(cs, 2)
            @test mass(cv, 1) ≈ mass(cs, 2)
            @test cell(cv) == cell(cs)
            @test length(cv) == 1
        end
    end

end
