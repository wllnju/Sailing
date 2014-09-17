//#####################################################################
// This file is part of code for Sailing project which was built in PhysBAM
// whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Parallel_Computation/MPI_SOLIDS.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include <Solids_Evolution/SOLIDS_PARAMETERS.h>
#include <Standard_Tests/STANDARD_OPTIONS.h>
#include "SAILING.h"
#include "Water/SAILING_WATER.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    RW rw=RW();STREAM_TYPE stream_type(rw); // gcc 3.3.2 workaround

    MPI_WORLD mpi_world(argc,argv);

    PARSE_ARGS parse_args;STANDARD_OPTIONS options(parse_args,&mpi_world);
    parse_args.Add_Double_Argument("-poisson",.5,"","poisson's ratio for test 24");
    parse_args.Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args.Add_Double_Argument("-stiffen_bending",10,"","stiffness multiplier for bending springs in various cloth tests");
    parse_args.Add_Double_Argument("-dampen_bending",1,"","damping multiplier for bending springs in various cloth tests");
    parse_args.Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
    parse_args.Add_Option_Argument("-noself","disable self-collisions");
    parse_args.Add_Option_Argument("-noglobalrepulsions","disable global repulsions");
    parse_args.Add_Option_Argument("-nopertimesteprepulsions","disable per time step repulsions");
    parse_args.Add_Integer_Argument("-repulsion_pair_update_frequency",INT_MAX,"How many time steps before repulsion pairs are recomputed");
    parse_args.Add_Integer_Argument("-topological_hierarchy_build_frequency",INT_MAX,"How many collision steps before a hierarchy topology rebuild");
    parse_args.Add_Integer_Argument("-side_panels",INT_MAX,"Cloth side panels");
    parse_args.Add_Integer_Argument("-cloth_triangles",INT_MAX,"Cloth number of triangles");
    parse_args.Add_Option_Argument("-noalt","don't use altitude springs");
    parse_args.Add_Option_Argument("-backward","use backward Euler evolution");
    parse_args.Add_Option_Argument("-printpairs","output interaction pairs");
    parse_args.Add_Option_Argument("-spectrum","run spectral analysis during timestepping");
    parse_args.Add_Option_Argument("-residuals","print residuals during timestepping");
    parse_args.Add_Option_Argument("-binding_springs","use binding springs for drift particles");
    parse_args.Add_Option_Argument("-velocity_prune","turn on velocity culling for collision pairs");
    parse_args.Add_Integer_Argument("-totalloops",0,"","set total loops");
    parse_args.Add_String_Argument("-example","","example name");
    parse_args.Add_Option_Argument("-noattractions","disable attractions on repulsion pair inversions");
    parse_args.Add_Double_Argument("-attractionthreshold",-.3,"","threshold for attractions of inverted repulsion pairs");
    parse_args.Add_Double_Argument("-clothcfl",4.,"Cloth CFL");
    // incompressible-specific arguments
    parse_args.Add_Option_Argument("-incomp","use INCOMPRESSIBLE_TESTS");
    parse_args.Add_Double_Argument("-hittime",1,"time required for splat");
    parse_args.Add_Double_Argument("-framerate",24,"frame rate for incompressible tests");
    parse_args.Add_Option_Argument("-hires","use high resolution objects");
    parse_args.Add_Integer_Argument("-cgincomp",20,"maximum number of CG iterations for incompressible pressure solves");
    parse_args.Add_Option_Argument("-project","Combine boundary one-rings with neighbors");
    parse_args.Add_Double_Argument("-cgsolids",1e-3,"CG tolerance for backward Euler");
    parse_args.Add_Double_Argument("-min_volume_recovery_time_scale",0,"Minimum dt over which entire volume may be recovered");
    parse_args.Add_Option_Argument("-noneumann","Do not apply Neumann boundary conditions");
    parse_args.Add_Double_Argument("-ground_friction",0,"Ground friction");
    // SIMPLE_HARD_BINDING options
    parse_args.Add_Integer_Argument("-subsamples",0,"subsamples");
    parse_args.Add_Double_Argument("-sphere_scale",1,"sphere scale");
    parse_args.Add_Double_Argument("-radius",1,"radius");
    parse_args.Add_Option_Argument("-dynamic","dynamic");
    parse_args.Add_Double_Argument("-binding_stiffness",100,"");
    parse_args.Add_Double_Argument("-binding_overdamping_fraction",1,"");
    parse_args.Add_Option_Argument("-backward_euler","use backward Euler update instead of default modified Newmark");
    parse_args.Add_Option_Argument("-backward_euler_position","do the modified Newmark position update step using backward Euler");
    // Hair options
    parse_args.Add_String_Argument("-hairsim","","the hair sime to run");
    parse_args.Add_String_Argument("-modelname","","the rigid model to bind to");
    parse_args.Add_String_Argument("-guide","","the guide hair sim to read from");
    // mocap flesh
    parse_args.Add_Integer_Argument("-steps",1,"steps per frame");
    parse_args.Add_Option_Argument("-water");
    parse_args.Add_Integer_Argument("-resolution",1);

    parse_args.Parse(argc,argv);
    options.Initialize_Logging();
    parse_args.Print_Arguments(argc,argv);

    //bool incompressible=parse_args.Get_Option_Value("-incomp");
    std::string example_name=parse_args.Get_String_Value("-example");

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID_3D<T> >* example;

    //example=new STANDARD_TESTS<T>(stream_type,options.Subexample(),parse_args);
    if(parse_args.Get_Option_Value("-water")) example = new SAILING_WATER<T>(stream_type,options.Subexample(),parse_args);
    else example = new SAILING<T>(stream_type,options.Subexample(),parse_args);

    if(mpi_world.initialized) example->solids_parameters.mpi_solids=new MPI_SOLIDS<VECTOR<T,3> >();
    options.Apply_Options(*example,example->solids_parameters.mpi_solids);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID_3D<T> > driver(*example);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
