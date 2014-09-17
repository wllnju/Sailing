//#####################################################################
// This file is part of code for Sailing project which was built in PhysBAM 
// whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class Sailing
//#####################################################################
//   1. equal_length cloth sail
//	 2. small_top cloth sail
//   3. small_bottom cloth sail
//#####################################################################
#ifndef __SAILING__
#define __SAILING__
#include <Articulated_Rigid_Bodies/POINT_JOINT.h>
#include <Articulated_Rigid_Bodies/ANGLE_JOINT.h>
#include <Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Articulated_Rigid_Bodies/CONSTRAINT_FUNCTION.h>
#include <Articulated_Rigid_Bodies/JOINT_FUNCTION_3D.h>
#include <Articulated_Rigid_Bodies/POINT_JOINT.h>
#include <Articulated_Rigid_Bodies/PRISMATIC_TWIST_JOINT.h>
#include <Articulated_Rigid_Bodies/RIGID_JOINT.h>
#include <Articulated_Rigid_Bodies/NORMAL_JOINT.h>
#include <Rigid_Bodies/RIGID_BODY_BASIC_FORCES.h>
#include <Rigid_Bodies/RIGID_BODY_INTERSECTIONS.h>
#include <Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Constitutive_Models/NEO_HOOKEAN.h>
#include <Constitutive_Models/ROTATED_LINEAR.h>
#include <Constitutive_Models/SPLINE_MODEL.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT.h>
#include <Deformable_Objects/LINEAR_BINDING.h>
#include <Deformable_Objects/RIGID_BODY_BINDING.h>
#include <Forces_And_Torques/BINDING_SPRINGS.h>
#include <Forces_And_Torques/ETHER_DRAG.h>
#include <Forces_And_Torques/FINITE_VOLUME.h>
#include <Forces_And_Torques/GRAVITY.h>
#include <Forces_And_Torques/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <Forces_And_Torques/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <Forces_And_Torques/LINEAR_SPRINGS.h>
#include <Forces_And_Torques/SEGMENT_BENDING_SPRINGS.h>
#include <Forces_And_Torques/TRIANGLE_BENDING_ELEMENTS.h>
#include <Forces_And_Torques/TRIANGLE_BENDING_SPRINGS.h>
#include <Forces_And_Torques/WIND_DRAG_3D.h>
#include <Fracture/EMBEDDED_MATERIAL_SURFACE.h>
#include <Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <Geometry/FREE_PARTICLES.h>
#include <Geometry/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/SPHERE.h>
#include <Geometry/TORUS.h>
#include <Grids/UNIFORM_GRID_ITERATOR_CELL_3D.h>
#include <Grids/UNIFORM_GRID_ITERATOR_FACE_3D.h>
#include <Grids/UNIFORM_GRID_ITERATOR_NODE_3D.h>
#include <Interpolation/INTERPOLATION_CURVE.h>
#include <Meshing/RED_GREEN_TRIANGLES.h>
#include <Random_Numbers/RANDOM_NUMBERS.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Solids_Evolution/QUASISTATIC_EVOLUTION.h>
#include <Solids_Evolution/SOLIDS_PARAMETERS.h>
#include <Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Utilities/DEBUG_PRINT.h>
#include <Utilities/PARSE_ARGS.h>
#include <boost/math/special_functions/asinh.hpp>
namespace PhysBAM{

using boost::math::asinh;

template<class T_input>
class SAILING:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID_3D<T_input> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef UNIFORM_GRID_ITERATOR_CELL_3D<T> CELL_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_FACE_3D<T> FACE_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_NODE_3D<T> NODE_ITERATOR;
    typedef typename RIGID_BODY_POLICY<TV>::RIGID_BODY T_RIGID_BODY;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID_3D<T> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;

    int test_number;    
    SOLIDS_STANDARD_TESTS<TV> tests;
    ARRAY<TV> deformable_body_rest_positions;
    
    int number_side_panels;
    T aspect_ratio,side_length;
    int constrained_particle,suspended_particle,drifting_particle;
    bool no_altitude_springs;
    T stiffness_multiplier,damping_multiplier,bending_stiffness_multiplier,bending_damping_multiplier;
	bool use_forces_for_drift;
	int cloth_triangles;

    T test_30_constrained_off;
    T test_30_friction_off;
    T test_30_wind_off;

    //#####################################################################
    // Constructor
    //#####################################################################    
    SAILING(const STREAM_TYPE stream_type,const int test_number,const PARSE_ARGS& parse_args)
        :BASE(stream_type,0,fluids_parameters.NONE),test_number(test_number),tests(*this),number_side_panels(40),aspect_ratio((T)1.0),side_length((T)4.0),
        constrained_particle(0),suspended_particle(0),drifting_particle(0)
    {
        LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
        output_directory=str(boost::format("./Test_%d")%test_number);
        frame_rate=24;
        T cloth_clamp_fraction=(T).03; // from curtain and ball

        int total_loops=parse_args.Get_Integer_Value("-totalloops");
        if(total_loops!=0) solids_parameters.total_collision_loops=total_loops;
        no_altitude_springs=parse_args.Get_Option_Value("-noalt");
        stiffness_multiplier=(T)parse_args.Get_Double_Value("-stiffen");
        damping_multiplier=(T)parse_args.Get_Double_Value("-dampen");
        bending_stiffness_multiplier=(T)parse_args.Get_Double_Value("-stiffen_bending");
        bending_damping_multiplier=(T)parse_args.Get_Double_Value("-dampen_bending");
        use_forces_for_drift=parse_args.Get_Option_Value("-binding_springs");
        T cloth_cfl=(T)parse_args.Get_Double_Value("-clothcfl");
        
        solids_parameters.use_rigid_deformable_contact=true;
        solids_parameters.deformable_object.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
        solids_parameters.verbose_super_fragments=true;
        solids_parameters.use_push_out=true;
        
                frame_rate=120;last_frame=(int)(20*frame_rate);                
                solids_parameters.throw_exception_on_backward_euler_failure=false;
                solids_parameters.cfl=cloth_cfl; // was 4
                solids_parameters.self_collision_friction_coefficient=(T)3.2;
                solids_parameters.cg_tolerance=(T)1e-2;
                if(parse_args.Is_Value_Set("-cgsolids")) solids_parameters.cg_tolerance=(T)parse_args.Get_Double_Value("-cgsolids");
                solids_parameters.cg_iterations=200;
                solids_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
                solids_parameters.perform_self_collision=false;//true;               
                    test_30_constrained_off=(T).3;
                    test_30_friction_off=(T)3;
                    test_30_wind_off=(T)1000;                
                PHYSBAM_DEBUG_PRINT("Basic settings",solids_parameters.cfl,solids_parameters.cg_tolerance,solids_parameters.cg_iterations);                

        if(parse_args.Get_Option_Value("-backward")){
            if(typeid(*solids_parameters.solids_evolution)!=typeid(NEWMARK_EVOLUTION<TV>)){
                LOG::cerr<<"refusing to replace non-Newmark evolution with -backward"<<std::endl;exit(1);}
            solids_parameters.cg_tolerance=(T)1e-2;
            solids_parameters.cfl*=10;
            solids_parameters.Set_Evolution(new BACKWARD_EULER_EVOLUTION<TV>(solids_parameters));
            output_directory+="_backward";}

        if(parse_args.Get_Option_Value("-noself")){
            if(!solids_parameters.perform_self_collision){LOG::cerr<<"-noself invalid for examples without self-collisions"<<std::endl;exit(1);}
            solids_parameters.perform_self_collision=false;
            output_directory+="_noself";}

        if(parse_args.Get_Option_Value("-noglobalrepulsions")){
            solids_parameters.perform_per_collision_step_repulsions=false;
            output_directory+="_noglobalrepulsions";}

        if(parse_args.Get_Option_Value("-nopertimesteprepulsions")){
            solids_parameters.perform_per_time_step_repulsions=false;
            output_directory+="_nopertimesteprepulsions";}

        if(parse_args.Get_Option_Value("-noattractions")){
            solids_parameters.perform_repulsion_pair_attractions=false;
            output_directory+="_noattractions";}

        if(parse_args.Is_Value_Set("-repulsion_pair_update_frequency")){
            solids_parameters.repulsion_pair_update_frequency=parse_args.Get_Integer_Value("-repulsion_pair_update_frequency");
            output_directory+=str(boost::format("_repulsionpairupdatefrequency=%d")%solids_parameters.repulsion_pair_update_frequency);}

        if(parse_args.Get_Option_Value("-velocity_prune")){
            //solids_parameters.collision_pair_velocity_pruning=true;
            output_directory+="_velocityprune";}

        if(parse_args.Is_Value_Set("-topological_hierarchy_build_frequency")){
            solids_parameters.topological_hierarchy_build_frequency=parse_args.Get_Integer_Value("-topological_hierarchy_build_frequency");
            output_directory+=str(boost::format("_topologicalhierarchybuildfrequency=%d")%solids_parameters.topological_hierarchy_build_frequency);}

        if(parse_args.Is_Value_Set("-side_panels")){
            number_side_panels=parse_args.Get_Integer_Value("-side_panels");
            output_directory+=str(boost::format("_sidepanels=%d")%number_side_panels);}
        cloth_triangles=2*number_side_panels*(int)(number_side_panels*aspect_ratio);

        if(parse_args.Is_Value_Set("-cloth_triangles")){
            if(parse_args.Is_Value_Set("-side_panels")) throw std::runtime_error("Cannot have both side panels and cloth triangles set");
            cloth_triangles=parse_args.Get_Integer_Value("-cloth_triangles");
            if(test_number==35) cloth_triangles/=100;
            number_side_panels=(int)ceil(sqrt((T).5*(T)cloth_triangles/aspect_ratio));
            int actual_cloth_triangles=2*((int)(aspect_ratio*number_side_panels))*number_side_panels;
            PHYSBAM_DEBUG_PRINT("cloth resolution",cloth_triangles,actual_cloth_triangles,number_side_panels);
            cloth_triangles=actual_cloth_triangles;
            if(test_number==35) cloth_triangles*=100;}

        solids_parameters.output_interaction_pairs=parse_args.Get_Option_Value("-printpairs");
        solids_parameters.spectral_analysis=parse_args.Get_Option_Value("-spectrum");
        solids_parameters.deformable_object.Print_Residuals(parse_args.Get_Option_Value("-residuals"));
        solids_parameters.repulsion_pair_attractions_threshold=(T)parse_args.Get_Double_Value("-attractionthreshold");
    }

    ~SAILING()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T dt,const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(RAW_ARRAY<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(RAW_ARRAY<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(RAW_ARRAY<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(RAW_ARRAY<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(RAW_ARRAY<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(RAW_ARRAY<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(RAW_ARRAY<FRAME<TV> > X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(RAW_ARRAY<TWIST<TV> > twist,const T velocity_time,const T current_position_time,const SUPER_FRAGMENT<TV>& super_fragment) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Fragments() PHYSBAM_OVERRIDE {}
    void Postprocess_Fragments() PHYSBAM_OVERRIDE {}
    void Postprocess_Super_Fragments(const LIST_ARRAY<PAIR<SUPER_FRAGMENT_ID,SUPER_FRAGMENT_ID> >& swap_pairs,const LIST_ARRAY<SUPER_FRAGMENT_ID>& rebuild,
        SUPER_FRAGMENT_ID old_max) PHYSBAM_OVERRIDE {}
    void Update_Super_Fragments(const LIST_ARRAY<PAIR<SUPER_FRAGMENT_ID,SUPER_FRAGMENT_ID> >& swap_pairs,const LIST_ARRAY<SUPER_FRAGMENT_ID>& rebuild,
        SUPER_FRAGMENT_ID old_max) PHYSBAM_OVERRIDE {}
    void Add_External_Super_Fragment_Connectivity(SPARSE_UNION_FIND<FRAGMENT_ID>& union_find) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    
    
void Initialize_Sail(bool small_top, bool small_bottom)
{
	//equal_length case
	if(!small_top && !small_bottom ){	
        tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,side_length,0),ROTATION<TV>((T)pi/2,TV((T)0,0,(T)1).Normalized()))));
    }
  	//small_top case
  	else if(small_top){
	   	TRIANGULATED_SURFACE<T>& cloth_panel=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,side_length,0),ROTATION<TV>(-(T)pi/2,TV((T)0,0,(T)1).Normalized()))));
	   	//remove unwanted triangular mesh
		if(small_top || small_bottom){  		
		    int N = number_side_panels;
		    int m;
		    for(int i = N-1; i >= 0; i--){
		    	m = N-i-1;//(N-i-1)/2 *2;
		    	for(int j = 0; j < m-1; j++){
		    		cloth_panel.mesh.elements.Remove_Index_Lazy((i+1)*N*2-j);
		    	}
		    	for(int j = m; j > 0; j--){
		    		cloth_panel.mesh.elements.Remove_Index_Lazy(i*N*2+j);
		    	}
		    }
		    cloth_panel.Discard_Valence_Zero_Particles_And_Renumber();
		}
	}
	//small_bottom case
	else if(small_bottom){	   		
	   	TRIANGULATED_SURFACE<T>& cloth_panel=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,side_length,0),ROTATION<TV>((T)pi/2,TV((T)0,0,(T)1).Normalized()))));
	   	//remove unwanted triangular mesh
		if(small_top || small_bottom){  		
		    int N = number_side_panels;
		    int m;
		    for(int i = N-1; i >= 0; i--){
		    	m = N-i-1;//(N-i-1)/2 *2;
		    	for(int j = 0; j < m-1; j++){
		    		cloth_panel.mesh.elements.Remove_Index_Lazy((i+1)*N*2-j);
		    	}
		    	for(int j = m; j > 0; j--){
		    		cloth_panel.mesh.elements.Remove_Index_Lazy(i*N*2+j);
		    	}
		    }
		    cloth_panel.Discard_Valence_Zero_Particles_And_Renumber();
		}
	}	
}    
    
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    bool automatically_add_to_collision_structures=false;//true;
    // deformable bodies
    DEFORMABLE_OBJECT<TV>& deformable_object=solids_parameters.deformable_object;
    //RIGID_BODY_PARTICLE<TV>& rigid_body_particles=deformable_object.rigid_body_particles;
    SOLIDS_PARTICLE<TV>& particles=deformable_object.particles;
    BINDING_LIST<TV>& binding_list=deformable_object.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=deformable_object.soft_bindings;
    
    switch(test_number){
        case 1: Initialize_Sail(false,false);break;//equal_length case
        case 2: Initialize_Sail(true,false);break;//small_top case
        case 3: Initialize_Sail(false,true);break;//small_bottom case
        
      default: PHYSBAM_FATAL_ERROR(str(boost::format("Unrecognized test number %d")%test_number));}
    
	  ARTICULATED_RIGID_BODY_3D<T>& arb=solids_parameters.deformable_object.articulated_rigid_body;
	  JOINT<TV>* joint;
	  PRISMATIC_TWIST_JOINT<TV>* joint_p;
            
    T sidey=side_length*aspect_ratio/2.0;
        
   	RIGID_BODY_3D<T>& bar_top=tests.Add_Analytic_Cylinder(side_length,0.1,10,4);bar_top.Frame().t=TV(0,side_length+sidey,0); 
		RIGID_BODY_3D<T>& bar_bottom=tests.Add_Analytic_Cylinder(side_length,0.1,10,4);bar_bottom.Frame().t=TV(0,side_length-sidey,0); bar_bottom.Set_Mass(0.1);
		RIGID_BODY_3D<T>& mast=tests.Add_Analytic_Box(TV(0.2,2*side_length,0.2));mast.Frame().t=TV((T)-0.1,side_length+0.1,0);  
        					
    tests.Bind_Particles_In_Rigid_Body(bar_top);
		tests.Bind_Particles_In_Rigid_Body(bar_bottom);		    
		    
    mast.Set_Name("parent");
    bar_top.Set_Name("child");
    joint=new RIGID_JOINT<TV>();((RIGID_JOINT<TV>*)joint)->Set_Prismatic_Component_Translation(TV((T)0,0,0));    
    arb.joint_mesh.Add_Articulation(mast.Id_Number(),bar_top.Id_Number(),joint);
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0.1,sidey-0.1,0)));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,0)));
          
    mast.Set_Name("parent");
		bar_bottom.Set_Name("child");
		joint_p=new PRISMATIC_TWIST_JOINT<TV>();
		arb.joint_mesh.Add_Articulation(mast.Id_Number(),bar_bottom.Id_Number(),joint_p);
		//joint_p->Set_Prismatic_Constraints(VECTOR<bool,3>(true,true,true),TV(0,0.01,0),TV(0,-0.1,0));
		joint_p->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0.1,-sidey-0.1,0)));
		joint_p->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,0)));
			
		RIGID_BODY_3D<T>& box=tests.Add_Analytic_Box(TV(4,0.2,4));
		box.Set_Mass(10);
		box.Frame().t=TV((T)-0.1,(T)0,0);
		box.Set_Name("parent");
		mast.Set_Name("child");
		joint=new RIGID_JOINT<TV>();((RIGID_JOINT<TV>*)joint)->Set_Prismatic_Component_Translation(TV((T)0,0,0));    
		arb.joint_mesh.Add_Articulation(box.Id_Number(),mast.Id_Number(),joint);
		joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,0.1,0)));
		joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-side_length,0)));
		
		//RIGID_BODY_3D<T>& ground=tests.Add_Analytic_Box(TV(40,2,20));ground.Frame().t=TV((T)16,(T)-1.1,-6);ground.is_static=true;
		RIGID_BODY_3D<T>& ground=tests.Add_Ground((T)0.2,(T)-.1);
		ground.Set_Coefficient_Of_Friction(0.1);
		box.Set_Coefficient_Of_Friction(0.1);
		ground.Set_Name("parent");
		box.Set_Name("child");
		T ground_angle_rad = 0;
		FRAME<TV> frame(TV(-sin(ground_angle_rad),cos(ground_angle_rad),0),ROTATION<TV>(ground_angle_rad,TV(0,0,1)));
		joint=new NORMAL_JOINT<TV>();
		arb.joint_mesh.Add_Articulation(ground.Id_Number(),box.Id_Number(),joint);
		FRAME<TV> J(frame*TV((T)0,-(T)1,(T)0),ROTATION<TV>::From_Rotated_Vector(TV::Axis_Vector(1),frame.r.Rotated_Axis(2)));
		joint->Set_Joint_To_Parent_Frame(ground.Frame().Inverse()*J);
		joint->Set_Joint_To_Child_Frame(box.Frame().Inverse()*J);
  
		
    // add structures and rigid bodies to collisions
    if(automatically_add_to_collision_structures) deformable_object.collisions.collision_structures.Append_Elements(deformable_object.structures.Raw_Pointers());
  solids_parameters.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_object.structures.Raw_Pointers());
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.deformable_object.rigid_body_particles.rigid_body_list);
    //solids_parameters.deformable_object.collisions.Use_Structure_Skip_Collision_Body();
    //solids_parameters.deformable_object.collisions.Use_Structure_Skip_Collision_Body();

    // correct number nodes
    for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.mass.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solids_parameters.deformable_object;
    //SOFT_BINDINGS<TV>& soft_bindings=deformable_object.soft_bindings;
    //SOLIDS_PARTICLE<TV>& particles=deformable_object.particles;

    Get_Initial_Data();

            deformable_object.Add_Force(new GRAVITY<TV>(deformable_object.particles,deformable_object.rigid_body_particles,true,true));
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            for(int i=1;i<=deformable_object.structures.m;i++) if(TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_object.structures(i).get())){
               // if(surface->mesh.elements.m != cloth_triangles) PHYSBAM_FATAL_ERROR();
                deformable_object.Add_Force(Create_Edge_Springs(*surface,linear_stiffness,linear_damping)); // were *2 and *10
                deformable_object.Add_Force(Create_Bending_Springs(*surface,bending_stiffness,bending_damping));
                WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(*surface);deformable_object.Add_Force(drag);
                drag->Use_Linear_Normal_Viscosity((T)5);drag->Use_Constant_Wind(0,TV((T)1,(T)0,(T)0));}

}
//#####################################################################
// Function Add_External_Fragment_Connectivity
//#####################################################################
void Add_External_Fragment_Connectivity(PARTICLE_CONNECTIVITY<TV>& particle_connectivity,LIST_ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    BASE::Read_Output_Files_Solids(frame);
    solids_parameters.deformable_object.Update_Fragments();
    /*int count=1;
    for(int i=1;i<=solids_parameters.deformable_object.binding_list.bindings.m;i++)
        if(RIGID_BODY_BINDING<T,TV>* binding=dynamic_cast<RIGID_BODY_BINDING<T,TV>*>(solids_parameters.deformable_object.binding_list.bindings(i))){
            binding->rigid_body_particles=&solids_parameters.deformable_object.rigid_body_particles;
            binding->rigid_body_particle_index=stored_bindings(count).x;
            binding->object_space_position=stored_bindings(count).y;
            count++;}*/
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    LOG::cout<<"Preprocess Frame "<<frame<<std::endl;    
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{/*
    if(LINEAR_SPRINGS<TV>* linear_springs=solids_parameters.deformable_object.template Find_Force<LINEAR_SPRINGS<TV>*>())
        linear_springs->Print_Deformation_Statistics();*/
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(RAW_ARRAY<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
	//Pointwise_Object_Velocity();

}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(RAW_ARRAY<TV> V,const T velocity_time,const T current_position_time,const SUPER_FRAGMENT<TV>& super_fragment) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(RAW_ARRAY<TV> X,const T time) PHYSBAM_OVERRIDE
{
    return;
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(RAW_ARRAY<TV> X,const T time,const SUPER_FRAGMENT<TV>& super_fragment) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const RIGID_BODY_ID id) PHYSBAM_OVERRIDE
{
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const RIGID_BODY_ID id) PHYSBAM_OVERRIDE
{
}
//#####################################################################
};
}
#endif
