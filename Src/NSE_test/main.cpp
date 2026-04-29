#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <string>

#include <FlowField.H>
#include <ProjectionWorkspace.H>


using namespace amrex;

void test_print()
{
    amrex::Print() << "This is the start of my first NSE solver." << "\n";
}

void writeStaggeredPlotFile(int step, amrex::Real time, const FlowField& state, const amrex::BoxArray ba, const amrex::DistributionMapping dm, const amrex::Geometry& geom, int n_cell, std::string plot_prefix)
{
    
    // building a multiFab with 3 components for plotting
   amrex::MultiFab plotFab(ba, dm, 5, 0);

    // 1. Use AMReX's native hardware-accelerated face-to-cell averaging
    #if AMREX_SPACEDIM == 1
        amrex::average_face_to_cellcenter(plotFab, 0, amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>{&state.getVel(0)});
    #elif AMREX_SPACEDIM == 2
        amrex::average_face_to_cellcenter(plotFab, 0, amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>{&state.getVel(0), &state.getVel(1)});
    #elif AMREX_SPACEDIM == 3
        amrex::average_face_to_cellcenter(plotFab, 0, amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>{&state.getVel(0), &state.getVel(1), &state.getVel(2)});
    #endif
    

    // 2. Pressure and TagRegion are already cell-centered, just copy them!
    amrex::MultiFab::Copy(plotFab, state.getPres(), 0, 3, 1, 0); 
    amrex::MultiFab::Copy(plotFab, state.getTagRegion(), 0, 4, 1, 0);

    // exporting the names of the MultiFabs
    amrex::Vector<std::string> varnames = {"pressure", "x_velocity", "x_velocity", "x_velocity", "Active_Box_Tag"};

    // writing a simple plotfile
    const std::string& plotfile_name = amrex::Concatenate(plot_prefix, step);
    amrex::Print() << "Writing plotfile to: " << plotfile_name << "\n";
    WriteSingleLevelPlotfile(plotfile_name, plotFab, varnames, geom, time, step);
    amrex::Print() << "Plotfile written to: " << plotfile_name << "\n";

}

void extendedMain()
{
    BL_PROFILE("extendedMain()");

    auto overall_start_time = amrex::second();

    // necessary ParmParse parameters
    int n_cell, max_grid_size, rk_order, plot_int, max_steps, n_ghost;
    amrex::Real t_start, t_stop, cfl, source_tag_thresh;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dom_lo, dom_hi;

    // pre-defined ParmParse parameters
    int n_comp = 1;
    bool write_plot = true;
    std::string plot_prefix = "./Results/plt";

    // reading inputs file
    amrex::ParmParse pp;

    // domain data
    pp.get("n_cell",n_cell);
    pp.get("max_grid_size",max_grid_size);
    pp.get("domain_lo", dom_lo);
    pp.get("domain_hi", dom_hi);
    pp.get("n_ghost", n_ghost);
    pp.query("n_comp", n_comp);


    // solver data
    pp.get("tagging_threshold", source_tag_thresh);

    // time stepping data
    pp.get("t_start", t_start);
    pp.get("t_stop", t_stop);
    pp.get("cfl", cfl);        
    pp.get("rk_order", rk_order);
    pp.get("max_steps", max_steps);

    // plotting data
    pp.query("write_plot", write_plot);
    pp.get("plot_int", plot_int);
    pp.query("plot_prefix", plot_prefix);


    // creating domain data objects
    amrex::IntVect dom_lo_iv(AMREX_D_DECL(0, 0, 0));
    amrex::IntVect dom_hi_iv(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
    amrex::Box domain(dom_lo_iv, dom_hi_iv);

    amrex::BoxArray ba(domain);
    ba.maxSize(max_grid_size);

    amrex::DistributionMapping dm(ba);

    amrex::RealBox real_box(dom_lo, dom_hi);
    amrex::Vector<int> is_periodic(AMREX_SPACEDIM, 0); // infinite domain using zero-grad BC
    amrex::Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic.data());

    // initializing solver objects
    FlowField state_n(geom, ba, dm, n_comp, n_ghost);
    ProjectionWorkspace workspace(geom, ba, dm, n_comp, n_ghost);

    // PENDING: write a function that performs initialization of state_n as per problem
    
    // performing time stepping
    amrex::Real dt;
    amrex::Real time = t_start;
    int step = 0;

    // plotting initial conditions
    if (write_plot)
    {
        BL_PROFILE("<IO> Initial Plot()");
        writeStaggeredPlotFile(step, time, state_n, ba, dm, geom, n_cell, plot_prefix);

    }

    while(time < t_stop && step < max_steps)
    {
        auto step_start_time = amrex::second();
        
        // PENDING: write a function to dynamically compute dt and satisfy CFL criterion
        dt = 0.0001;

        // advance time using RK for time, KEP Morinishi for space and LGF for pressure poisson
        workspace.advanceTimeStep(state_n, dt, rk_order);

        // update counters
        time += dt;
        step++;

        //  plot in specified intervals
        if (step % plot_int == 0 && write_plot)
        {
            BL_PROFILE("<IO> Interval Plot()");
            writeStaggeredPlotFile(step, time, state_n, ba, dm, geom, n_cell, plot_prefix);
        }

        // track duration of simulation
        auto step_stop_time = amrex::second();
        auto step_duration = step_stop_time - step_start_time;

        // print to terminal each timestep
        amrex::Print() << "Step: " << step << " | Time: " << time << " | dt: " << dt 
                       << " | WallTime: " << (step_duration) << "s | divU_max: " << workspace.divU_max_norm << "\n";
    }

    // overall code walltime tracking
    auto overall_end_time = amrex::second();
    auto elapsed_time = overall_end_time - overall_start_time;

    // making copies to track slowest and fastest processor
    amrex::Real max_time = elapsed_time;
    amrex::Real min_time = elapsed_time;

    // performing a reduction over all the processors to track the slowest and
    // fastest MPI rank
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::ReduceRealMax(max_time, IOProc);
    amrex::ParallelDescriptor::ReduceRealMin(min_time, IOProc);

    amrex::Print() << "Max compute time (Slowest Rank): " << max_time << " s\n"
                   << "Min compute time (Fastest Rank): " << min_time << " s\n"
                   << "Time spread (Load Imbalance)   : " << (max_time - min_time) << " s\n";
}   

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    test_print();
    extendedMain();

    amrex::Finalize();
    return 0;
}