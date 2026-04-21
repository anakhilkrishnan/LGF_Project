#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <string>

#include <MyFunctions.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    test_print();
    extendedMain();

    amrex::Finalize();
    return 0;
}

void test_print()
{
    amrex::Print() << "This is finally the start of my AMReX code." << "\n";
}

void extendedMain()
{
    // marking start to track runtimes
    auto start_time = amrex::second();

    // variables to be read from ParmParse
    int n_cell, max_grid_size;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> dom_lo, dom_hi;

    // setting a default plotfile prefix in case not specified in inputs
    std::string plot_prefix = "../Results/plt";

    // reading inputs file
    amrex::ParmParse pp;
    pp.get("n_cell",n_cell);
    pp.get("max_grid_size",max_grid_size);
    pp.get("domain_lo", dom_lo);
    pp.get("domain_hi", dom_hi);

    pp.query("plot_prefix", plot_prefix);

    // initializing parameters for MultiFabs
    int n_ghost = 1;
    int n_comp = 1;

    // creating source and target multifabs
    MeshData sourceFab = createMeshData(dom_lo, dom_hi, n_cell, max_grid_size, n_ghost, n_comp);
    MeshData targetFab = createMeshData(dom_lo, dom_hi, n_cell, max_grid_size, n_ghost, n_comp);

    // initializing multifabs
    // initializeMultiFab(sourceFab.mf, 1.0);
    initializeSourceMultiFab(sourceFab);
    initializeMultiFab(targetFab.mf, 0.0);

    auto compute_start_time = amrex::second();
    // perform addition of box values
    addEverySourceBox(sourceFab, targetFab);

    // marking end time and elapsed time
    auto compute_end_time = amrex::second();
    auto compute_time = compute_end_time - compute_start_time;
    amrex::Print() << "Time taken for computation: " << compute_time << "\n";
    
    {
        // adding profiling blocks for Tiny/Base profilers
        BL_PROFILE("<I/O> writingPlotfile");

        // writing a simple plotfile
        const std::string& plotfile_name = amrex::Concatenate(plot_prefix, n_cell);
        amrex::Print() << "Writing plotfile to: " << plotfile_name << "\n";
        WriteSingleLevelPlotfile(plotfile_name, targetFab.mf, {"phi"}, targetFab.geom, 0.0, 0);
        amrex::Print() << "Plotfile written to: " << plotfile_name << "\n";
    }
    

    auto end_time = amrex::second();
    auto elapsed_time = end_time - start_time;
    amrex::Print() << "Total program runtime: " << elapsed_time << "\n";
}
