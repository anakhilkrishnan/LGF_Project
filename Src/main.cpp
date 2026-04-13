#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

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
    // this is funny, idk why it paniqued, this is NOT needed anymore since
    // the use of addEverySourceBox()
    MFIter::allowMultipleMFIters(true);

    // marking start to track runtimes
    auto start_time = amrex::second();

    // variables to be read from ParmParse
    int n_cell, max_grid_size;

    // reading inputs file
    ParmParse pp;
    pp.get("n_cell",n_cell);
    pp.get("max_grid_size",max_grid_size);

    // initializing parameters for MultiFabs
    int n_ghost = 1;
    int n_comp = 1;

    // creating source and target multifabs
    MeshData sourceFab = createMeshData(n_cell, max_grid_size, n_ghost, n_comp);
    MeshData targetFab = createMeshData(n_cell, max_grid_size, n_ghost, n_comp);

    // initializing multifabs
    initializeMultiFab(sourceFab.mf, 1.0);
    initializeMultiFab(targetFab.mf, 0.0);

    // perform addition of box values
    addEverySourceBox(sourceFab.mf, targetFab.mf);

    // writing a simple plotfile
    std::string dir_name = "../Results";
    std::string plotfile_name = dir_name + "/plt00000";
    WriteSingleLevelPlotfile(plotfile_name, targetFab.mf, {"phi"}, targetFab.geom, 0.0, 0);
    amrex::Print() << "Plotfile written to: " << plotfile_name << "\n";

    // marking end time and elapsed time
    auto end_time = amrex::second();
    auto elapsed_time = end_time - start_time;
    amrex::Print() << "Time taken for computation: " << elapsed_time << "\n";
}
