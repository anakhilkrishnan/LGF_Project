#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <string>



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
    amrex::Print() << "This is the start of my first NSE solver." << "\n";
}

void extendedMain()
{
    
}