#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void test_print();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    test_print();

    amrex::Finalize();
    return 0;
}

void test_print()
{
    amrex::Print() << "This is finally the start of my AMReX code." << "\n";
}
