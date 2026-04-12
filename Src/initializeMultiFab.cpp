#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include <MyFunctions.H>

using namespace amrex;

void initializeMultiFab(MultiFab& phi, amrex::Real init_val)
{
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phiarr = phi.array(mfi);
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    phiarr(i,j,k) = init_val;
                });
    }
}