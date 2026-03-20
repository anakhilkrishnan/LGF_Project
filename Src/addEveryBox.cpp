#include <AMReX_MultiFab.H>

#include <MyFunctions.H>

using namespace amrex;

void addEveryBox(const MultiFab& source, const MultiFab& target, MultiFab& result) 
{
    // Create a result MultiFab with the same structure as target
    result.setVal(0.0);

    // Loop over the boxes that THIS MPI rank owns in the target
    for (MFIter mfi_t(target); mfi_t.isValid(); ++mfi_t)
    {
        const Box& bx_t = mfi_t.validbox();
        auto const& res_arr = result.array(mfi_t);
        
        // NESTED LOOP: Now loop over ALL boxes in the source
        // Note: This is computationally expensive and assumes 
        // all source boxes are accessible to this MPI rank!
        for (MFIter mfi_s(source); mfi_s.isValid(); ++mfi_s)
        {
            auto const& src_arr = source.array(mfi_s);
            const Box& bx_s = mfi_s.validbox();

            // ParallelFor over the TARGET box pixels
            amrex::ParallelFor(bx_t, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Here is the tricky part: how do you map indices 
                // from the source box to the target box? 
                // For this "silly" exercise, let's just add the value 
                // from the same (i,j,k) if it exists in the source box.
                
                if (bx_s.contains(i,j,k)) {
                    res_arr(i,j,k) += src_arr(i,j,k);
                }
            });
        }
    }
}