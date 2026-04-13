#include <AMReX_MultiFab.H>

#include <MyFunctions.H>
#include <ComputeLGF.H>

using namespace amrex;

void addEverySourceBox(const MultiFab& source, MultiFab& target) 
{
    
    // Read data from the source MultiFab and make it available to all processes
    ConsolidatedData consolSource = consolidateMultiFab(source);
    
    // Further work: Exporting data for computations on GPU
    // ConsolidatedDataToGPUPointers();

    // Export the consolidated data as pointers to the target MFIter
    int num_blocks = consolSource.metadata.size();
    const Real* data_ptr = consolSource.data.data();
    const FabMetaData* meta_ptr = consolSource.metadata.data();
    
    // Loop over target boxes in a separate MFIter
    for (amrex::MFIter mfi(target); mfi.isValid(); ++mfi)
    {
        const amrex::Box& targetbox = mfi.validbox();
        const Array4<Real>& phi = target.array(mfi);

        amrex::ParallelFor(targetbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {   
            amrex::Real total_contribution = 0.0;

            // Iterate over every source block gathered from all MPI processes
            for (int b = 0; b < num_blocks; ++b) 
            {
                const auto& block = meta_ptr[b];
                int idx = block.offset;
                
                // Unpack and sum every cell in the source block
                for (int sk = AMREX_D_PICK(0, 0, block.lo[2]); sk <= AMREX_D_PICK(0, 0, block.hi[2]); ++sk) 
                {
                    for (int sj = AMREX_D_PICK(0, block.lo[1], block.lo[1]); sj <= AMREX_D_PICK(0, block.hi[1], block.hi[1]); ++sj) 
                    {
                        for (int si = block.lo[0]; si <= block.hi[0]; ++si) 
                        {
                            total_contribution += data_ptr[idx++];
                        }
                    }
                }
            }
            phi(i, j, k) = total_contribution;
        });
    }
}