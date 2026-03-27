#include <AMReX_MultiFab.H>

#include <MyFunctions.H>

using namespace amrex;

void addEverySourceBox(const MultiFab& source, const MultiFab& target, MultiFab& result) 
{
    // Initialize the result MultiFab for += operations
    result.setVal(0.0);

    // Export the source data to a single big SourceBlock
    amrex::Vector<SourceBlock> full_source_data;

    for (amrex::MFIter mfi(source); mfi.isValid(); ++mfi) 
    {
        const auto& bx = mfi.validbox();
        // Logic for checking whether this source is to be taken or not goes here
        full_source_data.push_back({source.array(mfi), bx.smallEnd(), bx.bigEnd()});
    }

    // Transfer source data to the GPU. In case USE_OMP = TRUE, it falls back to
    // the typical CPU setup
    amrex::Gpu::DeviceVector<SourceBlock> gpu_full_source_data;
    gpu_full_source_data.assign(full_source_data.begin(), full_source_data.end());

    // Obtain pointers for kernels to be used inside ParallelFor
    const SourceBlock* src_lib_ptr = gpu_full_source_data.data();
    int num_blocks = gpu_full_source_data.size();

    // Loop over target boxes in a separate MFIter
    for (amrex::MFIter mfi(target); mfi.isValid(); ++mfi)
    {
        const amrex::Box& targetbox = mfi.validbox();
        const Array4<Real>& phi = result.array(mfi);

        amrex::ParallelFor(targetbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {   
            amrex::Real total_contribution = 0.0;

            // Iterate over every block in gpu_full_source_data
            for (int b = 0; b < num_blocks; ++b) 
            {
                auto const& block = src_lib_ptr[b];
                auto const& src_val = block.data;

                // Iterate over every cell in the block
                for (int sj = block.lo[1]; sj <= block.hi[1]; ++sj) 
                {
                    for (int si = block.lo[0]; si <= block.hi[0]; ++si) 
                    {
                        // Add all the values in each source box
                        total_contribution += src_val(si, sj, 0);
                    }
                }
            }
            phi(i, j, k) = total_contribution;
        });
    }
}