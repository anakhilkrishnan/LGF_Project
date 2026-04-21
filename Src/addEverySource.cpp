#include <AMReX_MultiFab.H>
#include <AMReX_BLProfiler.H>

#include <MyFunctions.H>
#include <ComputeLGF.H>

using namespace amrex;

void addEverySourceBox(const MeshData& sourceMeshData, MeshData& targetMeshData) 
{
    // adding profiling blocks for Tiny/Base profilers
    BL_PROFILE("<Compute> addEverySourceBox()");

    //extract MultiFabs and Geometries
    const amrex::MultiFab& source = sourceMeshData.mf;
    amrex::MultiFab& target = targetMeshData.mf;
    const amrex::Geometry& target_geom = targetMeshData.geom;

    GpuArray<amrex::Real, AMREX_SPACEDIM> dx = target_geom.CellSizeArray();
    GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo = target_geom.ProbLoArray();

    // Read data from the source MultiFab and make it available to all processes
    ConsolidatedData consolSource = consolidateMultiFab(source);
    
    // allocating space in VRAM for the source data and metadata
    amrex::Gpu::DeviceVector<Real> d_data(consolSource.data.size());
    amrex::Gpu::DeviceVector<FabMetaData> d_meta(consolSource.metadata.size());

    // copying the consolidated data from CPU to GPU
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, 
                     consolSource.data.begin(), consolSource.data.end(), 
                     d_data.begin());
                     
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, 
                     consolSource.metadata.begin(), consolSource.metadata.end(), 
                     d_meta.begin());

    // the above lines collapse into CPU use when compiled without USE_CUDA=TRUE

    // export the consolidated data as pointers to the target MFIter
    int num_blocks = consolSource.metadata.size();
    const Real* data_ptr = d_data.data();
    const FabMetaData* meta_ptr = d_meta.data();
    
    // Loop over target boxes in a separate MFIter
    for (amrex::MFIter mfi(target, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& targetbox = mfi.validbox();
        const Array4<Real>& phi = target.array(mfi);

        amrex::ParallelFor(targetbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {   
            // extract physical coordinates of target cell
            amrex::Real x_tar = prob_lo[0] + ((i + 0.5) * dx[0]);

            #if AMREX_SPACEDIM >= 2
                amrex::Real y_tar = prob_lo[1] + ((j + 0.5) * dx[1]);
            #else
                amrex::Real y_tar = 0.0;
            #endif

            #if AMREX_SPACEDIM == 3
                amrex::Real z_tar = prob_lo[2] + ((k + 0.5) * dx[2]);
            #else
                amrex::Real z_tar = 0.0;
            #endif

            amrex::Real total_contribution = 0.0;

            // Iterate over every source block gathered from all MPI processes
            for (int b = 0; b < num_blocks; ++b) 
            {
                const auto& block = meta_ptr[b];
                int idx = block.offset;
                
                // Unpack and sum every cell in the source block
                for (int sk = AMREX_D_PICK(0, 0, block.lo[2]); sk <= AMREX_D_PICK(0, 0, block.hi[2]); ++sk) 
                {
                    // conditionally extract physical coordinates of source cell
                    #if AMREX_SPACEDIM == 3
                        amrex::Real z_src = prob_lo[2] + ((sk + 0.5) * dx[2]);
                    #else
                        amrex::Real z_src = 0.0;
                    #endif

                    for (int sj = AMREX_D_PICK(0, block.lo[1], block.lo[1]); sj <= AMREX_D_PICK(0, block.hi[1], block.hi[1]); ++sj) 
                    {
                        #if AMREX_SPACEDIM >= 2
                            amrex::Real y_src = prob_lo[1] + ((sj + 0.5) * dx[1]);
                        #else
                            amrex::Real y_src = 0.0;
                        #endif

                        for (int si = block.lo[0]; si <= block.hi[0]; ++si) 
                        {
                            amrex::Real x_src = prob_lo[0] + ((si + 0.5) * dx[0]);
                            
                            // compute the LGF kernel for the current source-target cell pair
                            amrex::Real lgf = computeLGF(x_tar, y_tar, z_tar, x_src, y_src, z_src);
                            amrex::Real dvol = AMREX_D_TERM(dx[0], * dx[1], * dx[2]);
                            
                            // add the contribution of source cell based on lgf
                            total_contribution += (data_ptr[idx++] * lgf * dvol);
                        }
                    }
                }
            }
            phi(i, j, k) = total_contribution;
        });
    }

    // wait for the GPU to synchnorize
    amrex::Gpu::streamSynchronize();
}