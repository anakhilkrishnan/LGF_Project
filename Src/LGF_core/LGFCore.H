#ifndef LGFCORE_H_
#define LGFCORE_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Math.H>


// defining structs used within LGFcore
struct MeshData 
{
    amrex::MultiFab mf;
    amrex::Geometry geom;
};

struct FabMetaData
{
    int offset;
    amrex::IntVect lo;
    amrex::IntVect hi;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx;
};

struct ConsolidatedData 
{
    amrex::Vector<amrex::Real> data;
    amrex::Vector<FabMetaData> metadata;
};

// defining key functions within LGFcore
MeshData createMeshData(amrex::Array<amrex::Real, AMREX_SPACEDIM> phy_dom_lo, amrex::Array<amrex::Real, AMREX_SPACEDIM> phy_dom_hi, int n_cell, int max_grid_size, int n_ghost, int n_comp);
ConsolidatedData consolidateMeshData(const MeshData& phifab);
void addEverySourceBox(const MeshData& sourceMeshData, MeshData& targetMeshData);

// function for computing lgf kernel inside parallelFor
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real computeLGF(AMREX_D_DECL(amrex::Real x_tar, amrex::Real y_tar, amrex::Real z_tar),
                        AMREX_D_DECL(amrex::Real x_src, amrex::Real y_src, amrex::Real z_src))
    {
        // compute distance squared
        amrex::Real rsquare = AMREX_D_TERM((x_tar - x_src)*(x_tar - x_src),
                                           + (y_tar - y_src)*(y_tar - y_src),
                                           + (z_tar - z_src)*(z_tar - z_src));
        
        // avoiding self interaction
        if (rsquare < 1e-12)
        {
            return 0.0;
        }
            
        // compute and return the green's function value of the given pair
        // written as preprocessor directives to keep the function light on GPU
        #if AMREX_SPACEDIM == 1
            amrex::Real r = std::sqrt(rsquare);
            // 1D: G(x) = 0.5 * |x|
            return 0.5 * r;
            
        #elif AMREX_SPACEDIM == 2
            // 2D: G(r) = (1 / 2*pi) * ln(r)
            return (1.0 / (4.0 * amrex::Math::pi<amrex::Real>())) * std::log(rsquare);
            
        #elif AMREX_SPACEDIM == 3
            amrex::Real r = std::sqrt(rsquare);
            // 3D: G(r) = -1 / (4*pi*r)
            return -1.0 / (4.0 * amrex::Math::pi<amrex::Real>() * r);
            
        #endif        
    }
#endif