#ifndef LGFFUNCTIONS_H_
#define LGFFUNCTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

using namespace amrex;

// function for computing lgf kernel inside parallelFor
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real computeLGF(amrex::Real x_tar, amrex::Real y_tar, amrex::Real z_tar,
                        amrex::Real x_src, amrex::Real y_tar, amrex::Real z_tar)
    {
        // compute distance squared
        amrex::Real rsquare = AMREX_D_TERM((x_tar - x_src)*(x_tar - x_src),
                                            (y_tar - y_src)*(y_tar - y_src),
                                            (z_tar - z_src)*(z_tar - z_src));
        
        // avoiding self interaction
        if (rsquare < 1e-12)
        {
            return 0.0;
        }
            
        amrex::Real r = amrex::sqrt(rsquare);
        
        // compute and return the green's function value of the given pair
        // written as preprocessor directives to keep the function light on GPU
        #if AMREX_SPACEDIM == 1
            // 1D: G(x) = 0.5 * |x|
            return 0.5 * r;
            
        #elif AMREX_SPACEDIM == 2
            // 2D: G(r) = (1 / 2*pi) * ln(r)
            return (1.0 / (2.0 * M_PI)) * std::log(r);
            
        #elif AMREX_SPACEDIM == 3
            // 3D: G(r) = -1 / (4*pi*r)
            return -1.0 / (4.0 * M_PI * r);
            
        #endif        
    }
#endif