#ifndef SPATIALDISCRETIZATION_H_
#define SPATIALDISCRETIZATION_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_Math.H>

using namespace amrex;

// PENDING: Morinishi spatial discretization is expected here
template <int dim>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real morinishiUFlux(int i, int j, int k,
                            amrex::GpuArray<amrex::Array4<amrex::Real const>, AMREX_SPACEDIM> vel_arr,
                            amrex::Array4<amrex::Real const> const& p,
                            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                            amrex::Real Re)
                            {

                                // extract velocities from vel_arr
                                auto const& u = vel_arr[0];
                                #if AMREX_SPACEDIM >= 2
                                    auto const& v = vel_arr[1];
                                #endif
                                #if AMREX_SPACEDIM == 3
                                    auto const& w = vel_arr[2];
                                #endif

                                // precompute important divisions
                                amrex::Real Re_inv = 1/Re;

                                amrex::Real idx  = 1.0 / dx[0];
                                amrex::Real idx2 = idx * idx;
                                
                                #if AMREX_SPACEDIM >= 2
                                    amrex::Real idy  = 1.0 / dx[1];
                                    amrex::Real idy2 = idy * idy;
                                #endif

                                #if AMREX_SPACEDIM == 3
                                    amrex::Real idz  = 1.0 / dx[2];
                                    amrex::Real idz2 = idz * idz;
                                #endif

                                amrex::Real presterm = 0.0;
                                amrex::Real viscterm = 0.0;
                                amrex::Real convterm = 0.0;

                                if constexpr (dim == 0)
                                {
                                    // computing the pressure term
                                    presterm = idx * (p(i,j,k) - p(i-1,j,k));

                                    // computing the viscous term dimension by dimension
                                    amrex::Real d2udx2 = (u(i+1,j,k) - 2.0*u(i,j,k) + u(i-1,j,k)) * idx2;
                                    viscterm = d2udx2;

                                    #if AMREX_SPACEDIM >= 2
                                        d2udy2 = (u(i,j+1,k) - 2.0*u(i,j,k) + u(i,j-1,k)) * idy2;
                                        viscterm += d2udy2;
                                    #endif

                                    #if AMREX_SPACEDIM == 3
                                        d2udz2 = (u(i,j,k+1) - 2.0*u(i,j,k) + u(i,j,k-1)) * idz2;
                                        viscterm += d2udz2;
                                    #endif

                                    viscterm *= Re_inv;

                                    // computing the convective term with the following references
                                    // considering a box around the point of stored data
                                    // F/B -> +- z-direction
                                    // N/S -> +- y-direction
                                    // E/W -> +- x-direction

                                    amrex::Real u_E = 0.5 * (u(i+1,j,k) + u(i,j,k));
                                    amrex::Real u_W = 0.5 * (u(i,j,k) + u(i-1,j,k));
                                    
                                    amrex::Real conv_x = (u_E * u_E - u_W * u_W) * idx;
                                    convterm = conv_x;

                                    #if AMREX_SPACEDIM >= 2
                                        amrex::Real u_N = 0.5 * (u(i,j+1,k) + u(i,j,k));
                                        amrex::Real u_S = 0.5 * (u(i,j,k) + u(i,j-1,k));

                                        amrex::Real v_N = 0.5 * (v(i,j+1,k) + v(i-1,j+1,k));
                                        amrex::Real v_S = 0.5 * (u(i,j,k) + u(i-1,j,k));

                                        amrex::Real conv_y = (u_N * v_N - u_S * v_S) * idy;
                                        convterm += conv_y;
                                    #endif
                                    
                                    #if AMREX_SPACEDIM == 3
                                        amrex::Real u_F = 0.5 * (u(i,j,k+1) + u(i,j,k))
                                        amrex::Real u_B = 0.5 * (u(i,j,k) + u(i,j,k-1))

                                        amrex::Real w_F = 0.5 * (w(i,j,k+1) + w(i-1,j,k+1));
                                        amrex::Real w_B = 0.5 * (w(i,j,k) + w(i-1,j,k));

                                        amrex::Real conv_z = (u_T * w_T - u_B * w_B) * idz;
                                        convterm += conv_z;
                                    #endif
                                }
                                else if constexpr (dim == 1)
                                {
                                    // computing the pressure term
                                    presterm = idy * (p(i,j,k) - p(i,j-1,k));

                                    // computing the viscous term dimension by dimension
                                    amrex::Real d2vdx2 = (v(i+1,j,k) - 2.0*v(i,j,k) + v(i-1,j,k)) * idx2;
                                    viscterm = d2vdx2;

                                    #if AMREX_SPACEDIM >= 2
                                        amrex::Real d2vdy2 = (v(i,j+1,k) - 2.0*v(i,j,k) + v(i,j-1,k)) * idy2;
                                        viscterm += d2vdy2;
                                    #endif

                                    #if AMREX_SPACEDIM == 3
                                        amrex::Real d2vdz2 = (v(i,j,k+1) - 2.0*v(i,j,k) + v(i,j,k-1)) * idz2;
                                        viscterm += d2vdz2;
                                    #endif

                                    viscterm *= Re_inv;

                                    // computing the convective term with the following references
                                    // considering a box around the point of stored data
                                    // F/B -> +- z-direction
                                    // N/S -> +- y-direction
                                    // E/W -> +- x-direction

                                    amrex::Real u_E = 0.5 * (u(i+1,j,k) + u(i,j,k));
                                    amrex::Real u_W = 0.5 * (u(i,j,k) + u(i-1,j,k));
                                    
                                    amrex::Real conv_x = (u_E * u_E - u_W * u_W) * idx;
                                    convterm = conv_x;

                                    #if AMREX_SPACEDIM >= 2
                                        amrex::Real u_N = 0.5 * (u(i,j+1,k) + u(i,j,k));
                                        amrex::Real u_S = 0.5 * (u(i,j,k) + u(i,j-1,k));

                                        amrex::Real v_N = 0.5 * (v(i,j+1,k) + v(i-1,j+1,k));
                                        amrex::Real v_S = 0.5 * (u(i,j,k) + u(i-1,j,k));

                                        amrex::Real conv_y = (u_N * v_N - u_S * v_S) * idy;
                                        convterm += conv_y;
                                    #endif
                                    
                                    #if AMREX_SPACEDIM == 3
                                        amrex::Real u_F = 0.5 * (u(i,j,k+1) + u(i,j,k))
                                        amrex::Real u_B = 0.5 * (u(i,j,k) + u(i,j,k-1))

                                        amrex::Real w_F = 0.5 * (w(i,j,k+1) + w(i-1,j,k+1));
                                        amrex::Real w_B = 0.5 * (w(i,j,k) + w(i-1,j,k));

                                        amrex::Real conv_z = (u_T * w_T - u_B * w_B) * idz;
                                        convterm += conv_z;
                                    #endif
                                }
                                elseif constexpr (dim == 2)
                                {
                                    // computing the pressure term
                                    presterm = idz * (p(i,j,k) - p(i,j,k-1));

                                    // computing the viscous term dimension by dimension
                                    amrex::Real d2wdx2 = (w(i+1,j,k) - 2.0*w(i,j,k) + w(i-1,j,k)) * idx2;
                                    viscterm = d2wdx2;

                                    #if AMREX_SPACEDIM >= 2
                                        amrex::Real d2wdy2 = (w(i,j+1,k) - 2.0*w(i,j,k) + w(i,j-1,k)) * idy2;
                                        viscterm += d2wdy2;
                                    #endif

                                    #if AMREX_SPACEDIM == 3
                                        amrex::Real d2wdz2 = (w(i,j,k+1) - 2.0*w(i,j,k) + w(i,j,k-1)) * idz2;
                                        viscterm += d2wdz2;
                                    #endif

                                    viscterm *= Re_inv;

                                    // computing the convective term with the following references
                                    // considering a box around the point of stored data
                                    // F/B -> +- z-direction
                                    // N/S -> +- y-direction
                                    // E/W -> +- x-direction

                                    amrex::Real u_E = 0.5 * (u(i+1,j,k) + u(i,j,k));
                                    amrex::Real u_W = 0.5 * (u(i,j,k) + u(i-1,j,k));
                                    
                                    amrex::Real conv_x = (u_E * u_E - u_W * u_W) * idx;
                                    convterm = conv_x;

                                    #if AMREX_SPACEDIM >= 2
                                        amrex::Real u_N = 0.5 * (u(i,j+1,k) + u(i,j,k));
                                        amrex::Real u_S = 0.5 * (u(i,j,k) + u(i,j-1,k));

                                        amrex::Real v_N = 0.5 * (v(i,j+1,k) + v(i-1,j+1,k));
                                        amrex::Real v_S = 0.5 * (u(i,j,k) + u(i-1,j,k));

                                        amrex::Real conv_y = (u_N * v_N - u_S * v_S) * idy;
                                        convterm += conv_y;
                                    #endif
                                    
                                    #if AMREX_SPACEDIM == 3
                                        amrex::Real u_F = 0.5 * (u(i,j,k+1) + u(i,j,k))
                                        amrex::Real u_B = 0.5 * (u(i,j,k) + u(i,j,k-1))

                                        amrex::Real w_F = 0.5 * (w(i,j,k+1) + w(i-1,j,k+1));
                                        amrex::Real w_B = 0.5 * (w(i,j,k) + w(i-1,j,k));

                                        amrex::Real conv_z = (u_T * w_T - u_B * w_B) * idz;
                                        convterm += conv_z;
                                    #endif
                                }

                                return (viscterm - presterm - convterm);
                            }
#endif