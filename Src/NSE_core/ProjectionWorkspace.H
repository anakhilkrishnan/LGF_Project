#ifndef PROJECTIONWORKSPACE_H_
#define PROJECTIONWORKSPACE_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_Math.H>
#include <AMReX_Array.H>

#include <FlowField.H>

using namespace amrex;

class ProjectionWorkspace
{
    private:
        amrex::Geometry geom;
        amrex::BoxArray ba;
        amrex::DistributionMapping dm;

        amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> rhs_vel;
        amrex::MultiFab divU;
        amrex::MultiFab corr_pres;
        

    public: 
        amrex::Real divU_max_norm;
        ProjectionWorkspace(const amrex::Geometry& geom_in, const amrex::BoxArray& ba_in, const amrex::DistributionMapping& dm_in, const int n_comp, const int n_ghost);
        void computeConvectiveFluxes(const FlowField& state);
        void predictVelocity(FlowField& state_n, FlowField& stage, amrex::Real dt, amrex::Real alpha, amrex::Real beta, amrex::Real gamma);
        void computePressure(FlowField& stage);
        void correctVelocity(FlowField& stage);
        void advanceTimeStep(FlowField& state_n, amrex::Real dt, int rk_order);
};

#endif