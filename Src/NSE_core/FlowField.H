#ifndef FLOWFIELD_H_
#define FLOWFIELD_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>
#include <AMReX_Math.H>
#include <AMReX_Array.H>

using namespace amrex;

class FlowField {
private:
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> vel;
    amrex::MultiFab pres;

public:
    amrex::MultiFab tagRegion;
    FlowField(const amrex::Geometry& geom, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const int n_comp, const int n_ghost);
    FlowField(const FlowField& other);
    FlowField& operator=(const FlowField& other);

    // Dimension-agnostic getter! Pass 0 for X, 1 for Y, 2 for Z
    amrex::MultiFab& getVel(int idim);
    const amrex::MultiFab& getVel(int idim) const;
    amrex::MultiFab& getPres();
    const amrex::MultiFab& getPres() const;
    amrex::MultiFab& getTagRegion();
    const amrex::MultiFab& getTagRegion() const;

    void setBoundary(const amrex::Geometry& geom);
};

#endif