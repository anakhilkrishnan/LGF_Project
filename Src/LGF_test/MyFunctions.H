#ifndef MYFUNCTIONS_H_
#define MYFUNCTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>

#include <LGFCore.H>

using namespace amrex;

void test_print();
void extendedMain();

void initializeMultiFab(MultiFab& phi, amrex::Real init_val);
void initializeSourceMultiFab(MeshData& phi);
MeshData createMeshData(Array<Real, AMREX_SPACEDIM> phy_dom_lo, Array<Real, AMREX_SPACEDIM> phy_dom_hi, int n_cell, int max_grid_size, int n_ghost, int n_comp);

#endif