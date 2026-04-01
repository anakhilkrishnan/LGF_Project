#ifndef MYFUNCTIONS_H_
#define MYFUNCTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

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
};

struct ConsolidatedData 
{
    amrex::Vector<amrex::Real> data;
    amrex::Vector<FabMetaData> metadata;
};

using namespace amrex;

void test_print();
void extendedMain();

MeshData createMultiFab(int n_cell, int max_grid_size, int n_ghost, int n_comp);
void initializeMultiFab(MultiFab& phi);
ConsolidatedData consolidateMultiFab(const MultiFab& phifab);
void addEverySourceBox(const MultiFab& source, const MultiFab& target, MultiFab& result);

#endif