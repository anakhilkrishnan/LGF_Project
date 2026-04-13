#include <AMReX_MultiFab.H>

#include <MyFunctions.H>

using namespace amrex;

MeshData createMeshData(Array<Real, AMREX_SPACEDIM> phy_dom_lo, Array<Real, AMREX_SPACEDIM> phy_dom_hi, int n_cell, int max_grid_size, int n_ghost, int n_comp)
{

    BoxArray ba;
    Geometry geom;

    // Define the computational domain
    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
    Box domain(dom_lo, dom_hi);

    // Define the periodicity
    Vector<int> is_periodic(AMREX_SPACEDIM, 0); // 0 = not periodic

    // Initialize the boxarray "ba" from the single box "bx"
    ba.define(domain);
    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box(phy_dom_lo, phy_dom_hi);

    // This defines a Geometry object
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // Create multifab
    MultiFab mfans(ba, dm, n_comp, n_ghost);

    // Return reslt
    return {std::move(mfans), geom};
}
