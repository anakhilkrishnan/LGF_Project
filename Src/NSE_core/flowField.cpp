#include <FlowField.H>

FlowField::FlowField(const amrex::Geometry& geom, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const int n_comp, const int n_ghost)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // convert the box array to face centered
        amrex::BoxArray ba_face = amrex::convert(ba, amrex::IntVect::TheDimensionVector(idim));

        // declare the specific velocity component
        vel[idim].define(ba_face, dm, n_comp, n_ghost);

        // initialize velocities upon creation
        vel[idim].setVal(0.0);
    }

    // initialize pressure upon creation
    pres.define(ba, dm, n_comp, n_ghost);
    pres.setVal(0.0);

    // initialize tagging regions
    tagRegion.define(ba, dm, n_comp, n_ghost);
    tagRegion.setVal(0.0);

}

FlowField::FlowField(const FlowField& other)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
        vel[idim].define(other.vel[idim].boxArray(), 
                         other.vel[idim].DistributionMap(), 
                         other.vel[idim].nComp(), 
                         other.vel[idim].nGrow());

        amrex::MultiFab::Copy(vel[idim], other.vel[idim], 0, 0, 1, vel[idim].nGrow());
    }

    pres.define(other.pres.boxArray(), other. pres.DistributionMap(), other.pres.nComp(), other.pres.nGrow());
    amrex::MultiFab::Copy(pres, other.pres, 0, 0, 1, pres.nGrow());

    tagRegion.define(other.tagRegion.boxArray(), other.tagRegion.DistributionMap(), other.tagRegion.nComp(), other.tagRegion.nGrow());
    amrex::MultiFab::Copy(tagRegion, other.tagRegion, 0, 0, 1, tagRegion.nGrow());
}

FlowField& FlowField::operator=(const FlowField& other) 
{
    if (this != &other) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            amrex::MultiFab::Copy(vel[idim], other.vel[idim], 0, 0, 1, vel[idim].nGrow());
        }
        
        amrex::MultiFab::Copy(pres, other.pres, 0, 0, 1, pres.nGrow());
        amrex::MultiFab::Copy(tagRegion, other.tagRegion, 0, 0, 1, tagRegion.nGrow());
    }
    return *this;
}

amrex::MultiFab& FlowField::getVel(int idim) 
{
    return vel[idim]; 
}

const amrex::MultiFab& FlowField::getVel(int idim) const 
{
    return vel[idim]; 
}

amrex::MultiFab& FlowField::getPres() 
{
    return pres; 
}

const amrex::MultiFab& FlowField::getPres() const 
{
    return pres; 
}

amrex::MultiFab& FlowField::getTagRegion() 
{
    return tagRegion; 
}

const amrex::MultiFab& FlowField::getTagRegion() const 
{
    return tagRegion; 
}

void FlowField::setBoundary(const amrex::Geometry& geom)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        vel[idim].FillBoundary(geom.periodicity());
    }

    pres.FillBoundary(geom.periodicity());
}

