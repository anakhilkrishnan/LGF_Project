#include <AMReX_BLProfiler.H>

#include <ProjectionWorkspace.H>
#include <spatialDiscretization.H>
#include <RKCoefficients.H>
#include <LGFCore.H>

ProjectionWorkspace::ProjectionWorkspace(const amrex::Geometry& geom_in, const amrex::BoxArray& ba_in, const amrex::DistributionMapping& dm_in, const int n_comp, const int n_ghost)
    : geom(geom_in), ba(ba_in), dm(dm_in)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // convert the box array to face centered
        amrex::BoxArray ba_face = amrex::convert(ba, amrex::IntVect::TheDimensionVector(idim));

        // declare the specific velocity component
        rhs_vel[idim].define(ba_face, dm, n_comp, n_ghost);

        // initialize velocities upon creation
        rhs_vel[idim].setVal(0.0);
    }

    // initialize divU upon creation
    divU.define(ba, dm, n_comp, n_ghost);
    divU.setVal(0.0);

    // initialize corr_pres upon creation
    corr_pres.define(ba, dm, n_comp, n_ghost);
    corr_pres.setVal(0.0);

    divU_max_norm = 0.0;
}

void ProjectionWorkspace::computeConvectiveFluxes(const FlowField& state)
{
    // compute the right hand side which is of the form
    // 1/Re(laplacian(u)) - grad(P) - u.divergence(u)
    // all taken at the n^th timestep
    // discretized using a second order finite difference KEP scheme
    // as outlined in Morinishi et al.

}

void ProjectionWorkspace::predictVelocity(FlowField& state_n, FlowField& stage, amrex::Real dt, amrex::Real alpha, amrex::Real beta, amrex::Real gamma)
{

    // use the right hand side to predict velocity at the next time step,
    // before enforcing divergence free condition

}

void ProjectionWorkspace::computePressure(FlowField& state)
{

    // use the custom lgf solver to compute the pressure at the next time step
}

void ProjectionWorkspace::correctVelocity(FlowField& state)
{
    // use the updated pressure to correct velocity to a divergence free field
}

void ProjectionWorkspace::advanceTimeStep(FlowField& state_n, amrex::Real dt, int rk_order)
{

    BL_PROFILE("<Compute> advanceTimeStep()");
    FlowField stage = state_n;
    amrex::Vector<RKCoeffs> coeffs = getRKCoeffs(rk_order);

    for(int k = 0; k < rk_order; ++k)
    {
        // extracting RK coefficients
        amrex::Real alpha = coeffs[k].alp;
        amrex::Real beta = coeffs[k].bet;
        amrex::Real gamma = coeffs[k].gam;

        // compute and store fluxes in workspace
        computeConvectiveFluxes(stage);

        // compute predicted velocity without divergence free condition
        // store predicted velocity within stage
        // find divergence of predicted velocity, store in workspace
        predictVelocity(state_n, stage, dt, alpha, beta, gamma);

        // use custom LGF solver to find pressure correction
        // store in workspace 
        computePressure(stage);

        // use corr_pres to correct velocity field and store back into stage
        // update stage to reflect p_stage + corr_p
        correctVelocity(stage);

    }

    state_n = stage;
}