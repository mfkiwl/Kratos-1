//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//

/**
 *
 *                          WARNING! THIS FILE IS READ-ONLY
 *
 * This file has been auto-generated by the compressible navier stokes symbolic generator
 * located in the symbolic_generation directories of the FluidDynamicsApplication
 *
 * Any modifications to this file will be overwritten if and when that script is run again.
 *
 * In order to do any lasting changes, modify the template used by the script:
 * templates/compressible_navier_stokes_explicit_neumann_condition_template_2D2N.cpp
 * located in the symbolic_generation directories of the FluidDynamicsApplication.
 *
 * In order to change the formulation you will have to modify the script itself.
 */

// System includes


// External includes


// Project includes


// Application includes
#include "compressible_navier_stokes_explicit_neumann_condition.h"


namespace Kratos {

/**
 * Returns the integration method for computation of midpoint magnitudes.
 * Computation of RHS integration method is chosen in the symbolic generator.
 */
template<>
GeometryData::IntegrationMethod CompressibleNavierStokesExplicitNeumannCondition<2,2>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

template <>
void CompressibleNavierStokesExplicitNeumannCondition<2, 2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != DofSize) {
        rResult.resize(DofSize);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicitNeumannCondition<2, 2>::GetDofList(
    DofsVectorType& ConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (ConditionDofList.size() != DofSize) {
        ConditionDofList.resize(DofSize);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        ConditionDofList[local_index++] = r_geometry[i_node].pGetDof(DENSITY, den_pos);
        ConditionDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ConditionDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ConditionDofList[local_index++] = r_geometry[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template<>
BoundedVector<double, 8> CompressibleNavierStokesExplicitNeumannCondition<2,2>::CalculateRightHandSideInternal(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BoundedVector<double, BlockSize*NumNodes> rRightHandSideBoundedVector = ZeroVector(BlockSize*NumNodes);

    const auto data = ConditionData();

    rRightHandSideBoundedVector[0] = -0.78867513459481287*data.fluxes[0].density - 0.21132486540518713*data.fluxes[1].density;
    rRightHandSideBoundedVector[1] = -0.78867513459481287*data.fluxes[0].momentum(0) - 0.21132486540518713*data.fluxes[1].momentum(0);
    rRightHandSideBoundedVector[2] = -0.78867513459481287*data.fluxes[0].momentum(1) - 0.21132486540518713*data.fluxes[1].momentum(1);
    rRightHandSideBoundedVector[3] = -0.78867513459481287*data.fluxes[0].total_energy - 0.21132486540518713*data.fluxes[1].total_energy;
    rRightHandSideBoundedVector[4] = -0.21132486540518713*data.fluxes[0].density - 0.78867513459481287*data.fluxes[1].density;
    rRightHandSideBoundedVector[5] = -0.21132486540518713*data.fluxes[0].momentum(0) - 0.78867513459481287*data.fluxes[1].momentum(0);
    rRightHandSideBoundedVector[6] = -0.21132486540518713*data.fluxes[0].momentum(1) - 0.78867513459481287*data.fluxes[1].momentum(1);
    rRightHandSideBoundedVector[7] = -0.21132486540518713*data.fluxes[0].total_energy - 0.78867513459481287*data.fluxes[1].total_energy;


    rRightHandSideBoundedVector *= data.volume / ConditionDataStruct::NGauss;

    return rRightHandSideBoundedVector;
    KRATOS_CATCH("")
}


template class CompressibleNavierStokesExplicitNeumannCondition<2,2>;
using CompressibleNavierStokesExplicitNeumannCondition2D2N = CompressibleNavierStokesExplicitNeumannCondition<2,2>;

}