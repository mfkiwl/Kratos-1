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
 * templates/compressible_navier_stokes_explicit_condition_template_2D2N.cpp
 * located in the symbolic_generation directories of the FluidDynamicsApplication.
 *
 * In order to change the formulation you will have to modify the script itself.
 */

// System includes


// External includes


// Project includes


// Application includes
#include "compressible_navier_stokes_explicit_condition.h"


namespace Kratos {

/**
 * Returns the integration method for computation of midpoint magnitudes.
 * Computation of RHS integration method is chosen in the symbolic generator.
 */
template<>
GeometryData::IntegrationMethod CompressibleNavierStokesExplicitCondition<2,2>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

template <>
void CompressibleNavierStokesExplicitCondition<2, 2>::EquationIdVector(
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
void CompressibleNavierStokesExplicitCondition<2, 2>::GetDofList(
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
BoundedVector<double, 8> CompressibleNavierStokesExplicitCondition<2,2>::CalculateRightHandSideInternal(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BoundedVector<double, BlockSize*NumNodes> rRightHandSideBoundedVector = ZeroVector(BlockSize*NumNodes);

    const auto data = ConditionData(rCurrentProcessInfo);

    const double crRightHandSideBoundedVector0 = 0.21132486540518713*data.alpha_sc_nodes(0) + 0.78867513459481287*data.alpha_sc_nodes(1);
    const double crRightHandSideBoundedVector1 = 0.21132486540518713*data.U(0,1) + 0.78867513459481287*data.U(1,1);
    const double crRightHandSideBoundedVector2 = crRightHandSideBoundedVector0*data.gradients[1].density(0) + crRightHandSideBoundedVector1;
    const double crRightHandSideBoundedVector3 = 0.21132486540518713*data.unit_normal(0);
    const double crRightHandSideBoundedVector4 = 0.78867513459481287*data.alpha_sc_nodes(0) + 0.21132486540518713*data.alpha_sc_nodes(1);
    const double crRightHandSideBoundedVector5 = 0.78867513459481287*data.U(0,1) + 0.21132486540518713*data.U(1,1);
    const double crRightHandSideBoundedVector6 = crRightHandSideBoundedVector4*data.gradients[0].density(0) + crRightHandSideBoundedVector5;
    const double crRightHandSideBoundedVector7 = 0.78867513459481287*data.unit_normal(0);
    const double crRightHandSideBoundedVector8 = 0.21132486540518713*data.U(0,2) + 0.78867513459481287*data.U(1,2);
    const double crRightHandSideBoundedVector9 = crRightHandSideBoundedVector0*data.gradients[1].density(1) + crRightHandSideBoundedVector8;
    const double crRightHandSideBoundedVector10 = 0.21132486540518713*data.unit_normal(1);
    const double crRightHandSideBoundedVector11 = 0.78867513459481287*data.U(0,2) + 0.21132486540518713*data.U(1,2);
    const double crRightHandSideBoundedVector12 = crRightHandSideBoundedVector11 + crRightHandSideBoundedVector4*data.gradients[0].density(1);
    const double crRightHandSideBoundedVector13 = 0.78867513459481287*data.unit_normal(1);
    const double crRightHandSideBoundedVector14 = 0.21132486540518713*data.U(0,0) + 0.78867513459481287*data.U(1,0);
    const double crRightHandSideBoundedVector15 = 1.0/crRightHandSideBoundedVector14;
    const double crRightHandSideBoundedVector16 = crRightHandSideBoundedVector1*crRightHandSideBoundedVector15;
    const double crRightHandSideBoundedVector17 = 0.26794919243112275*data.U(0,0) + data.U(1,0);
    const double crRightHandSideBoundedVector18 = pow(crRightHandSideBoundedVector17, -2);
    const double crRightHandSideBoundedVector19 = 1.6076951545867364*crRightHandSideBoundedVector18;
    const double crRightHandSideBoundedVector20 = data.mu + 0.21132486540518713*data.mu_sc_nodes(0) + 0.78867513459481287*data.mu_sc_nodes(1);
    const double crRightHandSideBoundedVector21 = crRightHandSideBoundedVector20*(-crRightHandSideBoundedVector1*data.gradients[1].density(1) + crRightHandSideBoundedVector14*data.gradients[1].momentum(0,1) + crRightHandSideBoundedVector14*data.gradients[1].momentum(1,0) - crRightHandSideBoundedVector8*data.gradients[1].density(0));
    const double crRightHandSideBoundedVector22 = crRightHandSideBoundedVector19*crRightHandSideBoundedVector21;
    const double crRightHandSideBoundedVector23 = -crRightHandSideBoundedVector16*crRightHandSideBoundedVector8 + crRightHandSideBoundedVector22;
    const double crRightHandSideBoundedVector24 = 0.78867513459481287*data.U(0,0) + 0.21132486540518713*data.U(1,0);
    const double crRightHandSideBoundedVector25 = 1.0/crRightHandSideBoundedVector24;
    const double crRightHandSideBoundedVector26 = crRightHandSideBoundedVector25*crRightHandSideBoundedVector5;
    const double crRightHandSideBoundedVector27 = data.U(0,0) + 0.26794919243112275*data.U(1,0);
    const double crRightHandSideBoundedVector28 = pow(crRightHandSideBoundedVector27, -2);
    const double crRightHandSideBoundedVector29 = 1.6076951545867364*crRightHandSideBoundedVector28;
    const double crRightHandSideBoundedVector30 = data.mu + 0.78867513459481287*data.mu_sc_nodes(0) + 0.21132486540518713*data.mu_sc_nodes(1);
    const double crRightHandSideBoundedVector31 = crRightHandSideBoundedVector30*(-crRightHandSideBoundedVector11*data.gradients[0].density(0) + crRightHandSideBoundedVector24*data.gradients[0].momentum(0,1) + crRightHandSideBoundedVector24*data.gradients[0].momentum(1,0) - crRightHandSideBoundedVector5*data.gradients[0].density(1));
    const double crRightHandSideBoundedVector32 = crRightHandSideBoundedVector29*crRightHandSideBoundedVector31;
    const double crRightHandSideBoundedVector33 = -crRightHandSideBoundedVector11*crRightHandSideBoundedVector26 + crRightHandSideBoundedVector32;
    const double crRightHandSideBoundedVector34 = pow(0.26794919243112275*data.U(0,1) + data.U(1,1), 2);
    const double crRightHandSideBoundedVector35 = 0.62200846792814624*crRightHandSideBoundedVector15;
    const double crRightHandSideBoundedVector36 = crRightHandSideBoundedVector1*data.gradients[1].density(0);
    const double crRightHandSideBoundedVector37 = crRightHandSideBoundedVector14*data.gradients[1].momentum(0,0);
    const double crRightHandSideBoundedVector38 = crRightHandSideBoundedVector36 - crRightHandSideBoundedVector37;
    const double crRightHandSideBoundedVector39 = 3.2153903091734728*crRightHandSideBoundedVector20;
    const double crRightHandSideBoundedVector40 = crRightHandSideBoundedVector18*crRightHandSideBoundedVector39;
    const double crRightHandSideBoundedVector41 = data.gamma - 1;
    const double crRightHandSideBoundedVector42 = pow(0.26794919243112275*data.U(0,2) + data.U(1,2), 2);
    const double crRightHandSideBoundedVector43 = 0.21132486540518713*data.U(0,3) + 0.78867513459481287*data.U(1,3);
    const double crRightHandSideBoundedVector44 = crRightHandSideBoundedVector41*(-crRightHandSideBoundedVector15*(0.31100423396407312*crRightHandSideBoundedVector34 + 0.31100423396407312*crRightHandSideBoundedVector42) + crRightHandSideBoundedVector43);
    const double crRightHandSideBoundedVector45 = crRightHandSideBoundedVector8*data.gradients[1].density(1);
    const double crRightHandSideBoundedVector46 = crRightHandSideBoundedVector14*data.gradients[1].momentum(1,1);
    const double crRightHandSideBoundedVector47 = crRightHandSideBoundedVector45 - crRightHandSideBoundedVector46;
    const double crRightHandSideBoundedVector48 = 0.21132486540518713*data.beta_sc_nodes(0);
    const double crRightHandSideBoundedVector49 = 0.78867513459481287*data.beta_sc_nodes(1);
    const double crRightHandSideBoundedVector50 = 0.66666666666666663*data.mu;
    const double crRightHandSideBoundedVector51 = 0.14088324360345808*data.mu_sc_nodes(0);
    const double crRightHandSideBoundedVector52 = 0.52578342306320858*data.mu_sc_nodes(1);
    const double crRightHandSideBoundedVector53 = -crRightHandSideBoundedVector19*(crRightHandSideBoundedVector38 + crRightHandSideBoundedVector47)*(-crRightHandSideBoundedVector48 - crRightHandSideBoundedVector49 + crRightHandSideBoundedVector50 + crRightHandSideBoundedVector51 + crRightHandSideBoundedVector52) + crRightHandSideBoundedVector44;
    const double crRightHandSideBoundedVector54 = crRightHandSideBoundedVector34*crRightHandSideBoundedVector35 + crRightHandSideBoundedVector38*crRightHandSideBoundedVector40 + crRightHandSideBoundedVector53;
    const double crRightHandSideBoundedVector55 = pow(data.U(0,1) + 0.26794919243112275*data.U(1,1), 2);
    const double crRightHandSideBoundedVector56 = 0.62200846792814624*crRightHandSideBoundedVector25;
    const double crRightHandSideBoundedVector57 = crRightHandSideBoundedVector5*data.gradients[0].density(0);
    const double crRightHandSideBoundedVector58 = crRightHandSideBoundedVector24*data.gradients[0].momentum(0,0);
    const double crRightHandSideBoundedVector59 = crRightHandSideBoundedVector57 - crRightHandSideBoundedVector58;
    const double crRightHandSideBoundedVector60 = 3.2153903091734728*crRightHandSideBoundedVector30;
    const double crRightHandSideBoundedVector61 = crRightHandSideBoundedVector28*crRightHandSideBoundedVector60;
    const double crRightHandSideBoundedVector62 = pow(data.U(0,2) + 0.26794919243112275*data.U(1,2), 2);
    const double crRightHandSideBoundedVector63 = 0.78867513459481287*data.U(0,3) + 0.21132486540518713*data.U(1,3);
    const double crRightHandSideBoundedVector64 = crRightHandSideBoundedVector41*(-crRightHandSideBoundedVector25*(0.31100423396407312*crRightHandSideBoundedVector55 + 0.31100423396407312*crRightHandSideBoundedVector62) + crRightHandSideBoundedVector63);
    const double crRightHandSideBoundedVector65 = crRightHandSideBoundedVector11*data.gradients[0].density(1);
    const double crRightHandSideBoundedVector66 = crRightHandSideBoundedVector24*data.gradients[0].momentum(1,1);
    const double crRightHandSideBoundedVector67 = crRightHandSideBoundedVector65 - crRightHandSideBoundedVector66;
    const double crRightHandSideBoundedVector68 = 0.78867513459481287*data.beta_sc_nodes(0);
    const double crRightHandSideBoundedVector69 = 0.21132486540518713*data.beta_sc_nodes(1);
    const double crRightHandSideBoundedVector70 = 0.52578342306320858*data.mu_sc_nodes(0);
    const double crRightHandSideBoundedVector71 = 0.14088324360345808*data.mu_sc_nodes(1);
    const double crRightHandSideBoundedVector72 = -crRightHandSideBoundedVector29*(crRightHandSideBoundedVector59 + crRightHandSideBoundedVector67)*(crRightHandSideBoundedVector50 - crRightHandSideBoundedVector68 - crRightHandSideBoundedVector69 + crRightHandSideBoundedVector70 + crRightHandSideBoundedVector71) + crRightHandSideBoundedVector64;
    const double crRightHandSideBoundedVector73 = crRightHandSideBoundedVector55*crRightHandSideBoundedVector56 + crRightHandSideBoundedVector59*crRightHandSideBoundedVector61 + crRightHandSideBoundedVector72;
    const double crRightHandSideBoundedVector74 = crRightHandSideBoundedVector35*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector40*crRightHandSideBoundedVector47 + crRightHandSideBoundedVector53;
    const double crRightHandSideBoundedVector75 = crRightHandSideBoundedVector56*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector61*crRightHandSideBoundedVector67 + crRightHandSideBoundedVector72;
    const double crRightHandSideBoundedVector76 = crRightHandSideBoundedVector43 + crRightHandSideBoundedVector44;
    const double crRightHandSideBoundedVector77 = 1.6076951545867364*crRightHandSideBoundedVector8;
    const double crRightHandSideBoundedVector78 = 1.6076951545867364*crRightHandSideBoundedVector43;
    const double crRightHandSideBoundedVector79 = 1.6076951545867364*crRightHandSideBoundedVector1;
    const double crRightHandSideBoundedVector80 = 1.6076951545867364*crRightHandSideBoundedVector14;
    const double crRightHandSideBoundedVector81 = 1.267949192431123/crRightHandSideBoundedVector17;
    const double crRightHandSideBoundedVector82 = crRightHandSideBoundedVector81*data.gradients[1].density(0);
    const double crRightHandSideBoundedVector83 = 1.0/data.c_v;
    const double crRightHandSideBoundedVector84 = crRightHandSideBoundedVector18*crRightHandSideBoundedVector83*(0.21132486540518713*data.lamb_sc_nodes(0) + 0.78867513459481287*data.lamb_sc_nodes(1) + data.lambda);
    const double crRightHandSideBoundedVector85 = -crRightHandSideBoundedVector36 + crRightHandSideBoundedVector37;
    const double crRightHandSideBoundedVector86 = -crRightHandSideBoundedVector50;
    const double crRightHandSideBoundedVector87 = -crRightHandSideBoundedVector45 + crRightHandSideBoundedVector46;
    const double crRightHandSideBoundedVector88 = 1.6076951545867364*(crRightHandSideBoundedVector85 + crRightHandSideBoundedVector87)*(crRightHandSideBoundedVector48 + crRightHandSideBoundedVector49 - crRightHandSideBoundedVector51 - crRightHandSideBoundedVector52 + crRightHandSideBoundedVector86);
    const double crRightHandSideBoundedVector89 = crRightHandSideBoundedVector15*crRightHandSideBoundedVector18*crRightHandSideBoundedVector21*crRightHandSideBoundedVector77 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector18*(crRightHandSideBoundedVector39*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector88) - crRightHandSideBoundedVector16*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector84*(crRightHandSideBoundedVector34*crRightHandSideBoundedVector82 + crRightHandSideBoundedVector42*crRightHandSideBoundedVector82 - crRightHandSideBoundedVector77*data.gradients[1].momentum(0,1) - crRightHandSideBoundedVector78*data.gradients[1].density(0) - crRightHandSideBoundedVector79*data.gradients[1].momentum(0,0) + crRightHandSideBoundedVector80*data.gradients[1].total_energy(0));
    const double crRightHandSideBoundedVector90 = crRightHandSideBoundedVector63 + crRightHandSideBoundedVector64;
    const double crRightHandSideBoundedVector91 = 1.6076951545867364*crRightHandSideBoundedVector11;
    const double crRightHandSideBoundedVector92 = 1.6076951545867364*crRightHandSideBoundedVector63;
    const double crRightHandSideBoundedVector93 = 1.6076951545867364*crRightHandSideBoundedVector5;
    const double crRightHandSideBoundedVector94 = 1.6076951545867364*crRightHandSideBoundedVector24;
    const double crRightHandSideBoundedVector95 = 1.267949192431123/crRightHandSideBoundedVector27;
    const double crRightHandSideBoundedVector96 = crRightHandSideBoundedVector95*data.gradients[0].density(0);
    const double crRightHandSideBoundedVector97 = crRightHandSideBoundedVector28*crRightHandSideBoundedVector83*(0.78867513459481287*data.lamb_sc_nodes(0) + 0.21132486540518713*data.lamb_sc_nodes(1) + data.lambda);
    const double crRightHandSideBoundedVector98 = -crRightHandSideBoundedVector57 + crRightHandSideBoundedVector58;
    const double crRightHandSideBoundedVector99 = -crRightHandSideBoundedVector65 + crRightHandSideBoundedVector66;
    const double crRightHandSideBoundedVector100 = 1.6076951545867364*(crRightHandSideBoundedVector98 + crRightHandSideBoundedVector99)*(crRightHandSideBoundedVector68 + crRightHandSideBoundedVector69 - crRightHandSideBoundedVector70 - crRightHandSideBoundedVector71 + crRightHandSideBoundedVector86);
    const double crRightHandSideBoundedVector101 = crRightHandSideBoundedVector25*crRightHandSideBoundedVector28*crRightHandSideBoundedVector31*crRightHandSideBoundedVector91 + crRightHandSideBoundedVector26*crRightHandSideBoundedVector28*(crRightHandSideBoundedVector100 + crRightHandSideBoundedVector60*crRightHandSideBoundedVector98) - crRightHandSideBoundedVector26*crRightHandSideBoundedVector90 + crRightHandSideBoundedVector97*(crRightHandSideBoundedVector55*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector62*crRightHandSideBoundedVector96 - crRightHandSideBoundedVector91*data.gradients[0].momentum(0,1) - crRightHandSideBoundedVector92*data.gradients[0].density(0) - crRightHandSideBoundedVector93*data.gradients[0].momentum(0,0) + crRightHandSideBoundedVector94*data.gradients[0].total_energy(0));
    const double crRightHandSideBoundedVector102 = crRightHandSideBoundedVector15*crRightHandSideBoundedVector8;
    const double crRightHandSideBoundedVector103 = crRightHandSideBoundedVector81*data.gradients[1].density(1);
    const double crRightHandSideBoundedVector104 = crRightHandSideBoundedVector102*crRightHandSideBoundedVector18*(crRightHandSideBoundedVector39*crRightHandSideBoundedVector87 + crRightHandSideBoundedVector88) - crRightHandSideBoundedVector102*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector22 + crRightHandSideBoundedVector84*(crRightHandSideBoundedVector103*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector103*crRightHandSideBoundedVector42 - crRightHandSideBoundedVector77*data.gradients[1].momentum(1,1) - crRightHandSideBoundedVector78*data.gradients[1].density(1) - crRightHandSideBoundedVector79*data.gradients[1].momentum(1,0) + crRightHandSideBoundedVector80*data.gradients[1].total_energy(1));
    const double crRightHandSideBoundedVector105 = crRightHandSideBoundedVector11*crRightHandSideBoundedVector25;
    const double crRightHandSideBoundedVector106 = crRightHandSideBoundedVector95*data.gradients[0].density(1);
    const double crRightHandSideBoundedVector107 = crRightHandSideBoundedVector105*crRightHandSideBoundedVector28*(crRightHandSideBoundedVector100 + crRightHandSideBoundedVector60*crRightHandSideBoundedVector99) - crRightHandSideBoundedVector105*crRightHandSideBoundedVector90 + crRightHandSideBoundedVector26*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector97*(crRightHandSideBoundedVector106*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector106*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector91*data.gradients[0].momentum(1,1) - crRightHandSideBoundedVector92*data.gradients[0].density(1) - crRightHandSideBoundedVector93*data.gradients[0].momentum(1,0) + crRightHandSideBoundedVector94*data.gradients[0].total_energy(1));
    rRightHandSideBoundedVector[0] = -crRightHandSideBoundedVector10*crRightHandSideBoundedVector9 - crRightHandSideBoundedVector12*crRightHandSideBoundedVector13 - crRightHandSideBoundedVector2*crRightHandSideBoundedVector3 - crRightHandSideBoundedVector6*crRightHandSideBoundedVector7;
    rRightHandSideBoundedVector[1] = crRightHandSideBoundedVector10*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector33 - crRightHandSideBoundedVector3*crRightHandSideBoundedVector54 - crRightHandSideBoundedVector7*crRightHandSideBoundedVector73;
    rRightHandSideBoundedVector[2] = -crRightHandSideBoundedVector10*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector13*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector7;
    rRightHandSideBoundedVector[3] = crRightHandSideBoundedVector10*crRightHandSideBoundedVector104 + crRightHandSideBoundedVector101*crRightHandSideBoundedVector7 + crRightHandSideBoundedVector107*crRightHandSideBoundedVector13 + crRightHandSideBoundedVector3*crRightHandSideBoundedVector89;
    rRightHandSideBoundedVector[4] = -crRightHandSideBoundedVector10*crRightHandSideBoundedVector12 - crRightHandSideBoundedVector13*crRightHandSideBoundedVector9 - crRightHandSideBoundedVector2*crRightHandSideBoundedVector7 - crRightHandSideBoundedVector3*crRightHandSideBoundedVector6;
    rRightHandSideBoundedVector[5] = crRightHandSideBoundedVector10*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector23 - crRightHandSideBoundedVector3*crRightHandSideBoundedVector73 - crRightHandSideBoundedVector54*crRightHandSideBoundedVector7;
    rRightHandSideBoundedVector[6] = -crRightHandSideBoundedVector10*crRightHandSideBoundedVector75 - crRightHandSideBoundedVector13*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector7 + crRightHandSideBoundedVector3*crRightHandSideBoundedVector33;
    rRightHandSideBoundedVector[7] = crRightHandSideBoundedVector10*crRightHandSideBoundedVector107 + crRightHandSideBoundedVector101*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector13 + crRightHandSideBoundedVector7*crRightHandSideBoundedVector89;


    rRightHandSideBoundedVector *= data.volume / ConditionDataStruct::NGauss;

    return rRightHandSideBoundedVector;
    KRATOS_CATCH("")
}


template class CompressibleNavierStokesExplicitCondition<2,2>;
using CompressibleNavierStokesExplicitCondition2D2N = CompressibleNavierStokesExplicitCondition<2,2>;

}