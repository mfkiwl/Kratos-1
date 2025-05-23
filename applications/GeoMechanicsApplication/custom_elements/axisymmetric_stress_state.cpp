// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Marjan Fathian
//

#include "axisymmetric_stress_state.h"
#include "custom_utilities/element_utilities.hpp"
#include "includes/serializer.h"

namespace Kratos
{

Matrix AxisymmetricStressState::CalculateBMatrix(const Matrix&         rDN_DX,
                                                 const Vector&         rN,
                                                 const Geometry<Node>& rGeometry) const
{
    const auto radius = GeoElementUtilities::CalculateRadius(rN, rGeometry);

    const auto dimension       = rGeometry.WorkingSpaceDimension();
    const auto number_of_nodes = rGeometry.size();
    Matrix     result = ZeroMatrix(VOIGT_SIZE_2D_AXISYMMETRIC, dimension * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const IndexType index = dimension * i;

        result(INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X) = rDN_DX(i, INDEX_X);
        result(INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y) = rDN_DX(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_ZZ, index + INDEX_X) = rN[i] / radius;
        result(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X) = rDN_DX(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y) = rDN_DX(i, INDEX_X);
    }

    return result;
}

std::unique_ptr<StressStatePolicy> AxisymmetricStressState::Clone() const
{
    return std::make_unique<AxisymmetricStressState>();
}

Vector AxisymmetricStressState::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
{
    KRATOS_ERROR << "The calculation of Green Lagrange strain is not implemented for axisymmetric "
                    "configurations.\n";
}

const Vector& AxisymmetricStressState::GetVoigtVector() const { return VoigtVector2D; }

SizeType AxisymmetricStressState::GetVoigtSize() const { return GetVoigtSize2D(); }

SizeType AxisymmetricStressState::GetStressTensorSize() const { return GetStressTensorSize2D(); }

void AxisymmetricStressState::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void AxisymmetricStressState::load(Serializer&)
{
    // No data members to be loaded (yet)
}

} // namespace Kratos
