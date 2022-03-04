// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#include "includes/model_part.h"
#include "custom_processes/set_cylindrical_local_axes_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SetCylindricalLocalAxesProcess::SetCylindricalLocalAxesProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetCylindricalLocalAxesProcess::ExecuteInitialize()
{
    KRATOS_TRY
    const array_1d<double, 3>& r_generatrix_axis  = mThisParameters["cylindrical_generatrix_axis"].GetVector();
    const array_1d<double, 3>& r_generatrix_point = mThisParameters["cylindrical_generatrix_point"].GetVector();

    KRATOS_ERROR_IF(MathUtils<double>::Norm3(r_generatrix_axis) < std::numeric_limits<double>::epsilon()) << "The r_generatrix_axis has norm zero" << std::endl;

    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
        array_1d<double, 3> local_axis_1;
        const array_1d<double, 3> coords = rElement.GetGeometry().Center();
        const double c = -r_generatrix_axis[0] * coords[0] - r_generatrix_axis[1] * coords[1] - r_generatrix_axis[2] * coords[2];
        const double lambda = -(r_generatrix_axis[0] * r_generatrix_point[0] + r_generatrix_axis[1] * r_generatrix_point[1] + r_generatrix_axis[2] * r_generatrix_point[2] + c) / (std::pow(r_generatrix_axis[0], 2) + std::pow(r_generatrix_axis[1], 2) + std::pow(r_generatrix_axis[2], 2));
        array_1d<double, 3> intersection;
        noalias(intersection) = r_generatrix_point + lambda * r_generatrix_axis;

        noalias(local_axis_1) = coords - intersection;
        ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_1);
        rElement.SetValue(LOCAL_AXIS_1, local_axis_1);

        if (mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3) {
            array_1d<double, 3> local_axis_2;
            noalias(local_axis_2) = r_generatrix_axis;
            ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_2);
            rElement.SetValue(LOCAL_AXIS_2, local_axis_2);
        }
    });

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SetCylindricalLocalAxesProcess::ExecuteInitializeSolutionStep()
{
    if (mThisParameters["update_at_each_step"].GetBool()) {
        auto &r_process_info = mrThisModelPart.GetProcessInfo();
        block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
            std::vector <Matrix> F;
            // Here we update the local axis according to dx = F*dX
            rElement.CalculateOnIntegrationPoints(DEFORMATION_GRADIENT, F, r_process_info);

            auto local_axis_1 = rElement.GetValue(LOCAL_AXIS_1);
            auto local_axis_2 = rElement.GetValue(LOCAL_AXIS_2);

            local_axis_1 = prod(F[0], local_axis_1);
            local_axis_2 = prod(F[0], local_axis_2);

            ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_1);
            ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_2);

            rElement.SetValue(LOCAL_AXIS_1, local_axis_1);
            if (r_process_info[DOMAIN_SIZE] == 3)
                rElement.SetValue(LOCAL_AXIS_2, local_axis_2);

        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SetCylindricalLocalAxesProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY
    if (mThisParameters["update_at_each_step"].GetBool()) {
        const array_1d<double, 3>& r_generatrix_axis  = mThisParameters["cylindrical_generatrix_axis"].GetVector();
        const array_1d<double, 3>& r_generatrix_point = mThisParameters["cylindrical_generatrix_point"].GetVector();

        KRATOS_ERROR_IF(MathUtils<double>::Norm3(r_generatrix_axis) < std::numeric_limits<double>::epsilon()) << "The r_generatrix_axis has norm zero" << std::endl;

        block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
            array_1d<double, 3> local_axis_1;
            // now we compute the original center of the element
            const auto &r_geometry = rElement.GetGeometry();
            const std::size_t points_number = r_geometry.size();
            array_1d<double, 3> coords;
            for (int i = 1 ; i < points_number ; i++) {
                coords[0] += r_geometry[i].X0();
                coords[1] += r_geometry[i].Y0();
                coords[2] += r_geometry[i].Z0();
            }
            coords /= double(points_number);

            const double c = -r_generatrix_axis[0] * coords[0] - r_generatrix_axis[1] * coords[1] - r_generatrix_axis[2] * coords[2];
            const double lambda = -(r_generatrix_axis[0] * r_generatrix_point[0] + r_generatrix_axis[1] * r_generatrix_point[1] + r_generatrix_axis[2] * r_generatrix_point[2] + c) / (std::pow(r_generatrix_axis[0], 2) + std::pow(r_generatrix_axis[1], 2) + std::pow(r_generatrix_axis[2], 2));
            array_1d<double, 3> intersection;
            noalias(intersection) = r_generatrix_point + lambda * r_generatrix_axis;

            noalias(local_axis_1) = coords - intersection;
            ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_1);
            rElement.SetValue(LOCAL_AXIS_1, local_axis_1);

            if (mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3) {
                array_1d<double, 3> local_axis_2;
                noalias(local_axis_2) = r_generatrix_axis;
                ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_2);
                rElement.SetValue(LOCAL_AXIS_2, local_axis_2);
            }
        });
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetCylindricalLocalAxesProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "cylindrical_generatrix_axis"   : [0.0,0.0,1.0],
        "cylindrical_generatrix_point"  : [0.0,0.0,0.0],
        "update_at_each_step"           : false
    })");

    return default_parameters;
}

// class SetCylindricalLocalAxesProcess
} // namespace Kratos.
