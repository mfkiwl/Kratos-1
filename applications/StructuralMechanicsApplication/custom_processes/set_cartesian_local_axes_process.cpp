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
#include "custom_processes/set_cartesian_local_axes_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SetCartesianLocalAxesProcess::SetCartesianLocalAxesProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetCartesianLocalAxesProcess::ExecuteInitialize()
{
    KRATOS_TRY

    array_1d<double, 3> local_axis_1;

    if (mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3) {
        array_1d<double, 3> local_axis_2;
        const Matrix& cartesian_local_axes_matrix = mThisParameters["cartesian_local_axis"].GetMatrix();
        local_axis_1[0] = cartesian_local_axes_matrix(0, 0);
        local_axis_1[1] = cartesian_local_axes_matrix(0, 1);
        local_axis_1[2] = cartesian_local_axes_matrix(0, 2);
        local_axis_2[0] = cartesian_local_axes_matrix(1, 0);
        local_axis_2[1] = cartesian_local_axes_matrix(1, 1);
        local_axis_2[2] = cartesian_local_axes_matrix(1, 2);
        ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_1);
        ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_2);

        block_for_each(mrThisModelPart.Elements(), [&](Element& rElement) {
            rElement.SetValue(LOCAL_AXIS_1, local_axis_1);
            rElement.SetValue(LOCAL_AXIS_2, local_axis_2);
        });
    } else if (mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2) {
        const Vector& cartesian_local_axes_matrix = mThisParameters["cartesian_local_axis"].GetVector();
        local_axis_1[0] = cartesian_local_axes_matrix[0];
        local_axis_1[1] = cartesian_local_axes_matrix[1];
        local_axis_1[2] = cartesian_local_axes_matrix[2];
        ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(local_axis_1);
        block_for_each(mrThisModelPart.Elements(), [&](Element& rElement) {
            rElement.SetValue(LOCAL_AXIS_1, local_axis_1);
        });
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SetCartesianLocalAxesProcess::ExecuteInitializeSolutionStep()
{
    if (mThisParameters["update_at_each_step"].GetBool()) {
        auto &r_process_info = mrThisModelPart.GetProcessInfo();
        block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
            std::vector <Matrix> F;
            // Here we update the local axis according to dx = F*dX
            rElement.CalculateOnIntegrationPoints(DEFORMATION_GRADIENT, F, r_process_info);

            auto& r_local_axis_1 = rElement.GetValue(LOCAL_AXIS_1);
            r_local_axis_1 = prod(F[0], r_local_axis_1);
            ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(r_local_axis_1);

            if (r_process_info[DOMAIN_SIZE] == 3) {
                auto& r_local_axis_2 = rElement.GetValue(LOCAL_AXIS_2);
                r_local_axis_2 = prod(F[0], r_local_axis_2);
                ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(r_local_axis_2);
            }
        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SetCartesianLocalAxesProcess::ExecuteFinalizeSolutionStep()
{
    if (mThisParameters["update_at_each_step"].GetBool()) {
        // we have to reset them to be correctly modified in the
        // ExecuteInitializeSolutionStep
        ExecuteInitialize();
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetCartesianLocalAxesProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "cartesian_local_axis"          : [[1.0,0.0,0.0],[0.0,1.0,0.0]],
        "update_at_each_step"           : false
    })" );

    return default_parameters;
}

// class SetCartesianLocalAxesProcess
} // namespace Kratos.
