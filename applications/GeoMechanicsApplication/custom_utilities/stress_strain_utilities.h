// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    JMCarbonell,
//                   Vahid Galavi
//

#pragma once

/* Project includes */
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) StressStrainUtilities
{
public:
    static double CalculateVonMisesStress(const Vector& rStressVector);
    static double CalculateTrace(const Vector& rStressVector);
    static double CalculateMeanStress(const Vector& rStressVector);
    static double CalculateLodeAngle(const Vector& rStressVector);
    static double CalculateMohrCoulombShearCapacity(const Vector& rStressVector, double C, double Phi);
    static double CalculateMohrCoulombPressureCapacity(const Vector& rStressVector, double C, double Phi);
    static double CalculateVonMisesStrain(const Vector& rStrainVector);
    static Vector CalculateHenckyStrain(const Matrix& rDeformationGradient, size_t VoigtSize);
    static Matrix CalculateGreenLagrangeStrainTensor(const Matrix& rDeformationGradient);

private:
    static double CalculateQMohrCoulomb(const Vector& rStressVector, double C, double Phi);
    static double CalculateDenominator(const Vector& rStressVector, double Phi);
};

} // namespace Kratos