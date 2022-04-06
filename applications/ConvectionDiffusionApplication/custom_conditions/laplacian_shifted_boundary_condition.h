// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_LAPLACIAN_SHIFTED_BOUNDARY_CONDITION_H_INCLUDED)
#define  KRATOS_LAPLACIAN_SHIFTED_BOUNDARY_CONDITION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/condition.h"
#include "geometries/geometry.h"
#include "includes/variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) LaplacianShiftedBoundaryCondition: public Condition
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of LaplacianShiftedBoundaryCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LaplacianShiftedBoundaryCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    LaplacianShiftedBoundaryCondition(
        IndexType NewId,
        Geometry< Node<3> >::Pointer pGeometry);

    LaplacianShiftedBoundaryCondition(
        IndexType NewId,
        Geometry< Node<3> >::Pointer pGeometry,
        Properties::Pointer pProperties);

    LaplacianShiftedBoundaryCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes);

    /// Destructor.
    ~LaplacianShiftedBoundaryCondition() override;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& ConditionalDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

protected:

    ///@name Protected Life Cycle
    ///@{

    // Internal default constructor for serialization
    LaplacianShiftedBoundaryCondition();

    ///@}
    ///@name Protected Operations
    ///@{


    ///@}

private:

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    LaplacianShiftedBoundaryCondition& operator=(LaplacianShiftedBoundaryCondition const& rOther);

    /// Copy constructor.
    LaplacianShiftedBoundaryCondition(LaplacianShiftedBoundaryCondition const& rOther);

    ///@}

}; // Class LaplacianShiftedBoundaryCondition

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LAPLACIAN_SHIFTED_BOUNDARY_CONDITION_H_INCLUDED  defined