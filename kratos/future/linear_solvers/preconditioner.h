//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos::Future
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Preconditioner class.
/** Base class for preconditioners for linesr system solvers which defining
    standard interface for all the preconditioners derived from it.

    Considering a linear solver FooSolver with a method FooSolver::Solve.
    A typical code using this type of Preconditioners would be:

    \begin{verbatim}
     FooSolver::Solve(A,b,x,preconditioner)
     {
        preconditioner.Initialize(A,x,b);
    ...
    ...
    while(...) // Start iteration.
    {
            preconditioner.ApplyLeft(x);
        mult(a,x)
        preconditioner.ApplyRight(x)
    } // End iteration

    preconditioner.Finalize(A,x,b);
     }
     \end{verbatim}
*/
template<class TMatrixType, class TVectorType>
class Preconditioner
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Preconditioner
    KRATOS_CLASS_POINTER_DEFINITION(Preconditioner);

    using VectorType = TVectorType;

    using SparseMatrixType = TMatrixType;

    using DataType = typename SparseMatrixType::DataType;

    using DenseMatrixType = DenseMatrix<DataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Preconditioner() = default;

    /// Copy constructor.
    Preconditioner(const Preconditioner& Other) = default;

    /// Destructor.
    virtual ~Preconditioner() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Preconditioner& operator=(const Preconditioner& Other)
    {
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /** Preconditioner Initialize
    Initialize preconditioner for linear system rA*rX=rB
    @param rA  system matrix.
    @param rX Unknows vector
    @param rB Right side linear system of equations.
    */
    virtual void Initialize(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB)
    {
    }

    virtual void Initialize(
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB)
    {
        VectorType x(rX.size1());
        VectorType b(rB.size1());

        for (IndexType i = 0; i < rX.size1(); ++i) {
            x[i] = rX(i,0);
            b[i] = rB(i,0);
        }

        Initialize(rA, x, b);
    }

    /** This function is designed to be called every time the coefficients change in the system
    * that is, normally at the beginning of each solve.
    * For example if we are implementing a direct solver, this is the place to do the factorization
    * so that then the backward substitution can be performed effectively more than once
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    virtual void InitializeSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB)
    {
    }

    /** This function is designed to be called at the end of the solve step.
     * for example this is the place to remove any data that we do not want to save for later
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    virtual void FinalizeSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB)
    {
    }

    /** This function is designed to clean up all internal data in the solver.
     * Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    virtual void Clear()
    {
    }

    /** Some preconditioners may require a minimum degree of knowledge of the structure of the matrix. To make an example
    * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
    * another example is the automatic prescription of rotation null-space for smoothed-aggregation preconditioners
    * which require knowledge on the spatial position of the nodes associated to a given dof.
    * This function tells if the solver requires such data
    */
    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return false;
    }

    /** Some preconditioners may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation preconditioners
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part)
    {}

    virtual void Mult(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rY)
    {
        VectorType z = rX;
        ApplyRight(z);
        rA.SpMV(z, rY);
        ApplyLeft(rY);
    }

    virtual void TransposeMult(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rY)
    {
        VectorType z = rX;
        ApplyTransposeLeft(z);
        rA.TransposeSpMV(z, rY);
        ApplyTransposeRight(rY);
    }

    virtual VectorType& ApplyLeft(VectorType& rX)
    {
        return rX;
    }

    virtual VectorType& ApplyRight(VectorType& rX)
    {
        return rX;
    }

    /** Preconditioner transpose solver.
    Solving tranpose preconditioner system M^T*x=y, where m^T means transpose.
    @param rX  Unknows of preconditioner suystem
    */
    virtual VectorType& ApplyTransposeLeft(VectorType& rX)
    {
        return rX;
    }

    virtual VectorType& ApplyTransposeRight(VectorType& rX)
    {
        return rX;
    }

    virtual VectorType& ApplyInverseRight(VectorType& rX)
    {
        return rX;
    }

    /* The method Finalize is used to recover the value of rX.
       In principle, it is enough to multiply by the right preconditioner.
    See the diagoinal preconditioner for a nontrivial example. */
    virtual VectorType& Finalize(VectorType& rX)
    {
        return ApplyRight(rX);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Preconditioner";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Preconditioner";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; // Class Preconditioner

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<class TMatrixType, class TVectorType>
inline std::istream& operator >> (
    std::istream& IStream,
    Preconditioner<TMatrixType, TVectorType>& rThis)
    {
        return IStream;
    }

/// output stream function
template<class TMatrixType, class TVectorType>
inline std::ostream& operator << (
    std::ostream& OStream,
    const Preconditioner<TMatrixType, TVectorType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}

///@}

}  // namespace Kratos::Future.
