//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                     PNavas $
//   Last modified by:    $Co-Author:               LMonforte $
//   Date:                $Date:                 October 2017 $
//   Revision:            $Revision:                     -0.1 $
//
//
//   Implementation of the Fluid Saturated porous media in a u-w Formulation

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_updated_lagrangian_U_J_W_wP_DME_element.hpp"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   AxisymUpdatedLagrangianUJWwPDMEElement::AxisymUpdatedLagrangianUJWwPDMEElement()
      : AxisymUpdatedLagrangianUJWwPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPDMEElement::AxisymUpdatedLagrangianUJWwPDMEElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : AxisymUpdatedLagrangianUJWwPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPDMEElement::AxisymUpdatedLagrangianUJWwPDMEElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : AxisymUpdatedLagrangianUJWwPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPDMEElement::AxisymUpdatedLagrangianUJWwPDMEElement( AxisymUpdatedLagrangianUJWwPDMEElement const& rOther)
      :AxisymUpdatedLagrangianUJWwPElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPDMEElement&  AxisymUpdatedLagrangianUJWwPDMEElement::operator=(AxisymUpdatedLagrangianUJWwPDMEElement const& rOther)
   {
      AxisymUpdatedLagrangianUJWwPElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUJWwPDMEElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new AxisymUpdatedLagrangianUJWwPDMEElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUJWwPDMEElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      AxisymUpdatedLagrangianUJWwPDMEElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

      //-----------//

      NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


      if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
         NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

         if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }

      for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
         NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

      //-----------//

      if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
         NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

      for(unsigned int i=0; i<mDeformationGradientF0.size(); i++)
      {
         NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
      }

      NewElement.mDeterminantF0 = mDeterminantF0;

      NewElement.SetData(this->GetData());
      NewElement.SetFlags(this->GetFlags());

      return Element::Pointer( new AxisymUpdatedLagrangianUJWwPDMEElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJWwPDMEElement::~AxisymUpdatedLagrangianUJWwPDMEElement()
   {
   }

   // *********************************************************************************
   //         Calculate the Damping matrix part due to the stabilization
   void AxisymUpdatedLagrangianUJWwPDMEElement::CalculateAndAddDampingStabilizationMatrix( MatrixType & rDampingMatrix, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = 2*dimension + 2;

      double ElementSize = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         double aux = 0;
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            aux += rVariables.DN_DX(i, iDim);
         }
         ElementSize += fabs( aux); 
      }
      ElementSize *= sqrt( double(dimension) );
      ElementSize = 4.0/ ElementSize;

      ProcessInfo SomeProcessInfo; 
      std::vector< double> Mmodulus;
      this->CalculateOnIntegrationPoints( M_MODULUS, Mmodulus, SomeProcessInfo);
      double ConstrainedModulus = Mmodulus[0];
      if ( ConstrainedModulus < 1e-5)
      {
         const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
         const double& nu    = GetProperties()[POISSON_RATIO];
         ConstrainedModulus =  YoungModulus * ( 1.0-nu)/(1.0+nu) / (1.0-2.0*nu);
      }


      double density_mixture0 = GetProperties()[DENSITY];
      double WaterDensity =GetProperties().GetValue(DENSITY_WATER);
      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);

      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0; 
      double density_solid = (density_mixture0 - porosity0*WaterDensity) / ( 1.0 - porosity0);
      double CurrentDensity = ( 1.0 - porosity) * density_solid + porosity * WaterDensity;

      double StabFactor = GetProperties().GetValue( STABILIZATION_FACTOR_WP);
      double tau_factor=0;
      double alpha_factor=0;
      if (WaterDensity > 1e-3 && CurrentDensity > 1e-3){
         //   double WaveVelocity = sqrt((ConstrainedModulus+rWaterBulk)/CurrentDensity);
         double WaveVelocity = sqrt((ConstrainedModulus)/CurrentDensity);
         alpha_factor = StabFactor/CurrentDensity*pow(ElementSize/WaveVelocity,2);
         tau_factor = (CurrentDensity-WaterDensity)*porosity/WaterDensity/(porosity-1);
      }

      Matrix SmallMatrix = ZeroMatrix( number_of_nodes, number_of_nodes);

      if ( ( fabs(StabFactor) > 1.0e-9) && dimension==2)  {

         // assume that K is a unit tensor by k
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               for (unsigned p = 0; p < dimension; p++) {
                  SmallMatrix(i,j) += rVariables.DN_DX(i,p) * rVariables.DN_DX(j,p) * alpha_factor*(tau_factor-1) * rIntegrationWeight;
               	  //if (p==0)
                  //	SmallMatrix( i, j ) += rVariables.N[i] * rVariables.N[j] / rVariables.CurrentRadius /rVariables.CurrentRadius * alpha_factor*tau_factor*rIntegrationWeight;
               }
            }
         } 
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               rDampingMatrix( (i+1)*dofs_per_node-1, (j+1)*dofs_per_node-1) -=SmallMatrix(i,j);
            }
         }
      } 


      KRATOS_CATCH("")
   }

   // *********************************************************************************
   //         Calculate the MASS matrix part due to the stabilization   
   void AxisymUpdatedLagrangianUJWwPDMEElement::CalculateAndAddMassStabilizationMatrix( MatrixType & rMassMatrix, ElementDataType & rVariables, double & rIntegrationWeight)
   {

      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = 2*dimension + 2;

      double ElementSize = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         double aux = 0;
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            aux += rVariables.DN_DX(i, iDim);
         }
         ElementSize += fabs( aux); 
      }
      ElementSize *= sqrt( double(dimension) );
      ElementSize = 4.0/ ElementSize;

      ProcessInfo SomeProcessInfo; 
      std::vector< double> Mmodulus;
      this->CalculateOnIntegrationPoints( M_MODULUS, Mmodulus, SomeProcessInfo);
      double ConstrainedModulus = Mmodulus[0];
      if ( ConstrainedModulus < 1e-5)
      {
         const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
         const double& nu    = GetProperties()[POISSON_RATIO];
         ConstrainedModulus =  YoungModulus * ( 1.0-nu)/(1.0+nu) / (1.0-2.0*nu);
      }


      double density_mixture0 = GetProperties()[DENSITY];
      double WaterDensity =GetProperties().GetValue(DENSITY_WATER);
      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);

      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0; 
      double density_solid = (density_mixture0 - porosity0*WaterDensity) / ( 1.0 - porosity0);
      double CurrentDensity = ( 1.0 - porosity) * density_solid + porosity * WaterDensity;

      double StabFactor = GetProperties().GetValue( STABILIZATION_FACTOR_WP);
      double alpha_factor=0;
      double tau_factor=0;
      if (WaterDensity > 1e-3 && CurrentDensity > 1e-3){
         double WaveVelocity = sqrt((ConstrainedModulus)/CurrentDensity);
         alpha_factor = StabFactor/CurrentDensity*pow(ElementSize/WaveVelocity,2);
         tau_factor = (CurrentDensity-WaterDensity)*porosity/WaterDensity/(porosity-1);
      }

      const double & rPermeability = GetProperties()[PERMEABILITY_WATER];

      if ( ( fabs(StabFactor) > 1.0e-9) && dimension==2)  {

         Matrix Q = ZeroMatrix( number_of_nodes, dimension*number_of_nodes);
         unsigned int voigtSize = 4;
         Matrix m = ZeroMatrix( 1, voigtSize);
         for ( unsigned int i = 0; i < 3; i++)
            m(0,i) = 1.0;
         Matrix partial =  prod( m, rVariables.B );
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < dimension*number_of_nodes; j++ ) {
               Q(i,j) = rVariables.N[i] * partial(0,j) * rIntegrationWeight;
            }
         }

         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               for (unsigned int jDim = 0; jDim < dimension; jDim++) {
                  rMassMatrix( (i+1)*dofs_per_node -1, j*dofs_per_node + jDim + (dimension+1) ) -= alpha_factor*tau_factor/rPermeability*Q(i, j*dimension+jDim);
               }
            }
         }

      }


      KRATOS_CATCH("")
   }


   void AxisymUpdatedLagrangianUJWwPDMEElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUJWwPElement )
   }
   void AxisymUpdatedLagrangianUJWwPDMEElement::load( Serializer& rSerializer )

   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUJWwPElement )
   }

}  // END KRATOS NAMESPACE
