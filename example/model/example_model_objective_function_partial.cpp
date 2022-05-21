// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "model/bcl_model_objective_function_partial.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_partial_monitoring_dataset.cpp
  //!
  //! @author mendenjl
  //! @date Jul 22, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionPartial :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionPartial *Clone() const
    {
      return new ExampleModelObjectiveFunctionPartial( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // create 3 arrays with different largest columns
      float col_a[ 3] = { 6.0, 3.0, 1.0 };
      float col_b[ 3] = { 3.0, 6.0, 1.0 };
      float col_c[ 3] = { 1.0, 3.0, 6.0 };

      // make a matrix of experimental data
      linal::Matrix< float> experimental_data( 3, 3);
      experimental_data.ReplaceRow( 0, linal::Vector< float>( 3, col_a));
      experimental_data.ReplaceRow( 1, linal::Vector< float>( 3, col_b));
      experimental_data.ReplaceRow( 2, linal::Vector< float>( 3, col_c));
      model::FeatureDataSet< float> experimental_dataset( experimental_data);

      // create a matrix of predicted data
      linal::Matrix< float> predicted_data( 3, 3);
      predicted_data.ReplaceRow( 0, linal::Vector< float>( 3, col_c));
      predicted_data.ReplaceRow( 1, linal::Vector< float>( 3, col_a));
      predicted_data.ReplaceRow( 2, linal::Vector< float>( 3, col_b));
      model::FeatureDataSet< float> predicted_dataset( predicted_data);

      // MSD of col_a to col_b, first two elements
      const float msd_ab_1_2
      (
        linal::Vector< float>
        (
          linal::Vector< float>( 2, col_a) - linal::Vector< float>( 2, col_b)
        ).SquareNorm()
      );
      // MSD of col_a to col_b, second & third elements
      const float msd_ab_2_3
      (
        linal::Vector< float>
        (
          linal::Vector< float>( 2, col_a + 1) - linal::Vector< float>( 2, col_b + 1)
        ).SquareNorm()
      );
      // MSD of col_a to col_c, first two elements
      const float msd_ac_1_2
      (
        linal::Vector< float>
        (
          linal::Vector< float>( 2, col_a) - linal::Vector< float>( 2, col_c)
        ).SquareNorm()
      );
      // MSD of col_a to col_c, second & third elements
      const float msd_ac_2_3
      (
        linal::Vector< float>
        (
          linal::Vector< float>( 2, col_a + 1) - linal::Vector< float>( 2, col_c + 1)
        ).SquareNorm()
      );
      // MSD of col_b to col_c, first two elements
      const float msd_bc_1_2
      (
        linal::Vector< float>
        (
          linal::Vector< float>( 2, col_b) - linal::Vector< float>( 2, col_c)
        ).SquareNorm()
      );
      // MSD of col_b to col_c, second & third elements
      const float msd_bc_2_3
      (
        linal::Vector< float>
        (
          linal::Vector< float>( 2, col_b + 1) - linal::Vector< float>( 2, col_c + 1)
        ).SquareNorm()
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // the partial objective will essentially always be constructed from a label, except during static initialization,
      // so only construct with a label here

      // make a partial objective that takes the rmsd of the first two elements, weighted by 5
      util::Implementation< model::ObjectiveFunctionInterface> rmsd_1_to_2_weight_p5p1
      (
        "PartialObjective( function=RMSD, outputs=\"[0,1]\", weight=5.1)"
      );
      // make a partial objective that takes the rmsd of the second and third element, weighted by -2
      // as always, the spacing doesn't matter
      util::Implementation< model::ObjectiveFunctionInterface> rmsd_2_to_3_weight_neg2
      (
        "PartialObjective( function=RMSD,weight = -2.0, outputs=\"[1,2]\"  )"
      );

      // make sure the implementations could be created, otherwise the later checks would cause a crash
      BCL_ExampleAssert( rmsd_1_to_2_weight_p5p1.IsDefined(), true);
      BCL_ExampleAssert( rmsd_2_to_3_weight_neg2.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // check the directionality: for the positively weighted partial objective, smaller is better, just like rmsd
      BCL_ExampleCheck
      (
        rmsd_1_to_2_weight_p5p1->GetImprovementType(),
        opti::e_SmallerEqualIsBetter
      );
      // for the negatively weighted objective, larger is better
      BCL_ExampleCheck
      (
        rmsd_2_to_3_weight_neg2->GetImprovementType(),
        opti::e_LargerEqualIsBetter
      );

    ///////////////
    // operators //
    ///////////////

      // check the output of the function on identical datasets
      BCL_ExampleCheck
      (
        rmsd_1_to_2_weight_p5p1->operator()( experimental_dataset, experimental_dataset),
        0.0
      );
      BCL_ExampleCheck
      (
        rmsd_2_to_3_weight_neg2->operator()( experimental_dataset, experimental_dataset),
        0.0
      );

      // check the output of the function on disparate pairs of datasets
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          rmsd_1_to_2_weight_p5p1->operator()( experimental_dataset, predicted_dataset),
          float( 5.1) * math::Sqrt( float( msd_ab_1_2 + msd_ac_1_2 + msd_bc_1_2) / float( 3.0) / float( 2.0))
        ),
        true
      );
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          rmsd_2_to_3_weight_neg2->operator()( experimental_dataset, predicted_dataset),
          float( -2.0) * math::Sqrt( float( msd_ab_2_3 + msd_ac_2_3 + msd_bc_2_3) / float( 3.0) / float( 2.0))
        ),
        true
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionPartial

  const ExampleClass::EnumType ExampleModelObjectiveFunctionPartial::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionPartial())
  );

} // namespace bcl
