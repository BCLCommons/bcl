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
#include "model/bcl_model_objective_function_binary_operation.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_binary_operation_monitoring_dataset.cpp
  //!
  //! @author mendenjl
  //! @date Jul 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionBinaryOperation :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionBinaryOperation *Clone() const
    {
      return new ExampleModelObjectiveFunctionBinaryOperation( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // the binary objective will essentially always be constructed from a label, except during static initialization,
      // so only construct with a label here

      // string for the left hand side objective
      const std::string left_hand_side_objective( "Accuracy( cutoff=1.5)");
      const std::string operation( "-");
      const std::string right_hand_side_objective( "RMSD");
      // make a partial objective that takes the rmsd of the first two elements, weighted by 5
      util::Implementation< model::ObjectiveFunctionInterface> accuracy_minus_rmsd
      (
        "BinaryOperation( lhs=" + left_hand_side_objective + ", "
        + std::string( "op=") + operation + ", "
        + std::string( "rhs=") + right_hand_side_objective + ")"
      );
      util::Implementation< model::ObjectiveFunctionInterface> rmsd_minus_accuracy
      (
        "BinaryOperation( rhs=" + left_hand_side_objective + ", "
        + std::string( "op=") + operation + ", "
        + std::string( "lhs=") + right_hand_side_objective + ")"
      );
      util::Implementation< model::ObjectiveFunctionInterface> accuracy( left_hand_side_objective);
      util::Implementation< model::ObjectiveFunctionInterface> rmsd( right_hand_side_objective);

      // make sure the implementations could be created, otherwise the later checks would cause a crash
      BCL_ExampleAssert( accuracy_minus_rmsd.IsDefined(), true);
      BCL_ExampleAssert( accuracy.IsDefined(), true);
      BCL_ExampleAssert( rmsd.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // check the directionality: accuracy gets better with larger values, rmsd gets better with smaller values,
      // so accuracy minus rmsd gets better with larger values
      BCL_ExampleCheck
      (
        accuracy_minus_rmsd->GetImprovementType(),
        opti::e_LargerEqualIsBetter
      );
      // for the rmsd - accuracy, smaller values to be better
      BCL_ExampleCheck
      (
        rmsd_minus_accuracy->GetImprovementType(),
        opti::e_SmallerEqualIsBetter
      );

    ///////////////
    // operators //
    ///////////////

      // check the output of the function on identical datasets
      BCL_ExampleCheck
      (
        accuracy_minus_rmsd->operator()( experimental_dataset, experimental_dataset),
        1.0
      );
      BCL_ExampleCheck
      (
        rmsd_minus_accuracy->operator()( experimental_dataset, experimental_dataset),
        -1.0
      );

      // check the output of the function on disparate pairs of datasets
      BCL_ExampleCheck
      (
        accuracy_minus_rmsd->operator()( experimental_dataset, predicted_dataset),
        accuracy->operator()( experimental_dataset, predicted_dataset) -
        rmsd->operator()( experimental_dataset, predicted_dataset)
      );
      BCL_ExampleCheck
      (
        rmsd_minus_accuracy->operator()( experimental_dataset, predicted_dataset),
        rmsd->operator()( experimental_dataset, predicted_dataset) -
        accuracy->operator()( experimental_dataset, predicted_dataset)
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionBinaryOperation

  const ExampleClass::EnumType ExampleModelObjectiveFunctionBinaryOperation::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionBinaryOperation())
  );

} // namespace bcl
