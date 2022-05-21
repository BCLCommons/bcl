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
#include "model/bcl_model_dtree_roc_data_partition_function.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"
#include "math/bcl_math_template_instantiations.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_dtree_roc_data_partition_function.cpp
  //!
  //! @author lemmonwa, mendenjl
  //! @date Aug 26, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDtreeRocDataPartitionFunction :
    public ExampleInterface
  {
  public:

    ExampleModelDtreeRocDataPartitionFunction *Clone() const
    {
      return new ExampleModelDtreeRocDataPartitionFunction( *this);
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
      const size_t feature_size( 3), result_size( 1), number_features( 10);

      // create a training dataset
      linal::Matrix< float> features( number_features, feature_size), results( number_features, result_size);
      linal::Matrix< size_t> results_state( number_features, result_size);
      storage::Vector< model::FeatureResultAndState> feature_result_and_state( number_features);

      const double activity_cutoff( 0);
      for( size_t counter( 0); counter < number_features; ++counter)
      {
        features( counter, 0) = std::cos( double( counter));
        features( counter, 1) = std::sin( double( counter));
        features( counter, 2) = double( counter) + 2;
        results( counter, 0) = std::sin( double( counter));
        results_state( counter, 0) = results( counter, 0) <= activity_cutoff;
        feature_result_and_state( counter) =
          model::FeatureResultAndState
          (
            model::FeatureReference< float>( feature_size, features[ counter]),
            model::FeatureReference< float>( result_size, results[ counter]),
            model::FeatureReference< size_t>( result_size, results_state[ counter])
          );
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // check default constructor
      model::DtreeRocDataPartitionFunction default_function;

      // check constructor with parameters
      model::DtreeRocDataPartitionFunction function_parameter( 0.5, 1.0);

      // check clone
      util::ShPtr< model::DtreeRocDataPartitionFunction> function_clone( function_parameter.Clone());

    ///////////////
    // operators //
    ///////////////

      // make data partition function
      util::ShPtr< model::DtreeRocDataPartitionFunction> function( new model::DtreeRocDataPartitionFunction( 0.5, 1.0));

      //get partitioned data, which should have split on descriptor at index 1.
      model::DtreeBinaryPartition partition( function->operator ()( feature_result_and_state));

      BCL_ExampleCheck( partition.GetFeatureIndex(), size_t( 1));

    //////////////////////
    // input and output //
    //////////////////////

      // write bcl object
      WriteBCLObject( *function);
      // create default object
      model::DtreeRocDataPartitionFunction function_read;
      // read bcl object
      ReadBCLObject( function_read);
      BCL_ExampleIndirectCheck( function_read( feature_result_and_state).GetFeatureIndex(), size_t( 1), "I/O");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; // ExampleClass::ExampleModelDtreeRocDataPartitionFunction

  const ExampleClass::EnumType ExampleModelDtreeRocDataPartitionFunction::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDtreeRocDataPartitionFunction())
  );

} // namespace bcl
