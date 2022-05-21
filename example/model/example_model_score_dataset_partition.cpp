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
#include "model/bcl_model_score_dataset_partition.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_retrieve_data_set_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_score_dataset_partition.cpp
  //!
  //! @author mendenjl
  //! @date Apr 06, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelScoreDatasetPartition :
    public ExampleInterface
  {
  public:

    ExampleModelScoreDatasetPartition *Clone() const
    {
      return new ExampleModelScoreDatasetPartition( *this);
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
      // create objects using the implementation framework, which is the only way to set this objects members
      util::Implementation< model::ScoreDatasetInterface>
        info_gain( "Partition(partitioner=InformationGain,cutoff=0)");

      // ensure it is defined
      BCL_ExampleAssert
      (
        util::Implementation< model::ScoreDatasetInterface>
        (
          "Partition(partitioner=InformationGain,cutoff=0)"
        ).IsDefined(),
        true
      );

      // load in the example score dataset
      const std::string model_directory( AddExampleInputPathToFilename( e_Model, ""));
      util::Implementation< model::RetrieveDataSetBase>
        dataset_retriever( "File(filename=" + model_directory + "example_data_set_score.bcl)");

      util::ShPtr< descriptor::Dataset> dataset( dataset_retriever->GenerateDataSet());

      // Calculate the info_gains

      // the first column should get practically 0, since it has the same distribution for results above & below the cutoff
      BCL_ExampleCheckWithinAbsTolerance( info_gain->Score( *dataset)( 0), float( 0.0), float( 1.0e-5));

      const float expected_info_gain_col_b( 15.5348);
      // the second column should get 15.5348
      BCL_ExampleCheckWithinTolerance( info_gain->Score( *dataset)( 1), expected_info_gain_col_b, float( 1.0e-4));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelScoreDatasetPartition

  const ExampleClass::EnumType ExampleModelScoreDatasetPartition::s_Instance
  (
    GetExamples().AddEnum( ExampleModelScoreDatasetPartition())
  );

} // namespace bcl
