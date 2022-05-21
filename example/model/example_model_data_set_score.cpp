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
#include "model/bcl_model_data_set_score.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_data_set_score.cpp
  //!
  //! @author mendenjl
  //! @date Apr 04, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDataSetScore :
    public ExampleInterface
  {
  public:

    ExampleModelDataSetScore *Clone() const
    {
      return new ExampleModelDataSetScore( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // make sure that the default score is empty
      BCL_ExampleCheck( model::DataSetScore().GetFeatures().GetString(), model::FeatureLabelSet().GetString());
      BCL_ExampleCheck( model::DataSetScore().GetScores(), linal::Vector< float>());

      model::FeatureLabelSet labels;
      labels.PushBack( util::ObjectDataLabel( "Example"), 5);

      linal::Vector< float> scores( 5);
      scores( 0) = 1.2;
      scores( 1) = 2.4;
      scores( 2) = 3.6;
      scores( 3) = 4.8;
      scores( 4) = 6.0;

    /////////////////
    // data access //
    /////////////////

      // try set functions
      model::DataSetScore score;
      score.SetFeatures( labels);

      // make sure that the score has the correct values
      BCL_ExampleIndirectCheck( score.GetFeatures().GetString(), labels.GetString(), "SetFeatures");

      score.SetScores( scores);
      BCL_ExampleIndirectCheck( score.GetScores(), scores, "SetScores");

      // end
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDataSetScore

  const ExampleClass::EnumType ExampleModelDataSetScore::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDataSetScore())
  );

} // namespace bcl
