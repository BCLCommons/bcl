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
#include "fold/bcl_fold_default_scores.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_fixed_line_width_writer.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_default_scores.cpp
  //!
  //! @author weinerbe
  //! @date Oct 10, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldDefaultScores :
    public ExampleInterface
  {
  public:

    ExampleFoldDefaultScores *Clone() const
    {
      return new ExampleFoldDefaultScores( *this);
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
    ////////////////
    // operations //
    ////////////////

      // InitializeScores
      BCL_MessageStd
      (
        "Initial score enums: " + util::Format()( fold::GetScores().GetEnumCount())
      );
      fold::DefaultScores::GetInstance().InitializeScores();
      BCL_ExampleCheck( fold::GetScores().GetEnumCount() != 0, true);
      BCL_MessageStd
      (
        "Final score enums: " + util::Format()( fold::GetScores().GetEnumCount())
      );

      // ModifyScoreWeightSet
      fold::ScoreWeightSet weight_set;
      fold::DefaultScores::GetInstance().ModifyScoreWeightSet( weight_set);
      BCL_ExampleCheck( weight_set.GetWeightMap().IsEmpty(), false);
      BCL_MessageStd
      (
        "Added weights for " + util::Format()( weight_set.GetWeightMap().GetSize()) + " scores"
      );
      io::FixedLineWidthWriter flww;
      weight_set.WriteHelp( flww);
      BCL_Debug( flww.String());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldDefaultScores

  const ExampleClass::EnumType ExampleFoldDefaultScores::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldDefaultScores())
  );

} // namespace bcl
