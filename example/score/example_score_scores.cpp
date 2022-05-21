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
#include "score/bcl_score_scores.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_scores.h"

// external includes - sorted alphabetically

namespace bcl
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_scores.cpp
  //! @brief This example tests the implementation of the class score::Scores, which initializes all scores for
  //! evaluating protein models and adds them to the enumerate.
  //!
  //! @author fischea
  //! @date Aug 18, 2016
  //! @remarks status complete
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleScoreScores :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief clone function
    //! @return pointer to a new ExampleScoreScores
    ExampleScoreScores *Clone() const
    {
      return new ExampleScoreScores( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! @detail this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! get the single instance of this class
      score::Scores &scores( score::Scores::GetInstance());

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( scores.GetClassIdentifier(), ( GetStaticClassName< score::Scores>()));

    ////////////////
    // operations //
    ////////////////

      // initialize scoring functions and check if they were added to the enumerator
      BCL_MessageStd
      (
        "Initial score enums: " + util::Format()( fold::GetScores().GetEnumCount())
      );
      scores.Initialize();
      BCL_ExampleCheck( fold::GetScores().GetEnumCount() != 0, true);
      BCL_MessageStd
      (
        "Final score enums: " + util::Format()( fold::GetScores().GetEnumCount())
      );

      return 0;
    }

  }; // class ExampleScoreScores

  //! single instance of this class
  const ExampleClass::EnumType ExampleScoreScores::s_Instance
  (
     GetExamples().AddEnum( ExampleScoreScores())
  );

} // namespace bcl
