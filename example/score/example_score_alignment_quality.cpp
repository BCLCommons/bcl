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
#include "score/bcl_score_alignment_quality.h"

// includes from bcl - sorted alphabetically
#include "quality/bcl_quality_measures.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_alignment_quality.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author woetzen
  //! @date Mar 23, 2012
  //! @remarks status empty
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAlignmentQuality :
    public ExampleInterface
  {
  public:

    ExampleScoreAlignmentQuality *Clone() const
    {
      return new ExampleScoreAlignmentQuality( *this);
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

      // default constructor
      score::AlignmentQuality score_alignment;

    /////////////////
    // data access //
    /////////////////

      // set the measure
      score_alignment.SetMeasure( quality::GetMeasures().e_DME);

      // set the atoms
      score_alignment.SetAtoms( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA));

    ///////////////
    // operators //
    ///////////////

      // score an alignment

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( score_alignment);

      // read the object back in
      score::AlignmentQuality score_alignment_read;
      ReadBCLObject( score_alignment_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAlignmentQuality

  const ExampleClass::EnumType ExampleScoreAlignmentQuality::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAlignmentQuality())
  );

} // namespace bcl
