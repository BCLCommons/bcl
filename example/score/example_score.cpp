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
#include "score/bcl_score.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory_entry.h"
#include "score/bcl_score_aa_pair_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScore :
    public ExampleInterface
  {
  public:

    ExampleScore *Clone() const
    {
      return new ExampleScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // check if get namespace identifier is working
      BCL_MessageStd
      (
        "this is the bcl::score namespace identifier: " + score::GetNamespaceIdentifier()
      );
      BCL_Example_Check
      (
        score::GetNamespaceIdentifier() == "bcl::score",
        "namespace identifier for score is incorrect!"
      );

      // check if prepending of the histogram path is working
      const std::string file_name( score::Score::AddHistogramPath( score::AAPairDistance::GetDefaultHistogramFilename()));
      BCL_MessageStd
      (
        "this is a histogram filename with the prepended histogram path: " + file_name
      );

      BCL_Example_Check
      (
        io::DirectoryEntry( file_name).DoesExist(),
        "the file with prepended histogram path could not be found. Either file does not exist, given histogram path "
        "is wrong, or the working directory is incorrect: " + file_name
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScore

  const ExampleClass::EnumType ExampleScore::s_Instance
  (
    GetExamples().AddEnum( ExampleScore())
  );

} // namespace bcl

