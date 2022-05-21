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
#include "chemistry/bcl_chemistry_fragment_split_linear_fragments.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_split_linear_fragments.cpp
  //! @details this tests the implementation of ECFP-like fragment splitting of molecules
  //!
  //! @author geanesar
  //! @date Sep 22, 2015
  //! @remarks status complete 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentSplitLinearFragments :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentSplitLinearFragments *Clone() const
    { return new ExampleChemistryFragmentSplitLinearFragments( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      chemistry::FragmentSplitLinearFragments linear_frags;

    /////////////////
    // data access //
    /////////////////

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      chemistry::FragmentEnsemble linear_ensemble( linear_frags( ensemble.GetMolecules().FirstElement()));

      const std::string linear_split( AddExampleOutputPathToFilename( linear_ensemble, "split_linear.sdf"));

      io::OFStream output;

      BCL_ExampleMustOpenOutputFile( output, linear_split);
      linear_ensemble.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename_linear( linear_split + ".correct");

      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( linear_split, correct_filename_linear),
          true,
          "Linear fragments were generated correctly"
        )
      )
      {
        remove( linear_split.c_str());
      }

      return 0;
    }

    static const ExampleClass::EnumType ExampleChemistryFragmentSplitLinearFragments_Instance;

  }; //end ExampleChemistryFragmentSplitLinearFragments

  const ExampleClass::EnumType ExampleChemistryFragmentSplitLinearFragments::ExampleChemistryFragmentSplitLinearFragments_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentSplitLinearFragments())
  );

} // namespace bcl
