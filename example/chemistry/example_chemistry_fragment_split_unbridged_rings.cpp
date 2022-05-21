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
#include "chemistry/bcl_chemistry_fragment_split_unbridged_rings.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_split_unbridged_rings.cpp
  //!
  //! @author kothiwsk
  //! @date Oct 23, 2014
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentSplitUnbridgedRings :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentSplitUnbridgedRings *Clone() const
    { return new ExampleChemistryFragmentSplitUnbridgedRings( *this);}

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
    { /*

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      chemistry::FragmentSplitUnbridgedRings split_rings( true);

      chemistry::FragmentSplitUnbridgedRings chains( false);

    /////////////////
    // data access //
    /////////////////

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      chemistry::FragmentEnsemble ring_ensemble( split_rings( ensemble.GetMolecules().FirstElement()));
      chemistry::FragmentEnsemble chain_ensemble( chains( ensemble.GetMolecules().FirstElement()));

      const std::string rings_split( AddExampleOutputPathToFilename( ring_ensemble, "split_rings.sdf"));
      const std::string chain_split( AddExampleOutputPathToFilename( chain_ensemble, "split_chains.sdf"));

      io::OFStream output;

      BCL_ExampleMustOpenOutputFile( output, rings_split);
      ring_ensemble.WriteMDL( output);
      io::File::CloseClearFStream( output);

      BCL_ExampleMustOpenOutputFile( output, chain_split);
      chain_ensemble.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename_rings( rings_split + ".correct");
      std::string correct_filename_chains( chain_split + ".correct");

      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( rings_split, correct_filename_rings),
          true,
          "graph was highlighted correctly"
        )
        &&
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( chain_split, correct_filename_chains),
          true,
          "rings and chains isolatedy"
        )
      )
      {
        remove( rings_split.c_str());
        remove( chain_split.c_str());
      }

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

      */ return 0;
    }

    static const ExampleClass::EnumType ExampleChemistryFragmentSplitUnbridgedRings_Instance;

  }; //end ExampleChemistryFragmentSplitUnbridgedRings

  const ExampleClass::EnumType ExampleChemistryFragmentSplitUnbridgedRings::ExampleChemistryFragmentSplitUnbridgedRings_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentSplitUnbridgedRings())
  );

} // namespace bcl
