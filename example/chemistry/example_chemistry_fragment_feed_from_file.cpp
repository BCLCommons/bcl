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
#include "chemistry/bcl_chemistry_fragment_feed_from_file.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_feed_from_file.cpp
  //!
  //! @author mendenjl
  //! @date May 25, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentFeedFromFile :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentFeedFromFile *Clone() const
    {
      return new ExampleChemistryFragmentFeedFromFile( *this);
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
      // example sdf filename
      const std::string mglur5_filename( AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));

    //////////////////
    // construction //
    //////////////////

      // create the feed
      chemistry::FragmentFeedFromFile ffff;

    /////////////////
    // data access //
    /////////////////

      // check that the file can be opened
      BCL_ExampleCheck( ffff.Open( mglur5_filename), true);

      // check the size
      BCL_ExampleCheck( ffff.GetSize(), 5);

      // manually open the file and get the molecules out

      // initialize input stream
      io::IFStream read;
      // input stream for diazepam
      BCL_ExampleMustOpenInputFile( read, mglur5_filename);
      // create shptr on ensemble and save molecules from sdf file
      chemistry::FragmentEnsemble ensemble( read); // TODO change this, this relies on FragmentFeed
      io::File::CloseClearFStream( read);

    ////////////////
    // operations //
    ////////////////

      BCL_ExampleIndirectCheck( ensemble.GetSize(), 5, "FragmentEnsemble used for checking FragmentFeedFromFile read molecules appropriately");

      // now walk through the molecules loaded in and the feed
      size_t number_mols( 0);

      chemistry::FragmentEnsemble::const_iterator itr_ensemble( ensemble.Begin());
      sdf::MdlHandler handler;
      while( ffff.RetrieveNextMolecule( handler))
      {
        // make the handler into a fragment
        chemistry::FragmentComplete fragment( sdf::FragmentFactory::MakeFragment( handler));

        // check that the atom info is identical
        BCL_ExampleIndirectAssert
        (
          fragment.GetAtomInfo(),
          itr_ensemble->GetAtomInfo(),
          "ffff.RetrieveNextMolecule( handler)"
        );

        // check that the bond info is identical
        BCL_ExampleIndirectAssert
        (
          fragment.GetBondInfo(),
          itr_ensemble->GetBondInfo(),
          "ffff.RetrieveNextMolecule( handler)"
        );

        ++itr_ensemble;
        ++number_mols;
      }

      // check that the fragment feed read all five molecules
      BCL_ExampleIndirectCheck( number_mols, 5, "fragment feed read all five molecules");

      // test skip
      BCL_ExampleIndirectCheck( ffff.Open( mglur5_filename), true, "reopen the same feed");

      // skip the first molecule
      BCL_ExampleCheck( ffff.Skip(), true);

      // read the second
      ffff.RetrieveNextMolecule( handler);

      BCL_ExampleAssert( ensemble.IsEmpty(), false);

      // test atom info between the feed and the ensemble
      BCL_ExampleIndirectCheck
      (
        chemistry::FragmentComplete( sdf::FragmentFactory::MakeFragment( handler)).GetAtomInfo(),
        ( ++ensemble.GetMolecules().Begin())->GetAtomInfo(),
        "ffff.Skip()"
      );

      // try opening a non-existant file
      BCL_ExampleCheck( ffff.Open( "not-a-file"), false);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryFragmentFeedFromFile

  const ExampleClass::EnumType ExampleChemistryFragmentFeedFromFile::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentFeedFromFile())
  );

} // namespace bcl
