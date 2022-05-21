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
#include "chemistry/bcl_chemistry_fragment_feed_from_stream.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_feed_from_stream.cpp
  //!
  //! @author mendenjl
  //! @date May 20, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentFeedFromStream :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentFeedFromStream *Clone() const
    {
      return new ExampleChemistryFragmentFeedFromStream( *this);
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

      // input stream for diazepam
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, mglur5_filename);
      chemistry::FragmentEnsemble ensemble( read);
      io::File::CloseClearFStream( read);

      // initialize input stream
      // input stream for diazepam
      BCL_ExampleMustOpenInputFile( read, mglur5_filename);
      // create the feed
      chemistry::FragmentFeedFromStream fffs( read);

    /////////////////
    // data access //
    /////////////////

      // check that the file can be opened
      BCL_ExampleCheck( fffs.Open( mglur5_filename), true);

      // check the size
      BCL_ExampleCheck( fffs.GetSize(), 5);

      // manually open the file and get the molecules out

    ////////////////
    // operations //
    ////////////////

      // now walk through the molecules loaded in and the feed
      size_t number_mols( 0);
      chemistry::FragmentEnsemble::const_iterator itr_ensemble( ensemble.Begin());
      sdf::MdlHandler handler;
      while( fffs.RetrieveNextMolecule( handler))
      {
        // make the handler into a fragment
        chemistry::FragmentComplete fragment( sdf::FragmentFactory::MakeFragment( handler));

        // check that the atom info is identical
        BCL_ExampleIndirectAssert
        (
          fragment.GetAtomInfo(),
          itr_ensemble->GetAtomInfo(),
          "fffs.RetrieveNextMolecule( handler)"
        );

        // check that the bond info is identical
        BCL_ExampleIndirectAssert
        (
          fragment.GetBondInfo(),
          itr_ensemble->GetBondInfo(),
          "fffs.RetrieveNextMolecule( handler)"
        );

        ++itr_ensemble;
        ++number_mols;
      }

      // check that the fragment feed read all five molecules
      BCL_ExampleIndirectCheck( number_mols, 5, "fragment feed read all five molecules");

      // test skip
      BCL_ExampleIndirectCheck( fffs.Open( mglur5_filename), true, "reopen the same feed");

      // skip the first molecule
      BCL_ExampleCheck( fffs.Skip(), true);

      // read the second
      fffs.RetrieveNextMolecule( handler);

      // test atom info between the feed and the ensemble
      BCL_ExampleIndirectCheck
      (
        chemistry::FragmentComplete( sdf::FragmentFactory::MakeFragment( handler)).GetAtomInfo(),
        ( ++ensemble.GetMolecules().Begin())->GetAtomInfo(),
        "fffs.Skip()"
      );

      io::File::CloseClearFStream( read);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryFragmentFeedFromStream

  const ExampleClass::EnumType ExampleChemistryFragmentFeedFromStream::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentFeedFromStream())
  );

} // namespace bcl
