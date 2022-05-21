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
#include "descriptor/bcl_descriptor_atom_hbond_info.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_hbond_info.cpp
  //!
  //! @author kothiwsk, mendenjl, geanesar
  //! @date Dec 14, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomHBondInfo :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomHBondInfo *Clone() const
    {
      return new ExampleDescriptorAtomHBondInfo( *this);
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
      descriptor::AtomHBondInfo acceptors( descriptor::AtomHBondInfo::e_Acceptor);

      // default constructor
      descriptor::AtomHBondInfo donors( descriptor::AtomHBondInfo::e_Donor);

      // copy constructor
      descriptor::AtomHBondInfo acceptors_copy( acceptors);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( acceptors.GetAlias(), "Atom_HbondAcceptors");

      BCL_ExampleCheck( donors.GetAlias(), "Atom_HbondDonors");

      BCL_ExampleCheck( acceptors.GetString(), "Atom_HbondAcceptors");

    ///////////////
    // operators //
    ///////////////

      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input);
      // close stream
      io::File::CloseClearFStream( input);

      // make a miscellaneous property retriever for hbond_acceptors
      std::string hbond_acceptors_from_adriana_name( "HAcc");

      // make a miscellaneous property retriever for hbond_acceptors
      std::string hbond_donors_from_adriana_name( "HDon");

      // check that the number of hbond acceptors from adriana is the same as the number we get from the atom property
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        // check that the misc property gave back just a single value
        BCL_ExampleAssert( itr->IsPropertyStored( hbond_acceptors_from_adriana_name), true);
        BCL_ExampleAssert( itr->IsPropertyStored( hbond_donors_from_adriana_name), true);

        // check that the adriana-calculated property is the same as the bcl-calculated property
        BCL_ExampleAssert
        (
          acceptors.CollectValuesOnEachElementOfObject( *itr).Sum(),
          itr->GetMDLPropertyAsVector( hbond_acceptors_from_adriana_name)( 0)
        );

        BCL_ExampleAssert
        (
          donors.CollectValuesOnEachElementOfObject( *itr).Sum(),
          itr->GetMDLPropertyAsVector( hbond_donors_from_adriana_name)( 0)
        );
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomHBondInfo

  const ExampleClass::EnumType ExampleDescriptorAtomHBondInfo::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomHBondInfo())
  );

} // namespace bcl
