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
#include "chemistry/bcl_chemistry_molecule_ensemble.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_molecule_ensemble.cpp
  //!
  //! @author mendenjl
  //! @date Mar 08, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMoleculeEnsemble :
    public ExampleInterface
  {
  public:

    ExampleChemistryMoleculeEnsemble *Clone() const
    {
      return new ExampleChemistryMoleculeEnsemble( *this);
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
      chemistry::MoleculeEnsemble small_mols_ensemble;

      // copy constructor
      chemistry::MoleculeEnsemble ensemble_copy( small_mols_ensemble);

      // clone
//      chemistry::MoleculeEnsemble ensemble_clone( ensemble_copy.Clone());

    /////////////////
    // data access //
    /////////////////

      // check static class name
      BCL_ExampleCheck
      (
        GetStaticClassName< chemistry::MoleculeEnsemble>(),
        small_mols_ensemble.GetClassIdentifier()
      );

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));
      chemistry::MoleculeEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));

      // load a partial ensemble of 3 molecules
      chemistry::MoleculeEnsemble partial_ensemble( input, math::Range< size_t>( 2, 4));

      // check whether GetMolecules retrieves the right number of molecules
      BCL_ExampleCheck( ensemble.GetMolecules().GetSize(), 5);

      // check whether partial ensemble has right number of molecules
      BCL_ExampleCheck( partial_ensemble.GetMolecules().GetSize(), 3);

      // close the input stream
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));

      // pushback to ensemble
      ensemble.PushBack( sdf::Factory::MakeMolecule( sdf::MdlHandler( input)));

      // close the input stream
      io::File::CloseClearFStream( input);

      // check whether GetMolecules retrieves the right number of molecules
      BCL_ExampleIndirectCheck
      (
        ensemble.GetMolecules().GetSize(),
        6,
        "chemistry::MoleculeEnsemble::PushBack( ShPtr< SmallMolecule>)"
      );

// TODO
//      // create a property to count the number of hydrogens in the molecule to ensure that it is zero
//      chemistry::HasProperties hydrogen_count
//      (
//        util::ObjectDataLabel( "Sum(Atom_Mask(property=AtomicNumbers,mask=[1]))")
//      );
//
//      // check that only hydrogen was removed from each molecule
//      for
//      (
//        storage::List< chemistry::FragmentComplete>::const_iterator itr( partial_ensemble.Begin()),
//        itr_end( partial_ensemble.End());
//        itr != itr_end;
//        ++itr
//      )
//      {
//        BCL_ExampleCheck( hydrogen_count( *itr).FirstElement(), 0);
//        BCL_ExampleCheck( ( *itr).GetNumberAtoms() > 0, true);
//      }

      chemistry::MoleculeEnsemble combined_ensemble( partial_ensemble);

      // check that append function works
      combined_ensemble.Append( ensemble);
      BCL_ExampleIndirectCheck( combined_ensemble.GetSize(), 9, "Append - Append ensemble to partial ensemble");

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read write function
//      WriteBCLObject( ensemble);
//      chemistry::MoleculeEnsemble ensemble_read_in;
//      ReadBCLObject( ensemble_read_in);
//
//      // check whether read in ensemble has the right number of molecules
//      BCL_ExampleIndirectCheck( ensemble_read_in.GetMolecules().GetSize(), 6, "I/O");

//      {
//        BCL_ExampleMustOpenInputFile( input, "/home/mendenjl/fragment_search/csd_whole.sdf");
//        util::MemoryUsage::WriteCurrentMemoryUsageInfo( util::GetLogger());
//        util::Stopwatch timer( "time to load csd for old small molecule");
//        chemistry::MoleculeEnsemble ensemble( input);
//        io::File::CloseClearFStream( input);
//        BCL_Message
//        (
//          util::Message::e_Standard,
//          "finished reading ensemble with " + util::Format()( ensemble.GetSize()) + " molecules."
//        );
//        util::MemoryUsage::WriteCurrentMemoryUsageInfo( util::GetLogger());
//      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryMoleculeEnsemble

  const ExampleClass::EnumType ExampleChemistryMoleculeEnsemble::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMoleculeEnsemble())
  );

} // namespace bcl
