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
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_constitution_shared.cpp
  //!
  //! @author kothiwsk
  //! @date Dec 20, 2011
  //! @remarks status complete
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentConstitutionShared :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentConstitutionShared *Clone() const
    {
      return new ExampleChemistryFragmentConstitutionShared( *this);
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

      // create initializer for atom_constitution vector
      storage::Vector< sdf::AtomInfo> atom_types
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().H_S, chemistry::e_UnknownChirality),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_DiDiPiPi, chemistry::e_UnknownChirality),
          sdf::AtomInfo( chemistry::GetAtomTypes().N_Di2DiPiPi, chemistry::e_UnknownChirality)
        )
      );

      // create constitutional bonds
      storage::Vector< sdf::BondInfo> constitution_bonds
      (
        storage::Vector< sdf::BondInfo>::Create
        (
          sdf::BondInfo( 0, 1, chemistry::GetConstitutionalBondTypes().e_NonConjugatedSingleBond),
          sdf::BondInfo( 1, 2, chemistry::GetConstitutionalBondTypes().e_ConjugatedTripleBond)
        )
      );

      // create atom_constitutional atom_vector
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> atom_vector( atom_types, constitution_bonds);

      // default constructor
      chemistry::FragmentConstitutionShared constitution_default;

      // construct constitution from atom_vector
      chemistry::FragmentConstitutionShared constitution_a( atom_vector);

      //  read in molecule and construct constitution from conformation interface
      io::IFStream input_sdf;
      const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into small_mol_conformation
      chemistry::FragmentComplete small_mol_conformation;
      small_mol_conformation
        = sdf::FragmentFactory::MakeFragment( input_sdf);
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

      chemistry::FragmentConstitutionShared constitution_b( small_mol_conformation);

    /////////////////
    // data access //
    /////////////////

      // check GetAtomsIterator for constitution_a
      size_t atom_index( 0);
      for
      (
        iterate::Generic< const chemistry::AtomConstitutionalInterface>
        itr( constitution_a.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        if( itr->GetAtomType() != atom_vector( atom_index).GetAtomType())
        {
          break;
        }
        ++atom_index;
      }

      // atom index should be equal to number of atoms
      BCL_ExampleIndirectCheck
      (
        atom_index,
        atom_vector.GetSize(),
        "constitutional iterator succeeds to iterate over underlying atom vector"
      );

      // check GetBonds for constitution_a
      BCL_ExampleCheck( constitution_a.GetBondInfo(), atom_vector.GetBondInfo());

      // check GetBonds for constitution_b
      BCL_ExampleCheck( constitution_b.GetNumberBonds(), 6);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

//      BCL_ExampleIndirectCheck
//      (
//        TestBCLObjectOutputDiffers( constitution_a, chemistry::FragmentConstitutionShared()),
//        true,
//        "I/O"
//      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryFragmentConstitutionShared

  const ExampleClass::EnumType ExampleChemistryFragmentConstitutionShared::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentConstitutionShared())
  );

} // namespace bcl
