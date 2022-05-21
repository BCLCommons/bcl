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
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_configuration_shared.cpp
  //! @details Tests ChemistryFragmentConfigurationShared class which contains small molecule configuration data
  //!
  //! @author kothiwsk
  //! @date
  //! @remarks status complete
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentConfigurationShared :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentConfigurationShared *Clone() const
    {
      return new ExampleChemistryFragmentConfigurationShared( *this);
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

      // create AtomVector of constitution
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> constitution_vector
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().H_S, chemistry::e_UnknownChirality),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_DiDiPiPi, chemistry::e_UnknownChirality),
          sdf::AtomInfo( chemistry::GetAtomTypes().N_Di2DiPiPi, chemistry::e_UnknownChirality)
        ),
        storage::Vector< sdf::BondInfo>::Create
        (
          sdf::BondInfo( 0, 1, chemistry::GetConstitutionalBondTypes().e_NonConjugatedSingleBond),
          sdf::BondInfo( 1, 2, chemistry::GetConstitutionalBondTypes().e_ConjugatedTripleBond)
        )
      );

      // create atom_configuration vector using intializer and bonds
      chemistry::AtomVector< chemistry::AtomConfigurationalShared> atom_vector
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().H_S, chemistry::e_NonChiral),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_DiDiPiPi, chemistry::e_NonChiral),
          sdf::AtomInfo( chemistry::GetAtomTypes().N_Di2DiPiPi, chemistry::e_NonChiral)
        ),
        storage::Vector< sdf::BondInfo>::Create
        (
          sdf::BondInfo( 0, 1, chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond),
          sdf::BondInfo( 1, 2, chemistry::GetConfigurationalBondTypes().e_ConjugatedTripleBond)
        )
      );
      atom_vector.LinkToLayer
      (
        iterate::Generic< const chemistry::AtomConstitutionalInterface>
        (
          constitution_vector.Begin(),
          constitution_vector.End()
        )
      );

      // create FragmentConstitutionShared using atom_constitution atom vector
      util::ShPtr< chemistry::FragmentConstitutionShared> fragment_constitution
      (
        new chemistry::FragmentConstitutionShared( constitution_vector)
      );

      // default constructor
      chemistry::FragmentConfigurationShared configuration_default;

      // construct FragmentConfigurationShared
      chemistry::FragmentConfigurationShared configuration_a( fragment_constitution, atom_vector);

      // read in a molecule and construct from conformation interface
      io::IFStream input_sdf;
      const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into FragmentComplete
      chemistry::FragmentComplete small_mol_conformation;
      small_mol_conformation
        = sdf::FragmentFactory::MakeFragment( input_sdf);
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

      // construct configuration from conformation
      chemistry::FragmentConfigurationShared configuration_b( small_mol_conformation);

    /////////////////
    // data access //
    /////////////////

      // check GetAtomsIterator for configuration_a
      size_t atom_index( 0);
      for
      (
        iterate::Generic< const chemistry::AtomConfigurationalInterface>
        itr( configuration_a.GetAtomsIterator());
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
        "Configurational iterator succeeds to iterate over underlying atom vector"
      );

      // check GetBonds for configuration_a
      BCL_ExampleCheck( configuration_a.GetNumberBonds(), atom_vector.GetNumberBonds());

      // check GetAtomsIterator for configuration_b
      atom_index = 0;
      for
      (
        iterate::Generic< const chemistry::AtomConfigurationalInterface>
        itr( configuration_b.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        if( itr->GetAtomType() != chemistry::GetAtomTypes().C_TeTeTeTe)
        {
          break;
        }
        ++atom_index;
      }

      // atom_index should be equal to number of atoms
      BCL_ExampleIndirectCheck
      (
        atom_index,
        6,
        "Configurational iterator succeeds to iterate over underlying atom vector"
      );

      // check GetBonds for configuration_b
      BCL_ExampleCheck( configuration_b.GetNumberBonds(), 6);

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
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistryFragmentConfigurationShared::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentConfigurationShared())
  );

} // namespace bcl
