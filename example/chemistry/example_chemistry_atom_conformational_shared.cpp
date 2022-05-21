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
#include "chemistry/bcl_chemistry_atom_conformational_shared.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_conformational_shared.cpp
  //!
  //! @author kothiwsk
  //! @date Dec 05, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomConformationalShared :
    public ExampleInterface
  {
  public:

    ExampleChemistryAtomConformationalShared *Clone() const
    {
      return new ExampleChemistryAtomConformationalShared( *this);
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

      // construct atom_conformation using AtomVector
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> constitution_vector
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_UnknownChirality),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_UnknownChirality)
        ),
        storage::Vector< sdf::BondInfo>
        (
          1,
          sdf::BondInfo( 0, 1, chemistry::GetConstitutionalBondTypes().e_ConjugatedDoubleBond)
        )
      );

      util::ShPtr< chemistry::FragmentConstitutionShared> constitution
      (
        new chemistry::FragmentConstitutionShared( constitution_vector)
      );

      // construct atom_conformation using AtomVector
      chemistry::AtomVector< chemistry::AtomConfigurationalShared> configuration_vector
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_SChirality),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_SChirality)
        ),
        storage::Vector< sdf::BondInfo>
        (
          1,
          sdf::BondInfo( 0, 1, chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E)
        )
      );
      configuration_vector.LinkToLayer( constitution->GetAtomsIterator());

      util::ShPtr< chemistry::FragmentConfigurationShared> configuration
      (
          new chemistry::FragmentConfigurationShared( constitution, configuration_vector)
      );

      // construct atom_conformation using AtomVector
      chemistry::AtomVector< chemistry::AtomConformationalShared> atom_vector
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_SChirality, linal::Vector3D( double( 1.0))),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_SChirality, linal::Vector3D( double( 0.0)))
        ),
        storage::Vector< sdf::BondInfo>
        (
          1,
          sdf::BondInfo( 0, 1, chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E)
        )
      );

      atom_vector.LinkToLayer( configuration->GetAtomsIterator());

    /////////////////
    // data access //
    /////////////////

      // test GetChirality
      BCL_ExampleCheck( atom_vector.First().GetChirality(), chemistry::e_SChirality);

      // test GetAtomType
      BCL_ExampleCheck( atom_vector.First().GetAtomType(), chemistry::GetAtomTypes().C_TeTeTeTe);

      // test GetElementType
      BCL_ExampleCheck( atom_vector.First().GetElementType(), chemistry::GetElementTypes().e_Carbon);

      // test GetCharge
      BCL_ExampleCheck( atom_vector.First().GetCharge(), 0);

      // test GetBonds
      BCL_ExampleCheck( atom_vector.First().GetBonds().GetSize(), 1);

      // test GetBonds
      BCL_ExampleCheck
      (
        atom_vector.First().GetBonds().FirstElement().GetBondType(),
        chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E
      );

      // test GetBonds
      BCL_ExampleCheck
      (
        atom_vector( 1).GetBonds().FirstElement().GetBondType(),
        chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E
      );

      // test GetBonds
      BCL_ExampleCheck( atom_vector( 1).GetBonds().GetSize(), 1);

      // test GetNumberValenceBonds
      BCL_ExampleCheck( atom_vector.First().GetNumberValenceBonds(), 3);

      // test GetNumberElectronsInValenceBonds
      BCL_ExampleCheck( atom_vector.First().GetNumberElectronsInValenceBonds(), 2);

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

  }; //end ExampleChemistryAtomConformationalShared

  const ExampleClass::EnumType ExampleChemistryAtomConformationalShared::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomConformationalShared())
  );

} // namespace bcl
