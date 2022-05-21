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
#include "sdf/bcl_sdf_mdl_property.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_mdl_property.cpp
  //!
  //! @author mendenjl
  //! @date Feb 24, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfMdlProperty :
    public ExampleInterface
  {
  public:

    ExampleSdfMdlProperty *Clone() const
    {
      return new ExampleSdfMdlProperty( *this);
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
      sdf::MdlProperty entry_type;

      // create some mdl atoms for which properties will be generated
      storage::Vector< sdf::AtomInfo> mdl_atoms;
      mdl_atoms.PushBack
      (
        sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_SChirality)
      );
      mdl_atoms.PushBack
      (
        sdf::AtomInfo( chemistry::GetAtomTypes().C_TrTrTrPi, chemistry::e_NonChiral)
      );
      mdl_atoms.PushBack
      (
        sdf::AtomInfo( chemistry::GetAtomTypes().N_Di2DiPi2Pi, chemistry::e_NonChiral)
      );
      mdl_atoms.PushBack
      (
        sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_RChirality)
      );
      mdl_atoms.PushBack
      (
        sdf::AtomInfo( chemistry::GetAtomTypes().C_TrTrPi, chemistry::e_NonChiral)
      );

      // create some bonds
      // for this example, the atom indices are irrelevant; the only information used from teh bond
      // is the double bond isometry
      storage::Vector< sdf::BondInfo> mdl_bonds;
      mdl_bonds.PushBack
      (
        sdf::BondInfo( 0, 1, chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond)
      );
      mdl_bonds.PushBack
      (
        sdf::BondInfo( 0, 3, chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond)
      );
      mdl_bonds.PushBack
      (
        sdf::BondInfo( 1, 2, chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E)
      );
      mdl_bonds.PushBack
      (
        sdf::BondInfo( 1, 2, chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_Z)
      );
      mdl_bonds.PushBack
      (
        sdf::BondInfo( 1, 2, chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_X)
      );
      mdl_bonds.PushBack
      (
        sdf::BondInfo( 1, 2, chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond)
      );

      // create each mdl property
      sdf::MdlProperty mdl_charge( sdf::MdlProperty::e_Charge, mdl_atoms, mdl_bonds);
      sdf::MdlProperty mdl_atom_types( sdf::MdlProperty::e_BclAtomType, mdl_atoms, mdl_bonds);
      sdf::MdlProperty mdl_double_bond_isometries( sdf::MdlProperty::e_BclDoubleBondIsometry, mdl_atoms, mdl_bonds);
      sdf::MdlProperty mdl_chirality( sdf::MdlProperty::e_BclChirality, mdl_atoms, mdl_bonds);

    /////////////////
    // data access //
    /////////////////

      // test get property
      BCL_ExampleCheck( mdl_charge.GetProperty(), sdf::MdlProperty::e_Charge);

      // test get string on each property
      const std::string expected_charge_string( "M  CHG  2   3  -1   5   1");
      const std::string expected_atom_types( "M  BCL ATM C_TeTeTeTe C_TrTrTrPi N_Di2DiPi2Pi C_TeTeTeTe C_TrTrPi");
      const std::string expected_bond_isometries( "M  BCL DBI E Z X -");
      const std::string expected_chiralities( "M  BCL CHI   1   S   4   R");
      BCL_ExampleCheck( mdl_charge.GetString(), expected_charge_string);
      BCL_ExampleCheck( mdl_atom_types.GetString(), expected_atom_types);
      BCL_ExampleCheck( mdl_double_bond_isometries.GetString(), expected_bond_isometries);
      BCL_ExampleCheck( mdl_chirality.GetString(), expected_chiralities);

      // now iterate over the mdl atoms, change all the atom types to basic types, all the chiralities to unknown,
      // all charges to 0
      for
      (
        storage::Vector< sdf::AtomInfo>::iterator itr( mdl_atoms.Begin()), itr_end( mdl_atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->SetAtomType( chemistry::GetAtomTypes().GetAtomType( itr->GetAtomType()->GetElementType(), 0));
        if( itr->GetChirality() != chemistry::e_NonChiral)
        {
          itr->SetChirality( chemistry::e_UnknownChirality);
        }
      }

      // also iterate over bonds and change all isometries to unknown
      for
      (
        storage::Vector< sdf::BondInfo>::iterator itr( mdl_bonds.Begin()), itr_end( mdl_bonds.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->SetIsometry( chemistry::e_UnknownIsometry);
      }

      sdf::MdlProperty mdl_charge_unk( sdf::MdlProperty::e_Charge, mdl_atoms, mdl_bonds);
      sdf::MdlProperty mdl_atom_types_unk( sdf::MdlProperty::e_BclAtomType, mdl_atoms, mdl_bonds);
      sdf::MdlProperty mdl_double_bond_isometries_unk( sdf::MdlProperty::e_BclDoubleBondIsometry, mdl_atoms, mdl_bonds);
      sdf::MdlProperty mdl_chirality_unk( sdf::MdlProperty::e_BclChirality, mdl_atoms, mdl_bonds);

      BCL_ExampleCheck( mdl_charge_unk.GetString(), "");
      BCL_ExampleCheck( mdl_atom_types_unk.GetString(), "M  BCL ATM Carbon_0 Carbon_0 Nitrogen_0 Carbon_0 Carbon_0");
      BCL_ExampleCheck( mdl_double_bond_isometries_unk.GetString(), "M  BCL DBI X X X X");
      BCL_ExampleCheck( mdl_chirality_unk.GetString(), "");

    ////////////////
    // operations //
    ////////////////

      // now use each property to set the associated atom / bond
      mdl_charge.ApplyProperty( mdl_atoms, mdl_bonds);
      BCL_ExampleIndirectCheck
      (
        sdf::MdlProperty( sdf::MdlProperty::e_Charge, mdl_atoms, mdl_bonds).GetString(),
        expected_charge_string,
        "mdl_charge.ApplyProperty( mdl_atoms, mdl_bonds)"
      );

      mdl_atom_types.ApplyProperty( mdl_atoms, mdl_bonds);
      BCL_ExampleIndirectCheck
      (
        sdf::MdlProperty( sdf::MdlProperty::e_BclAtomType, mdl_atoms, mdl_bonds).GetString(),
        expected_atom_types,
        "mdl_atom_types.ApplyProperty( mdl_atoms, mdl_bonds)"
      );

      mdl_chirality.ApplyProperty( mdl_atoms, mdl_bonds);
      BCL_ExampleIndirectCheck
      (
        sdf::MdlProperty( sdf::MdlProperty::e_BclChirality, mdl_atoms, mdl_bonds).GetString(),
        expected_chiralities,
        "mdl_chirality.ApplyProperty( mdl_atoms, mdl_bonds)"
      );

      mdl_double_bond_isometries.ApplyProperty( mdl_atoms, mdl_bonds);
      BCL_ExampleIndirectCheck
      (
        sdf::MdlProperty( sdf::MdlProperty::e_BclDoubleBondIsometry, mdl_atoms, mdl_bonds).GetString(),
        expected_bond_isometries,
        "mdl_double_bond_isometries.ApplyProperty( mdl_atoms, mdl_bonds)"
      );

      // now test constructor from string
      BCL_ExampleCheck( sdf::MdlProperty( expected_atom_types).GetString(), expected_atom_types);
      BCL_ExampleCheck( sdf::MdlProperty( expected_chiralities).GetString(), expected_chiralities);
      BCL_ExampleCheck( sdf::MdlProperty( expected_bond_isometries).GetString(), expected_bond_isometries);
      BCL_ExampleCheck( sdf::MdlProperty( expected_charge_string).GetString(), expected_charge_string);

      // also test that the property type is recognized from the string
      BCL_ExampleCheck( sdf::MdlProperty( expected_atom_types).GetProperty(), sdf::MdlProperty::e_BclAtomType);
      BCL_ExampleCheck( sdf::MdlProperty( expected_chiralities).GetProperty(), sdf::MdlProperty::e_BclChirality);
      BCL_ExampleCheck( sdf::MdlProperty( expected_bond_isometries).GetProperty(), sdf::MdlProperty::e_BclDoubleBondIsometry);
      BCL_ExampleCheck( sdf::MdlProperty( expected_charge_string).GetProperty(), sdf::MdlProperty::e_Charge);
      BCL_ExampleCheck( sdf::MdlProperty( "M  END").GetProperty(), sdf::MdlProperty::e_BlockTerminator);

    //////////////////////
    // helper functions //
    //////////////////////

      // test io
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( mdl_charge, sdf::MdlProperty()), true);
      BCL_ExampleCheck( TestBCLObjectOutputDiffers( mdl_charge, sdf::MdlProperty()), true);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSdfMdlProperty

  const ExampleClass::EnumType ExampleSdfMdlProperty::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfMdlProperty())
  );

} // namespace bcl
