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
#include "chemistry/bcl_chemistry_has_properties.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_has_properties.cpp
  //! @details Tests ChemistryHasProperties class which contains small molecule configuration data
  //!
  //! @author kothiwsk
  //! @date
  //! @remarks status empty
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryHasProperties :
    public ExampleInterface
  {
  public:

    ExampleChemistryHasProperties *Clone() const
    {
      return new ExampleChemistryHasProperties( *this);
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

//      // test default constructor
//      chemistry::HasProperties small_mol_configuration_a;
//
//      // test constructor given given PairFragments list, ComplexFragments list, and Constitution
//
//      // PairFragments list
//      util::ShPtr< chemistry::AtomWithConfiguration> first_atom
//      (
//        new chemistry::AtomWithConfiguration
//        (
//          util::ShPtr< chemistry::Atom>( new chemistry::Atom( chemistry::GetAtomTypes().C_TrTrTrPi)),
//          chemistry::e_R
//        )
//      );
//      util::ShPtr< chemistry::AtomWithConfiguration> second_atom
//      (
//        new chemistry::AtomWithConfiguration
//        (
//          util::ShPtr< chemistry::Atom>( new chemistry::Atom( chemistry::GetAtomTypes().C_TrTrTrPi)),
//          chemistry::e_S
//        )
//      );
//
//      util::ShPtrVector< chemistry::AtomWithConfiguration> atoms_with_configuration;
//
//      atoms_with_configuration.PushBack( first_atom);
//      atoms_with_configuration.PushBack( second_atom);
//
//      // Constitution
//      chemistry::SmallMoleculeConstitution constitution;
//      util::ShPtr< chemistry::SmallMoleculeConstitution> sp_constitution
//      (
//        new chemistry::SmallMoleculeConstitution( constitution)
//      );
//
//      // now call the constructor
//      chemistry::HasProperties small_mol_configuration_b
//      (
//        atoms_with_configuration,
//        sp_constitution
//      );
//
//      // test copy constructor
//      chemistry::HasProperties small_mol_configuration_c
//      (
//        small_mol_configuration_b
//      );
//
//      // test Clone
//      util::ShPtr< chemistry::HasProperties> sp_small_mol_configuration_d
//      (
//        small_mol_configuration_c.Clone()
//      );
//
//      //Test GetAtomsWithConfiguration function
//      BCL_ExampleCheck( sp_small_mol_configuration_d.IsDefined(), true);
//
//    /////////////////
//    // data access //
//    /////////////////
//
//      //Test GetAtomsWithConfiguration function
//      BCL_ExampleCheck( small_mol_configuration_b.GetAtomsWithConfiguration(), atoms_with_configuration);
//
//      //Test GetConstitution function
//      BCL_ExampleCheck( small_mol_configuration_b.GetConstitution(), sp_constitution);
//
//    ///////////////
//    // operators //
//    ///////////////
//
//    ////////////////
//    // operations //
//    ////////////////
//
//    //////////////////////
//    // input and output //
//    //////////////////////
//
//      // test the stream operators
//      BCL_Message
//      (
//        util::Message::e_Standard,
//        "Outputting small_mol_configuration_c: " + util::Format()( small_mol_configuration_c)
//      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistryHasProperties::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryHasProperties())
  );

} // namespace bcl
