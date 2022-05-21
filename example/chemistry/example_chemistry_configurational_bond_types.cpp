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
#include "chemistry/bcl_chemistry_configurational_bond_types.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_configurational_bond_types.cpp
  //!
  //! @author mendenjl
  //! @date Dec 02, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConfigurationalBondTypes :
    public ExampleInterface
  {
  public:

    ExampleChemistryConfigurationalBondTypes *Clone() const
    {
      return new ExampleChemistryConfigurationalBondTypes( *this);
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

    /////////////////
    // data access //
    /////////////////

      // # of electrons
      BCL_ExampleCheck( chemistry::GetConfigurationalBondTypes().e_AromaticDoubleBond->GetNumberOfElectrons(), 4);
      BCL_ExampleCheck( chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond->GetNumberOfElectrons(), 2);

      // sdf file id of the basic type
      BCL_ExampleCheck( chemistry::GetConfigurationalBondTypes().e_AromaticDoubleBond->GetSDFileID(), 2);
      BCL_ExampleCheck( chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond->GetSDFileID(), 1);

      // alternative sdf file id (for conjugated types)
      BCL_ExampleCheck( chemistry::GetConfigurationalBondTypes().e_AromaticDoubleBond->GetSDAltFileID(), 7);
      BCL_ExampleCheck( chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond->GetSDAltFileID(), 1);

      // conjugation
      BCL_ExampleCheck( chemistry::GetConfigurationalBondTypes().e_AromaticDoubleBond->GetConjugation(), chemistry::ConstitutionalBondTypeData::e_Aromatic);
      BCL_ExampleCheck( chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond->GetConjugation(), chemistry::ConstitutionalBondTypeData::e_Nonconjugated);

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
  }; //end ExampleChemistryConfigurationalBondTypes

  const ExampleClass::EnumType ExampleChemistryConfigurationalBondTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConfigurationalBondTypes())
  );
} // namespace bcl
