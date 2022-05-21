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
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_constitutional_bond_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_configurational_bond_type_data.cpp
  //!
  //! @author mendenjl
  //! @date Dec 02, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConfigurationalBondTypeData :
    public ExampleInterface
  {
  public:

    ExampleChemistryConfigurationalBondTypeData *Clone() const
    {
      return new ExampleChemistryConfigurationalBondTypeData( *this);
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

      const chemistry::ConfigurationalBondTypeData aromatic_double_bond( chemistry::GetConstitutionalBondTypes().e_AromaticDoubleBond, chemistry::e_EIsometry);
      const chemistry::ConfigurationalBondTypeData basic_single_bond( chemistry::GetConstitutionalBondTypes().e_NonConjugatedSingleBond, chemistry::e_NonIsometric);

    /////////////////
    // data access //
    /////////////////

      // # of electrons
      BCL_ExampleCheck( aromatic_double_bond.GetNumberOfElectrons(), 4);
      BCL_ExampleCheck( basic_single_bond.GetNumberOfElectrons(), 2);

      // sdf file id of the basic type
      BCL_ExampleCheck( aromatic_double_bond.GetSDFileID(), 2);
      BCL_ExampleCheck( basic_single_bond.GetSDFileID(), 1);

      // alternative sdf file id (for conjugated types)
      BCL_ExampleCheck( aromatic_double_bond.GetSDAltFileID(), 7);
      BCL_ExampleCheck( basic_single_bond.GetSDAltFileID(), 1);

      // conjugation
      BCL_ExampleCheck( aromatic_double_bond.GetConjugation(), chemistry::ConstitutionalBondTypeData::e_Aromatic);
      BCL_ExampleCheck( basic_single_bond.GetConjugation(), chemistry::ConstitutionalBondTypeData::e_Nonconjugated);

      // isometry
      BCL_ExampleCheck( aromatic_double_bond.GetIsometry(), chemistry::e_EIsometry);
      BCL_ExampleCheck( basic_single_bond.GetIsometry(), chemistry::e_NonIsometric);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // Checks Read and Write
      BCL_ExampleIndirectCheck
      (
        ExampleInterface::TestBCLObjectIOForSymmetry( aromatic_double_bond, chemistry::ConfigurationalBondTypeData()),
        true,
        "Chemistry::ConfigurationalBondTypeData I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleChemistryConfigurationalBondTypeData

  const ExampleClass::EnumType ExampleChemistryConfigurationalBondTypeData::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConfigurationalBondTypeData())
  );
} // namespace bcl
