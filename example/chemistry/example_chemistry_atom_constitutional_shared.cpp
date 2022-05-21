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
#include "chemistry/bcl_chemistry_atom_constitutional_shared.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_constitutional_shared.cpp
  //!
  //! @author kothiwsk
  //! @date Dec 05, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomConstitutionalShared :
    public ExampleInterface
  {
  public:

    ExampleChemistryAtomConstitutionalShared *Clone() const
    {
      return new ExampleChemistryAtomConstitutionalShared( *this);
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

      // test construction using AtomVector
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> atom_vector
      (
        storage::Vector< sdf::AtomInfo>( 2, sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_UnknownChirality)),
        storage::Vector< sdf::BondInfo>
        (
          1,
          sdf::BondInfo( 0, 1, chemistry::GetConstitutionalBondTypes().e_NonConjugatedSingleBond)
        )
      );

    /////////////////
    // data access //
    /////////////////

      // test GetAtomType
      BCL_ExampleCheck( atom_vector.First().GetAtomType(), chemistry::GetAtomTypes().C_TeTeTeTe);

      // test GetElementType
      BCL_ExampleCheck( atom_vector.First().GetElementType(), chemistry::GetElementTypes().e_Carbon);

      // test GetCharge
      BCL_ExampleCheck( atom_vector.First().GetCharge(), 0);

      // test GetBonds
      BCL_ExampleCheck( atom_vector.First().GetBonds().GetSize(), 1);

      // test GetBonds
      BCL_ExampleCheck( atom_vector( 1).GetBonds().GetSize(), 1);

      // test GetBonds
      BCL_ExampleCheck
      (
        atom_vector.First().GetBonds().FirstElement().GetBondType(),
        chemistry::GetConstitutionalBondTypes().e_NonConjugatedSingleBond
      );

      // test GetBonds
      BCL_ExampleCheck
      (
        atom_vector( 1).GetBonds().FirstElement().GetBondType(),
        chemistry::GetConstitutionalBondTypes().e_NonConjugatedSingleBond
      );

      // test GetNumberValenceBonds
      BCL_ExampleCheck( atom_vector.First().GetNumberValenceBonds(), 3);

      // test GetNumberElectronsInValenceBonds
      BCL_ExampleCheck( atom_vector.First().GetNumberElectronsInValenceBonds(), 3);

      // test GetValenceBonds
      BCL_ExampleCheck( atom_vector.First().GetValenceBonds().FirstElement(), chemistry::GetConstitutionalBondTypes().e_NonConjugatedSingleBond);

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

  }; //end ExampleChemistryAtomConstitutionalShared

  const ExampleClass::EnumType ExampleChemistryAtomConstitutionalShared::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomConstitutionalShared())
  );

} // namespace bcl
