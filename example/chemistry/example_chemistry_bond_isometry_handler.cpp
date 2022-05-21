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
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_molecule_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_bond_isometry_handler.cpp
  //!
  //! @author mendenjl
  //! @date May 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryBondIsometryHandler :
    public ExampleInterface
  {
  public:

    ExampleChemistryBondIsometryHandler *Clone() const
    {
      return new ExampleChemistryBondIsometryHandler( *this);
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
    ////////////////
    // operations //
    ////////////////

      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_E.sdf"));

      // load in the conformation of 1,3 pentadiene, the trans conformation
      chemistry::MoleculeComplete pentadiene_E_conformation( sdf::Factory::MakeMolecule( input));
      io::File::CloseClearFStream( input);

      // there should be one isometric bond
      BCL_ExampleIndirectCheck
      (
        pentadiene_E_conformation.CountNonValenceBondsWithProperty
        (
          chemistry::ConfigurationalBondTypeData::e_IsIsometric,
          size_t( 1)
        ),
        1,
        "loading conformation automatically determines E/Z isometry"
      );

      // there should be one bond with E isometry
      BCL_ExampleIndirectCheck
      (
        pentadiene_E_conformation.CountNonValenceBondsWithProperty
        (
          chemistry::ConfigurationalBondTypeData::e_Isometry,
          size_t( chemistry::e_EIsometry)
        ),
        1,
        "loading conformation automatically determines E/Z isometry"
      );

      // there should not be any bonds with Z isometry
      BCL_ExampleIndirectCheck
      (
        pentadiene_E_conformation.CountNonValenceBondsWithProperty
        (
          chemistry::ConfigurationalBondTypeData::e_Isometry,
          size_t( chemistry::e_ZIsometry)
        ),
        0,
        "loading conformation automatically determines E/Z isometry"
      );

      // all other bonds should be non-isometric
      BCL_ExampleIndirectCheck
      (
        pentadiene_E_conformation.CountNonValenceBondsWithProperty
        (
          chemistry::ConfigurationalBondTypeData::e_Isometry,
          size_t( chemistry::e_NonIsometric)
        ),
        pentadiene_E_conformation.GetNumberBonds() - 1,
        "loading conformation automatically determines E/Z isometry"
      );

      // now check the z stereoisomer
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_Z.sdf"));

      // load in the conformation of 1,3 pentadiene, the trans conformation
      chemistry::MoleculeComplete pentadiene_Z_conformation( sdf::Factory::MakeMolecule( input));
      io::File::CloseClearFStream( input);

      // there should be one isometric bond
      BCL_ExampleIndirectCheck
      (
        pentadiene_Z_conformation.CountNonValenceBondsWithProperty
        (
          chemistry::ConfigurationalBondTypeData::e_IsIsometric,
          size_t( 1)
        ),
        1,
        "loading conformation automatically determines E/Z isometry"
      );

      // there should be no bonds with E isometry
      BCL_ExampleIndirectCheck
      (
        pentadiene_Z_conformation.CountNonValenceBondsWithProperty
        (
          chemistry::ConfigurationalBondTypeData::e_Isometry,
          size_t( chemistry::e_EIsometry)
        ),
        0,
        "loading conformation automatically determines E/Z isometry"
      );

      // there should 1 bond with Z isometry
      BCL_ExampleIndirectCheck
      (
        pentadiene_Z_conformation.CountNonValenceBondsWithProperty
        (
          chemistry::ConfigurationalBondTypeData::e_Isometry,
          size_t( chemistry::e_ZIsometry)
        ),
        1,
        "loading conformation automatically determines Z isometry"
      );

      // all other bonds should be non-isometric
      BCL_ExampleIndirectCheck
      (
        pentadiene_Z_conformation.CountNonValenceBondsWithProperty
        (
          chemistry::ConfigurationalBondTypeData::e_Isometry,
          size_t( chemistry::e_NonIsometric)
        ),
        pentadiene_Z_conformation.GetNumberBonds() - 1,
        "loading conformation automatically determines E/Z isometry"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryBondIsometryHandler

  const ExampleClass::EnumType ExampleChemistryBondIsometryHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryBondIsometryHandler())
  );

} // namespace bcl
