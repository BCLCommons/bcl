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
#include "chemistry/bcl_chemistry_atom_complete.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_complete.cpp
  //!
  //! @author kothiwsk
  //! @date Jan 11, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomComplete :
    public ExampleInterface
  {
  public:

    ExampleChemistryAtomComplete *Clone() const
    {
      return new ExampleChemistryAtomComplete( *this);
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

      // construct initializer for AtomComplete
      storage::Vector< sdf::AtomInfo> complete_intializer
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_NonChiral, linal::Vector3D( double( 0.0))),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_NonChiral, linal::Vector3D( double( -1.0))),
          sdf::AtomInfo( chemistry::GetAtomTypes().N_Di2DiPiPi, chemistry::e_NonChiral, linal::Vector3D( double( 1.0)))
        )
      );

      // construct bond connectivities
      storage::Vector< sdf::BondInfo> atom_bonds
      (
        storage::Vector< sdf::BondInfo>::Create
        (
          sdf::BondInfo( 0, 1, chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond),
          sdf::BondInfo( 0, 2, chemistry::GetConfigurationalBondTypes().e_ConjugatedTripleBond)
        )
      );

      // test construction using AtomVector
      chemistry::AtomVector< chemistry::AtomComplete> atom_vector( complete_intializer, atom_bonds);

    /////////////////
    // data access //
    /////////////////

      // test GetChirality
      BCL_ExampleCheck( atom_vector( 1).GetChirality(), chemistry::e_NonChiral);

      // test AtomType
      BCL_ExampleCheck( atom_vector( 1).GetAtomType(), chemistry::GetAtomTypes().C_TeTeTeTe);

      // test GetBonds
      BCL_ExampleCheck( atom_vector( 1).GetBonds().GetSize(), 1);

      // test GetPosition
      BCL_ExampleCheck( atom_vector( 1).GetPosition(), linal::Vector3D( double( -1.0)));

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

  }; //end ExampleChemistryAtomComplete

  const ExampleClass::EnumType ExampleChemistryAtomComplete::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomComplete())
  );

} // namespace bcl
