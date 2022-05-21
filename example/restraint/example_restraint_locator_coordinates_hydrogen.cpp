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
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_locator_coordinates_hydrogen.cpp
  //!
  //! @author weinerbe
  //! @date May 31, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintLocatorCoordinatesHydrogen :
    public ExampleInterface
  {
  public:

    ExampleRestraintLocatorCoordinatesHydrogen *Clone() const
    {
      return new ExampleRestraintLocatorCoordinatesHydrogen( *this);
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
      // get a protein model
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      restraint::LocatorCoordinatesHydrogen def_construct;
      BCL_ExampleIndirectCheck
      (
        def_construct.GetChainID() == 'A' &&
          def_construct.GetSeqID() == util::GetUndefined< int>() &&
          def_construct.GetAtomType() == biol::GetAtomTypes().e_Undefined,
        true,
        "default constructor"
      );

      // test constructor from data members
      const char chain_id( 'A');
      const int seq_id( 26);
      const biol::AtomType &atom_type( biol::GetAtomTypes().HG13);
      const restraint::LocatorCoordinatesHydrogen hg_26( chain_id, seq_id, atom_type);

    /////////////////
    // data access //
    /////////////////

      // test GetChainID
      BCL_ExampleCheck( hg_26.GetChainID(), chain_id);

      // test GetSeqID
      BCL_ExampleCheck( hg_26.GetSeqID(), seq_id);

      // test GetAtomType
      BCL_ExampleCheck( hg_26.GetAtomType(), atom_type);

      // test GetAAType
      BCL_ExampleCheck( hg_26.GetAAType(), biol::GetAtomTypes().e_Undefined);

    ///////////////
    // operators //
    ///////////////

      // test Locate
      const linal::Vector3D cb_coords( 32.115, 25.287, 13.367);
      BCL_ExampleCheckWithinTolerance( hg_26.Locate( protein_model), cb_coords, 0.001);

      // test LocateAtomCopy
      const biol::Atom located_atom( hg_26.LocateAtomCopy( protein_model));
      BCL_ExampleCheck( hg_26.LocateAtomCopy( protein_model).GetType(), atom_type);
      BCL_ExampleCheckWithinTolerance
      (
        hg_26.LocateAtomCopy( protein_model).GetCoordinates(),
        cb_coords,
        0.001
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( hg_26);
      ReadBCLObject( def_construct);
      BCL_ExampleIndirectCheck
      (
        def_construct.GetChainID() == chain_id &&
          def_construct.GetSeqID() == seq_id &&
          def_construct.GetAtomType() == atom_type,
        true,
        "read and write"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintLocatorCoordinatesHydrogen

  const ExampleClass::EnumType ExampleRestraintLocatorCoordinatesHydrogen::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintLocatorCoordinatesHydrogen())
  );

} // namespace bcl
