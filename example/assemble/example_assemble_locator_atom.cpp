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
#include "assemble/bcl_assemble_locator_atom.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_locator_atom.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleLocatorAtom :
    public ExampleInterface
  {
  public:

    ExampleAssembleLocatorAtom *Clone() const
    {
      return new ExampleAssembleLocatorAtom( *this);
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

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::LocatorAtom def_const;

      // test construct from chain ID and AA identifiers
      assemble::LocatorAtom atom_known( 'B', 52, biol::GetAtomTypes().CB);
      BCL_MessageStd( "test construct from chain id and amino acid and atom identifiers");
      BCL_Example_Check
      (
        atom_known.GetAtomType() == biol::GetAtomTypes().CB &&
        atom_known.GetLocatorAA().GetLocatorChain().GetChainID() == 'B' && atom_known.GetLocatorAA().GetAAID() == 52,
        "AtomID should be " + util::Format()( biol::GetAtomTypes().CB) + " but is " +
        util::Format()( atom_known.GetAtomType()) + " and ChainID should be \"B\" but is " +
        util::Format()( atom_known.GetLocatorAA().GetLocatorChain().GetChainID()) +
        " and residue id should be 52 but is " + util::Format()( atom_known.GetLocatorAA().GetAAID())
      );

      // test copy constructor
      assemble::LocatorAtom copy_constr( atom_known);
      BCL_MessageStd( "test copy constructor");
      BCL_Example_Check
      (
        copy_constr.GetLocatorAA().GetAAID() == atom_known.GetLocatorAA().GetAAID() &&
        copy_constr.GetLocatorAA().GetLocatorChain().GetChainID() ==
          atom_known.GetLocatorAA().GetLocatorChain().GetChainID() &&
        copy_constr.GetAtomType() == atom_known.GetAtomType(),
        "AtomID should be " + util::Format()( atom_known.GetAtomType()) +
        " but is " + util::Format()( copy_constr.GetAtomType()) +
        " and ChainID should be " + util::Format()( atom_known.GetLocatorAA().GetLocatorChain().GetChainID()) +
        " but is " + util::Format()( copy_constr.GetLocatorAA().GetLocatorChain().GetChainID()) +
        " and residue id should be " + util::Format()( atom_known.GetLocatorAA().GetAAID()) +
        " but is " + util::Format()( copy_constr.GetLocatorAA().GetAAID())
      );

      // test clone constructor
      util::ShPtr< assemble::LocatorAtom> clone_constr( copy_constr.Clone());
      BCL_Example_Check
      (
        clone_constr->GetLocatorAA().GetAAID() == copy_constr.GetLocatorAA().GetAAID() &&
        clone_constr->GetLocatorAA().GetLocatorChain().GetChainID() ==
          copy_constr.GetLocatorAA().GetLocatorChain().GetChainID() &&
        clone_constr->GetAtomType() == copy_constr.GetAtomType(),
        "AtomID should be " + util::Format()( atom_known.GetAtomType()) +
        " but is " + util::Format()( copy_constr.GetAtomType()) +
        " and ChainID should be " + util::Format()( atom_known.GetLocatorAA().GetLocatorChain().GetChainID()) +
        " but is " + util::Format()( copy_constr.GetLocatorAA().GetLocatorChain().GetChainID()) +
        " and residue id should be " + util::Format()( atom_known.GetLocatorAA().GetAAID()) +
        " but is " + util::Format()( copy_constr.GetLocatorAA().GetAAID())
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::LocatorAtom>(), clone_constr->GetClassIdentifier());

      // check GetIdentification
      BCL_ExampleCheck( "  B  52  \"UsePDBID\"  0  \"CB\"", clone_constr->GetIdentification());

      // check GetPymolName
      BCL_ExampleCheck( "B_52_CB", clone_constr->GetPymolName());

      // test SetAtomID and GetAtomID functions
      BCL_MessageStd( "test SetAtomID and GetAtomID functions");
      def_const.SetAtomType( biol::GetAtomTypes().CB);
      BCL_ExampleIndirectCheck( def_const.GetAtomType(), biol::GetAtomTypes().CB, "SetAtomType( GetAtomTypes().CB)");

      // test SetLocatorAA and GetLocatorAA functions
      def_const.SetLocatorAA( assemble::LocatorAA( 'B', 52));
      BCL_Example_Check
      (
        def_const.GetLocatorAA().GetLocatorChain().GetChainID() == 'B' && def_const.GetLocatorAA().GetAAID() == 52,
        "ChainID should be \"B\" but is " + util::Format()( def_const.GetLocatorAA().GetLocatorChain().GetChainID()) +
        " and residue id should be 52 but is " + util::Format()( def_const.GetLocatorAA().GetAAID())
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test Locate atom function
      BCL_MessageStd( "test GetAtom function");
      BCL_ExampleCheck( copy_constr.LocateAtom( protein_model)->GetType(), biol::GetAtomTypes().CB);

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( copy_constr);
      assemble::LocatorAtom read_locator_atom;
      ReadBCLObject( read_locator_atom);

      // make sure that the written locator and the read-in locator are the same
      BCL_ExampleIndirectAssert
      (
        copy_constr.GetLocatorAA().GetString(),
        read_locator_atom.GetLocatorAA().GetString(),
        "I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // operator <
      BCL_ExampleCheck( *clone_constr < assemble::LocatorAtom( 'B', 1, biol::GetAtomTypes().CA), false);
      BCL_ExampleCheck( *clone_constr < assemble::LocatorAtom( 'B', 100, biol::GetAtomTypes().CA), true);
      BCL_ExampleCheck( *clone_constr < assemble::LocatorAtom( 'A', 100, biol::GetAtomTypes().CA), false);
      BCL_ExampleCheck( *clone_constr < assemble::LocatorAtom( 'C', 1, biol::GetAtomTypes().CA), true);
      BCL_ExampleCheck( *clone_constr < assemble::LocatorAtom( 'A', 100, biol::GetAtomTypes().CA), false);
      BCL_ExampleCheck( *clone_constr < assemble::LocatorAtom( 'B', 52, biol::GetAtomTypes().CA), false);
      BCL_ExampleCheck( *clone_constr < assemble::LocatorAtom( 'B', 52, biol::GetAtomTypes().CG), true);

      // operator !=
      BCL_ExampleCheck( *clone_constr != assemble::LocatorAtom( 'B', 52, biol::GetAtomTypes().CB), false);
      BCL_ExampleCheck( *clone_constr != assemble::LocatorAtom( 'B', 51, biol::GetAtomTypes().CB), true);
      BCL_ExampleCheck( *clone_constr != assemble::LocatorAtom( 'A', 52, biol::GetAtomTypes().CB), true);
      BCL_ExampleCheck( *clone_constr != assemble::LocatorAtom( 'B', 52, biol::GetAtomTypes().CG), true);

      // operator ==
      BCL_ExampleCheck( *clone_constr == assemble::LocatorAtom( 'B', 52, biol::GetAtomTypes().CB), true);
      BCL_ExampleCheck( *clone_constr == assemble::LocatorAtom( 'B', 51, biol::GetAtomTypes().CB), false);

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleLocatorAtom

  const ExampleClass::EnumType ExampleAssembleLocatorAtom::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleLocatorAtom())
  );

} // namespace bcl
