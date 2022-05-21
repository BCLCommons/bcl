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
#include "chemistry/bcl_chemistry_possible_atom_types_for_atom.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_possible_atom_types_for_atom.cpp
  //!
  //! @author mendenjl
  //! @date August 01, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryPossibleAtomTypesForAtom :
    public ExampleInterface
  {
  public:

    ExampleChemistryPossibleAtomTypesForAtom *Clone() const
    {
      return new ExampleChemistryPossibleAtomTypesForAtom( *this);
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

      // test default constructor
      chemistry::PossibleAtomTypesForAtom possible_atom_types_def;
      chemistry::PossibleAtomTypesForAtom possible_atom_types_c_trtrtrpi;
      chemistry::PossibleAtomTypesForAtom possible_atom_types_uncharged_carbons;

      // add C_TrTrTrPi
      possible_atom_types_c_trtrtrpi.AddAtomType( chemistry::GetAtomTypes().C_TrTrTrPi);

      // add C_TeTeTeTe, C_TrTrTrPi, C_DiDiPiPi, C_SPPP to possible_atom_types_uncharged_carbons
      possible_atom_types_uncharged_carbons.AddAtomType( chemistry::GetAtomTypes().C_TeTeTeTe);
      possible_atom_types_uncharged_carbons.AddAtomType( chemistry::GetAtomTypes().C_TrTrTrPi);
      possible_atom_types_uncharged_carbons.AddAtomType( chemistry::GetAtomTypes().C_DiDiPiPi);
      possible_atom_types_uncharged_carbons.AddAtomType( chemistry::GetAtomTypes().C_SPPP);

      // test clone
      util::ShPtr< chemistry::PossibleAtomTypesForAtom> possible_atom_types_sp( possible_atom_types_def.Clone());
      BCL_ExampleCheck( possible_atom_types_sp.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // test name
      BCL_ExampleCheck( possible_atom_types_def.GetClassIdentifier(), "bcl::chemistry::PossibleAtomTypesForAtom");

      // default constructor should not have any possible types
      BCL_ExampleCheck( possible_atom_types_def.GetNumberPossibleTypes(), 0);
      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.GetNumberPossibleTypes(), 1);
      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.GetNumberPossibleTypes(), 4);

      // default constructor has no defined types, so the most stable type is undefined
      BCL_ExampleCheck( possible_atom_types_def.GetMostStableType(), chemistry::GetAtomTypes().e_Undefined);
      BCL_ExampleCheck( possible_atom_types_def.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP), false);
      BCL_ExampleCheck( possible_atom_types_def.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP2), false);
      BCL_ExampleCheck( possible_atom_types_def.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP3), false);
      BCL_ExampleCheck( possible_atom_types_def.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_Unhybridized), false);
      BCL_ExampleCheck( possible_atom_types_def.GetNumberPossibleTypes(), 0);
      BCL_ExampleCheck( possible_atom_types_def.CouldBeConjugated(), false);
      BCL_ExampleCheck( possible_atom_types_def.MustBeConjugated(), false);
      BCL_ExampleCheck( possible_atom_types_def.GetMaxElectronsParticipatingInPiSystem(), 0);

      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.GetMostStableType(), chemistry::GetAtomTypes().C_TrTrTrPi);
      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP), false);
      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP2), true);
      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP3), false);
      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_Unhybridized), false);
      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.CouldBeConjugated(), true);
      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.MustBeConjugated(), true);

      // C_TrTrTrPi has 1 electron participating in pi systems
      BCL_ExampleCheck( possible_atom_types_c_trtrtrpi.GetMaxElectronsParticipatingInPiSystem(), 1);

      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.GetMostStableType(), chemistry::GetAtomTypes().C_TeTeTeTe);
      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP), true);
      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP2), true);
      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP3), true);
      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_Unhybridized), true);
      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.CouldBeConjugated(), true);
      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.MustBeConjugated(), false);

      // C_SPPP has two electrons participating in pi-systems
      BCL_ExampleCheck( possible_atom_types_uncharged_carbons.GetMaxElectronsParticipatingInPiSystem(), 2);

    //////////////////////
    // helper functions //
    //////////////////////

      const chemistry::ElementType chlorine( chemistry::GetElementTypes().e_Chlorine);
      const chemistry::ElementType carbon(   chemistry::GetElementTypes().e_Carbon);
      const chemistry::ElementType boron(   chemistry::GetElementTypes().e_Boron);

      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( chlorine, 1, 1, 0, false).GetMostStableType(),
        chemistry::GetAtomTypes().Cl_S2P2P2P,
        "FinalizeUnhybridized"
      );
      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( chlorine, 1, 1, 0, false).GetNumberPossibleTypes(),
        1,
        "FinalizeUnhybridized"
      );

      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( carbon, 3, 3, 0, false).GetMostStableType(),
        chemistry::GetAtomTypes().C_TrTrTr,
        "Internal Finalize"
      );
      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( carbon, 3, 3, 0, false).GetNumberPossibleTypes(),
        1,
        "Internal Finalize"
      );
      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( boron, 3, 3, 0, false).GetMostStableType(),
        chemistry::GetAtomTypes().B_TrTrTr,
        "Internal Finalize"
      );

      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( carbon, 4, 4, 0, false).GetMostStableType(),
        chemistry::GetAtomTypes().C_TeTeTeTe,
        "Internal Finalize"
      );

      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( carbon, 4, 4, 0, true).GetMostStableType().IsDefined(),
        false,
        "Finalize Aromatic with type that should not exist in an aromatic ring"
      );

      // nitrogen and oxygen should require additional calls to finalize to determine their type
      // except in cases w/ 4 bonds or 4 e- in bonds
      const chemistry::ElementType nitrogen( chemistry::GetElementTypes().e_Nitrogen);
      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( nitrogen, 3, 3, 0, false).GetNumberPossibleTypes(),
        2,
        "Nitrogen with three bonds could be trigonal or tetrahedral"
      );
      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( nitrogen, 3, 2, 0, false).GetNumberPossibleTypes(),
        2,
        "Nitrogen with two bonds, one unsaturated could be trigonal or tetrahedral"
      );
      BCL_ExampleIndirectCheck
      (
        chemistry::PossibleAtomTypesForAtom( nitrogen, 3, 2, 0, false).CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP2)
        && chemistry::PossibleAtomTypesForAtom( nitrogen, 3, 2, 0, false).CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP)
        && !chemistry::PossibleAtomTypesForAtom( nitrogen, 3, 2, 0, false).CouldHaveHybridization( chemistry::GetHybridOrbitalTypes().e_SP3),
        true,
        "Nitrogen with two bonds, one unsaturated could be trigonal or digonal"
      );

      // Write out the atom typing scheme
      std::stringstream test_output;
      chemistry::PossibleAtomTypesForAtom::WriteDetailedScheme( test_output);

      const std::string filename
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "atom_typing_scheme.txt")
      );
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, filename);
      if( !BCL_ExampleIndirectCheck( io::File::StreamsMatch( input, test_output), true, "atom typing scheme"))
      {
        // atom typing scheme differed, write out the new one; which can easily be compared against its value in the
        // base revision
        io::OFStream output;
        BCL_ExampleMustOpenOutputFile( output, filename);
        chemistry::PossibleAtomTypesForAtom::WriteDetailedScheme( output);
        io::File::CloseClearFStream( output);
      }

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
        ExampleInterface::TestBCLObjectIOForSymmetry
        (
          possible_atom_types_c_trtrtrpi,
          chemistry::PossibleAtomTypesForAtom()
        ),
        true,
        "Chemistry::PossibleAtomTypesForAtom read and/or write are broken"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleChemistryPossibleAtomTypesForAtom

  const ExampleClass::EnumType ExampleChemistryPossibleAtomTypesForAtom::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryPossibleAtomTypesForAtom())
  );
} // namespace bcl
