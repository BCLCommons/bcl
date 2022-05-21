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
#include "chemistry/bcl_chemistry_atom_vector.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_vector.cpp
  //!
  //! @author kothiwsk
  //! @date 01/06/2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomVector :
    public ExampleInterface
  {
  public:

    ExampleChemistryAtomVector *Clone() const
    {
      return new ExampleChemistryAtomVector( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create initializer for  atom type for the molecule HCN ( hydrogen cyanide)
      storage::Vector< sdf::AtomInfo> initialize
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().H_S, chemistry::e_NonChiral),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_DiDiPiPi, chemistry::e_NonChiral),
          sdf::AtomInfo( chemistry::GetAtomTypes().N_Di2DiPiPi, chemistry::e_NonChiral)
        )
      );

      // create the bond triplet for HCN
      storage::Vector< sdf::BondInfo> bonds
      (
        storage::Vector< sdf::BondInfo>::Create
        (
          sdf::BondInfo( 0, 1, chemistry::GetConstitutionalBondTypes().e_NonConjugatedSingleBond),
          sdf::BondInfo( 1, 2, chemistry::GetConstitutionalBondTypes().e_ConjugatedTripleBond)
        )
      );

      // default constructor
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> default_atom_vector;

      // constructor from HCN initalizer and bonds
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> HCN_atom_vector( initialize, bonds);

      // copy constructor from HCN
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> copy_atom_vector( HCN_atom_vector);

    /////////////////
    // data access //
    /////////////////

      // test GetSize
      BCL_ExampleCheck( HCN_atom_vector.GetSize(), initialize.GetSize());

      // test operator()
      BCL_ExampleCheck( HCN_atom_vector( 1).GetAtomType(), initialize( 1).GetAtomType());

      // test first()
      BCL_ExampleCheck( HCN_atom_vector.First().GetAtomType(), initialize.FirstElement().GetAtomType());

      // test last()
      BCL_ExampleCheck( HCN_atom_vector.Last().GetAtomType(), initialize.LastElement().GetAtomType());

      // test the copied atom vector

      BCL_ExampleCheck( copy_atom_vector.GetSize(), initialize.GetSize());
      BCL_ExampleCheck( copy_atom_vector( 1).GetAtomType(), initialize( 1).GetAtomType());
      BCL_ExampleCheck( copy_atom_vector.First().GetAtomType(), initialize.FirstElement().GetAtomType());
      BCL_ExampleCheck( copy_atom_vector.Last().GetAtomType(), initialize.LastElement().GetAtomType());

      // test that copy doesnt point to the same object that original points to
      BCL_ExampleCheck( HCN_atom_vector.Begin() != copy_atom_vector.Begin(), true);

      // check bidirectionality
      BCL_ExampleCheck( HCN_atom_vector( 1).GetBonds().GetSize(), 2);
      BCL_ExampleCheck( HCN_atom_vector( 2).GetBonds().GetSize(), 1);
      BCL_ExampleCheck
      (
        HCN_atom_vector( 2).GetBonds().FirstElement().GetTargetAtom().GetAtomType(),
        HCN_atom_vector( 1).GetAtomType()
      );
    ///////////////
    // operators //
    ///////////////

      // test operator =
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> copy1_atom_vector;
      copy1_atom_vector = HCN_atom_vector;
      BCL_ExampleCheck( copy1_atom_vector.GetSize(), initialize.GetSize());
      BCL_ExampleCheck( copy1_atom_vector( 1).GetAtomType(), initialize( 1).GetAtomType());
      BCL_ExampleCheck( copy1_atom_vector.First().GetAtomType(), initialize.FirstElement().GetAtomType());
      BCL_ExampleCheck( copy1_atom_vector.Last().GetAtomType(), initialize.LastElement().GetAtomType());

      copy1_atom_vector = default_atom_vector;

      BCL_ExampleCheck( HCN_atom_vector.GetBondInfo(), bonds);
      BCL_ExampleCheck( HCN_atom_vector.GetBondInfo(), copy_atom_vector.GetBondInfo());

      BCL_ExampleCheck( HCN_atom_vector.IsEmpty(), false);
      BCL_ExampleCheck( default_atom_vector.IsEmpty(), true);

    ////////////////
    // operations //
    ////////////////

      // test reorder
      storage::Vector< size_t> reordered_indices( storage::Vector< size_t>::Create( 2, 0, 1));
      chemistry::AtomVector< chemistry::AtomConstitutionalShared> NCH_atom_vector( HCN_atom_vector);
      NCH_atom_vector.Reorder( reordered_indices);
      BCL_ExampleCheck( NCH_atom_vector.GetSize(), initialize.GetSize());
      for( size_t i( 0); i < 3; ++i)
      {
        BCL_ExampleIndirectCheck
        (
          NCH_atom_vector( i).GetAtomType(),
          HCN_atom_vector( reordered_indices( i)).GetAtomType(),
          "reorder"
        );
        BCL_ExampleIndirectCheck
        (
          NCH_atom_vector( i).GetBonds().GetSize(),
          HCN_atom_vector( reordered_indices( i)).GetBonds().GetSize(),
          "reorder"
        );
        BCL_ExampleIndirectCheck
        (
          NCH_atom_vector( i).GetBonds().FirstElement().GetTargetAtom().GetAtomType(),
          HCN_atom_vector( reordered_indices( i)).GetBonds().LastElement().GetTargetAtom().GetAtomType(),
          "reorder"
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryAtomVector

  const ExampleClass::EnumType ExampleChemistryAtomVector::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomVector())
  );

} // namespace bcl

