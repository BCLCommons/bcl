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
#include "biol/bcl_biol_aa_base.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_back_bone.h"
#include "biol/bcl_biol_aa_complete.h"
#include "biol/bcl_biol_rotamer.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_base.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAABase :
    public ExampleInterface
  {
  public:

    ExampleBiolAABase *Clone() const
    {
      return new ExampleBiolAABase( *this);
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
      // construct the backbone amino acids
      biol::AABackBone backbone_a( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().MET, 1, 1, 'X', 'A')));
      biol::AABackBone backbone_b( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().GLN, 2, 2, 'X', 'A')));
      biol::AABackBone backbone_c( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ILE, 3, 3, 'X', 'A')));

      // construct the atoms necessary for building the backbone
      const biol::Atom atom_a_N ( linal::Vector3D( 27.343, 24.294, 2.683), biol::GetAtomTypes().N );
      const biol::Atom atom_a_CA( linal::Vector3D( 26.381, 25.361, 2.894), biol::GetAtomTypes().CA);
      const biol::Atom atom_a_C ( linal::Vector3D( 26.997, 26.557, 3.583), biol::GetAtomTypes().C );
      const biol::Atom atom_a_O ( linal::Vector3D( 27.973, 26.439, 4.351), biol::GetAtomTypes().O );
      const biol::Atom atom_a_CB( linal::Vector3D( 25.112, 24.879, 3.647), biol::GetAtomTypes().CB);
      const biol::Atom atom_b_N ( linal::Vector3D( 26.410, 27.694, 3.332), biol::GetAtomTypes().N );
      const biol::Atom atom_b_CA( linal::Vector3D( 26.865, 28.934, 3.898), biol::GetAtomTypes().CA);
      const biol::Atom atom_b_C ( linal::Vector3D( 26.136, 29.272, 5.204), biol::GetAtomTypes().C );
      const biol::Atom atom_b_O ( linal::Vector3D( 24.911, 29.191, 5.289), biol::GetAtomTypes().O );
      const biol::Atom atom_b_CB( linal::Vector3D( 26.705, 30.075, 2.869), biol::GetAtomTypes().CB);
      const biol::Atom atom_c_N ( linal::Vector3D( 26.906, 29.647, 6.221), biol::GetAtomTypes().N );
      const biol::Atom atom_c_CA( linal::Vector3D( 26.337, 30.097, 7.481), biol::GetAtomTypes().CA);
      const biol::Atom atom_c_C ( linal::Vector3D( 26.951, 31.407, 7.906), biol::GetAtomTypes().C );
      const biol::Atom atom_c_O ( linal::Vector3D( 28.012, 31.794, 7.418), biol::GetAtomTypes().O );
      const biol::Atom atom_c_CB( linal::Vector3D( 26.384, 29.049, 8.615), biol::GetAtomTypes().CB);

      // construct atom lists and set atoms for the corresponding amino acid
      util::SiPtrVector< const biol::Atom> atoms_a
      (
        util::SiPtrVector< const biol::Atom>::Create
        (
          atom_a_N,
          atom_a_CA,
          atom_a_C,
          atom_a_O,
          atom_a_CB
        )
      );
      backbone_a.SetAtoms( atoms_a);

      util::SiPtrVector< const biol::Atom> atoms_b
      (
        util::SiPtrVector< const biol::Atom>::Create
        (
          atom_b_N,
          atom_b_CA,
          atom_b_C,
          atom_b_O,
          atom_b_CB
        )
      );
      backbone_b.SetAtoms( atoms_b);

      util::SiPtrVector< const biol::Atom> atoms_c
      (
        util::SiPtrVector< const biol::Atom>::Create
        (
          atom_c_N,
          atom_c_CA,
          atom_c_C,
          atom_c_O,
          atom_c_CB
        )
      );
      backbone_c.SetAtoms( atoms_c);

      // construct a new aa complete
      biol::AAComplete aa_complete( backbone_a);
      const biol::Atom atom_CG( linal::Vector3D( 25.341, 24.685,  5.142), biol::GetAtomTypes().CG);
      const biol::Atom atom_SD( linal::Vector3D( 23.944, 23.887,  5.986), biol::GetAtomTypes().SD);
      const biol::Atom atom_CE( linal::Vector3D( 24.546, 23.890,  7.694), biol::GetAtomTypes().CE);
      aa_complete.SetAtom( atom_CG);
      aa_complete.SetAtom( atom_SD);
      aa_complete.SetAtom( atom_CE);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a reference from an AABackBone
      BCL_MessageStd( "testing constructors from an AABackbone's AABase base");
      biol::AABase &aa_base_a( backbone_a);
      biol::AABase &aa_base_b( backbone_b);
      biol::AABase &aa_base_c( backbone_c);

      // create a aa base from another aa base
      BCL_MessageStd( "testing Clone()")
      util::ShPtr< biol::AABase> sp_aa_base_a( aa_base_a.Clone());

      // use the Empty() constructor
      BCL_MessageStd( "testing Empty() constructor")
      util::ShPtr< biol::AABase> sp_aa_base_empty_a( aa_base_a.Empty( util::ShPtr< biol::AAData>( new biol::AAData())));
      util::ShPtr< biol::AABase> sp_aa_base_empty_b( aa_base_a.Empty( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().VAL))));

    /////////////////
    // data access //
    /////////////////

      // check the class identifier
      BCL_ExampleCheck( GetStaticClassName( aa_base_a), "bcl::biol::AABase");

      // class identifier
      BCL_ExampleCheck( aa_base_a.GetClassIdentifier(), backbone_a.GetClassIdentifier());

      // check GetType()
      BCL_ExampleCheck( aa_base_a.GetType(), biol::GetAATypes().MET);

      // check GetSeqID()
      BCL_ExampleCheck( aa_base_b.GetSeqID(), 2);

      // check GetPdbID()
      BCL_ExampleCheck( aa_base_b.GetPdbID(), 2);

      // check GetPdbICode()
      BCL_ExampleCheck( aa_base_b.GetPdbICode(), 'X');

      // check GetChainID()
      BCL_ExampleCheck( aa_base_b.GetChainID(), 'A');

      // check CalculatePhi()
      const double expected_phi( -1.62432);
      BCL_ExampleCheckWithinTolerance
      (
        expected_phi,
        aa_base_b.CalculatePhi( aa_base_a.GetAtom( biol::GetAtomTypes().C)),
        0.001
      );
      // check CalculatePsi()
      const double expected_psi( 2.3137);
      BCL_ExampleCheckWithinTolerance
      (
        expected_psi,
        aa_base_b.CalculatePsi( aa_base_c.GetAtom( biol::GetAtomTypes().N)),
        0.001
      );

      // CalculateSideChainDihedralAngles
      {
        const storage::Vector< double> expected_dihedral
        (
          storage::Vector< double>::Create( 1.33168622, -2.99015, -3.10577)
        );
        BCL_MessageDbg
        (
          "expected angles are " + util::Format()( expected_dihedral) +
          " calculated angles are " + util::Format()( aa_complete.CalculateSideChainDihedralAngles())
        );
        BCL_ExampleCheck
        (
          math::EqualWithinTolerance
          (
            aa_complete.CalculateSideChainDihedralAngles().GetAngle( biol::ChiAngle::e_One, math::Angle::e_Radian),
            expected_dihedral( 0), 0.0001
          ), true
        );
        BCL_ExampleCheck
        (
          math::EqualWithinTolerance
          (
            aa_complete.CalculateSideChainDihedralAngles().GetAngle( biol::ChiAngle::e_Two, math::Angle::e_Radian),
            expected_dihedral( 1), 0.0001
          ), true
        );
        BCL_ExampleCheck
        (
          math::EqualWithinTolerance
          (
            aa_complete.CalculateSideChainDihedralAngles().GetAngle( biol::ChiAngle::e_Three, math::Angle::e_Radian),
            expected_dihedral( 2), 0.0001
          ), true
        );
      }
      {
        BCL_ExampleCheck( backbone_a.CalculateSideChainDihedralAngles().GetSize(), 3);
      }

    ///////////////
    // operators //
    ///////////////

      // check operator = (AA_BASE)
      // first create a new empty AABase
      util::ShPtr< biol::AABase> sp_aa_base_empty_c( aa_base_b.Empty( util::ShPtr< biol::AAData>( new biol::AAData())));
      *sp_aa_base_empty_c = aa_base_a;
      BCL_ExampleCheck( sp_aa_base_empty_c->GetType(), aa_base_a.GetType());
      BCL_ExampleCheck( sp_aa_base_empty_c->GetSeqID(), aa_base_a.GetSeqID());
      BCL_ExampleCheck( sp_aa_base_empty_c->GetPdbID(), aa_base_a.GetPdbID());

    ////////////////
    // operations //
    ////////////////

      // check CalculateCenterOfMassOfSideChain()
      // first check with an AAComplete
      const linal::Vector3D expected_center_of_mass_sc_a( 24.7357, 24.3353, 5.61725);
      BCL_ExampleCheckWithinTolerance
      (
        expected_center_of_mass_sc_a,
        aa_complete.CalculateCenterOfMassOfSideChain(),
        0.001
      );
      // now check with just an AABackbone, this should just return the coordinates of the CB atom
      const linal::Vector3D expected_center_of_mass_sc_b( 25.112, 24.879, 3.647);
      BCL_ExampleCheckWithinTolerance
      (
        expected_center_of_mass_sc_b,
        aa_base_a.CalculateCenterOfMassOfSideChain(),
        0.001
      );

    //////////////////////
    // input and output //
    //////////////////////

      // check GetIdentification()
      const std::string expected_identification( "    1 M MET U");
      BCL_ExampleCheck( aa_base_a.GetIdentification(), expected_identification);

      // write the object
      WriteBCLObject( aa_base_b, ".aa_base");

      // check the read object
      ReadBCLObject( *sp_aa_base_empty_b, ".aa_base");
      BCL_ExampleCheck( aa_base_b.GetType(), sp_aa_base_empty_b->GetType());
      BCL_ExampleCheck( aa_base_b.GetSeqID(), sp_aa_base_empty_b->GetSeqID());
      BCL_ExampleCheck( aa_base_b.GetPdbID(), sp_aa_base_empty_b->GetPdbID());

      // check WriteFasta()
      std::stringstream sstream;
      aa_base_a.WriteFasta( sstream);
      const std::string expected_fasta_output( "M");
      BCL_ExampleCheck( sstream.str(), expected_fasta_output);

      // check ReadFasta
      sp_aa_base_empty_a->ReadFasta( sstream, 1);
      BCL_ExampleCheck( sp_aa_base_empty_a->GetType(), biol::GetAATypes().MET);
      BCL_ExampleCheck( sp_aa_base_empty_a->GetSeqID(), 1);
      BCL_ExampleCheck( sp_aa_base_empty_a->GetPdbID(), 1);
      sstream.clear();

    //////////////////////
    // helper functions //
    //////////////////////

      // check operator ==
      BCL_ExampleCheck( aa_base_a == aa_base_a, true);

      // check FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B)
      const double expected_first_sc_distance( 5.49011);
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          expected_first_sc_distance, biol::FirstSidechainAtomDistance( aa_base_a, aa_base_b)
        ), true
      );

      // check SequenceSeparation( AMINO_ACID_A, AMINO_ACID_B)
      const size_t expected_sequence_separation( 1);
      BCL_ExampleCheck( biol::SequenceSeparation( aa_base_a, aa_base_c), expected_sequence_separation);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAABase

  const ExampleClass::EnumType ExampleBiolAABase::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAABase())
  );

} // namespace bcl

