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
#include "biol/bcl_biol_aa_sequence_flexibility.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_sequence_flexibility.cpp
  //!
  //! @author karakam
  //! @date Jan 22, 2011
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAASequenceFlexibility :
    public ExampleInterface
  {
  public:

    ExampleBiolAASequenceFlexibility *Clone() const
    {
      return new ExampleBiolAASequenceFlexibility( *this);
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

    //! @brief convenience function to write phi and psi values to a string
    //! @param PHI_PSI phi psi values
    //! @return string that has the phi psi values
    std::string WritePhiPsiToString( const storage::VectorND< 2, double> &PHI_PSI) const
    {
      std::string temp
      (
        "phi: " + util::Format()( PHI_PSI.First()) + " ( " + util::Format()( math::Angle::Degree( PHI_PSI.First())) +
        ")\tpsi: " + util::Format()( PHI_PSI.Second()) + " ( " + util::Format()( math::Angle::Degree( PHI_PSI.Second())) + ")"
      );

      return temp;
    }

    int Run() const
    {
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // locate helix 23-34
      const util::ShPtr< assemble::SSE> helix_23_34( Proteins::GetSSE( pdb_filename, 'A', 23, 34));
      const util::ShPtr< assemble::SSE> strand_64_72( Proteins::GetSSE( pdb_filename, 'A', 64, 72));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      biol::AASequenceFlexibility flexibility;

    ////////////////
    // operations //
    ////////////////

      // make a copy of the helix 23_34
      util::ShPtr< assemble::SSE> sp_helix( helix_23_34->Clone());

      // calculate phi, psi for residue 23_34
      const storage::VectorND< 2, double> phi_psi_29( sp_helix->CalculatePhiPsi( 29));
      BCL_MessageStd( "residue 29 native\t\t" + WritePhiPsiToString( phi_psi_29));

      // test CalculatePhiPsiChange function
      BCL_MessageStd( "Checking CalculatePhiPsiChange()");
      // construct the expected values
      const storage::VectorND< 2, double> phi_psi_29_new( math::Angle::Radian( -70.0), math::Angle::Radian( -35.0));
      const storage::VectorND< 2, double> phi_psi_29_change
      (
        phi_psi_29_new.First() - phi_psi_29.First(),
        phi_psi_29_new.Second() - phi_psi_29.Second()
      );
      BCL_MessageStd
      (
        "residue 29 change expected\t" + WritePhiPsiToString( phi_psi_29_change)
      );

      // test calculate change function
      const storage::VectorND< 2, double> phi_psi_29_calc_change
      (
        flexibility.CalculatePhiPsiChange( *sp_helix, 29, phi_psi_29_new)
      );

      BCL_MessageStd
      (
        "residue 29 change calc\t" + WritePhiPsiToString( phi_psi_29_calc_change)
      );

      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( phi_psi_29_change.First(), phi_psi_29_calc_change.First()) &&
        math::EqualWithinTolerance( phi_psi_29_change.Second(), phi_psi_29_calc_change.Second()),
        true,
        "The CalculatePhiPsiChange function failed"
      );

      // test GetNumberDifferentPhiPsi()
      BCL_MessageStd( "Checking GetNumberDifferentPhiPsi() with 0 changes");
      BCL_ExampleCheck( flexibility.GetNumberDifferentPhiPsi( *sp_helix, *sp_helix), 0);

      // test set phi psi
      BCL_MessageStd( "Checking SetPhiPsi()");
      flexibility.SetPhiPsi( *sp_helix, 29, phi_psi_29_new, flexibility.e_NTerminal);
      storage::VectorND< 2, double> phi_psi_29_after( sp_helix->CalculatePhiPsi( 29));
      BCL_MessageStd( "residue 29 after set\t\t" + WritePhiPsiToString( phi_psi_29_after));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( phi_psi_29_new.First(), phi_psi_29_after.First()) &&
        math::EqualWithinTolerance( phi_psi_29_new.Second(), phi_psi_29_after.Second()),
        true,
        "SetPhiPsi function on sp_helix failed"
      );

      // test GetNumberDifferentPhiPsi()
      BCL_MessageStd( "Checking GetNumberDifferentPhiPsi() with 1 change");
      BCL_ExampleCheck( flexibility.GetNumberDifferentPhiPsi( *helix_23_34, *sp_helix), 1);
      Proteins::WriteSSEToPDB( sp_helix, AddExampleOutputPathToFilename( flexibility, "phi_psi_helix_a.pdb"));

      // make a copy of the helix again
      sp_helix = util::ShPtr< assemble::SSE>( helix_23_34->Clone());

      // test change phi psi
      BCL_MessageStd( "Checking ChangePhiPsi()");
      flexibility.ChangePhiPsi( *sp_helix, 29, phi_psi_29_change, flexibility.e_CTerminal);
      phi_psi_29_after = storage::VectorND< 2, double>( sp_helix->CalculatePhiPsi( 29));
      BCL_MessageStd( "residue 29 after change\t\t" + WritePhiPsiToString( phi_psi_29_after));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( phi_psi_29_new.First(), phi_psi_29_after.First()) &&
        math::EqualWithinTolerance( phi_psi_29_new.Second(), phi_psi_29_after.Second()),
        true,
        "ChangePhiPsi function on sp_helix failed"
      );
      Proteins::WriteSSEToPDB( sp_helix, AddExampleOutputPathToFilename( flexibility, "phi_psi_helix_b.pdb"));

      // do the same with strand now
      // calculate phi, psi for residue 23_34
      const storage::VectorND< 2, double> phi_psi_67( strand_64_72->CalculatePhiPsi( 67));
      BCL_MessageStd( "residue 67 native\t\t" + WritePhiPsiToString( phi_psi_67));

      // construct the expected values
      const storage::VectorND< 2, double> phi_psi_67_new( math::Angle::Radian( 150.0), math::Angle::Radian( 110.0));
      const storage::VectorND< 2, double> phi_psi_67_change
      (
        phi_psi_67_new.First() - phi_psi_67.First(),
        phi_psi_67_new.Second() - phi_psi_67.Second()
      );

      // make a copy of the helix 23_34
      util::ShPtr< assemble::SSE> sp_strand( strand_64_72->Clone());

      // test set phi psi
      BCL_MessageStd( "Checking SetPhiPsi()");
      flexibility.SetPhiPsi( *sp_strand, 67, phi_psi_67_new, flexibility.e_Bidirectional);
      storage::VectorND< 2, double> phi_psi_67_after( sp_strand->CalculatePhiPsi( 67));
      BCL_MessageStd( "residue 67 after set\t\t" + WritePhiPsiToString( phi_psi_67_after));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( phi_psi_67_new.First(), phi_psi_67_after.First()) &&
        math::EqualWithinTolerance( phi_psi_67_new.Second(), phi_psi_67_after.Second()),
        true,
        "SetPhiPsi function on sp_strand failed"
      );
      Proteins::WriteSSEToPDB
      (
        sp_strand, AddExampleOutputPathToFilename( flexibility, "phi_psi_strand_a.pdb")
      );

      // make a copy of the helix again
      sp_strand = util::ShPtr< assemble::SSE>( strand_64_72->Clone());

      // test change phi psi
      BCL_MessageStd( "Checking ChangePhiPsi()");
      flexibility.ChangePhiPsi( *sp_strand, 67, phi_psi_67_change, flexibility.e_CTerminal);
      phi_psi_67_after = storage::VectorND< 2, double>( sp_strand->CalculatePhiPsi( 67));
      BCL_MessageStd( "residue 67 after change\t" + WritePhiPsiToString( phi_psi_67_after));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( phi_psi_67_new.First(), phi_psi_67_after.First()) &&
        math::EqualWithinTolerance( phi_psi_67_new.Second(), phi_psi_67_after.Second()),
        true,
        "ChangePhiPsi function on sp_strand failed"
      );
      Proteins::WriteSSEToPDB
      (
        sp_strand, AddExampleOutputPathToFilename( flexibility, "phi_psi_strand_b.pdb")
      );

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAASequenceFlexibility

  const ExampleClass::EnumType ExampleBiolAASequenceFlexibility::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAASequenceFlexibility())
  );

} // namespace bcl
