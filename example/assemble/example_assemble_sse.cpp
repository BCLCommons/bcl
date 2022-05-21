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
#include "assemble/bcl_assemble_sse.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSE :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSE *Clone() const
    {
      return new ExampleAssembleSSE( *this);
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
      BCL_MessageStd( "read pdbfile: " + pdb_filename);

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");

      // create protein model made out of coil, strand, and helix
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::SSE def_construct;
      BCL_Example_Check
      (
        def_construct.GetType() == biol::GetSSTypes().e_Undefined && def_construct.GetSize() == 0,
        "Default constructor failed"
      );

      // test constructor from SSType
      BCL_MessageStd( "testing constructor from SSType");
      assemble::SSE helix_no_sequence( biol::GetSSTypes().HELIX);
      BCL_Example_Check
      (
        helix_no_sequence.GetType() == biol::GetSSTypes().HELIX,
        "sstype of helix should be helix - but is: " + helix_no_sequence.GetType().GetName()
      );

      // test constructor from AASequence and SSType
      biol::AASequence helix_sequence( protein_model.GetSequences().FirstElement()->SubSequence( 47, 23));
      BCL_MessageStd( "testing constructor from AASequence and SSType");
      assemble::SSE helix_with_atoms( helix_sequence, biol::GetSSTypes().HELIX);
      BCL_Example_Check
      (
        helix_with_atoms.GetType() == biol::GetSSTypes().HELIX && helix_with_atoms.GetSize() == 23,
        "failed to construct properly - the SSE is: " + util::Format()( helix_with_atoms)
      );
      BCL_MessageStd( "Created SSE: " + helix_with_atoms.GetIdentification());

      // test copy constructor
      assemble::SSE copy_sse( helix_with_atoms);
      BCL_Example_Check
      (
        copy_sse == helix_with_atoms,
        "copy construct gave: " + util::Format()( copy_sse) + " instead of: " + util::Format()( helix_with_atoms)
      );

      // test Clone function
      util::ShPtr< assemble::SSE> clone_sse( helix_with_atoms.Clone());
      BCL_Example_Check
      (
        *clone_sse == helix_with_atoms,
        "Clone gave: " + util::Format()( *clone_sse) + " instead of: " + util::Format()( helix_with_atoms)
      );

      // test HardCopy function
      util::SiPtr< assemble::SSE> hard_copy_sse( helix_with_atoms.HardCopy());
      BCL_Example_Check
      (
        *hard_copy_sse == helix_with_atoms,
        "HardCopy gave: " + util::Format()( *hard_copy_sse) + " instead of: " + util::Format()( helix_with_atoms)
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::assemble::SSE");
      BCL_Example_Check
      (
        GetStaticClassName< assemble::SSE>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< assemble::SSE>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< assemble::SSE>() == clone_sse->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_sse->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // check GetGeometries
      BCL_Example_Check
      (
        helix_with_atoms.GetGeometries().GetSize() == 19,
        "GetGeometries should have 19 fragments but has: " + util::Format()( copy_sse.GetGeometries().GetSize())
      );

      // check GetMainAxis
      const double known_length( 34.5568);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.GetMainAxis().GetLength(), known_length),
        "GetMainAxis has a length of : " + util::Format()( helix_with_atoms.GetMainAxis().GetLength()) +
        " but should be: " + util::Format()( known_length)
      );

      // check GetExtent
      const double known_extent( 4.24);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.GetExtent( coord::GetAxes().e_X), known_extent),
        "GetExtent gives : " + util::Format()( helix_with_atoms.GetExtent( coord::GetAxes().e_X)) +
        " but should be: " + util::Format()( known_extent)
      );

      // check GetRadialExtent
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.GetRadialExtent(), known_extent),
        "GetRadialExtent gives : " + util::Format()( helix_with_atoms.GetRadialExtent()) + " but should be: " +
        util::Format()( known_extent)
      );

      // check GetOrientation
      const linal::Vector3D known_orientation_x_axis( 0.489166, 0.85493, 0.172659);
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          helix_with_atoms.GetOrientation().GetAxis( coord::GetAxes().e_X), known_orientation_x_axis
        ),
        "GetOrientation has X-axis : " +
        util::Format()( helix_with_atoms.GetOrientation().GetAxis( coord::GetAxes().e_X)) + " but should be: " +
        util::Format()( known_orientation_x_axis)
      );

      // check GetCenter
      const linal::Vector3D known_center( 16.2606, 34.742, -22.8355);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.GetCenter(), known_center),
        "GetCenter gives : " + util::Format()( helix_with_atoms.GetCenter()) + " but should be: " +
        util::Format()( known_center)
      );

      // check GetAxis
      const linal::Vector3D known_axis( 0.151004, -0.277985, 0.948642);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.GetAxis( coord::GetAxes().e_Z), known_axis),
        "GetAxis gives : " + util::Format()( helix_with_atoms.GetAxis( coord::GetAxes().e_Z)) + " but should be: " +
        util::Format()( known_axis)
      );

      // check GetRotation
      const linal::Vector3D known_rotation_x_axis( 0.489166, -0.859019, 0.151004);
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          helix_with_atoms.GetRotation().GetAxis( coord::GetAxes().e_X), known_rotation_x_axis
        ),
        "GetRotation has X-axis : " +
        util::Format()( helix_with_atoms.GetRotation().GetAxis( coord::GetAxes().e_X)) + " but should be: " +
        util::Format()( known_rotation_x_axis)
      );

      // check EndOfZ
      const linal::Vector3D known_end_of_z( 18.8698, 29.9389, -6.44451);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.EndOfZ(), known_end_of_z),
        "EndOfZ gives : " + util::Format()( helix_with_atoms.EndOfZ()) + " but should be: " +
        util::Format()( known_end_of_z)
      );

      // check BeginOfZ
      const linal::Vector3D known_begin_of_z( 13.6515, 39.5452, -39.2266);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.BeginOfZ(), known_begin_of_z),
        "BeginOfZ gives : " + util::Format()( helix_with_atoms.BeginOfZ()) + " but should be: " +
        util::Format()( known_begin_of_z)
      );

      // check SetType
      assemble::SSE new_strand( helix_with_atoms);
      new_strand.SetType( biol::GetSSTypes().STRAND);
      BCL_Example_Check
      (
        new_strand.GetType() == biol::GetSSTypes().STRAND,
        "SetType gave : " + util::Format()( new_strand.GetType()) + " but should be: " +
        util::Format()( biol::GetSSTypes().STRAND)
      );

      // check IsDefined
      BCL_Example_Check
      (
        helix_with_atoms.IsDefined() && !def_construct.IsDefined(),
        "IsDefined failed to give correct definition for " + util::Format()( helix_with_atoms) + " and " +
        util::Format()( def_construct)
      );

      // the sse is derived from aa sequence and now its secondary structure type
      assemble::SSE strand( biol::GetSSTypes().HELIX);

      BCL_MessageStd( "the sstype of the strand: " + strand.GetType().GetName());
      BCL_Example_Check
      (
        strand.GetType() == biol::GetSSTypes().HELIX,
        "sstype of strand should be helix - but is: " + strand.GetType().GetName()
      );

      // set the type to the correct strand type
      strand.SetType( biol::GetSSTypes().STRAND);

      BCL_MessageStd( "the sstype of the strand: " + strand.GetType().GetName());
      BCL_Example_Check
      (
        strand.GetType() == biol::GetSSTypes().STRAND,
        "sstype of strand should be strand - but is: " + strand.GetType().GetName()
      );

    ////////////////
    // operations //
    ////////////////

      // check Translate
      const linal::Vector3D translate( 1.0, 1.0, 1.0);
      helix_with_atoms.Translate( translate);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.GetCenter(), known_center + translate),
        "GetCenter (after translate) gives : " + util::Format()( helix_with_atoms.GetCenter()) + " but should be: " +
        util::Format()( known_center + translate)
      );

      // check Transform
      math::RotationMatrix3D rotation( helix_with_atoms.GetAxis( coord::GetAxes().e_Z), math::g_Pi);
      math::TransformationMatrix3D transform;
      transform( -known_center - translate);
      transform( rotation);
      transform( known_center);
      helix_with_atoms.Transform( transform);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( helix_with_atoms.GetCenter(), known_center),
        "GetCenter (after transform) gives : " + util::Format()( helix_with_atoms.GetCenter()) + " but should be: " +
        util::Format()( known_center)
      );

      // check Rotate
      math::RotationMatrix3D inverse_rotation( helix_with_atoms.GetAxis( coord::GetAxes().e_Z), math::g_Pi);
      helix_with_atoms.Translate( -known_center);
      helix_with_atoms.Rotate( inverse_rotation);
      helix_with_atoms.Translate( known_center);
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          helix_with_atoms.GetOrientation().GetAxis( coord::GetAxes().e_X), known_orientation_x_axis
        ),
        "GetOrientation (after rotate) gives : " + util::Format()( helix_with_atoms.GetOrientation()) +
        " but should be: " + util::Format()( known_orientation_x_axis)
      );

      // all aasequence functions are available
      BCL_MessageStd( "read sequence from 1eco_.fasta with aa backbone amino acids");
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.fasta"));
      strand = assemble::SSE
      (
        biol::AASequenceFactory::BuildSequenceFromFASTA( read, biol::GetAAClasses().e_AABackBone), biol::GetSSTypes().STRAND
      );
      io::File::CloseClearFStream( read);

      BCL_MessageStd( util::Format()( strand.GetSize()) + " residues are in the strand");
      BCL_Example_Check
      (
        strand.GetSize() == 122,
        "122 residues should have been read, but: " + util::Format()( strand.GetSize()) + " are in the strand"
      );

      // the backbone amino acids do not have coordinates yet, so we set them to an ideal conformation
      strand.SetToIdealConformationAtOrigin();

      // create SiPtrVector "sse_vector" and initialize with the helix SSEs of "protein_model"
      util::SiPtrVector< const assemble::SSE> sse_vector( protein_model.GetSSEs( biol::GetSSTypes().HELIX));

      // create SSE "helix" and initialize with the first SSE of "sse_vector"
      assemble::SSE helix( **sse_vector.Begin());

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename_b( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");

      // create protein model made out of coil, strand, and helix
      assemble::ProteinModel protein_model_b
      (
        Proteins::GetModel( pdb_filename_b, biol::GetAAClasses().e_AABackBone)
      );

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename_c( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // create protein model made out of coil, strand, and helix
      assemble::ProteinModel protein_model_c( Proteins::GetModel( pdb_filename_c));

      // store the first helix and first strand in 1UBI
      util::ShPtr< assemble::SSE> first_helix( ( *protein_model_c.GetSSEs( biol::GetSSTypes().HELIX).Begin())->Clone());
      util::ShPtr< assemble::SSE> first_strand( ( *protein_model_c.GetSSEs( biol::GetSSTypes().STRAND).Begin())->Clone());
      first_helix->SetToIdealConformationAtOrigin();
      first_strand->SetToIdealConformationAtOrigin();

      storage::List< storage::Pair< double, double> > degrees;
      storage::List< storage::Pair< double, double> > radians;
      degrees.PushBack( storage::Pair< double, double>( 0, 0));
      for( int i( -10); i != 10; i += 5)
      {
        degrees.PushBack( storage::Pair< double, double>( i, 0));
        degrees.PushBack( storage::Pair< double, double>( 0, i));
        degrees.PushBack( storage::Pair< double, double>( i, i));
      }

    //////////////////////////////////
    // Prepend/Append functionality //
    //////////////////////////////////

      // create the locators that finds helix 23-34 and strand 64-72
      assemble::LocatorSSE helix_23_34_locator( 'A', 23, 34);
      assemble::LocatorSSE strand_64_72_locator( 'A', 64, 72);

      // check append prepend functions
      BCL_MessageStd( "Checking Append and Prepend functions");

      // now locate the SSEs from the protein model
      util::SiPtr< const assemble::SSE> sp_helix_23_34( helix_23_34_locator.Locate( protein_model_b));
      util::SiPtr< const assemble::SSE> sp_strand_64_72( strand_64_72_locator.Locate( protein_model_b));

      // make sure the sse was correctly located
      BCL_ExampleIndirectAssert
      (
        sp_helix_23_34.IsDefined(),
        true,
        "Helix 23-34 was not located in 1ubi pdb " + util::Format()( helix_23_34_locator)
      );

      // make sure the sse was correctly located
      BCL_ExampleIndirectAssert
      (
        sp_strand_64_72.IsDefined(),
        true,
        "Strand 64-72 was not located in 1ubi pdb " + util::Format()( strand_64_72_locator)
      );

      // make a copy of the helix and strand
      util::ShPtr< assemble::SSE> helix_23_34( sp_helix_23_34->Clone());
      util::ShPtr< assemble::SSE> strand_64_72( sp_strand_64_72->Clone());

      // make a copy of the ShPtr to the sequence of chain A
      util::ShPtr< biol::AASequence> sp_sequence( protein_model_b.GetChain( 'A')->GetSequence());

      // create subsequences to prepend and append
      util::ShPtr< biol::AABase> helix_aa_to_prepend( sp_sequence->GetAA( 21));
      util::ShPtr< biol::AABase> helix_aa_to_append( sp_sequence->GetAA( 34));
      biol::AASequence helix_seq_to_prepend( sp_sequence->SubSequence( 18, 4));
      biol::AASequence helix_seq_to_append( sp_sequence->SubSequence( 34, 4));
      util::ShPtr< biol::AABase> strand_aa_to_prepend( sp_sequence->GetAA( 62));
      util::ShPtr< biol::AABase> strand_aa_to_append( sp_sequence->GetAA( 72));
      biol::AASequence strand_seq_to_prepend( sp_sequence->SubSequence( 60, 3));
      biol::AASequence strand_seq_to_append( sp_sequence->SubSequence( 72, 3));

      // prepend amino acid
      BCL_MessageStd( "Prepending amino acid to SSEs");
      // make a copy of the helix and strand
      util::ShPtr< assemble::SSE> helix_22_34( helix_23_34->Clone());
      helix_22_34->Prepend( *helix_aa_to_prepend);
      util::ShPtr< assemble::SSE> strand_63_72( strand_64_72->Clone());
      strand_63_72->Prepend( *strand_aa_to_prepend);
      // initialize chain
      assemble::Chain this_chain_extended( sp_sequence);
      // insert this fragment into chain
      this_chain_extended.Insert( helix_22_34);
      this_chain_extended.Insert( strand_63_72);
      // output the sequence information
      BCL_MessageStd( "helix_22_34 " + helix_22_34->GetIdentification());
      BCL_MessageStd( "strand_63_72 " + strand_63_72->GetIdentification());
      // open the write stream and write the pdb
      std::string extended_sse_filename( AddExampleOutputPathToFilename( assemble::SSE(), "1ubi.pdb.extended_a"));
      Proteins::WriteChainToPDB( this_chain_extended, extended_sse_filename);

      // append amino acid
      BCL_MessageStd( "Appending amino acid to SSEs");
      // make a copy of the helix and strand
      util::ShPtr< assemble::SSE> helix_23_35( helix_23_34->Clone());
      helix_23_35->Append( *helix_aa_to_append);
      util::ShPtr< assemble::SSE> strand_64_73( strand_64_72->Clone());
      strand_64_73->Append( *strand_aa_to_append);
      this_chain_extended.ReplaceWithOverlapping( helix_23_35);
      this_chain_extended.ReplaceWithOverlapping( strand_64_73);
      BCL_MessageStd( "helix_23_35 " + helix_23_35->GetIdentification());
      BCL_MessageStd( "strand_64_73 " + strand_64_73->GetIdentification());
      extended_sse_filename = AddExampleOutputPathToFilename( assemble::SSE(), "1ubi.pdb.extended_b");
      Proteins::WriteChainToPDB( this_chain_extended, extended_sse_filename);

      // prepend sequence
      BCL_MessageStd( "Prepending sequence to SSEs");
      // make a copy of the helix and strand
      util::ShPtr< assemble::SSE> helix_19_34( helix_23_34->Clone());
      helix_19_34->PrependSequence( helix_seq_to_prepend);
      util::ShPtr< assemble::SSE> strand_61_72( strand_64_72->Clone());
      strand_61_72->PrependSequence( strand_seq_to_prepend);
      this_chain_extended.ReplaceWithOverlapping( helix_19_34);
      this_chain_extended.ReplaceWithOverlapping( strand_61_72);
      BCL_MessageStd( "helix_19_34 " + helix_19_34->GetIdentification());
      BCL_MessageStd( "strand_61_72 " + strand_61_72->GetIdentification());
      extended_sse_filename = AddExampleOutputPathToFilename( assemble::SSE(), "1ubi.pdb.extended_c");
      Proteins::WriteChainToPDB( this_chain_extended, extended_sse_filename);

      // append sequence
      BCL_MessageStd( "Appending sequence to SSEs");
      // make a copy of the helix and strand
      util::ShPtr< assemble::SSE> helix_23_38( helix_23_34->Clone());
      helix_23_38->AppendSequence( helix_seq_to_append);
      util::ShPtr< assemble::SSE> strand_64_75( strand_64_72->Clone());
      strand_64_75->AppendSequence( strand_seq_to_append);
      this_chain_extended.ReplaceWithOverlapping( helix_23_38);
      this_chain_extended.ReplaceWithOverlapping( strand_64_75);
      BCL_MessageStd( "helix_23_38 " + helix_23_38->GetIdentification());
      BCL_MessageStd( "strand_64_75 " + strand_64_75->GetIdentification());
      extended_sse_filename = AddExampleOutputPathToFilename( assemble::SSE(), "1ubi.pdb.extended_d");
      Proteins::WriteChainToPDB( this_chain_extended, extended_sse_filename);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "Testing write and read functionalities");
      WriteBCLObject( *sp_helix_23_34);
      assemble::SSE read_sse;
      ReadBCLObject( read_sse);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSE

  const ExampleClass::EnumType ExampleAssembleSSE::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSE())
  );

} // namespace bcl
