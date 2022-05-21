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
#include "biol/bcl_biol_aa_sequence.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_sequence.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAASequence :
    public ExampleInterface
  {
  public:

    ExampleBiolAASequence *Clone() const
    { return new ExampleBiolAASequence( *this);}

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
      // declare to AASequences from AA and AACaCb
      biol::AASequence seq2( biol::GetAAClasses().e_AACaCb, 10, 'A', "> fun1");
      io::IFStream read;
      sspred::MethodHandler method_handler;

      // declare MethodHandler

      // read seq1 from fasta
      BCL_MessageStd( "read fasta: " + AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      BCL_MessageStd( "the fasta sequence is: ");
      seq1.WriteFasta( util::GetLogger());

      // read psipred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.psipred_ss2").c_str());
      method_handler.ReadPredictionsForAASequence( read, seq1, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.rdb6"));
      method_handler.ReadPredictionsForAASequence( read, seq1, sspred::GetMethods().e_SAM);
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq1);
      io::File::CloseClearFStream( read);

      // write fasta to util::GetLogger()
      seq1.WriteFasta( util::GetLogger());
      seq2.WriteFasta( util::GetLogger());
      // set seq1 as subsequence of original sequence
      seq1 = seq1.SubSequence( 10, 8);

      // output fastaTOP, psipred in sam format, blast
      seq1.WriteFasta( util::GetLogger());
      method_handler.WritePredictionsForAASequence( util::GetLogger(), seq1, sspred::GetMethods().e_PSIPRED);
      biol::BlastProfileHandler::WriteProfileForAASequence( util::GetLogger(), seq1);

      //read sequence from pdb for chopping
      BCL_MessageStd( "read pdb file: " + AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

      //instantiate sequences
      BCL_MessageStd( "building sequences from pdb chains");
      util::ShPtrVector< biol::AASequence> sequences
      (
        pdb::Factory( biol::GetAAClasses().e_AACaCb).AASequencesFromPDB( pdb)
      );

      BCL_MessageStd( "sequences are built");

      BCL_MessageStd
      (
        "the number of residues in the first chain: " + util::Format()( sequences.FirstElement()->GetSize())
      );
      sequences.FirstElement()->WriteFasta( util::GetLogger());

      //chop the first sequence
      for( size_t i( 23); i != 15; --i)
      {
        BCL_MessageStd( sequences.FirstElement()->SubSequence( 0, i).GetSequenceIdentification());

        // chop the sequence and store it
        util::ShPtrVector< biol::AASequence>
          chopped_sequences( sequences.FirstElement()->SubSequence( 0, i).ChopSequence( 5));
        // iterate over chopped sequences
        size_t ctr( 1);
        for
        (
          util::ShPtrVector< biol::AASequence>::const_iterator seq_itr( chopped_sequences.Begin()),
            seq_itr_end( chopped_sequences.End());
          seq_itr != seq_itr_end; ++seq_itr, ++ctr
        )
        {
          BCL_MessageStd
          (
            "chopped seq # " + util::Format()( ctr) + " " + ( *seq_itr)->GetSequenceIdentification()
          );
        }
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleBiolAASequence

  const ExampleClass::EnumType ExampleBiolAASequence::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAASequence())
  );

} // namespace bcl
