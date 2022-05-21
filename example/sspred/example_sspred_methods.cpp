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
#include "sspred/bcl_sspred_methods.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "io/bcl_io_file.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sspred_methods.cpp
  //!
  //! @author karakam
  //! @date Nov 10, 2008
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredMethods :
    public ExampleInterface
  {
  public:

    ExampleSspredMethods *Clone() const
    {
      return new ExampleSspredMethods( *this);
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
      // declare ifstream and MethodHandler
      io::IFStream read;
      sspred::MethodHandler method_handler;

      // read seq from fasta
      BCL_MessageStd( "read fasta: " + util::Format()( AddExampleInputPathToFilename( e_Biology, "1eco_.fasta")));
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.fasta"));
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // write fasta to util::GetLogger()
      seq.WriteFasta( util::GetLogger());

      // read blast profile
      BCL_MessageStd( "reading BLAST");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq);
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_MessageStd( "reading PSIPRED");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.psipred_ss2"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read jufo ss prediction
      BCL_MessageStd( "reading JUFO");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.jufo"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_JUFO);
      io::File::CloseClearFStream( read);

      // read sam ss prediction
      BCL_MessageStd( "reading SAM");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.rdb6"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_SAM);
      io::File::CloseClearFStream( read);

      // read conpred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_CONPRED.txt"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_CONPRED);
      io::File::CloseClearFStream( read);

      // read tmmod ss prediction
      BCL_MessageStd( "reading TMMOD");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_TMMOD.txt"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_TMMOD);
      io::File::CloseClearFStream( read);

      // read B2TMPRED ss prediction
      BCL_MessageStd( "reading B2TMPRED");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_B2TMPRED.txt"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_B2TMPRED);
      io::File::CloseClearFStream( read);

      // create a subsequence
      const biol::AASequence sub_seq( seq.SubSequence( 0, 5));

      // iterate over residues in the subsequence
      for
      (
        biol::AASequence::const_iterator aa_itr( sub_seq.Begin()), aa_itr_end( sub_seq.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // write out first residue so it prints out m_data
        ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_PSIPRED)->WriteThreeStatePredictions( util::GetLogger());
        ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_JUFO)->WriteThreeStatePredictions( util::GetLogger());
        ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_SAM)->WriteThreeStatePredictions( util::GetLogger());
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredMethods

  const ExampleClass::EnumType ExampleSspredMethods::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredMethods())
  );

} // namespace bcl
