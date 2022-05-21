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
#include "align/bcl_align_word_generator_subsequences.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_handler_pir.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_assignment_blosum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_word_generator_subsequences.cpp
  //!
  //! @author heinzes1
  //! @date Mar 24, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignWordGeneratorSubsequences :
    public ExampleInterface
  {
  public:

    ExampleAlignWordGeneratorSubsequences *Clone() const
    {
      return new ExampleAlignWordGeneratorSubsequences( *this);
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
      // prepare everything: read in a sequence to test
      io::IFStream read;

      // read seq_a from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      util::ShPtr< biol::AASequence> seq( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // create alignment and alignment_list
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment( new align::AlignmentLeaf< biol::AABase>( seq));

      align::HandlerPIR< biol::AABase> handler;
      handler.WriteAlignment( util::GetLogger(), *alignment);

      // create assignment score
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( double( 1.0) * score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -0.5, -0.5, -0.5, -0.5
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // use default constructor
      align::WordGeneratorSubsequences< biol::AABase> generator;

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // generate words
      BCL_MessageStd( "2: Generate()")
      const size_t word_length( 3);
      const size_t exp_words_size( alignment->GetSize() - word_length + 1);
      storage::Vector< align::AlignmentWord< biol::AABase> > words( generator.Generate( alignment, word_length));

      // print query words only when message level at least verbose
      if( util::GetMessenger().GetCurrentMessageLevel() > util::Message::e_Standard)
      {
        WriteAllQueryWords( words, assign_score);
      }

      BCL_ExampleIndirectCheck
      (
        words.GetSize(),
        exp_words_size,
        "2: Size==" + util::Format()( words.GetSize()) + "!=" + util::Format()( exp_words_size)
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    //! @brief write a single word
    void WriteAllQueryWords
    (
      const storage::Vector< align::AlignmentWord< biol::AABase> > &WORDS,
      const score::AssignmentWithGap< biol::AABase> &SCORE
    ) const
    {
      for
      (
        storage::Vector< align::AlignmentWord< biol::AABase> >::const_iterator
          itr( WORDS.Begin()),
          itr_end( WORDS.End());
        itr != itr_end;
        ++itr
      )
      {
        util::GetLogger() << "word=" << itr->ToString() << "|score=" << itr->ScoreSelf( SCORE) << '\n';
      }
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignWordGeneratorSubsequences

  const ExampleClass::EnumType ExampleAlignWordGeneratorSubsequences::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignWordGeneratorSubsequences())
  );
  
} // namespace bcl
