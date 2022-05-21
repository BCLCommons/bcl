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
#include "align/bcl_align_aligner_merge.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_aligner_merge.cpp
  //!
  //! @author heinzes1
  //! @date 2010/11/01
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignerMerge :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignerMerge *Clone() const
    {
      return new ExampleAlignAlignerMerge( *this);
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

      // read seq from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      util::ShPtr< biol::AASequence> seq
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read))
      );
      io::File::CloseClearFStream( read);

      // create alignment and alignment_list
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment( new align::AlignmentLeaf< biol::AABase>( seq));
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignment_list( 3, alignment);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      align::AlignerMerge< biol::AABase> aligner;

      // test Clone()
      util::ShPtr< align::AlignerMerge< biol::AABase> > aligner_clone( aligner.Clone());

    /////////////////
    // data access //
    /////////////////

      // no need to test SetScoringFunction, its empty

    ////////////////
    // operations //
    ////////////////

      // test AlignPair()
      align::AlignmentNode< biol::AABase> result_alignment( aligner.AlignPair( alignment, alignment).First());
      BCL_MessageStd( "1: aligner.AlignPair()");
      BCL_ExampleIndirectCheck
      (
        result_alignment.GetDepth(),
        2,
        "1a: result_alignment.GetDepth==" + util::Format()( result_alignment.GetDepth()) + "!=2"
      );
      BCL_ExampleIndirectCheck
      (
        result_alignment.GetSize(),
        alignment->GetSize(),
        "1b: result_alignment.GetSize==" + util::Format()( result_alignment.GetSize()) + "!="
          + util::Format()( alignment->GetSize())
      );

      // test AlignMultiple() and cloned aligner
      result_alignment = aligner_clone->AlignMultiple( alignment_list).First();
      BCL_MessageStd( "2: aligner_clone.AlignMultiple()");
      BCL_ExampleIndirectCheck
      (
        result_alignment.GetDepth(),
        3,
        "2a: result_alignment.GetDepth==" + util::Format()( result_alignment.GetDepth()) + "!=3"
      );
      BCL_ExampleIndirectCheck
      (
        result_alignment.GetSize(),
        alignment->GetSize(),
        "2b: result_alignment.GetSize==" + util::Format()( result_alignment.GetSize()) + "!="
          + util::Format()( alignment->GetSize())
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignEngineMerge

  const ExampleClass::EnumType ExampleAlignAlignerMerge::s_Instance( GetExamples().AddEnum( ExampleAlignAlignerMerge()));
  
} // namespace bcl
