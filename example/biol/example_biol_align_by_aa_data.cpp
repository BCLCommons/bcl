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
#include "biol/bcl_biol_align_by_aa_data.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "align/bcl_align_handler_pir.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_align_by_aa_data.cpp
  //!
  //! @author karakam
  //! @date Jan 13, 2011
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAlignByAAData :
    public ExampleInterface
  {
  public:

    ExampleBiolAlignByAAData *Clone() const
    {
      return new ExampleBiolAlignByAAData( *this);
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
      // initialize
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      const std::string pdb_filename_a( AddExampleInputPathToFilename( e_Biology, "1PY6.pdb"));
      const std::string pdb_filename_b( AddExampleInputPathToFilename( e_Biology, "2KSY.pdb"));

      // read in pairwise alignment
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1PY6_2KSY_A.pir"));
      const align::HandlerPIR< biol::AABase> handler;
      util::ShPtr< align::AlignmentNode< biol::AABase> > sp_alignment( handler.ReadAlignment( read, biol::AASequence()));
      io::File::CloseClearFStream( read);

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> sse_min_size;
      sse_min_size[ biol::GetSSTypes().HELIX] = 9;
      sse_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, sse_min_size));
      assemble::ProteinModel model_a( Proteins::GetModel( pdb_filename_a, biol::GetAAClasses().e_AABackBone, sse_min_size));
      assemble::ProteinModel model_b( Proteins::GetModel( pdb_filename_b, biol::GetAAClasses().e_AABackBone, sse_min_size));

      // lengths
      const size_t nr_aas( 76);
      const size_t nr_aas_in_sses( 42);

      // alignment from the SSEs
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > align_sses
      (
        assemble::Quality::CreateAlignmentFromProteinModelSSEs( model)
      );
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > align_sses_a
      (
        assemble::Quality::CreateAlignmentFromProteinModelSSEs( model_a)
      );
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > align_sses_b
      (
        assemble::Quality::CreateAlignmentFromProteinModelSSEs( model_b)
      );

      // alignment from the sequences
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > align_seq
      (
        assemble::Quality::CreateAlignmentFromProteinModelSequences( model)
      );

      // check the size
      BCL_ExampleCheck( align_sses[ 'A']->GetSize(), nr_aas_in_sses);

      // check the size
      BCL_ExampleCheck( align_seq[ 'A']->GetSize(), nr_aas);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      biol::AlignByAAData alignment_engine;

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // do the alignment
      const storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment_score_pair
      (
        alignment_engine.AlignPair( align_seq[ 'A'], align_sses[ 'A'])
      );

      // output the alignment
      BCL_MessageStd( "Alignment of SSE and sequence residues by AAData");
      handler.WriteAlignment( util::GetLogger(), alignment_score_pair.First());

      // check the size
      BCL_ExampleCheck( alignment_score_pair.First().GetSize(), nr_aas);

      // now remove the gaps
      align::AlignmentNode< biol::AABase> alignment_no_gaps
      (
        assemble::Quality::RemoveAssignmentsWithGapsFromAlignment( alignment_score_pair.First())
      );

      // output the alignment
      BCL_MessageStd( "Alignment of SSE and sequence residues by AAData with gaps removed");
      handler.WriteAlignment( util::GetLogger(), alignment_no_gaps);

      // check the size
      BCL_ExampleCheck( alignment_no_gaps.GetSize(), nr_aas_in_sses);

      // align two sequences using a third, template alignment
      const align::AlignmentNode< biol::AABase> final_alignment
      (
        alignment_engine.AlignPair( align_sses_a[ 'A'], align_sses_b[ 'A'], sp_alignment).First()
      );
      BCL_MessageStd( "Alignment of two sequences using a third, template alignment");
      handler.WriteAlignment( util::GetLogger(), final_alignment);

      // check the size
      BCL_ExampleCheck( final_alignment.GetSize(), sp_alignment->GetSize());

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAlignByAAData

  const ExampleClass::EnumType ExampleBiolAlignByAAData::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAlignByAAData())
  );

} // namespace bcl
