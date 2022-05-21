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
#include "pdb/bcl_pdb_factory.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_factory.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbFactory :
    public ExampleInterface
  {
  public:

    ExamplePdbFactory *Clone() const
    { return new ExamplePdbFactory( *this);}

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
      //read pdb file
      io::IFStream read;
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      BCL_MessageStd( "read pdb file: " + pdb_filename);
      BCL_ExampleMustOpenInputFile( read, pdb_filename);

      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

    //////////////
    // backbone //
    //////////////

      pdb::Factory factory( biol::GetAAClasses().e_AABackBone);

      //instantiate sequences
      BCL_MessageStd( "building sequences from pdb chains using backbone");
      util::ShPtrVector< biol::AASequence> sequences( factory.AASequencesFromPDB( pdb));

      //number of sequences
      BCL_MessageStd( "pdb has " + util::Format()( sequences.GetSize()) + " chains.");

      //all chains with their identifiers and the number of residues
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator chain_itr( sequences.Begin()),
          chain_itr_end( sequences.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        BCL_MessageStd
        (
          "chain with identifier '" + util::Format()( ( *chain_itr)->GetChainID()) +
          "' has " + util::Format()( ( *chain_itr)->GetSize()) + " residues."
        );

        BCL_MessageStd( "sequence for this chain is: " + ( *chain_itr)->Sequence());
      }

      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel protein_model( factory.ProteinModelFromPDB( pdb, ssetype_min_size));

      //write models which have only the missing loops and only the BackBone coordinates
      io::OFStream write;
      BCL_MessageStd( "write model_withoutloop_nonideal_backbone.pdb");
      std::string outfile( AddExampleOutputPathToFilename( factory, "model_withoutloop_nonideal_backbone.pdb"));
      BCL_ExampleMustOpenOutputFile( write, outfile);
      factory.WriteModelToPDB( protein_model, write);
      io::File::CloseClearFStream( write);

      // output of all secondary structure element types and starting and ending residues
      // of each model and number of helices and strands
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( protein_model.GetChains().Begin()),
          chain_itr_end( protein_model.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        BCL_MessageStd
        (
          "this is the number of helices and strands in the chain " + util::Format()( ( *chain_itr)->GetChainID()) +
          ": " + util::Format()( ( *chain_itr)->GetNumberSSE( biol::GetSSTypes().HELIX)) + " " +
          util::Format()( ( *chain_itr)->GetNumberSSE( biol::GetSSTypes().STRAND))
        );

        BCL_MessageStd
        (
          "these are beginning and ending amino acid seqids according to the secondary structure elements"
        );
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          BCL_MessageStd( util::Format()( ( *sse_itr)->GetType()));

          BCL_MessageStd
          (
            util::Format()( ( *sse_itr)->GetFirstAA()->GetSeqID()) + " " +
            util::Format()( ( *sse_itr)->GetFirstAA()->GetSeqID())
          );
        }
      }

    //////////
    // cacb //
    //////////

      pdb::Factory factory_cacb( biol::GetAAClasses().e_AACaCb);

      //instantiate sequences
      BCL_MessageStd( "building sequences from pdb chains using cacb");
      util::ShPtrVector< biol::AASequence> sequences_cacb( factory_cacb.AASequencesFromPDB( pdb));

      //number of sequences
      BCL_MessageStd( "pdb has " + util::Format()( sequences_cacb.GetSize()) + " chains.");

      //all chains with their identifiers and the number of residues
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator chain_itr( sequences_cacb.Begin()),
          chain_itr_end( sequences_cacb.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        BCL_MessageStd
        (
          "chain with identifier '" + util::Format()( ( *chain_itr)->GetChainID()) +
            "' has " + util::Format()( ( *chain_itr)->GetSize()) + " residues."
        );

        BCL_MessageStd( "sequence for this chain is: " + ( *chain_itr)->Sequence());
      }

      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      assemble::ProteinModel protein_model_cacb( factory_cacb.ProteinModelFromPDB( pdb, ssetype_min_size));

      //write models which have only the missing loops and only the BackBone coordinates
      BCL_MessageStd( "write model_withoutloop_nonideal_cacb.pdb");
      outfile = AddExampleOutputPathToFilename( factory, "model_withoutloop_nonideal_cacb.pdb");
      BCL_ExampleMustOpenOutputFile( write, outfile);
      factory_cacb.WriteModelToPDB( protein_model_cacb, write);
      io::File::CloseClearFStream( write);

      // output of all secondary structure element types and starting and ending residues
      // of each model and number of helices and strands
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( protein_model_cacb.GetChains().Begin()),
          chain_itr_end( protein_model_cacb.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        BCL_MessageStd
        (
          "this is the number of helices and strands in the chain " + util::Format()( ( *chain_itr)->GetChainID()) +
          ": " + util::Format()( ( *chain_itr)->GetNumberSSE( biol::GetSSTypes().HELIX)) + " " +
          util::Format()( ( *chain_itr)->GetNumberSSE( biol::GetSSTypes().STRAND))
        );

        BCL_MessageStd
        (
          "these are beginning and ending amino acid seqids according to the secondary structure elements"
        );
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          BCL_MessageStd( util::Format()( ( *sse_itr)->GetType()));

          BCL_MessageStd
          (
            util::Format()( ( *sse_itr)->GetFirstAA()->GetSeqID()) + " " +
            util::Format()( ( *sse_itr)->GetFirstAA()->GetSeqID())
          );
        }
      }

      return 0;

    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbFactory

  const ExampleClass::EnumType ExamplePdbFactory::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbFactory())
  );

} // namespace bcl
