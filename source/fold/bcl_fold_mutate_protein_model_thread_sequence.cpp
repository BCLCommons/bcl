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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "fold/bcl_fold_mutate_protein_model_thread_sequence.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_compare.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelThreadSequence::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelThreadSequence())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelThreadSequence::MutateProteinModelThreadSequence() :
      m_Alignment(),
      m_PrintStructuralProblemResidues(),
      m_Scheme()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param CHAIN_ALIGNMENT for each chain, the alignment that should be used
    //! @param SCHEME the scheme of the mutate
    MutateProteinModelThreadSequence::MutateProteinModelThreadSequence
    (
      const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &CHAIN_ALIGNMENT,
      const std::string &SCHEME,
      const bool PRINT_STRUCTURAL_PROBLEM_AAS
    ) :
      m_Alignment( CHAIN_ALIGNMENT),
      m_PrintStructuralProblemResidues( PRINT_STRUCTURAL_PROBLEM_AAS),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelThreadSequence
    MutateProteinModelThreadSequence *MutateProteinModelThreadSequence::Clone() const
    {
      return new MutateProteinModelThreadSequence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelThreadSequence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateProteinModelThreadSequence::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel>
    MutateProteinModelThreadSequence::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // copy the model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // iterate through the alignments for each chain
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
          align_itr( m_Alignment.Begin()), align_itr_end( m_Alignment.End());
        align_itr != align_itr_end;
        ++align_itr
      )
      {
        // the current chain the current alignment corresponds to
        const char chain_id( align_itr->first);

        // get the chain out of the protein model
        const util::ShPtr< assemble::Chain> current_chain( new_model->GetChain( chain_id));

        // true if the chain does not exist in the protein model
        if( !current_chain.IsDefined())
        {
          // go to next alignment and chain
          continue;
        }

        // thread the sequence onto the chain according to the current alignment
        const util::ShPtr< assemble::Chain> new_chain
        (
          ThreadSequenceOntoChain( *current_chain, *align_itr->second)
        );

        // replace the new chain in the protein model
        const util::ShPtr< assemble::ProteinModel> current_model( UpdateChainInModel( *new_model, new_chain));

        // set the new model to the current model with the current chain just threaded
        new_model = current_model;

        if( m_PrintStructuralProblemResidues)
        {
          util::GetLogger() << OutputStructurallyProbematicSequenceRegions( *current_chain, *align_itr->second)
            << "\n";
        }
      }

      // construct mutate result
      const math::MutateResult< assemble::ProteinModel> mutate_result( new_model, *this);

      // return mutate result
      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelThreadSequence::Read( std::istream &ISTREAM)
    {
      // read members
      BCL_Exit( "cannot read Alignments", 1);
      io::Serialize::Read( m_Alignment, ISTREAM);
      io::Serialize::Read( m_PrintStructuralProblemResidues, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelThreadSequence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Alignment, OSTREAM, INDENT);
      io::Serialize::Write( m_PrintStructuralProblemResidues, OSTREAM, INDENT);
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief uses a sequence alignment to assign coordinates from a chain onto a sequence with unknown structure
    //! @param CHAIN the coordinates that will be assigned to sequence of unknown structure
    //! @param ALIGNMENT the sequence alignment indicating which residues to assign which coordinates to
    //! @return shptr to chain which has the second sequence from the alignment but coordinates from the CHAIN param
    util::ShPtr< assemble::Chain> MutateProteinModelThreadSequence::ThreadSequenceOntoChain
    (
      const assemble::Chain &CHAIN, const align::AlignmentInterface< biol::AABase> &ALIGNMENT
    )
    {
      // make sure the alignment has two sequences
      BCL_Assert
      (
        ALIGNMENT.GetDepth() == 2, "the alignment needs to have 2 sequences but has "
        + util::Format()( ALIGNMENT.GetDepth())
      );

      // get the current chain sequence id
      const std::string chain_seq_id( CHAIN.GetSequence()->GetSequenceId());

      // get the sequence ids of the alignment
      const storage::List< std::string> alignment_seq_ids( ALIGNMENT.GetSequenceIds());

      // known coordinates for iterating through
      const biol::AASequence &coord_seq( *CHAIN.GetSequence());
      biol::AASequence::const_iterator seq_itr( coord_seq.Begin()), seq_itr_end( coord_seq.End());

      // will hold the sequence from the second in the alignment but coordinates from the CHAIN
      util::ShPtr< biol::AASequence> new_sequence( new biol::AASequence());

      // the assignments from the alignment
      const util::ShPtrList< align::Assignment< biol::AABase> > &assignments( ALIGNMENT.GetAssignments());

      // iterate through the assignments
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
          assign_itr( assignments.Begin()), assign_itr_end( assignments.End());
        assign_itr != assign_itr_end && seq_itr != seq_itr_end;
        ++assign_itr
      )
      {
        // make sure the current assignment is defined and get a reference to it
        BCL_Assert( assign_itr->IsDefined(), "assignment is not defined");
        const align::Assignment< biol::AABase> &current_assign( **assign_itr);

        // get the members from the assignment
        const util::SiPtrList< const biol::AABase> &members( current_assign.GetMembers());

        // should only be two members in the assignment
        BCL_Assert( members.GetSize() == 2, "members size is not 2 but " + util::Format()( members.GetSize()));

        // sequence with known coordinates (from CHAIN) is the first sequence
        const util::SiPtr< const biol::AABase> &member_start_seq( members.FirstElement());

        // sequence with unknown coordinates but is the desired sequence is the second sequence
        const util::SiPtr< const biol::AABase> &member_new_seq( members.LastElement());

        // true if both residues in the alignment are defined
        if( member_start_seq.IsDefined() && member_new_seq.IsDefined())
        {
          // make new aa data from the new sequence
          const util::ShPtr< biol::AAData> new_data( member_new_seq->GetData()->Clone());

          // make new residue from the new data
          util::ShPtr< biol::AABase> new_aa( ( *( *seq_itr)->GetAAClass())->Empty( new_data));

          //  assert matching aa types from the sequence of known coordinates and from the chain and alignment
          BCL_Assert
          (
            member_start_seq->GetType() == ( *seq_itr)->GetType(),
            "types differ for : template sequence from alignment : " + member_start_seq->GetIdentification() +
            " versus coordinate sequence from pdb :" + ( *seq_itr)->GetIdentification()
          );

          BCL_MessageDbg( "adding new residue " + new_aa->GetIdentification());

          // get the atoms from the sequence of known structure
          const util::SiPtrVector< const biol::Atom> &atoms( ( *seq_itr)->GetAtoms());

          // set the atoms and coordinates of the new residue
          new_aa->SetAtoms( atoms);

          // add the residue to the new sequence
          new_sequence->PushBack( new_aa);

          // go to next residue with known coordinates
          ++seq_itr;
        }
        // aa of desired sequence is not defined (gap in desired sequence)
        else if( member_start_seq.IsDefined() && !member_new_seq.IsDefined())
        {
          // go to next residue with known coordinates
          ++seq_itr;

          // go to next residues
          continue;
        }
        // aa of known coordinates is not defined (gap in sequence of known coordinates)
        else if( !member_start_seq.IsDefined() && member_new_seq.IsDefined())
        {
          // make new aa data
          const util::ShPtr< biol::AAData> new_data( member_new_seq->GetData()->Clone());

          // make new residue - won't have any coordinates
          const util::ShPtr< biol::AABase> new_aa( ( *( *seq_itr)->GetAAClass())->Empty( new_data));

          BCL_MessageDbg( "adding new residue " + new_aa->GetIdentification());

          // add new residue to the new sequence
          new_sequence->PushBack( new_aa);
        }
        else
        {
          BCL_Exit( "both the starting sequence and desired sequence are undefined in alignment", 1);
        }
      }

      // set the chain of the new sequence
      new_sequence->SetChainID( coord_seq.GetChainID());

      BCL_MessageDbg( "SetChainID to |" + util::Format()( new_sequence->GetChainID()) + "|");

      // make vector of sses with the entire new sequence as a coil sse
      const util::ShPtrVector< assemble::SSE> sses
      (
        1, util::ShPtr< assemble::SSE>( new assemble::SSE( *new_sequence, biol::GetSSTypes().COIL))
      );

      // create a new chain
      const util::ShPtr< assemble::Chain> new_chain( new assemble::Chain( new_sequence, sses));

      // return the new chain
      return new_chain;
    }

    //! @brief replaces a chain in a model
    //! @param MODEL the model in which the chain will be replaced
    //! @param CHAIN the chain which will replace existing chain in the provided model
    //! @return shptr to a new protein model with the desired chain replaced
    util::ShPtr< assemble::ProteinModel>
    MutateProteinModelThreadSequence::UpdateChainInModel
    (
      const assemble::ProteinModel &MODEL, const util::ShPtr< assemble::Chain> &CHAIN
    )
    {
      // make new model
      util::ShPtr< assemble::ProteinModel> current_model( new assemble::ProteinModel());

      // set the data of the model
      current_model->SetProteinModelData( MODEL.GetProteinModelData());

      // insert the chains into the current model
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( MODEL.GetChains().Begin()), chain_itr_end( MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // true if the current chain is not the one that should be replaced
        if( ( *chain_itr)->GetChainID() != CHAIN->GetChainID())
        {
          BCL_MessageDbg( "inserting chain " + util::Format()( ( *chain_itr)->GetChainID()));
          current_model->Insert( *chain_itr);
        }
        else //< need to insert the new chain
        {
          current_model->Insert( CHAIN);
        }
      }

      // return the updated model
      return current_model;
    }

    //! @brief uses a sequence alignment to and chain to determine which stretches of unknown structure need rebuilt
    //!        This could be due e.g. to not having coordinates from the template structure, or if there are gaps in
    //!        the sequence of unknown structure (since this will likely lead to a chain break).
    //! @param CHAIN the coordinates that will be assigned to sequence of unknown structure
    //! @param ALIGNMENT the sequence alignment indicating which residues to assign which coordinates to
    //! @return string which has the regions of the sequence of unknown structure which most likely need rebuilt
    std::string MutateProteinModelThreadSequence::OutputStructurallyProbematicSequenceRegions
    (
      const assemble::Chain &CHAIN, const align::AlignmentInterface< biol::AABase> &ALIGNMENT
    )
    {
      // make sure the alignment has two sequences
      BCL_Assert
      (
        ALIGNMENT.GetDepth() == 2, "the alignment needs to have 2 sequences but has "
        + util::Format()( ALIGNMENT.GetDepth())
      );
      // known coordinates for iterating through
      const biol::AASequence &coord_seq( *CHAIN.GetSequence());
      biol::AASequence::const_iterator seq_itr( coord_seq.Begin()), seq_itr_end( coord_seq.End());

      // the assignments from the alignment
      const util::ShPtrList< align::Assignment< biol::AABase> > &assignments( ALIGNMENT.GetAssignments());

      storage::Set< util::SiPtr< const biol::AABase>, biol::AALessThanSeqID> problem_resis;

      util::SiPtr< const biol::AABase> last_defined_new_seq;

      // iterate through the assignments
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
          assign_itr( assignments.Begin()), assign_itr_end( assignments.End());
        assign_itr != assign_itr_end && seq_itr != seq_itr_end;
        ++assign_itr
      )
      {
        // make sure the current assignment is defined and get a reference to it
        BCL_Assert( assign_itr->IsDefined(), "assignment is not defined");
        const align::Assignment< biol::AABase> &current_assign( **assign_itr);

        // get the members from the assignment
        const util::SiPtrList< const biol::AABase> &members( current_assign.GetMembers());

        // should only be two members in the assignment
        BCL_Assert( members.GetSize() == 2, "members size is not 2 but " + util::Format()( members.GetSize()));

        // sequence with known coordinates (from CHAIN) is the first sequence
        const util::SiPtr< const biol::AABase> &member_start_seq( members.FirstElement());

        // sequence with unknown coordinates but is the desired sequence is the second sequence
        const util::SiPtr< const biol::AABase> &member_new_seq( members.LastElement());

        if( member_new_seq.IsDefined())
        {
          last_defined_new_seq = member_new_seq;
        }

        // true if both residues in the alignment are defined
        if( member_start_seq.IsDefined() && member_new_seq.IsDefined())
        {
          // go to next residue with known coordinates
          ++seq_itr;
        }
        // aa of desired sequence is not defined (gap in desired sequence)
        else if( member_start_seq.IsDefined() && !member_new_seq.IsDefined())
        {
          if( last_defined_new_seq.IsDefined())
          {
//              BCL_Assert
//              (
//                problem_resis.Insert( last_defined_new_seq).second, "could not insert " + last_defined_new_seq->GetIdentification()
//              );

            problem_resis.Insert( last_defined_new_seq);
          }

          // go to next residue with known coordinates
          ++seq_itr;

          // go to next residues
          continue;
        }
        // aa of known coordinates is not defined (gap in sequence of known coordinates)
        else if( !member_start_seq.IsDefined() && member_new_seq.IsDefined())
        {
          BCL_Assert( member_new_seq.IsDefined(), "resi is not defined");
          BCL_Assert
          (
            problem_resis.Insert( member_new_seq).second, "could not insert " + member_new_seq->GetIdentification()
          );
        }
        else
        {
          BCL_Exit( "both the starting sequence and desired sequence are undefined in alignment", 1);
        }
      }

      std::string problem_residues;
      util::SiPtr< const biol::AABase> region_start_aa( *problem_resis.Begin());
      util::SiPtr< const biol::AABase> region_previous_aa( region_start_aa);
      int previous_seq_id( region_start_aa->GetSeqID() - 1);
      for
      (
        storage::Set< util::SiPtr< const biol::AABase>, biol::AALessThanSeqID>::const_iterator
          resi_itr( problem_resis.Begin()), resi_itr_end( problem_resis.End());
        resi_itr != resi_itr_end;
        ++resi_itr
      )
      {
        if( ( *resi_itr)->GetSeqID() == previous_seq_id + 1)
        {
          previous_seq_id = ( *resi_itr)->GetSeqID();
          region_previous_aa = ( *resi_itr);
        }
        else
        {
          BCL_Assert( resi_itr->IsDefined(), "resi is not defined");
          problem_residues += region_start_aa->GetIdentification() + " to " + region_previous_aa->GetIdentification() + "\n";
          region_start_aa = ( *resi_itr);
          region_previous_aa = ( *resi_itr);
          previous_seq_id = ( *resi_itr)->GetSeqID();
        }
      }

      problem_residues += region_start_aa->GetIdentification() + " to " + region_previous_aa->GetIdentification() + "\n";

      return problem_residues;
    }

    //! @brief determines if a residue has any defined coordinates
    //! @param RESI the residue that will be checked to see if it has any defined coordinates
    //! @return boolean true if there are any defined atom coordinates in RESI - false otherwise
    bool MutateProteinModelThreadSequence::HasAnyDefinedCoordinates( const biol::AABase &RESI)
    {
      const util::SiPtrVector< const linal::Vector3D> &coords( RESI.GetAtomCoordinates());

      bool has_any_defined_coords( false);
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( coords.Begin()), coord_itr_end( coords.End());
        coord_itr != coord_itr_end;
        ++coord_itr
      )
      {
        if( coord_itr->IsDefined())
        {
          has_any_defined_coords = true;
        }
      }

      return has_any_defined_coords;
    }

  } // namespace fold
} // namespace bcl
