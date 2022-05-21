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
#include "fold/bcl_fold_mutate_protein_model_sse_pair_clash.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"

// external includes - sorted alphabetically

// Define this to visualize all clash removals
//#define DEBUG_BCL_FOLD_MUTATE_CLASH_RESOLVER

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEPairClash::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSEPairClash())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEPairClash::MutateProteinModelSSEPairClash( const size_t &MAX_SSES_TO_MOVE) :
      m_Scheme
      (
        "sse_clash_remover_"
        + ( MAX_SSES_TO_MOVE == util::GetUndefinedSize_t() ? std::string( "all") : util::Format()( MAX_SSES_TO_MOVE))
      ),
      m_MaxSSEsToMove( MAX_SSES_TO_MOVE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSEPairClash
    MutateProteinModelSSEPairClash *MutateProteinModelSSEPairClash::Clone() const
    {
      return new MutateProteinModelSSEPairClash( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelSSEPairClash::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEPairClash::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      static score::AAPairHiResClash s_clash;
      static assemble::VoxelGridAA s_voxel_grid( s_clash.GetDistanceCutoff());
      double prev_clash_score( s_clash( PROTEIN_MODEL));

      static size_t n_times_called( 0);

      ++n_times_called;
      if( prev_clash_score < 0.05)
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      storage::Vector< storage::Pair< util::SiPtr< const assemble::SSE>, linal::Vector3D> >
      sse_moves
      (
        s_voxel_grid.GetMinSSEMoveIDsToRemoveClashes
        (
          PROTEIN_MODEL.GetSSEs(),
          PROTEIN_MODEL.GetAminoAcids(),
          false
        )
      );

      // if no pair is available
      if( sse_moves.IsEmpty())
      {
        // warn user
        BCL_MessageVrb( "unable to find a pair of clashing SSEs!");

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // make copy of proteinmodel
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      sse_moves.Shuffle();

      if( m_MaxSSEsToMove < sse_moves.GetSize())
      {
        sse_moves.Resize( m_MaxSSEsToMove);
      }
      // iterate over clashing SSEs
      util::ShPtr< assemble::ProteinModel> new_model_b;
      const size_t max_rounds( 20);
      size_t rnd( 0);
      double total_move_dist( 0.0);
      size_t n_sses_moved( 0);

#ifdef DEBUG_BCL_FOLD_MUTATE_CLASH_RESOLVER
      io::OFStream out;
      const std::string folder( "/tmp/clash_itr_" + util::Format()( n_times_called) + "/");
      pdb::Factory factory;
      io::Directory::MkDir( folder);
      io::File::MustOpenOFStream( out, folder + "itr0.pdb");
      factory.WriteModelToPDB( PROTEIN_MODEL, out);
      io::File::CloseClearFStream( out);
#endif
      while( !sse_moves.IsEmpty() && ++rnd <= max_rounds)
      {
        for( auto itr( sse_moves.Begin()), itr_end( sse_moves.End()); itr != itr_end; ++itr)
        {
          // calculate the sse_packing
          BCL_MessageVrb( "mutated sse " + itr->First()->GetIdentification() + " " + util::Format()( itr->Second().ToString()));
          util::ShPtr< assemble::SSE> mutated_sse( itr->First()->Clone());

          // apply transformation
          total_move_dist += itr->Second().Norm();
          mutated_sse->Translate( -itr->Second());

          // replace mutated sse
          new_model->Replace( mutated_sse);
          ++n_sses_moved;
#ifdef DEBUG_BCL_FOLD_MUTATE_CLASH_RESOLVER
          io::File::MustOpenOFStream( out, folder + "itr" + util::Format()( rnd) + "_sse" + util::Format()( n_sses_moved) + ".pdb");
          factory.WriteModelToPDB( *new_model, out);
          io::File::CloseClearFStream( out);
#endif
        }
        new_model_b = util::ShPtr< assemble::ProteinModel>( new_model->Clone());
        sse_moves = s_voxel_grid.GetMinSSEMoveIDsToRemoveClashes
                    (
                      new_model_b->GetSSEs(),
                      new_model_b->GetAminoAcids(),
                      false
                    );
      }
      #ifdef DEBUG_BCL_FOLD_MUTATE_CLASH_RESOLVER
      const double new_clash_score( s_clash( *new_model));
      if( sse_moves.GetSize() && new_clash_score > 0.5)
      {
        BCL_MessageVrb
        (
          "Iteration " + util::Format()( n_times_called) +
          " Gave up removing clashes... remaining clash score "
          + util::Format()( prev_clash_score)
          + ". Remaining moves: " + util::Format()( sse_moves.GetSize())
          + " initial score: + " + util::Format()( prev_clash_score)
          + " final score: " + util::Format()( new_clash_score)
          + " # moved "
          + util::Format()( n_sses_moved) + " distance moved: " + util::Format()( total_move_dist)
        );
      }
      else
      {
        BCL_MessageStd
        (
          "Iteration " + util::Format()( n_times_called) +
          " Fixed clashes in " + util::Format()( rnd) + " rounds of clash removal. # moved "
          + util::Format()( n_sses_moved) + " distance moved: " + util::Format()( total_move_dist)
          + " remaining score: " + util::Format()( new_clash_score)
        );
      }
      #endif
      // return
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSSEPairClash::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelSSEPairClash::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
