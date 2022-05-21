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
#include "fold/bcl_fold_mutate_protein_model_sse_pair_fix_loop_closure.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"
#include "score/bcl_score_loop_closure.h"

// external includes - sorted alphabetically

// Define this to visualize all clash removals
//#define DEBUG_BCL_FOLD_MUTATE_LOOP_CLOSURE_RESOLVER

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEPairFixLoopClosure::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSEPairFixLoopClosure())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEPairFixLoopClosure::MutateProteinModelSSEPairFixLoopClosure( const size_t &MAX_SSES_TO_MOVE) :
      m_Scheme
      (
        "directed_loop_closure_"
        + ( MAX_SSES_TO_MOVE == util::GetUndefinedSize_t() ? std::string( "all") : util::Format()( MAX_SSES_TO_MOVE))
      ),
      m_MaxSSEsToMove( MAX_SSES_TO_MOVE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSEPairFixLoopClosure
    MutateProteinModelSSEPairFixLoopClosure *MutateProteinModelSSEPairFixLoopClosure::Clone() const
    {
      return new MutateProteinModelSSEPairFixLoopClosure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelSSEPairFixLoopClosure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the weighted average distance between two sses based on the length of the coil between them
    double MutateProteinModelSSEPairFixLoopClosure::GetTypicalLoopLength( const size_t &N_RESIDUES) const
    {
      // R^2 = 0.986
      return 2.445 + 4.84 * std::log( double( N_RESIDUES));
    }

    //! @brief get the weighted average distance between two sses based on the length of the coil between them
    double MutateProteinModelSSEPairFixLoopClosure::GetMaxLoopLength( const size_t &N_RESIDUES) const
    {
      return 2.1136 + 2.5609 * double( N_RESIDUES);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEPairFixLoopClosure::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      static score::LoopClosure s_closure( size_t( 1), 0.0, 1.0);
      double prev_score( s_closure( PROTEIN_MODEL));

      // make copy of proteinmodel
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      static size_t n_times_called( 0);

      ++n_times_called;
      if( prev_score < 0.05)
      {
        return math::MutateResult< assemble::ProteinModel>( new_model, *this);
      }

#ifdef DEBUG_BCL_FOLD_MUTATE_LOOP_CLOSURE_RESOLVER
      io::OFStream out;
      const std::string folder( "/tmp/loop_itr_" + util::Format()( n_times_called) + "/");
      pdb::Factory factory;
      io::Directory::MkDir( folder);
      io::File::MustOpenOFStream( out, folder + "start.pdb");
      factory.WriteModelToPDB( PROTEIN_MODEL, out);
      io::File::CloseClearFStream( out);
#endif
      // Find all sses that violate the loop closure score on either side
      for
      (
        auto chain_itr( new_model->GetChains().Begin()), chain_itr_end( new_model->GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        assemble::Chain &chain( **chain_itr);
        // need at least two sses in the chain
        if( chain.GetNumberSSEs() < 2)
        {
          continue;
        }
        {
          auto sse_itr( chain.GetData().Begin()), sse_itr_b( ++chain.GetData().Begin());
          const storage::Pair< size_t, double> seq_euclid_dist_first
          (
            score::LoopClosure::SequenceAndEuclideanDistanceWithExclusion( **sse_itr, **sse_itr_b, 1)
          );
          // for the first two SSEs, check whether the N-terminal SSE is too far from the second.
          if( util::IsDefined( seq_euclid_dist_first.First()) && util::IsDefined( seq_euclid_dist_first.Second()))
          {
            const double max_length( GetMaxLoopLength( seq_euclid_dist_first.First()));
            if( seq_euclid_dist_first.Second() > max_length)
            {
              const assemble::SSE &sseb( **sse_itr_b);
              // too far, so move N-terminus to fix
              util::ShPtr< assemble::SSE> new_sse( ( *sse_itr)->Clone());
              const double typical_length( GetTypicalLoopLength( seq_euclid_dist_first.First()));
              new_sse->Translate( sseb.BeginOfZ() + linal::Vector3D().SetRandomTranslation( typical_length) - new_sse->EndOfZ());
              chain.Replace( new_sse);
            }
          }
        }
        if( chain.GetNumberSSEs() == size_t( 2))
        {
          continue;
        }
        {
          auto sse_itr_last( --chain.GetData().End()), sse_itr_pre_last( ----chain.GetData().End());
          const storage::Pair< size_t, double> seq_euclid_dist_last
          (
            score::LoopClosure::SequenceAndEuclideanDistanceWithExclusion( **sse_itr_pre_last, **sse_itr_last, 1)
          );
          // for the first two SSEs, check whether the C-terminal SSE is too far from the previous.
          if( util::IsDefined( seq_euclid_dist_last.First()) && util::IsDefined( seq_euclid_dist_last.Second()))
          {
            const double max_length( GetMaxLoopLength( seq_euclid_dist_last.First()));
            if( seq_euclid_dist_last.Second() > max_length)
            {
              const assemble::SSE &ssea( **sse_itr_pre_last);
              // too far, so move C-terminus to fix
              util::ShPtr< assemble::SSE> new_sse( ( *sse_itr_last)->Clone());
              const double typical_length( GetTypicalLoopLength( seq_euclid_dist_last.First()));
              new_sse->Translate( ssea.EndOfZ() + linal::Vector3D().SetRandomTranslation( typical_length) - new_sse->BeginOfZ());
              chain.Replace( new_sse);
            }
          }
        }
        if( chain.GetNumberSSEs() == size_t( 3))
        {
          continue;
        }

        // iterate through the sses of "CHAIN" and look for chains whose loops are too far away
        bool did_something( true);
        while( did_something)
        {
          did_something = false;
          for
          (
            auto sse_itr_b( ++chain.GetData().Begin()),
              sse_itr_c( ++++chain.GetData().Begin()),
              sse_itr_last( --chain.GetData().End());
            sse_itr_b != sse_itr_last;
            sse_itr_b = sse_itr_c++
          )
          {
            auto sse_itr( sse_itr_b);
            --sse_itr;
            const storage::Pair< size_t, double> seq_euclid_dist_left
            (
              score::LoopClosure::SequenceAndEuclideanDistanceWithExclusion( **sse_itr, **sse_itr_b, 1)
            );

            // loop is undefined
            if( !util::IsDefined( seq_euclid_dist_left.First()) || !util::IsDefined( seq_euclid_dist_left.Second()))
            {
              continue;
            }

            const storage::Pair< size_t, double> seq_euclid_dist_right
            (
              score::LoopClosure::SequenceAndEuclideanDistanceWithExclusion( **sse_itr_b, **sse_itr_c, 1)
            );

            // loop is undefined
            if( !util::IsDefined( seq_euclid_dist_right.First()) || !util::IsDefined( seq_euclid_dist_right.Second()))
            {
              continue;
            }

            const double max_left( GetMaxLoopLength( seq_euclid_dist_left.First()));
            const double max_right( GetMaxLoopLength( seq_euclid_dist_right.First()));
            if( seq_euclid_dist_left.Second() <= max_left && seq_euclid_dist_right.Second() <= max_right)
            {
              continue;
            }
            util::ShPtr< assemble::SSE> new_sse( ( *sse_itr_b)->Clone());
            bool is_fixed( false);

            const assemble::SSE &ssea( **sse_itr);
            const assemble::SSE &sseb( *new_sse);
            const assemble::SSE &ssec( **sse_itr_c);
            const double pref_left( GetTypicalLoopLength( seq_euclid_dist_left.First()));
            const double pref_right( GetTypicalLoopLength( seq_euclid_dist_right.First()));
            const coord::LineSegment3D original_line( sseb.GetMainAxis());
            for( size_t i( 0); i < size_t( 10); ++i)
            {
              const math::TransformationMatrix3D trans
              (
                sseb.GetMainAxis(),
                coord::LineSegment3D
                (
                  ssea.EndOfZ() + linal::Vector3D().SetRandomTranslation( random::GetGlobalRandom().Boolean() ? pref_left : max_left),
                  ssec.BeginOfZ() + linal::Vector3D().SetRandomTranslation( random::GetGlobalRandom().Boolean() ? pref_right : max_right)
                )
              );
              new_sse->Transform( trans);
              const storage::Pair< size_t, double> seq_euclid_dist_left_updated
              (
                score::LoopClosure::SequenceAndEuclideanDistanceWithExclusion( **sse_itr, *new_sse, 1)
              );
              const storage::Pair< size_t, double> seq_euclid_dist_right_updated
              (
                score::LoopClosure::SequenceAndEuclideanDistanceWithExclusion( *new_sse, **sse_itr_c, 1)
              );
              if( seq_euclid_dist_right_updated.Second() <= max_right && seq_euclid_dist_left_updated.Second() <= max_left)
              {
                is_fixed = true;
                BCL_MessageVrb( "Resolved loop_closure on " + util::Format()( i) + "th iteration");
              }
            }
            if( !is_fixed)
            {
              const math::TransformationMatrix3D trans
              (
                sseb.GetMainAxis(),
                original_line
              );
              new_sse->Transform( trans);
            }
            else
            {
              did_something = true;
              chain.Replace( new_sse);
            }
          }
        }
      }
      new_model->GetChangeSignal().Emit( *new_model);

#ifdef DEBUG_BCL_FOLD_MUTATE_LOOP_CLOSURE_RESOLVER
      io::File::MustOpenOFStream( out, folder + "end.pdb");
      factory.WriteModelToPDB( *new_model, out);
      io::File::CloseClearFStream( out);
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
    std::istream &MutateProteinModelSSEPairFixLoopClosure::Read( std::istream &ISTREAM)
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
    std::ostream &MutateProteinModelSSEPairFixLoopClosure::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
