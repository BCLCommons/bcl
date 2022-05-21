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
#include "score/bcl_score_protein_model_gap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelGap::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelGap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from default values
    //! @param SCHEME scheme of the score
    ProteinModelGap::ProteinModelGap( const std::string &SCHEME) :
      m_Scheme( SCHEME),
      m_Score()
    {
    }

    //! @brief returns a pointer to a new ProteinModelGap
    //! @return pointer to a new ProteinModelGap
    ProteinModelGap *ProteinModelGap::Clone() const
    {
      return new ProteinModelGap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &ProteinModelGap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the scheme of this score
    //! @return the scheme of this score
    const std::string &ProteinModelGap::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief returns the default scheme of this score
    //! @return the default scheme of this score
    const std::string &ProteinModelGap::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "gap");
      return s_default_scheme;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProteinModelGap::GetAlias() const
    {
      static const std::string s_name( "ProteinModelGap");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelGap::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Scoring function to penalize gaps in protein models.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Score a protein model regarding gaps.
    //! @param PROTEIN_MODEL protein model for which to compute the gap score
    //! @return gap score of the protein model
    double ProteinModelGap::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // sum up gap score for all chains in the model
      double score( 0.0);
      const util::SiPtrVector< const biol::AASequence> chains( PROTEIN_MODEL.GetSequences());
      for( auto chain_it( chains.Begin()), chain_it_end( chains.End()); chain_it != chain_it_end; ++chain_it)
      {
        score += ComputeGapScore( **chain_it);
      }

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Compute the gap score for the given sequence.
    //! @param SEQUENCE Sequence for which to compute the gap score.
    //! @return Gap score for the given sequence.
    double ProteinModelGap::ComputeGapScore( const biol::AASequence &SEQUENCE) const
    {
      double score( 0.0);
      const util::ShPtrVector< biol::AABase> aas( SEQUENCE.GetData());
      util::ShPtr< biol::AABase> def_n;
      util::ShPtr< biol::AABase> def_c;
      for( auto res_it( aas.Begin()), res_it_end( aas.End()); res_it != res_it_end; ++res_it)
      {
        const biol::AABase &current_aa( **res_it);
        if( current_aa.HasDefinedCoordinates())
        {
          def_n = *res_it;
        }
        else if( def_n.IsDefined() && def_n->GetSeqID() < current_aa.GetSeqID())
        {
          for( auto res_find_it( res_it + 1); res_find_it != res_it_end; ++res_find_it)
          {
            const util::ShPtr< biol::AABase> sp_c_cand( *res_find_it);
            if( sp_c_cand->HasDefinedCoordinates())
            {
              def_c = *res_find_it;
              res_it = res_find_it - 1;
              const size_t sequence_distance( def_c->GetSeqID() - def_n->GetSeqID() - 1);
              const double euclidean_distance( ( def_c->GetCA().GetCoordinates() - def_n->GetCA().GetCoordinates()).Norm());
              const storage::Pair< size_t, double> gap_parameters( sequence_distance, euclidean_distance);
              // score += m_Score.Score( gap_parameters);
              score += ComputeGapScore( *def_n, *def_c);
              break;
            }
          }
        }

      }

      return score;
    }

    //! @brief Compute the gap score for the given endpoints of the gap.
    //! @param AA_N n-terminal endpoint of the gap
    //! @param AA_C c-terminal endpoint of the gap
    //! @return Gap score for the given sequence.
    double ProteinModelGap::ComputeGapScore( const biol::AABase &AA_N, const biol::AABase &AA_C) const
    {
      const size_t seq_dist( AA_C.GetSeqID() - AA_N.GetSeqID() - 1);
      const double euclidean_dist( ( AA_C.GetCenter() - AA_N.GetCenter()).Norm());
      const double score( euclidean_dist / seq_dist);

      return score;
    }

  } // namespace score
} // namespace bcl
