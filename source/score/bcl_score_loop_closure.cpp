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
#include "score/bcl_score_loop_closure.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "math/bcl_math_linear_function.h"
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
    const util::SiPtr< const util::ObjectInterface> LoopClosure::s_Instance
    (
      util::Enumerated< math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double> >::AddInstance
      (
        new LoopClosure()
      )
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief get the name of the object
    //! @return the name of the object
    const std::string &LoopClosure::GetAlias() const
    {
      static const std::string s_name( "LoopClosure");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopClosure::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores distance to form loops between two SSEs");
      serializer.AddInitializer
      (
        "number excluded residues",
        "number of excluded residues from each end of loop",
        io::Serialization::GetAgent( &m_NrExcludedResidues)
       );
      serializer.AddInitializer
      (
        "sigmoid width",
        "width of the sigmoidal function to be used",
        io::Serialization::GetAgent( &m_SigmoidWidth)
      );
      serializer.AddInitializer
      (
        "fraction allowed distance",
        "fraction of allowed distance to use",
        io::Serialization::GetAgent( &m_FractionAllowedDistance)
      );
      serializer.AddInitializer
      (
        "exclude coil",
        "whether or not coil sses should be used in scoring",
        io::Serialization::GetAgent( &m_ExcludeCoil)
       );

      return serializer;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopClosure::LoopClosure() :
      m_NrExcludedResidues(),
      m_SigmoidWidth(),
      m_FractionAllowedDistance(),
      m_ExcludeCoil()
    {
    }

    //! @brief constructor from a sigmoid width
    //! @param NR_EXCLUDED_RESIDUES Number of excluded residues from each end of loop
    //! @param SIGMOID_WIDTH width of sigmoid function
    //! @param FRACTION_ALLOWED_DISTANCE the allowed distance is multiplied with the fraction before calculating the penalty
    //! @param EXCLUDE_COIL true if coil sses should not be considered in scoring
    LoopClosure::LoopClosure
    (
      const size_t NR_EXCLUDED_RESIDUES,
      const double SIGMOID_WIDTH,
      const double FRACTION_ALLOWED_DISTANCE,
      const bool EXCLUDE_COIL
    ) :
      m_NrExcludedResidues( NR_EXCLUDED_RESIDUES),
      m_SigmoidWidth( SIGMOID_WIDTH),
      m_FractionAllowedDistance( FRACTION_ALLOWED_DISTANCE),
      m_ExcludeCoil( EXCLUDE_COIL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LoopClosure
    LoopClosure *LoopClosure::Clone() const
    {
      return new LoopClosure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LoopClosure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate Sequence distance between aas and Euclidean distance between last and first amino acids of SSEs
    //! if the NR_EXCLUDED_RESIDUES is >= size of sse, it will be reduced, so that at least 1 amino acid remains
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @param NR_EXCLUDED_RESIDUES number of residues to be excluded from each end
    //! @return a pair which first member is the length of the loop and aminoacids and the second is the Euclidean distance
    storage::Pair< size_t, double> LoopClosure::SequenceAndEuclideanDistanceWithExclusion
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B,
      const size_t NR_EXCLUDED_RESIDUES
    )
    {
      // sequences of different chains
      if( SSE_A.GetChainID() != SSE_B.GetChainID())
      {
        BCL_MessageVrb( "passing two different chains");
        return storage::Pair< size_t, double>( util::GetUndefined< size_t>(), util::GetUndefined< double>());
      }

      // check if given sses are not empty
      if( SSE_A.GetSize() == 0 || SSE_B.GetSize() == 0)
      {
        BCL_MessageVrb
        (
          "the passed SSEs are empty, returning undefined" +
          SSE_A.GetIdentification() + " and " + SSE_B.GetIdentification()
        );
        return storage::Pair< size_t, double>( util::GetUndefined< size_t>(), util::GetUndefined< double>());
      }

      // order pair, so that first's Sequence first AA has the smaller Sequence ID
      storage::VectorND< 2, util::SiPtr< const assemble::SSE> > ordered_pair( SSE_A, SSE_B);

      // if the first SSE comes after the second SSE
      if( ordered_pair.First()->GetFirstAA()->GetSeqID() > ordered_pair.Second()->GetFirstAA()->GetSeqID())
      {
        std::swap( ordered_pair.First(), ordered_pair.Second());
      }

      const size_t nr_excluded_res1( std::min( NR_EXCLUDED_RESIDUES, ordered_pair.First()->GetSize() - 1));
      const size_t nr_excluded_res2( std::min( NR_EXCLUDED_RESIDUES, ordered_pair.Second()->GetSize() - 1));

      // make sure the SSEs are long enough
      if( nr_excluded_res1 != NR_EXCLUDED_RESIDUES || nr_excluded_res2 != NR_EXCLUDED_RESIDUES)
      {
        BCL_MessageVrb
        (
          "the SSEs are too short to cut by " + util::Format()( NR_EXCLUDED_RESIDUES) +
          " into " + SSE_A.GetIdentification() + " and " + SSE_B.GetIdentification() +
          " less residues are used to cut in"
        );
      }

      // construct the pair to hold sequence and euclidean distances
      storage::Pair< size_t, double> seq_euc_distance;

      // calculate the sequence distance
      seq_euc_distance.First() = biol::CalculateSequenceDistance( *ordered_pair.Second(), *ordered_pair.First())
        + nr_excluded_res1 + nr_excluded_res2;

      // now calculate the Euclidean distance
      seq_euc_distance.Second() = biol::GetPeptideBondLength
        (
          *ordered_pair.First()->GetAA( ordered_pair.First()->GetSize() - nr_excluded_res1 - 1),
          *ordered_pair.Second()->GetAA( nr_excluded_res2)
        );

      BCL_MessageDbg
      (
        "seq_euc_distance is " + util::Format()( seq_euc_distance) + " for " +
        ordered_pair.First()->GetIdentification() + " and " + ordered_pair.Second()->GetIdentification()
      );

      // end
      return seq_euc_distance;
    }

    //! @brief score the loop
    //! @param SEQ_EUC_DISTANCE pair of sequence and Euclidean distance
    //! @return the loop closure score
    double LoopClosure::ScoreLoop( const storage::Pair< size_t, double> &SEQ_EUC_DISTANCE) const
    {
      // initialize static linear function
      static const double s_slope( 2.5609);
      static const double s_intercept( 2.1136);
      static const math::LinearFunction s_linear_function( s_slope, s_intercept);

      // calculate the different to maximum allowed distance
      const double difference( SEQ_EUC_DISTANCE.Second() - m_FractionAllowedDistance * s_linear_function( SEQ_EUC_DISTANCE.First()));

      // if within the range
      if( difference <= 0.0)
      {
        return 0.0;
      }
      // if larger than the sigmoid width
      else if( difference >= m_SigmoidWidth)
      {
        if( !m_SigmoidWidth)
        {
          return difference;
        }
        return 1.0;
      }
      else
      {
        // apply the sigmoidal function to get a value between 0 and sigmoid width
        return math::WeightBetweenZeroAndPi( ( ( m_SigmoidWidth - difference) / m_SigmoidWidth) * math::g_Pi);
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief gives a penalty for two sses based on the difference between the observed distance and the distance the
    //! missing amino acids could bridge
    //! @param SSE_A the first sse
    //! @param SSE_B the second sse
    //! @return the score for two sses based on the possibility to close the loop between SSE_A and SSE_B
    double LoopClosure::operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      // if coils are excluded, return without scoring
      if( m_ExcludeCoil && ( SSE_A.GetType() == biol::GetSSTypes().COIL || SSE_B.GetType() == biol::GetSSTypes().COIL))
      {
        return 0.0;
      }

      // sum the score
      const storage::Pair< size_t, double> seq_euclid_dist( SequenceAndEuclideanDistanceWithExclusion( SSE_A, SSE_B, m_NrExcludedResidues));

      // loop is undefined
      if( !util::IsDefined( seq_euclid_dist.First()) || !util::IsDefined( seq_euclid_dist.Second()))
      {
        return 0.0;
      }

      const double score( ScoreLoop( seq_euclid_dist));
      if( !util::IsDefined( score))
      {
        return 0.0;
      }

      // return score
      return score;
    }

    //! @brief operator that scores the chain
    //! @param CHAIN the chain for which all neighbor scores are calculated
    //! @return score
    double LoopClosure::operator()( const assemble::Chain &CHAIN) const
    {
      // sum of all scores
      double score( 0.0);

      // need at least two sses in the chain
      if( CHAIN.GetNumberSSEs() < 2)
      {
        return score;
      }

      // iterate through the sses of "CHAIN" to calculate their score sum
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr_a( CHAIN.GetData().Begin()),
          sse_itr_b( ++CHAIN.GetData().Begin()),
          sse_itr_end( CHAIN.GetData().End());
        sse_itr_b != sse_itr_end;
        ++sse_itr_a, ++sse_itr_b
      )
      {
        score += this->operator ()( **sse_itr_a, **sse_itr_b);
      }

      // return the score sum
      return score;
    }

    //! @brief operator that scores the Protein model
    //! @param PROTEIN_MODEL the protein model for which all neighbor scores are calculated
    //! @return score
    double LoopClosure::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // score sum
      double score( 0.0);

      // iterate through the chains of PROTEIN_MODEL
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        score += operator()( **chain_itr);
      }

      // return score
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the Scheme and the Function value for the ARGUMENT to the STREAM
    //! @param SSE_A the first sse
    //! @param SSE_B the second sse
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    LoopClosure::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE_A, const assemble::SSE &SSE_B,
      std::ostream &OSTREAM
    ) const
    {
      // form the sse pair and calculate the sequence euclidean distances
      const storage::Pair< size_t, double> seq_euc_distance
      (
        SequenceAndEuclideanDistanceWithExclusion( SSE_A, SSE_B, m_NrExcludedResidues)
      );

      // write the detailed values and score
      OSTREAM << SSE_A.GetIdentification()    << '\t' << SSE_B.GetIdentification() << '\t'
              << seq_euc_distance.First()     << '\t' << seq_euc_distance.Second() << '\t'
              << ScoreLoop( seq_euc_distance) << '\n';

      // end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
