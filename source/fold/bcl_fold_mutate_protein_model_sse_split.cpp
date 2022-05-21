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
#include "fold/bcl_fold_mutate_protein_model_sse_split.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_placement_sse_short_loop.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! @brief the default scheme for this class
    const std::string &MutateProteinModelSSESplit::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "mutate_split_sse");
      return s_default_scheme;
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSESplit::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new MutateProteinModelSSESplit
        (
          sspred::GetMethods().e_Undefined,
          util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> >(),
          storage::Map< biol::SSType, size_t>(),
          math::Range< size_t>( 0, 0)
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from locator and mutate
    //! @param SS_METHOD the method to use to locate the largest drop in the ss prediction for the located sse
    //! @param SP_LOCATE_SSE locator of sse from domain
    //! @param MIN_SSE_SIZES min sse sizes for the resultant SSEs
    //! @param SPLIT_COIL_LENGTH_RANGE range of length of the coil to be inserted after the split
    //! @param SCHEME the scheme
    MutateProteinModelSSESplit::MutateProteinModelSSESplit
    (
      const sspred::Method &SS_METHOD,
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &SP_LOCATE_SSE,
      const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES,
      const math::Range< size_t> &SPLIT_COIL_LENGTH_RANGE,
      const std::string &SCHEME
    ) :
      m_Method( SS_METHOD),
      m_SSELocator( SP_LOCATE_SSE),
      m_MinSSESizes( MIN_SSE_SIZES),
      m_SplitCoilLengthRange( SPLIT_COIL_LENGTH_RANGE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSESplit
    MutateProteinModelSSESplit *MutateProteinModelSSESplit::Clone() const
    {
      return new MutateProteinModelSSESplit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSESplit::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateProteinModelSSESplit::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that mutates a single SSE in a given ProteinModel by splitting a single sse
    //! @param PROTEIN_MODEL ProteinModel from which a single SSE will be splitted
    //! @return MutateResult that results from mutating to the PROTEIN_MODEL
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSESplit::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize static undefined pointer to ProteinModel
      static const util::ShPtr< assemble::ProteinModel> sp_undefined_model;

      // locate a sse from the pool
      const util::SiPtr< const assemble::SSE> located_sse( m_SSELocator->Locate( PROTEIN_MODEL));

      // was locating successful
      if( !located_sse.IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( sp_undefined_model, *this);
      }

      // length of sse to split
      const size_t seq_length( located_sse->GetSize());

      // length of coil
      const size_t coil_length( random::GetGlobalRandom().Random( m_SplitCoilLengthRange));

      // find the min SSE size, set it to 0 if not containted in the min sse sizes map for that type
      storage::Map< biol::SSType, size_t>::const_iterator sse_size_itr( m_MinSSESizes.Find( located_sse->GetType()));
      const size_t min_size( sse_size_itr == m_MinSSESizes.End() ? 0 : sse_size_itr->second);

      // need enough residues to split
      if( seq_length < coil_length + 2 * min_size)
      {
        return math::MutateResult< assemble::ProteinModel>( sp_undefined_model, *this);
      }

      // initialize left and right coil lengths
      const size_t left_coil_length( coil_length / 2);
      const size_t right_coil_length( coil_length - left_coil_length);

      // find the position to split
      const util::SiPtr< const biol::AABase> sp_aa
      (
        FindSplitPosition( *located_sse, m_Method, min_size + left_coil_length, min_size + right_coil_length)
      );

      // no drop split pos could be identified
      if( !sp_aa.IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( sp_undefined_model, *this);
      }

      // determine the split positions for the left and the right SSE
      const size_t split_pos_left
      (
        sp_aa->GetSeqID() - located_sse->GetFirstAA()->GetSeqID() - left_coil_length
      );
      const size_t split_pos_right
      (
        sp_aa->GetSeqID() - located_sse->GetFirstAA()->GetSeqID() + right_coil_length
      );

      // construct two individual SSEs
      util::ShPtr< assemble::SSE> sp_sse_left
      (
        new assemble::SSE( located_sse->SubSequence( 0, split_pos_left), located_sse->GetType())
      );
      util::ShPtr< assemble::SSE> sp_sse_right
      (
        new assemble::SSE
        (
          located_sse->SubSequence( split_pos_right, seq_length - split_pos_right), located_sse->GetType()
        )
      );

      // decide on which one to keep intact in the same place as the original SSE
      bool keep_left_intact( random::GetGlobalRandom().Boolean());

      // if loop of 0 then always add to top, otherwise more chance to add side by side the longer the loop is
      const double add_to_top_probability
      (
        coil_length == 0 ? 1.0 : std::min( 1.0, 2.0 / coil_length)
      );
      // construct the placement object
      const PlacementSSEShortLoop placer( coil_length, add_to_top_probability);

      // determine which SSE to keep intact and which one to move
      const util::ShPtr< assemble::SSE> sp_sse_intact( keep_left_intact ? sp_sse_left : sp_sse_right);
      util::ShPtr< assemble::SSE> sp_sse_placed( keep_left_intact ? sp_sse_right : sp_sse_left);

      // calculate the placement
      const storage::Pair< math::TransformationMatrix3D, bool> placement( placer.Place( *sp_sse_placed, *sp_sse_intact));

      BCL_MessageVrb
      (
        "split sse: " + located_sse->GetIdentification() + " at : " +
        sp_aa->GetIdentification() +
        " into " + sp_sse_left->GetIdentification() + " and " + sp_sse_right->GetIdentification()
      );

      // make sure it is valid
      if( !placement.Second())
      {
        BCL_MessageStd
        (
          "Unable to calculate a placement for the split SSEs for SSE: " + located_sse->GetIdentification()
        );
      }
      else
      {
        // construct transformation initialized with the inverse that will move the selected sse to origin
        math::TransformationMatrix3D trans( sp_sse_placed->GetOrientation());
        trans.Invert();
        trans( placement.First());

        // apply the transformation to the split SSE to be moved
        sp_sse_placed->Transform( trans);
      }

      // clone the model
      util::ShPtr< assemble::ProteinModel> sp_new_model( PROTEIN_MODEL.Clone());

      // replace with overlapping for the moved SSE
      BCL_Assert
      (
        sp_new_model->ReplaceWithOverlapping( sp_sse_placed),
        "Replacement of the following SSE with overlapping one failed: " + sp_sse_placed->GetIdentification()
      );

      // insert the intact SSE and assert insertion
      BCL_Assert
      (
        sp_new_model->Insert( sp_sse_intact), "Insertion of following SSE failed: " + sp_sse_intact->GetIdentification()
      );

      // end
      return math::MutateResult< assemble::ProteinModel>( sp_new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSSESplit::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Method         , ISTREAM);
      io::Serialize::Read( m_SSELocator     , ISTREAM);
      io::Serialize::Read( m_MinSSESizes    , ISTREAM);
      io::Serialize::Read( m_SplitCoilLengthRange, ISTREAM);
      io::Serialize::Read( m_Scheme         , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelSSESplit::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Method         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSELocator     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinSSESizes    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SplitCoilLengthRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme         , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief find the amino acid right to the highest difference for the given sspred::Method and SSType and report the difference
    //! @param ELEMENT the secondary structure element
    //! @param SSMETHOD the ssprediction method
    //! @param NR_RESIDUES_TO_EXCLUDE number of residues to exclude from both sides
    //! @return Pair of difference and the amino acid to the left
    storage::Pair< util::SiPtr< const biol::AABase>, double>
    MutateProteinModelSSESplit::FindLargestSSPredDifferencePosition
    (
      const assemble::SSE &ELEMENT,
      const sspred::Method &SSMETHOD,
      const size_t NR_RESIDUES_TO_EXCLUDE
    )
    {
      // initialize the pair
      storage::Pair< util::SiPtr< const biol::AABase>, double> aa_diff( util::SiPtr< const biol::AABase>(), 0.0);

      // sequence too short
      if( ELEMENT.GetSize() < 2 * NR_RESIDUES_TO_EXCLUDE)
      {
        return aa_diff;
      }

      // sstype
      const biol::SSType sstype( ELEMENT.GetType());

      // iterator for amino acids
      assemble::SSE::const_iterator itr( ELEMENT.Begin() + NR_RESIDUES_TO_EXCLUDE);
      const assemble::SSE::const_iterator itr_end( ELEMENT.End() - NR_RESIDUES_TO_EXCLUDE);

      // first aa
      aa_diff.First() = *itr;
      double left_prediction( ( *itr)->GetSSPrediction( SSMETHOD)->GetThreeStatePrediction()( sstype));

      // next aa
      ++itr;

      // iterate over all aas
      while( itr != itr_end)
      {
        const double prediction( ( *itr)->GetSSPrediction( SSMETHOD)->GetThreeStatePrediction()( sstype));
        const double difference( prediction - left_prediction);

        // update the highest difference if necessary
        if( math::Absolute( difference) > math::Absolute( aa_diff.Second()))
        {
          aa_diff.First() = *itr;
          aa_diff.Second() = difference;
        }

        // go to next
        left_prediction = prediction;
        ++itr;
      }

      // if largest difference is 0, do not return a position
      if( aa_diff.Second() == 0.0)
      {
        aa_diff.First() = util::SiPtr< const biol::AABase>();
      }

      // end
      return aa_diff;
    }

    //! @brief find the amino acid right to the highest difference for the given sspred::Method and SSType and report the difference
    //! @param ELEMENT the secondary structure element
    //! @param SSMETHOD the ssprediction method
    //! @param NR_RESIDUES_TO_EXCLUDE_LEFT number of residues to exclude from the left side
    //! @param NR_RESIDUES_TO_EXCLUDE_RIGHT number of residues to exclude from the right side
    //! @return Pair of difference and the amino acid to the left
    util::SiPtr< const biol::AABase>
    MutateProteinModelSSESplit::FindSplitPosition
    (
      const assemble::SSE &ELEMENT,
      const sspred::Method &SSMETHOD,
      const size_t NR_RESIDUES_TO_EXCLUDE_LEFT,
      const size_t NR_RESIDUES_TO_EXCLUDE_RIGHT
    )
    {
      // initialize return value
      util::SiPtr< const biol::AABase> sp_aa;

      // sequence too short
      if( ELEMENT.GetSize() < NR_RESIDUES_TO_EXCLUDE_LEFT + NR_RESIDUES_TO_EXCLUDE_RIGHT)
      {
        return sp_aa;
      }

      // sstype
      const biol::SSType sstype( ELEMENT.GetType());

      // iterator for amino acids
      assemble::SSE::const_iterator itr( ELEMENT.Begin() + NR_RESIDUES_TO_EXCLUDE_LEFT);
      const assemble::SSE::const_iterator itr_end( ELEMENT.End() - NR_RESIDUES_TO_EXCLUDE_RIGHT);

      math::RunningAverageSD< double> stats_mean_sd;

      // iterate over all amino acid
      for( ; itr != itr_end; ++itr)
      {
        // make sure there is a prediction
        if( !( *itr)->GetSSPrediction( SSMETHOD).IsDefined())
        {
          return sp_aa;
        }

        const double this_prob( ( *itr)->GetSSPrediction( SSMETHOD)->GetThreeStatePrediction()( sstype));
        stats_mean_sd += this_prob;
      }

      // calculate the cutoff
      const double cutoff( stats_mean_sd.GetAverage() - stats_mean_sd.GetStandardDeviation());

      // store sum
      double sum( 0.0);

      // iterate again
      for( itr = ELEMENT.Begin() + NR_RESIDUES_TO_EXCLUDE_LEFT; itr != itr_end; ++itr)
      {
        const double this_prob( ( *itr)->GetSSPrediction( SSMETHOD)->GetThreeStatePrediction()( sstype));

        // if larger than cutoff skip
        if( this_prob > cutoff)
        {
          continue;
        }

        // sum
        sum += cutoff - this_prob;
      }

      // get a random number [0,sum]
      const double random( random::GetGlobalRandom().Random( sum));

      // initialize running sum
      double running_sum( 0.0);

      // re-initialize iterator
      itr = ELEMENT.Begin() + NR_RESIDUES_TO_EXCLUDE_LEFT;

      while( itr != itr_end)
      {
        // get sspred probability
        const double this_prob( ( *itr)->GetSSPrediction( SSMETHOD)->GetThreeStatePrediction()( sstype));

        // if larger than cutoff skip
        if( this_prob > cutoff)
        {
          ++itr;
          continue;
        }

        // update running sum
        running_sum += cutoff - this_prob;

        // if in correct sum interval identified by the random number, then break
        if( running_sum >= random)
        {
          sp_aa = *itr;
          break;
        }
        ++itr;
      }

      // end
      return sp_aa;
    }

  } // namespace fold
} // namespace bcl
