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
#include "score/bcl_score_data_set_pairwise_bipolar.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "math/bcl_math_statistics.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseBipolar::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseBipolar())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseBipolar::GetDefaultScheme()
    {
      static const std::string s_scheme( "bipolar");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DataSetPairwiseBipolar::DataSetPairwiseBipolar( const std::string &SCHEME) :
      m_Data(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking pool
    //! @param SSE_POOL the pool the sse definitions will come from
    //! @param SCHEME optional scheme
    DataSetPairwiseBipolar::DataSetPairwiseBipolar( const assemble::SSEPool &SSE_POOL, const std::string &SCHEME) :
      m_Data(),
      m_Scheme( SCHEME)
    {
      m_Data = CalculateResiduePositionWeight( SSE_POOL);
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseBipolar
    DataSetPairwiseBipolar *DataSetPairwiseBipolar::Clone() const
    {
      return new DataSetPairwiseBipolar( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseBipolar::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseBipolar::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

      //! @brief calculate the score of a data set
      //! @param DATA data set to be scored
      //! @return the score of the current data set
      double DataSetPairwiseBipolar::operator()( const restraint::DataSetPairwise &DATA) const
      {
        if( DATA.IsEmpty())
        {
          return 0;
        }

        double score( 0);

        // to keep track of the number of restraint data points that are in each sse
        storage::Map< util::ShPtr< assemble::SSE>, storage::Pair< size_t, double> > sse_restraintcount;

        storage::Map
        <
          util::ShPtr< assemble::SSE>,
          storage::Map< util::ShPtr< assemble::SSE>, storage::Vector< double> >
        > sse_connectionweights;

        // iterate over the data pairs to score them
        for
        (
          restraint::DataSetPairwise::const_iterator itr( DATA.Begin()), itr_end( DATA.End());
          itr != itr_end;
          ++itr
        )
        {
          storage::Map
          <
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
            storage::Pair< double, util::ShPtr< assemble::SSE> >,
            assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
          >::const_iterator resi_itr_a( m_Data.Find( itr->First()));

          storage::Map
          <
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
            storage::Pair< double, util::ShPtr< assemble::SSE> >,
            assemble::LocatorAtomCoordinatesInterface::PtrLessThan
          >::const_iterator resi_itr_b( m_Data.Find( itr->Second()));

          // if either residue could not found in the map
          if( resi_itr_a == m_Data.End() || resi_itr_b == m_Data.End())
          {
            ++score;
            continue;
          }

          // get the weights
          double weight_a( resi_itr_a->second.First());
          double weight_b( resi_itr_b->second.First());

          // optimal score here is zero, since that would mean the two residues are on the same position along their sse
          // relative to the membrane
          double current_score( math::Absolute( weight_a - weight_b));

          // want the residues to be away from the middle of their sses
          current_score += ( 1 - math::Absolute( weight_a - 0.5) - math::Absolute( weight_b - 0.5));

          // score goes from 0 to one at this point with 0 being best
          current_score /= double( 2);

          // add score to total score
          score += current_score;

          // add connection information
          sse_connectionweights[ resi_itr_a->second.Second()][ resi_itr_b->second.Second()].PushBack( weight_a + weight_b);
          sse_connectionweights[ resi_itr_b->second.Second()][ resi_itr_a->second.Second()].PushBack( weight_a + weight_b);

          // add count for the sses the two residues are in
          { // residue a
            storage::Map< util::ShPtr< assemble::SSE>, storage::Pair< size_t, double> >::const_iterator sse_itr
            (
              sse_restraintcount.Find( resi_itr_a->second.Second())
            );
            if( sse_itr != sse_restraintcount.End())
            {
              ++sse_restraintcount[ resi_itr_a->second.Second()].First();
              sse_restraintcount[ resi_itr_a->second.Second()].Second() += weight_a;
            }
            else
            {
              sse_restraintcount[ resi_itr_a->second.Second()].First() = 1;
              sse_restraintcount[ resi_itr_a->second.Second()].Second() = weight_a;
            }
          }
          { // residue b
            storage::Map< util::ShPtr< assemble::SSE>, storage::Pair< size_t, double> >::const_iterator sse_itr
            (
              sse_restraintcount.Find( resi_itr_b->second.Second())
            );
            if( sse_itr != sse_restraintcount.End())
            {
              ++sse_restraintcount[ resi_itr_b->second.Second()].First();
              sse_restraintcount[ resi_itr_b->second.Second()].Second() += weight_b;
            }
            else
            {
              sse_restraintcount[ resi_itr_b->second.Second()].First() = 1;
              sse_restraintcount[ resi_itr_b->second.Second()].Second() = weight_b;
            }
          }
        }

        score /= double( DATA.GetSize());

        // subtract one so best score will be -1 taking into account residues being towards end of sses and being on
        // the same side of the membrane. Does not take into account wanting residues within an sse to be on opposite
        // sides
        --score;

        double opposite_ends_score( 0);
        // iterate over the map of sses and the counts of datapoints and total residue weight in each sse
        for
        (
          storage::Map< util::ShPtr< assemble::SSE>, storage::Pair< size_t, double> >::const_iterator
            sse_itr( sse_restraintcount.Begin()), sse_itr_end( sse_restraintcount.End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          // count of data points in the current sse
          const size_t count( sse_itr->second.First());

          // sum of the residue weights of data points that are in the sse
          const double weight( sse_itr->second.Second());

          // calculate average weight
          const double average_weight( weight / double( count));

          // since weights of residues within sses go from 0 to 1, and we want residues to be on opposite sides of the
          // sse, optimal average weight is 0.5 i.e. every other data point has weight zero and every other data point
          // has weight 1
          // multiply by two so that it goes between 0 and 1
          opposite_ends_score += ( math::Absolute( average_weight - 0.5));
        }
        // this makes it go from 0 to 1 with 0 being best and 1 being worst
        opposite_ends_score = 2 * opposite_ends_score / double( sse_restraintcount.GetSize());

        // subtract one so score goes from -1 to 0
        --opposite_ends_score;

        // iterate over the connection information
        double total_connection_score( 0);
        size_t pair_count( 0);
        for
        (
          storage::Map
          <
            util::ShPtr< assemble::SSE>,
            storage::Map< util::ShPtr< assemble::SSE>, storage::Vector< double> >
          >::const_iterator sse_itr( sse_connectionweights.Begin()), sse_itr_end( sse_connectionweights.End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          const storage::Map< util::ShPtr< assemble::SSE>, storage::Vector< double> > &map_b( sse_itr->second);
          for
          (
            storage::Map< util::ShPtr< assemble::SSE>, storage::Vector< double> >::const_iterator
              sse_itr_b( map_b.Begin()), sse_itr_b_end( map_b.End());
            sse_itr_b != sse_itr_b_end;
            ++sse_itr_b
          )
          {
            const storage::Vector< double> &weights( sse_itr_b->second);
            const double weight_sum( math::Statistics::Sum( weights.Begin(), weights.End()));

            // divide the weightsum by the number of residues in it (i.e. 2 * the number of restraints)
            // optimal score is 0.5 so subtract the average restraint weight from 0.5
            // 0.5 is optimal since that would mean if you have two restraints involving four residues,
            // they would be on opposite sides of the sse (0 - 0, 1 - 1; 0 + 0 + 1 + 1 = 2; 2 / 2*2 = 0.5)
            const double current_score( math::Absolute( ( weight_sum / ( 2.0 * double( weights.GetSize()))) - 0.5));
            total_connection_score += current_score;
            ++pair_count;
          }
        }
        // since everything was entered in the map symmetrically by sse (i.e. twice) the scores are effectively
        // multiplied by 2, making them go from 0 to 1 with zero being best
        // divide by half the pair count and it will go between 0 and 1
        total_connection_score /= ( double( pair_count) / 2.0);
        // subtract 1 so it goes between -1 and 0
        --total_connection_score;

        // average scores
        score = ( 4 * score + opposite_ends_score + total_connection_score) / 6.0;

        return score;
      }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseBipolar::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseBipolar::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates the weight of every residue in helical sses according to its position relative to the ends
    //!        The residue weights go from 0 to 1 along an sse. The end of the sse that has a 0 and the end that has
    //!        a 1 doesn't matter, as long as it is consistent for all sses (which it is). The weights are calculated
    //!        such that the ends of sses that are on the same side of the membrane have the same weight, and the
    //!        weight of residues along sses is the same for each sse as it goes from one side of the membrane to the
    //!        other.
    //! @param SSE_POOL the pool of sses that will be used to calculate the weights
    //! @return map that has for every residue in helical sses the weight associated according to its position relative
    //!         to the end of the sse
    storage::Map
    <
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
      assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
    > DataSetPairwiseBipolar::CalculateResiduePositionWeight( const assemble::SSEPool &SSE_POOL)
    {
      // the sses from the pool
      storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sses
      (
        SSE_POOL.GetRandomNonOverlappingSet()
      );

      // map to hold the residues and their weights
      storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
        assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
      > residue_weights;

      // keep track of if the nterminus should be weighted as 0 or 1, since it is assumed the sses alternatively
      // go through the membrane every other n-terminus of an sse will have a 0, and every other will be weight 1.
      size_t n_terminus( 0);

      // iterate through the sses
      for
      (
        storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( sses.Begin()), sse_itr_end( sses.End());
          sse_itr != sse_itr_end;
          ++sse_itr
      )
      {
        // true if the sse is not helix, go to next sse
        if( ( *sse_itr)->GetType() != biol::GetSSTypes().HELIX)
        {
          continue;
        }

        // for keeping track of what sse each aa is in
        const util::ShPtr< assemble::SSE> sse_copy( ( *sse_itr)->Clone());

        // iterate over the residues of the sse to assign their weights
        for
        (
          biol::AASequence::const_iterator aa_itr( ( *sse_itr)->Begin()), aa_itr_end( ( *sse_itr)->End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {

          util::ShPtr< assemble::LocatorAtomCoordinatesInterface> current_resi
          (
            new restraint::LocatorCoordinatesFirstSideChainAtom( ( *aa_itr)->GetChainID(), ( *aa_itr)->GetSeqID())
          );

          const bool nterminus_bool( n_terminus % 2);
          const double weight( CalculateWeight( **sse_itr, **aa_itr, nterminus_bool));

          bool success
          (
            residue_weights.Insert
            (
              std::pair
              <
                util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
                storage::Pair< double, util::ShPtr< assemble::SSE> >
              >
              (
                current_resi, storage::Pair< double, util::ShPtr< assemble::SSE> >( weight, sse_copy)
              )
            ).second
          );

          BCL_Assert( success, "was not able to insert residue " + current_resi->GetIdentification());
        }

        ++n_terminus;
      }

      return residue_weights;
    }

    //! @brief calculates the weight of the given residue according to its position away from the end of the given sse
    //! @param SSE the sse whose ends will be used to calculate the weight
    //! @param AA_BASE the residue of interest whose weight will be calculated
    //! @param N_TERMINUS indicates true if nterminus of the sse should have the largest weight - otherwise cterminus
    //! @return double which is the weight of the residue
    double DataSetPairwiseBipolar::CalculateWeight( const assemble::SSE &SSE, const biol::AABase &AA_BASE, const bool N_TERMINUS)
    {
      double aa_resi_num( AA_BASE.GetSeqID());
      double weight( double( aa_resi_num - SSE.GetFirstAA()->GetSeqID()) / double( SSE.GetSize() - 1));
      if( N_TERMINUS)
      {
        weight = double( 1) - weight;
      }

      return weight;
    }

  } // namespace score
} // namespace bcl
