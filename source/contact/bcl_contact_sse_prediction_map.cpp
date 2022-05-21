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
#include "contact/bcl_contact_sse_prediction_map.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPredictionMap::SSEPredictionMap() :
      storage::ObjectNDHashMap< 4, biol::AAData, double>(),
      m_PredictionMap( new PredictionMap())
    {
    }

    //! @brief constructor from a prediction map and ShPtrVector of SSEs
    //! @param PREDICTION_MAP PredictionMap that contains contact predictions for amino acids
    //! @param SSE_POOL ShPtrVector of SSEs of interest
    SSEPredictionMap::SSEPredictionMap
    (
      const util::ShPtr< PredictionMap> &PREDICTION_MAP,
      const util::ShPtrVector< assemble::SSE> &SSE_POOL
    ) :
      storage::ObjectNDHashMap< 4, biol::AAData, double>(),
      m_PredictionMap( PREDICTION_MAP)
    {
      FillMap( SSE_POOL);
    }

    //! @brief virtual copy constructor
    SSEPredictionMap *SSEPredictionMap::Clone() const
    {
      return new SSEPredictionMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEPredictionMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns the prediction stored for given SSEs
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return prediction stored for SSE_A and SSE_B
    double SSEPredictionMap::GetPrediction( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      storage::VectorND< 4, util::SiPtr< const biol::AAData> > keys
      (
        util::SiPtr< const biol::AAData>( SSE_A.GetFirstAA()->GetData()),
        util::SiPtr< const biol::AAData>( SSE_A.GetLastAA()->GetData()),
        util::SiPtr< const biol::AAData>( SSE_B.GetFirstAA()->GetData()),
        util::SiPtr< const biol::AAData>( SSE_B.GetLastAA()->GetData())
      );

      storage::HashMap< size_t, double>::const_iterator itr = Find( keys);
      if( itr == End())
      {
        return 0;
      }
      return itr->second;
    }

    //! @brief iterates over every SSE pair for a given ShPtrVector of SSEs and calculates the probabilities
    //! @param SSE_POOL ShPtrVector of SSEs of interest
    void SSEPredictionMap::FillMap( const util::ShPtrVector< assemble::SSE> &SSE_POOL)
    {
      // iterate over every sse
      for
      (
        util::ShPtrVector< assemble::SSE>::const_iterator sse_itr_a( SSE_POOL.Begin()), sse_itr_end( SSE_POOL.End());
        sse_itr_a != sse_itr_end;
        ++sse_itr_a
      )
      {
        // versus every other sse
        for
        (
          util::ShPtrVector< assemble::SSE>::const_iterator sse_itr_b( SSE_POOL.Begin());
          sse_itr_b != sse_itr_end;
          ++sse_itr_b
        )
        {
          // if they are not the same SSE or do not overlap
          if( sse_itr_a != sse_itr_b && !biol::DoOverlap( **sse_itr_a, **sse_itr_b))
          {
            // value to store the sum of all predictions of every possible residue pairing between t two sses
            double sum( CalculateContactProbability( **sse_itr_a, **sse_itr_b));

            storage::VectorND< 4, util::SiPtr< const biol::AAData> > keys
            (
              util::SiPtr< const biol::AAData>( ( *sse_itr_a)->GetFirstAA()->GetData()),
              util::SiPtr< const biol::AAData>( ( *sse_itr_a)->GetLastAA()->GetData()),
              util::SiPtr< const biol::AAData>( ( *sse_itr_b)->GetFirstAA()->GetData()),
              util::SiPtr< const biol::AAData>( ( *sse_itr_b)->GetLastAA()->GetData())
            );

            // insert this value with corresponding keys into the ObjectNDHashMap
            Insert( keys, sum);
          }
        }
      }
    }

    //! @brief calculates the contact probability given two SSEs
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return contact probability for given pair of SSEs
    double SSEPredictionMap::CalculateContactProbability( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      // value to store the sum of all predictions of every possible residue pairing between t two sses
      double sum( 0.0);

      // value to store number of residue pairs for which no predictions existed
      size_t no_predictions( 0);

      // contact type of these two sses
      Type contact_type( GetTypes().TypeFromSSTypes( SSE_A, SSE_B));

      // iterate over every residue in SSE_A
      for
      (
        biol::AASequence::const_iterator aa_itr_a( SSE_A.GetData().Begin()),
          aa_itr_a_end( SSE_A.GetData().End());
        aa_itr_a != aa_itr_a_end;
        ++aa_itr_a
      )
      {
        // iterate over every residue in SSE_B
        for
        (
          biol::AASequence::const_iterator aa_itr_b( SSE_B.GetData().Begin()),
            aa_itr_b_end( SSE_B.GetData().End());
          aa_itr_b != aa_itr_b_end;
          ++aa_itr_b
        )
        {
          storage::VectorND< 2, util::SiPtr< const biol::AAData> > key
          (
            util::SiPtr< const biol::AAData>( ( *aa_itr_a)->GetData()),
            util::SiPtr< const biol::AAData>( ( *aa_itr_b)->GetData())
          );
          double prediction( ( m_PredictionMap->GetPredictions( key)).First()( contact_type));

          if( !util::IsDefined( prediction))
          {
            no_predictions++;
          }
          else
          {
            sum += prediction;
          }
        }
      }

      // return summed probability normalized with length of SSEs
      return sum / ( ( SSE_A.GetData().GetSize() * SSE_B.GetData().GetSize()) - no_predictions);
    }

  } // namespace contact
} // namespace bcl
