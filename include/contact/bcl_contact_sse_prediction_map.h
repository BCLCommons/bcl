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

#ifndef BCL_CONTACT_SSE_PREDICTION_MAP_H_
#define BCL_CONTACT_SSE_PREDICTION_MAP_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_contact_prediction_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPredictionMap
    //! @brief Provides contact probabilities for given SSE pairs
    //! @details  Calculates contact probabilities for each SSE pair by summing over contact predictions stored in
    //! contact map for each residue pair from these SSE pairs. The values are calculated once and then retrived from
    //! the stored map.
    //!
    //! @see @link example_contact_sse_prediction_map.cpp @endlink
    //! @author karakam
    //! @date 30.5.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPredictionMap :
      public storage::ObjectNDHashMap< 4, biol::AAData, double>
    {

    private:

    //////////
    // data //
    //////////

      //! PredictionMap for AA's
      util::ShPtr< PredictionMap> m_PredictionMap;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEPredictionMap();

      //! @brief constructor from a prediction map and ShPtrVector of SSEs
      //! @param PREDICTION_MAP PredictionMap that contains contact predictions for amino acids
      //! @param SSE_POOL ShPtrVector of SSEs of interest
      SSEPredictionMap
      (
        const util::ShPtr< PredictionMap> &PREDICTION_MAP,
        const util::ShPtrVector< assemble::SSE> &SSE_POOL
      );

      //! @brief virtual copy constructor
      SSEPredictionMap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the prediction stored for a given SSE pair
      //! @param SSE_PAIR SSE Pair of interest
      //! @return the prediction stored for the given SSE pair
      double GetPrediction( const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > &SSE_PAIR) const
      {
        return GetPrediction( *SSE_PAIR( 0), *SSE_PAIR( 1));
      }

      //! @brief returns the prediction stored for given SSEs
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @return prediction stored for SSE_A and SSE_B
      double GetPrediction( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

    private:

      //! @brief iterates over every SSE pair for a given ShPtrVector of SSEs and calculates the probabilities
      //! @param SSE_POOL ShPtrVector of SSEs of interest
      void FillMap( const util::ShPtrVector< assemble::SSE> &SSE_POOL);

      //! @brief calculates the contact probability given two SSEs
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @return contact probability for given pair of SSEs
      double CalculateContactProbability( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

    }; //class SSEPredictionMap

  } // namespace contact
} // namespace bcl

#endif //BCL_CONTACT_SSE_PREDICTION_MAP_H_
