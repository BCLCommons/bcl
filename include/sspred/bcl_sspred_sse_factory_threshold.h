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

#ifndef BCL_SSPRED_SSE_FACTORY_THRESHOLD_H_
#define BCL_SSPRED_SSE_FACTORY_THRESHOLD_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_methods.h"
#include "assemble/bcl_assemble_sse_factory_interface.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_ss_types.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEFactoryThreshold
    //! @brief This is an SSEFactory class derived from SSEFactoryInterface
    //! @details It uses provided SSPredictions to generate sses from a given sequence
    //!
    //! @see @link example_sspred_sse_factory_threshold.cpp @endlink
    //! @author karakam, linders
    //! @date 22.02.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEFactoryThreshold :
      public assemble::SSEFactoryInterface
    {

    private:

    //////////
    // data //
    //////////

      //! ssmethod to be used to generate sses from sequence
      Method m_SSMethod;

      //! threshold value to classify ss predictions as SSE or not
      storage::Map< biol::SSType, double> m_SSEThresholds;

      //! boolean to indicate whether to add systematically extended copies of sses in pool
      bool m_ExtendSSEs;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default map for sse thresholds
      //! @return default map for sse thresholds
      static const storage::Map< biol::SSType, double> &GetDefaultThresholdsMap();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEFactoryThreshold();

      //! @brief constructor from specified SSMethod
      //! @param SS_METHOD set of ssmethod to be used for generating SSEs from AASequence
      //! @param SSE_THRESHOLDS threshold values to classify ss predictions as helix or strand or not
      //! @param EXTEND_SSES whether or not to add systematically extended copies of sses in pool
      SSEFactoryThreshold
      (
        const Method &SS_METHOD,
        const storage::Map< biol::SSType, double> &SSE_THRESHOLDS = GetDefaultThresholdsMap(),
        const bool EXTEND_SSES = false
      );

      //! @brief virtual copy constructor
      SSEFactoryThreshold *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief set the ss method used
      //! @param SS_METHOD sspred method
      void SetMethod( const Method &SS_METHOD)
      {
        m_SSMethod = SS_METHOD;
      }

      //! @brief returns SSMethod being used
      //! @return SSMethod being used
      const Method &GetMethod() const
      {
        return m_SSMethod;
      }

      //! @brief set the thresholds to use
      //! @brief SSTYPE_THRESHOLDS the thresholds to use for the sstypes desired
      void SetThresholds( const storage::Map< biol::SSType, double> &SSTYPE_THRESHOLDS)
      {
        m_SSEThresholds = SSTYPE_THRESHOLDS;
      }

      //! @brief returns SSE Threshold being used
      //! @return SSE Threshold being used
      const storage::Map< biol::SSType, double> &SSEGetThresholds() const
      {
        return m_SSEThresholds;
      }

      //! @brief set to extend sses
      //! @param EXTEND_SSES bool to extend or not extend sses
      void SetExtendSSEs( const bool EXTEND_SSES)
      {
        m_ExtendSSEs = EXTEND_SSES;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that returns a set of SSEs for the given AASequence
      //! @param SEQUENCE AASequence from which the SSEPool is going to be built
      //! @return SSEPool built from provided SEQUENCE
      assemble::SSEPool
      operator()( const biol::AASequence &SEQUENCE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

      //! @brief identifies the regions with most likely probability to be in SSE and returns them in a ssepool
      //! @param SEQUENCE AASequence from which SSEs will be identified
      //! @param METHOD the prediction method considered
      //! @param SS_TYPE SSType of interest
      //! @param THRESHOLD the prediction threshold, that the prediction needs to be above to be considered for that sstype
      //! @return SSEPool composed of SSEs generated from SEQUENCE
      static assemble::SSEPool IdentifySSERegions
      (
        const biol::AASequence &SEQUENCE,
        const Method &METHOD,
        const biol::SSType &SS_TYPE,
        const double THRESHOLD
      );

      //! @brief elongates the SSE region of type SS_TYPE and beginning at AA_ITR
      //! @param AA_ITR iterator to the beginning of a new SSE region
      //! @param AA_ITR_END end iterator for the sequence
      //! @param SS_TYPE SSE type for which the window is going to be maximized
      //! @param METHOD the prediction method considered
      //! @param THRESHOLD the prediction threshold, that the prediction needs to be above to be extended
      //! @return iterator to the end of the region which SSE can be elongated
      static biol::AASequence::const_iterator
      ElongateSSERegion
      (
        const biol::AASequence::const_iterator &AA_ITR,
        const biol::AASequence::const_iterator &AA_ITR_END,
        const biol::SSType &SS_TYPE,
        const Method &METHOD,
        const double THRESHOLD
      );

      //! @brief systematically extends the SSEs in the given POOL and inserts them into POOL
      //! @param POOL SSEPool for which SSEs will be extended
      //! @param SEQUENCE AASequence on which the extensions will be based
      //! @return SSEPool which has SSEs in the original pool and the extended SSEs
      assemble::SSEPool SystematicallyExtendSSEs( assemble::SSEPool &POOL, const biol::AASequence &SEQUENCE) const;

    }; // class SSEFactoryThreshold

  } // namespace sspred
} // namespace bcl

#endif //BCL_SSPRED_SSE_FACTORY_THRESHOLD_H_
