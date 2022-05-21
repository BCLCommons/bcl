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

#ifndef BCL_SCORESTAT_PROTEIN_MODEL_SSE_TRIPLET_CHIRALITY_H_
#define BCL_SCORESTAT_PROTEIN_MODEL_SSE_TRIPLET_CHIRALITY_H_

// include the namespace header
#include "bcl_scorestat.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "math/bcl_math_histogram.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector.h"

namespace bcl
{
  namespace scorestat
  {
    /////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelSSETripletChirality
    //! @brief extracts triplet chirality log-odds from a set of proteins
    //!        for example, finding the number of times that Helix-Strand-Helix occurs with the Strand on the LHS vs RHS
    //!        of the plane formed by the helices
    //!
    //! @see @link example_scorestat_protein_model_sse_triplet_chirality.cpp @endlink
    //! @author mendenjl
    //! @date Mar 09, 2017
    //////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelSSETripletChirality :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      double m_InteractionDistance;     //!< Maximum distance between "contacting" SC atoms between SSEs
      size_t m_MinAtomsInContactHelix;  //!< Minimum contacting SC atoms between SSEs to count the SSEs as interacting
      size_t m_MinAtomsInContactStrand; //!< Minimum contacting SC atoms between SSEs to count the SSEs as interacting
      bool   m_OutputPerProtein;        //!< Output per-protein counts of each triplet type
      bool   m_ConsiderWhichSSEsInContact; //!< Whether to consider which SSE pair, if any, is not in contact

    public:

      //! Total Number of propensities computed by this class
      //! 8 (Helix/Strand triplets) X
      //! 8 (anti-parallel/parallel triplets) X
      //! 4 (in-contact/not-in-contact triplet combos with at least two sses in contact)
      //! 2 (adjacent/non-adjacent SSEs)
      //! 2 (LHS/RHS)
      enum
      {
        s_NumberHashesContactSplit = 1024,
        s_NumberHashesNoContactSplit = 256
      };

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, takes whether to consider which SSES are in contact
      explicit ProteinModelSSETripletChirality( const bool CONSIDER_WHICH_SSES_IN_CONTACT = false) :
        m_InteractionDistance( 7),
        m_MinAtomsInContactHelix( 1),
        m_MinAtomsInContactStrand( 1),
        m_OutputPerProtein( false),
        m_ConsiderWhichSSEsInContact( CONSIDER_WHICH_SSES_IN_CONTACT)
      {
      }

      //! @brief virtual copy constructor
      ProteinModelSSETripletChirality *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as constant reference to std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief get the minimum atoms in contact
      //! @return minimum atoms in contact
      const size_t &GetMinAtomsInContactHelix() const
      {
        return m_MinAtomsInContactHelix;
      }

      //! @brief get the minimum atoms in contact
      //! @return minimum atoms in contact
      const size_t &GetMinAtomsInContactStrand() const
      {
        return m_MinAtomsInContactStrand;
      }

      //! @brief get the minimum interaction distance between the SSEs
      //! @return minimum interaction distance between the SSEs
      const double &GetInteractionDistance() const
      {
        return m_InteractionDistance;
      }

      //! @brief get the expected number of hashes
      //! @return the expected number of unique hash objects
      size_t GetNumberHashes() const
      {
        return size_t( m_ConsiderWhichSSEsInContact ? s_NumberHashesContactSplit : s_NumberHashesNoContactSplit);
      }

      //! @brief hash a given packing. The given string is intended to be fast to hash, not necessarily easily readable
      std::string GetPackingTripletHash
      (
        const biol::SSType &TYPE_A,
        const biol::SSType &TYPE_B,
        const biol::SSType &TYPE_C,
        const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AB,
        const assemble::SSEGeometryPacking::OrientationEnum &PACKING_BC,
        const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AC,
        const size_t &COUNTS_PACK_A,
        const size_t &COUNTS_PACK_B,
        const size_t &COUNTS_PACK_C,
        const bool   &ADJACENT,
        const bool   &RHS
      ) const;

      //! @brief hash a given packing into a size_t for speed. The size_t will have the range 0 - 223 ( 7 * 8 * 4)
      size_t GetPackingTripletNumber
      (
        const biol::SSType &TYPE_A,
        const biol::SSType &TYPE_B,
        const biol::SSType &TYPE_C,
        const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AB,
        const assemble::SSEGeometryPacking::OrientationEnum &PACKING_BC,
        const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AC,
        const size_t &COUNTS_PACK_A,
        const size_t &COUNTS_PACK_B,
        const size_t &COUNTS_PACK_C,
        const bool   &ADJACENT,
        const bool   &RHS
      ) const;

      //! @brief convert a packing triplet number to a string
      std::string GetPackingTripletString( const size_t &HASH_NUMBER) const;

      //! @brief helper function to cache orientational information
      //! @param CACHE the cache to retrieve/store data in
      //! @param SSE_A_ID index of SSE_A
      //! @param SSE_B_ID index of SSE_B
      //! @param SSE_A, SSE_B the two sses of interest
      //! @return the orientation
      const assemble::SSEGeometryPacking::OrientationEnum &GetCacheOrientation
      (
        storage::Vector< storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> > &CACHE,
        const size_t &SSE_A_ID,
        const size_t &SSE_B_ID,
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B
      ) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the protein ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief convert a packing triplet number to a string
      linal::Vector< double> ComputeExpectedCounts
      (
        const linal::Matrix< double> &ADJ_CONTACT_SEEN, // Rows - ss triplet-type, Cols - contact type
        const linal::Matrix< double> &DIST_CONTACT_SEEN, // Rows - ss triplet-type, Cols - contact type
        const linal::Matrix< double> &PARALLEL_PROBS, // Rows - ss pair-type, Cols - SSE distance (1,2,3+)
        const linal::Vector< size_t> &SSE_TRIPLET_COUNTS_ADJACENT, // size - 7, sse triplet type, lowest value indicates first sse
        const linal::Vector< size_t> &SSE_TRIPLET_COUNTS // size - 7, sse triplet type, lowest value indicates first sse, 0 - helix, 1 - strand,
      ) const;
    }; // end class LoopDistanceStatistics
  } // namespace scorestat
} // namespace bcl

#endif // BCL_SCORESTAT_PROTEIN_MODEL_SSE_TRIPLET_CHIRALITY_H_
