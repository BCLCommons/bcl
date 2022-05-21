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

#ifndef BCL_SCORE_SSE_PAIRS_FRAGMENTS_H_
#define BCL_SCORE_SSE_PAIRS_FRAGMENTS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_sse_pack_interface.h"
#include "assemble/bcl_assemble_sse_geometry_packing_list_pickers.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPairsFragments
    //! @brief is scoring function that sums up scores for fragments of a given SSE Pair
    //! @details SSEPairsFragments generates all possible pairings between fragments of two given SSEs and calls the provided
    //! scoring function and returns the summed up result
    //!
    //! @author karakam
    //! @date 04/02/09
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPairsFragments :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      //! SSEGeometryPacker to use to calculate the packing between SSEPairsFragments
      assemble::SSEGeometryPackingListPicker m_Packer;

      //! ShPtr to function interface which will be called and summed over for all sse pair fragments
      util::Implementation< SSEPackInterface> m_SPScorePair;

      //! boolean to determine whether the sum should be normalized by total number of pairs
      bool m_NormalizeOverNumberOfPairs;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEPairsFragments();

      //! @brief construct from a Function Interface and a normalization flag boolean
      //! @param PACKING_PICKER Function that will be used to calculate the packings between fragments
      //! @param PACKING_FUNCTION Function that will be used to calculate scores for all packings
      //! @param NORMALIZE_OVER_NUMBER_OF_PAIRS boolean flag to determine whether to normalize over number of pairs,
      SSEPairsFragments
      (
        const assemble::SSEGeometryPackingListPicker &PACKING_PICKER,
        const SSEPackInterface &PACKING_FUNCTION,
        const bool NORMALIZE_OVER_NUMBER_OF_PAIRS
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new SSEPairsFragments instance copied from this one
      SSEPairsFragments *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the given pair of SSEs are valid
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @return whether the given pair of SSEs are valid
      bool AreValidSSEs( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

      //! @brief f to sum up all pairwise scores for the given SSE Pair
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @return summed up score for the given SSE Pair
      double operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @param INDENT indentation
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write the Scheme and the Function value for the ARGUMENT to the STREAM
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @param OSTREAM std::sotream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B,
        std::ostream &OSTREAM
      ) const;

    }; // class SSEPairsFragments

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_SSE_PAIRS_FRAGMENTS_H_
