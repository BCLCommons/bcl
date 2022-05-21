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

#ifndef BCL_SCORE_SSE_PAIR_GAP_H_
#define BCL_SCORE_SSE_PAIR_GAP_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPairGap
    //! @brief returns the number of amino acids between two SSEs as a penalty.
    //! @details Two SSEs that are neighbors within a chain should be passed. The operator returns the sequence
    //! separation between those two SSEs as a penalty, helping in completing loops in protein models.
    //!
    //! @see @link example_score_sse_pair_gap.cpp @endlink
    //! @author alexanns, woetzen
    //! @date Nov 7, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPairGap :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      //! scheme to be used in outputting schemes
      std::string m_Scheme;

      //! penalty if both ends are not coil
      double m_NonCoilPenalty;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a scheme
      //! @param NON_COIL_PENALTY additional penalty to add, if the ss type is not COIL
      //! @param SCHEME Scheme to be used
      SSEPairGap
      (
        const double NON_COIL_PENALTY,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new SSEPairGap
      SSEPairGap *Clone() const;

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief gives the score for two sses based on how close they are to being connected, if they are adjecent
      //! @param SSE_A the first sse whose closeness to being peptide bonded to the second sse will be scored
      //! @param SSE_B the second sse whose closeness to being peptide bonded to the first sse will be scored
      //! @return the score for two sses based on how close SSE_A and SSE_B are to being connected by a peptide bond
      double operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SSEPairGap

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_SSE_PAIR_GAP_H_
