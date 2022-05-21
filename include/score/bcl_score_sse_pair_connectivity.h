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

#ifndef BCL_SCORE_SSE_PAIR_CONNECTIVITY_H_
#define BCL_SCORE_SSE_PAIR_CONNECTIVITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPairConnectivity
    //! @brief Scores the ability to connect two SSEs in Euclidean space given the sequence separation between two SSEs.
    //! @details The sequence separation and Euclidean distance between the peptide bond atoms of both SSEs are
    //! calculated. The theoretical distance that could be "bridged" by the sequence separation and a residue extension
    //! with a peptide bond is subtracted from the Euclidean distance. If this result is smaller 1, it is a preferable
    //! distance, otherwise, the square of that difference is returned as a penalty score.
    //!
    //! @see @link example_score_sse_pair_connectivity.cpp @endlink
    //! @author alexanns, woetzen
    //! @date Nov 06, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPairConnectivity :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      //! scheme to be used in outputting schemes
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief GetExtendedResidueLength gives the length a completely extended residue can cover
      //! @return the length a completely extended residue can cover
      static const double &GetExtendedResidueLength();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a scheme
      //! @param SCHEME Scheme to be used
      SSEPairConnectivity
      (
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new SSEPairConnectivity
      SSEPairConnectivity *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

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

    public:

      //! @brief write the Scheme and the Function value for the ARGUMENT to the STREAM
      //! @param SSE_A the first sse whose closeness to being peptide bonded to the second sse will be scored
      //! @param SSE_B the second sse whose closeness to being peptide bonded to the first sse will be scored
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::SSE &SSE_A, const assemble::SSE &SSE_B,
        std::ostream &OSTREAM
      ) const;

    }; // class SSEPairConnectivity

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_SSE_PAIR_CONNECTIVITY_H_
