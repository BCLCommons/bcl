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

#ifndef BCL_SCORE_SSE_PAIR_CONTACT_H_
#define BCL_SCORE_SSE_PAIR_CONTACT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "contact/bcl_contact.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPairContact
    //! @brief This is a Function derived class for scoring contacts between secondary structure element pairs
    //!
    //! @see @link example_score_sse_pair_contact.cpp @endlink
    //! @author karakam, woetzen
    //! @date 06.06.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPairContact :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      //! scoring functions
      util::ShPtrVector< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> > m_ScoringFunctions;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEPairContact();

      //! @brief constructor with a prediction map
      //! @param SP_PREDICTION_MAP ShPtr to PredictionMap to be used
      SSEPairContact( const util::ShPtr< contact::PredictionMap> &SP_PREDICTION_MAP);

      //! @brief virtual copy constructor
      //! @return pointer to a new SSEPairContact copied from this one
      SSEPairContact *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the score for the given SSE pair
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @return the score for the given SSE pair
      double operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

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
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B,
        std::ostream &OSTREAM
      ) const;

    }; //class SSEPairContact

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_SSE_PAIR_CONTACT_H_
