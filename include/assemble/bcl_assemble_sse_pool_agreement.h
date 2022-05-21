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

#ifndef BCL_ASSEMBLE_SSE_POOL_AGREEMENT_H_
#define BCL_ASSEMBLE_SSE_POOL_AGREEMENT_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPoolAgreement
    //! @brief evaluate the agreement between two given pools
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_assemble_sse_pool_agreement.cpp @endlink
    //! @author woetzen
    //! @date Jul 15, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPoolAgreement :
      public math::BinaryFunctionInterfaceSerializable< SSEPool, SSEPool, double>
    {

    private:

    //////////
    // data //
    //////////

      //! scheme of the function
      std::string m_Scheme;

      //! log of the differences
      bool m_Log;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief default scheme as string
      //! @return string for default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param SCHEME scheme for function
      //! @param LOG before adding differences, take log
      SSEPoolAgreement( const bool LOG = true, const std::string &SCHEME = GetDefaultScheme());

      //! @brief Clone function
      //! @return pointer to new SSEPoolAgreement
      SSEPoolAgreement *Clone() const;

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

      //! @brief virtual operator taking two ARGUMENTs and returning a t_ResultType object
      //! @param POOL the pool to be evaluated
      //! @param POOL_TEMPLATE the template to evaluate against
      //! @return double that expresses the agreement of the pool to the template
      double AgreementToTemplate( const SSEPool &POOL, const SSEPool &POOL_TEMPLATE) const;

      //! @brief Q3 score for two SSEPools
      //! @param POOL_A the pool a
      //! @param POOL_B the pool b
      //! @return number_correct_sse/total_number_res * 100%
      double Q3Score( const SSEPool &POOL_A, const SSEPool &POOL_B) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking two pools and returning the symmetric agreement as the sum of the agreements
      //!        if each pool is treated as template to the other
      //! @param POOL_A the pool a
      //! @param POOL_B the pool b
      //! @return double that expresses the agreement between two pools
      double operator()( const SSEPool &POOL_A, const SSEPool &POOL_B) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief calculate the overlap between two sequences
      //! @param SEQUENCE the sequence of interest
      //! @param SEQUENCE_TEMPLATE the reference sequence
      //! @return pair where first is the overlap on the left (first seq id difference) and the right (last seqid difference),
      //!         where the sign indicates if the SEQUENCE is shorter (negative) or longer (positive) then the template, sum
      //!         of those two is the overall length difference
      static storage::VectorND< 2, int> Overlap
      (
        const biol::AASequence &SEQUENCE,
        const biol::AASequence &SEQUENCE_TEMPLATE
      );

    }; // class SSEPoolAgreement

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_POOL_AGREEMENT_H_ 
