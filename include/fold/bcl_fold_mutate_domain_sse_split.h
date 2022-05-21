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

#ifndef BCL_FOLD_MUTATE_DOMAIN_SSE_SPLIT_H_
#define BCL_FOLD_MUTATE_DOMAIN_SSE_SPLIT_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDomainSSESplit
    //! @brief splits each sse in a domain at a point biased towards the middle by a gaussian distribution
    //! @details Every sse in the domain will be split into two sses. The place of splitting is randomly determined but
    //!          biased towards the middle. The bias is based on gaussian distribution. The bias can be adjusted by
    //!          changing the standard deviation. Small standard deviation will greatly favor selecting the middle
    //!          of each sse as the split point. A large standard deviation will make it more equally probable to put
    //!          the cut point at any location in the sse.
    //!
    //! @see @link example_fold_mutate_domain_sse_split.cpp @endlink
    //! @author alexanns
    //! @date Sep 3, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDomainSSESplit :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

        //! the scheme for this mutate
        std::string m_Scheme;

        //! the standard deviation of the gaussian distribution that biases splitting towards the middle
        double m_StandardDeviation;

        //! the random number generator that should be used
        const random::DistributionInterface &m_RandomNumberGenerator;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateDomainSSESplit();

      //! @brief default constructor
      //! @param SCHEME the scheme for this mutate
      //! @param STANDARD_DEVIATION standard deviation of gaussian distribution that biases splitting towards the middle
      //! @param RNG the random number generator that should be used
      MutateDomainSSESplit
      (
        const std::string &SCHEME, const double STANDARD_DEVIATION,
        const random::DistributionInterface &RNG = random::GetGlobalRandom()
      );

      //! @brief Clone function
      //! @return pointer to new MutateDomainSSESplit
      MutateDomainSSESplit *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that takes a domain
      //! @param DOMAIN domain which will be mutated
      //! @return MutateResult with the mutated domain
      math::MutateResult< assemble::Domain> operator()( const assemble::Domain &THIS_DOMAIN) const;

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

    private:

      //! @brief splits an sse at a random position determined by the random number generator
      //! @param CURRENT_SSE the sse which will be split
      //! @return vectorND 2 with 2 shared pointers to sses which were created from the split sse
      storage::VectorND< 2, util::ShPtr< assemble::SSE> >
      SplitSSE( const assemble::SSE &CURRENT_SSE) const;

      size_t GetSplitPoint( const size_t SSE_LENGTH) const;

    }; // class MutateDomainSSESplit

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_DOMAIN_SSE_SPLIT_H_ 
