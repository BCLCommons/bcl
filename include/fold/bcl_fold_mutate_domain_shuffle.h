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

#ifndef BCL_FOLD_MUTATE_DOMAIN_SHUFFLE_H_
#define BCL_FOLD_MUTATE_DOMAIN_SHUFFLE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDomainShuffle
    //! @brief shuffles the locations and the orientations of SSEs in a domain
    //! @details This class works on a selected domain by picking a random number of SSE pairs to shuffle with each other within
    //! the domain
    //!
    //! @see @link example_fold_mutate_domain_shuffle.cpp @endlink
    //! @author karakam
    //! @date Mar 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDomainShuffle :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

      //! Maximum number of swaps within one mutate
      size_t m_MaxNumberSwaps;

      //! whether to bend the swapped SSE to match the phi/psi angles of the previous SSE at that location
      bool m_Bend;

      //! scheme
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a max number swaps and scheme
      //! @param MAX_NUMBER_SWAPS Maximum number of swaps within one mutate
      //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
      //! @param SCHEME Scheme to be used
      MutateDomainShuffle
      (
        const size_t MAX_NUMBER_SWAPS = 2,
        const bool BEND = false,
        const std::string &SCHEME = GetStaticClassName< MutateDomainShuffle>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateDomainShuffle
      MutateDomainShuffle *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get maximum number of swaps in one mutate
      //! @return maximum number of swaps in one mutate
      size_t GetMaxNumberSwaps() const
      {
        return m_MaxNumberSwaps;
      }

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that takes a Domain and return a mutated Domain
      //! @param THIS_DOMAIN Domain which will be mutated
      //! @return MutateResult with the mutated Domain
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

    }; // class MutateDomainShuffle

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_DOMAIN_SHUFFLE_H_ 
