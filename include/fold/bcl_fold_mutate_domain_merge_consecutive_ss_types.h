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

#ifndef BCL_FOLD_MUTATE_DOMAIN_MERGE_CONSECUTIVE_SS_TYPES_H_
#define BCL_FOLD_MUTATE_DOMAIN_MERGE_CONSECUTIVE_SS_TYPES_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "biol/bcl_biol_ss_types.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDomainMergeConsecutiveSSTypes
    //! @brief Mutate class to merge consecutive SSEs of a specified type in a domain
    //! @details mutates the SSE by iterating over all sses and merging consecutive ones into a new SSE and replacing
    //!          the old ones with the new one
    //!
    //! @see @link example_fold_mutate_domain_merge_consecutive_ss_types.cpp @endlink
    //! @author bitterd
    //! @date Sep 4, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDomainMergeConsecutiveSSTypes :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

        //! Scheme of the class
        std::string m_Scheme;

        //! SSE Types that should be merged
        biol::SSType m_SSType;

    public:

      //! single instance of that class
      static util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateDomainMergeConsecutiveSSTypes();

      //! @brief construct from SSType
      MutateDomainMergeConsecutiveSSTypes( const biol::SSType &SSTYPE);

      //! @brief Clone function
      //! @return pointer to new MutateDomainMergeConsecutiveSSTypes
      MutateDomainMergeConsecutiveSSTypes *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns scheme
      //! @return the scheme of the object
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      math::MutateResult< bcl::assemble::Domain> operator()( const assemble::Domain &DOMAINS) const;

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

    }; // class MutateDomainMergeConsecutiveSSTypes

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_DOMAIN_MERGE_CONSECUTIVE_SS_TYPES_H_ 
