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

#ifndef BCL_ASSEMBLE_SSE_POOL_MUTATE_SSE_H_
#define BCL_ASSEMBLE_SSE_POOL_MUTATE_SSE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPoolMutateSSE
    //! @brief mutate a single sse in a pool of sses
    //! @details  locate an sse in the pool, and replaces it with its mutated counterpart
    //!
    //! @see @link example_assemble_sse_pool_mutate_sse.cpp @endlink
    //! @author woetzen
    //! @date Jun 17, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPoolMutateSSE :
      public math::MutateInterface< SSEPool>
    {

    private:

    //////////
    // data //
    //////////

      //! locator for an sse
      util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > m_LocateSSE;

      //! mutate a single sse that was picked
      util::ShPtr< math::MutateInterface< SSE> > m_MutateSSE;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from a sse pick and mutate
      //! @param SP_LOCATE_SSE locator of sse from domain
      //! @param SP_SSE_MUTATE mutate for an sse
      SSEPoolMutateSSE
      (
        const util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > &SP_LOCATE_SSE,
        const util::ShPtr< math::MutateInterface< SSE> > &SP_SSE_MUTATE
      );

      //! @brief Clone function
      //! @return pointer to new SSEPoolMutateSSE
      SSEPoolMutateSSE *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that mutates a pool by mutating a single sse
      //! @param SSE_POOL
      //! @return MutateResult that results from mutating to the SSE_POOL
      math::MutateResult< SSEPool> operator()( const SSEPool &SSE_POOL) const;

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

    }; // class SSEPoolMutateSSE

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_POOL_MUTATE_SSE_H_ 
