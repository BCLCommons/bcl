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

#ifndef BCL_ASSEMBLE_LOCATOR_DOMAIN_SSE_POOL_OVERLAPPING_H_
#define BCL_ASSEMBLE_LOCATOR_DOMAIN_SSE_POOL_OVERLAPPING_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_pool.h"
#include "find/bcl_find_locator_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorDomainSSEPoolOverlapping
    //! @brief Locates a domain in a protein model that has SSEs from the model that overlap with SSEs from an SSEPool
    //! @details
    //!
    //! @see @link example_assemble_locator_domain_sse_pool_overlapping.cpp @endlink
    //! @author alexanns
    //! @date Jun 15, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorDomainSSEPoolOverlapping :
      public find::LocatorInterface< util::ShPtr< Domain>, ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! the pool that will be used as the basis for the domain
      SSEPool m_Pool;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorDomainSSEPoolOverlapping();

      //! @brief constructor taking parameters
      //! @param POOL the pool that will be used as the basis to locate the domain
      LocatorDomainSSEPoolOverlapping( const SSEPool &POOL);

      //! @brief Clone function
      //! @return pointer to new LocatorDomainSSEPoolOverlapping
      LocatorDomainSSEPoolOverlapping *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns reference to pool
      //! @return gives const reference to sse pool
      const SSEPool &GetPool() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief locate function to locate domain from protein model
      //! @param MODEL the model from which the domain will be located
      //! @return shptr to the located domain
      util::ShPtr< Domain> Locate( const ProteinModel &MODEL) const;

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

    }; // class LocatorDomainSSEPoolOverlapping

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_LOCATOR_DOMAIN_SSE_POOL_OVERLAPPING_H_
