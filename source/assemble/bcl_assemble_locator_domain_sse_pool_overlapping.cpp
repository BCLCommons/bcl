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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_locator_domain_sse_pool_overlapping.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorDomainSSEPoolOverlapping::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorDomainSSEPoolOverlapping())
    );
  
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorDomainSSEPoolOverlapping::LocatorDomainSSEPoolOverlapping() :
      m_Pool()
    {
    }

    //! @brief constructor taking parameters
    //! @param POOL the pool that will be used as the basis to locate the domain
    LocatorDomainSSEPoolOverlapping::LocatorDomainSSEPoolOverlapping( const SSEPool &POOL) :
      m_Pool( POOL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorDomainSSEPoolOverlapping
    LocatorDomainSSEPoolOverlapping *LocatorDomainSSEPoolOverlapping::Clone() const
    {
      return new LocatorDomainSSEPoolOverlapping( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorDomainSSEPoolOverlapping::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns reference to pool
    //! @return gives const reference to sse pool
    const SSEPool &LocatorDomainSSEPoolOverlapping::GetPool() const
    {
      return m_Pool;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief locate function to locate domain from protein model
    //! @param LOCATE_DOMAIN the model from which the domain will be located
    //! @return shptr to the located domain
    util::ShPtr< Domain> LocatorDomainSSEPoolOverlapping::Locate( const ProteinModel &MODEL) const
    {
      util::ShPtr< Domain> new_domain( new Domain());

      // get non overlapping set of sses from sse pool
      const util::SiPtrVector< const SSE> sses( m_Pool.GetSSEs());

      // get the sses in the protein model that overlap with them
      util::SiPtrList< const SSE> overlapping_sses
      (
        MODEL.GetOverlappingSSEs( util::SiPtrList< const SSE>( sses.Begin(), sses.End()))
      );

      // add the overlapping sses into the domain
      for
      (
        util::SiPtrList< const SSE>::const_iterator
        overlapping_sses_itr( overlapping_sses.Begin()), overlapping_sses_itr_end( overlapping_sses.End());
          overlapping_sses_itr != overlapping_sses_itr_end;
        ++overlapping_sses_itr
      )
      {
        // clone the sse and put it into the domain, any SSE repeats won't insert
        new_domain->Insert( util::ShPtr< SSE>( ( *overlapping_sses_itr)->Clone()));
      }

      // return the domain
      return new_domain;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorDomainSSEPoolOverlapping::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Pool, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorDomainSSEPoolOverlapping::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Pool, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble
} // namespace bcl
