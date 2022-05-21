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
#include "assemble/bcl_assemble_sse_pool_mutate_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPoolMutateSSE::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new SSEPoolMutateSSE
        (
          util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> >(),
          util::ShPtr< math::MutateInterface< SSE> >()
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a sse pick and mutate
    //! @param SP_LOCATE_SSE locator of sse from domain
    //! @param SP_SSE_MUTATE mutate for an sse
    SSEPoolMutateSSE::SSEPoolMutateSSE
    (
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > &SP_LOCATE_SSE,
      const util::ShPtr< math::MutateInterface< SSE> > &SP_SSE_MUTATE
    ) :
      m_LocateSSE( SP_LOCATE_SSE),
      m_MutateSSE( SP_SSE_MUTATE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPoolMutateSSE
    SSEPoolMutateSSE *SSEPoolMutateSSE::Clone() const
    {
      return new SSEPoolMutateSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPoolMutateSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &SSEPoolMutateSSE::GetScheme() const
    {
      return m_MutateSSE->GetScheme();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that mutates a pool by mutating a single sse
    //! @param SSE_POOL
    //! @return MutateResult that results from mutating to the SSE_POOL
    math::MutateResult< SSEPool> SSEPoolMutateSSE::operator()( const SSEPool &SSE_POOL) const
    {
      // pick a sse from the pool
      const util::SiPtr< const SSE> located_sse( m_LocateSSE->Locate( SSE_POOL));

      // was picking successful
      if( !located_sse.IsDefined())
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // mutate the sse
      math::MutateResult< SSE> result( m_MutateSSE->operator ()( *located_sse));

      // was mutating successful
      if( !result.GetArgument().IsDefined())
      {
        return math::MutateResult< SSEPool>( util::ShPtr< SSEPool>(), *this);
      }

      // replace with new sse
      util::SiPtrList< const SSE> pool_sses( SSE_POOL.Begin(), SSE_POOL.End());
      util::SiPtrList< const SSE>::iterator new_end
      (
        std::remove_if
        (
          pool_sses.Begin(), pool_sses.End(),
          SSECompare( *located_sse)
        )
      );

      util::ShPtr< SSEPool> sp_pool( new SSEPool( util::SiPtrList< const SSE>( pool_sses.Begin(), new_end), false));
      sp_pool->Insert( result.GetArgument());

      // end
      return math::MutateResult< SSEPool>( sp_pool, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPoolMutateSSE::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LocateSSE, ISTREAM);
      io::Serialize::Read( m_MutateSSE, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPoolMutateSSE::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LocateSSE, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MutateSSE, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
