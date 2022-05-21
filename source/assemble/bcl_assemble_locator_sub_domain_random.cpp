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
#include "assemble/bcl_assemble_locator_sub_domain_random.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorSubDomainRandom::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorSubDomainRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorSubDomainRandom::LocatorSubDomainRandom() :
      m_SizeRange( 0, 0),
      m_LocateConsecutive( false),
      m_UseTopologyOrder( false)
    {
    }

    //! @brief construct from a range for number of SSEs in sub-domain
    //! @param SIZE_RANGE range for number of SSEs
    LocatorSubDomainRandom::LocatorSubDomainRandom( const math::Range< size_t> SIZE_RANGE) :
      m_SizeRange( SIZE_RANGE),
      m_LocateConsecutive( false),
      m_UseTopologyOrder( false)
    {
    }

    //! @brief construct from a range for number of SSEs in sub-domain and whether they should be located consecutively
    //! @param SIZE_RANGE range for number of SSEs
    //! @param LOCATE_CONSECUTIVE boolean to whether to locate consecutive SSEs
    //! @param USE_TOPOLOGY_ORDER boolean to order SSEs by the topology order
    LocatorSubDomainRandom::LocatorSubDomainRandom
    (
      const math::Range< size_t> &SIZE_RANGE,
      bool LOCATE_CONSECUTIVE,
      bool USE_TOPOLOGY_ORDER
    ) :
      m_SizeRange( SIZE_RANGE),
      m_LocateConsecutive( LOCATE_CONSECUTIVE),
      m_UseTopologyOrder( USE_TOPOLOGY_ORDER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorSubDomainRandom
    LocatorSubDomainRandom *LocatorSubDomainRandom::Clone() const
    {
      return new LocatorSubDomainRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorSubDomainRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locate and return a random sub-domain from the given domain
    //! @param SP_DOMAIN ShPtr to the domain of interest
    //! @return ShPtr to a random sub-domain located within the given SP_DOMAIN
    util::ShPtr< Domain> LocatorSubDomainRandom::Locate( const util::ShPtr< Domain> &SP_DOMAIN) const
    {
      // make sure given domain is defined
      BCL_Assert( SP_DOMAIN.IsDefined(), "The given domain is not defined");

      // determine start and end indices for the strands to sort
      const size_t nr_sses( SP_DOMAIN->GetNumberSSEs());

      // determine the subset size
      size_t subset_size( 0);

      // if the max of the range is below the nr_sses
      if( m_SizeRange.GetMax() <= nr_sses)
      {
        subset_size = random::GetGlobalRandom().SizeT( m_SizeRange);
      }
      // if it's smaller
      else
      {
        // determine subset size
        subset_size =
          random::GetGlobalRandom().SizeT
          (
            math::Range< size_t>( ( m_SizeRange.GetMin() <= nr_sses ? m_SizeRange.GetMin() : nr_sses), nr_sses)
          );
      }

      // create a new domain
      util::ShPtr< Domain> sp_domain( new Domain());

      // get all the SSEs in the model
      util::SiPtrVector< const SSE> all_sses( SP_DOMAIN->GetSSEs());

      // initialize vector to hold the selected SSEs
      util::SiPtrVector< const SSE> selected_sses;

      // if consecutive is required
      if( m_LocateConsecutive)
      {
        // if by topology order
        if( m_UseTopologyOrder)
        {
          // make sure topology is defined
          BCL_Assert
          (
            SP_DOMAIN->GetTopology().IsDefined(),
            "given domain has undefined topology, although ordering by topology order is requested"
          );

          // get the elements vector and cast them back to SSEs
          all_sses = util::SiPtrVector< const SSE>( SP_DOMAIN->GetTopology()->GetElements());
          BCL_Assert( all_sses.IsDefined(), "The dynamic cast of topology SSEs failed!");
        }

        // otherwise all_sses is sorted by the sequence order
        // now determine a start index
        const size_t start_index
        (
          random::GetGlobalRandom().SizeT
          (
            math::Range< size_t>( 0, all_sses.GetSize() - subset_size)
          )
        );

        // get the sub-siptrvector and construct the subset
        selected_sses = all_sses.SubSiPtrVector( start_index, subset_size);

      }
      // if no consecutive selection is required
      else
      {
        // while all_sses still has more than subset_size elements
        while( all_sses.GetSize() > subset_size)
        {
          // remove a random element
          all_sses.RemoveRandomElement();
        }

        // update the selected_sses
        selected_sses = all_sses;
      }

      // now that we know which SSEs should be in it and in which order they should be added
      // we can go ahead and construct the domains SSEs
      // first iterate over the selected SSEs
      for
      (
        util::SiPtrVector< const SSE>::const_iterator
          sse_itr( selected_sses.Begin()), sse_itr_end( selected_sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // find SSE
        util::ShPtr< SSE> sp_sse( SP_DOMAIN->FindSSE( **sse_itr));
        BCL_Assert( sp_sse.IsDefined(), "could not find a corresponding SSE from domain!");

        // insert it into the new domain
        sp_domain->Insert( sp_sse);
      }

      // if the topology is not defined
      if( SP_DOMAIN->GetTopology().IsDefined())
      {
        // set topology for the new domain
        sp_domain->SetTopology( SP_DOMAIN->GetTopology()->GetSubTopology( selected_sses));
      }

      // end
      return sp_domain;

    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorSubDomainRandom::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SizeRange, ISTREAM);
      io::Serialize::Read( m_LocateConsecutive, ISTREAM);
      io::Serialize::Read( m_UseTopologyOrder, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorSubDomainRandom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SizeRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LocateConsecutive, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseTopologyOrder, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
