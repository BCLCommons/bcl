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
#include "assemble/bcl_assemble_locator_domain_random.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorDomainRandom::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorDomainRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorDomainRandom::LocatorDomainRandom() :
      m_DomainSizeRange( 1, 1),
      m_SSType( biol::GetSSTypes().COIL)
    {
    }

    //! @brief constructor from a domain size range and SSTypes
    //! @param DOMAIN_SIZE_RANGE min and max sizes of the domain to be collected
    //! @param SS_TYPE SSType to be collected
    LocatorDomainRandom::LocatorDomainRandom
    (
      const math::Range< size_t> &DOMAIN_SIZE_RANGE,
      const biol::SSType &SS_TYPE
    ) :
      m_DomainSizeRange( DOMAIN_SIZE_RANGE),
      m_SSType( SS_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorDomainRandom
    LocatorDomainRandom *LocatorDomainRandom::Clone() const
    {
      return new LocatorDomainRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorDomainRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the domain size range
    //! @return the domain size range
    const math::Range< size_t> &LocatorDomainRandom::GetDomainSizeRange() const
    {
      return m_DomainSizeRange;
    }

    //! @brief returns the sstype
    //! @return the sstype
    const biol::SSType &LocatorDomainRandom::GetSSType() const
    {
      return m_SSType;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locates a random domain and returns it
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return randomly located domain
    util::ShPtr< Domain> LocatorDomainRandom::Locate( const ProteinModel &PROTEIN_MODEL) const
    {
      // collect the SSEs from the given model as a domain
      util::ShPtr< Domain> sp_domain( new Domain( PROTEIN_MODEL.GetSSEsAsDomain( m_SSType)));

      // if there are not enough SSEs
      if( sp_domain->GetNumberSSEs() < m_DomainSizeRange.GetMin())
      {
        // warn user and return
        BCL_MessageVrb
        (
          "The given model does not have enough SSEs " + util::Format()( sp_domain->GetNumberSSEs()) + " vs " +
          util::Format()( m_DomainSizeRange.GetMin())
        )
        return util::ShPtr< Domain>();
      }

      // determine the size of the domain
      const size_t domain_size
      (
        random::GetGlobalRandom().SizeT
        (
          math::Range< size_t>
          (
            m_DomainSizeRange.GetMin(),
            std::min( sp_domain->GetNumberSSEs(), m_DomainSizeRange.GetMax())
          )
        )
      );

      // while not reached
      while( !sp_domain->IsEmpty() && sp_domain->GetNumberSSEs() != size_t( domain_size))
      {
        // get a random iterator
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator random_itr
        (
          random::GetGlobalRandom().Iterator
          (
            sp_domain->GetData().Begin(), sp_domain->GetData().End(), sp_domain->GetNumberSSEs()
          )
        );

        // remove the element from the domain
        sp_domain->Remove( **random_itr);
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
    std::istream &LocatorDomainRandom::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DomainSizeRange, ISTREAM);
      io::Serialize::Read( m_SSType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorDomainRandom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DomainSizeRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
