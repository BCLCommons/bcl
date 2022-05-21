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
#include "assemble/bcl_assemble_locator_domain_specified.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_sse.h"
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
    const util::SiPtr< const util::ObjectInterface> LocatorDomainSpecified::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorDomainSpecified())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorDomainSpecified::LocatorDomainSpecified() :
      m_Locators()
    {
    }

    //! @brief constructor taking member variable
    //! @param LOCATORS locators to specify the sses that make up the domain
    LocatorDomainSpecified::LocatorDomainSpecified
    (
      const util::ShPtrList
      <
        find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
      > &LOCATORS
    ) :
      m_Locators( LOCATORS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorDomainSpecified
    LocatorDomainSpecified *LocatorDomainSpecified::Clone() const
    {
      return new LocatorDomainSpecified( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorDomainSpecified::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief writes a pymol formatted script to stream that highlights this domain
    //! @param OSTREAM the stream to which the script will be written to
    //! @param PYMOL_NAME the name of the selection should be in pymol
    //! @return ostream that the script was written to
    std::ostream &LocatorDomainSpecified::WritePymolDomainFile
    (
      std::ostream &OSTREAM, const std::string &PYMOL_NAME
    ) const
    {
      OSTREAM << "select " << PYMOL_NAME << ", ";
      // iterate through the locators
      for
      (
        util::ShPtrList
        <
          find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
        >::const_iterator
          itr( m_Locators.Begin()), itr_end( m_Locators.End());
        itr != itr_end;
        ++itr
      )
      {
        // try to cast the pointer to LocatorSSE
        const util::ShPtr< LocatorSSE> sse_locator( *itr);

        // true if the sse pointer is not defined
        if( !sse_locator.IsDefined())
        {
          // skip to next sse
          BCL_MessageDbg( "locator cannot be cast to a LocatorSSE and won't be written");
          continue;
        }

        // true if not at first locator
        if( itr != m_Locators.Begin())
        {
          // print plus
          OSTREAM << " + ";
        }

        // write the information in pymol format for the current sse
        OSTREAM << " chain " << sse_locator->GetChainID() <<
          " and resi " << sse_locator->GetSSEID().First() << "-" << sse_locator->GetSSEID().Second();
      }

      return OSTREAM;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief locate function to locate domain from protein model
    //! @param MODEL the model from which the domain will be located
    //! @return shptr to the located domain
    util::ShPtr< Domain> LocatorDomainSpecified::Locate( const ProteinModel &MODEL) const
    {
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> sses;

      // iterate through the locators
      for
      (
        util::ShPtrList
        <
          find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
        >::const_iterator
          locator_itr( m_Locators.Begin()), locator_itr_end( m_Locators.End());
        locator_itr != locator_itr_end;
        ++locator_itr
      )
      {
        // locate the sse
        const util::SiPtr< const SSE> located_sse( ( *locator_itr)->Locate( MODEL));

        // true if the sse could not be located
        if( !located_sse.IsDefined())
        {
          // go to next
          BCL_MessageDbg( "could not locate sse " + util::Format()( **locator_itr));
          continue;
        }

        // clone the located sse
        const util::ShPtr< SSE> new_sse( located_sse->Clone());

        // insert the sse into the set of sses
        BCL_Assert( sses.Insert( new_sse).second, "could not insert sse " + new_sse->GetIdentification());
      }

      // create domain from the sses
      const util::ShPtr< Domain> domain( new Domain( sses));

      return domain;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorDomainSpecified::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Locators, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorDomainSpecified::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Locators, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble

} // namespace bcl
