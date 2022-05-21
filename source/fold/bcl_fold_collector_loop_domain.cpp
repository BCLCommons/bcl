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
#include "fold/bcl_fold_collector_loop_domain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "fold/bcl_fold_loop_domain.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorLoopDomain::s_Instance
    (
      util::Enumerated< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> >::AddInstance
      (
        new CollectorLoopDomain()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorLoopDomain::CollectorLoopDomain() :
      m_CollectUnclosedOnly( false),
      m_LoopClosureThreshold()
    {
    }

    //! @brief constructor from variables
    //! @param COLLECT_UNCLOSED_ONLY boolean to collect only the unclosed ones
    //! @param LOOP_CLOSURE_THRESHOLD Distance sum threshold for identifying closed loops
    CollectorLoopDomain::CollectorLoopDomain
    (
      const bool COLLECT_UNCLOSED_ONLY,
      const double LOOP_CLOSURE_THRESHOLD
    ) :
      m_CollectUnclosedOnly( COLLECT_UNCLOSED_ONLY),
      m_LoopClosureThreshold( LOOP_CLOSURE_THRESHOLD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorLoopDomain
    CollectorLoopDomain *CollectorLoopDomain::Clone() const
    {
      return new CollectorLoopDomain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorLoopDomain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &CollectorLoopDomain::GetAlias() const
    {
      static const std::string s_name( "CollectorLoopDomain");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorLoopDomain::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects loop domains in a protein model.");
      serializer.AddInitializer
      (
        "only unclosed",
        "collect only unclosed loop domains",
        io::Serialization::GetAgent( &m_CollectUnclosedOnly)
      );
      serializer.AddInitializer
      (
        "threshold",
        "loop closure threshold",
        io::Serialization::GetAgent( &m_LoopClosureThreshold)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Collect loop domains from given protein model
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return loop domains collected from given protein model
    util::ShPtrList< LoopDomain> CollectorLoopDomain::Collect( const assemble::DomainInterface &PROTEIN_MODEL) const
    {
      // try to cast the given DomainInterface to ProteinModel
      util::SiPtr< const assemble::ProteinModel> sp_model( &PROTEIN_MODEL);
      BCL_Assert( sp_model.IsDefined(), "The cast from DomainInterface to ProteinModel failed!");

      // first get the loop domain locators from protein model
      util::ShPtr< util::ShPtrList< LocatorLoopDomain> > sp_locators
      (
        sp_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_LoopDomainLocators)
      );

      // initialize list of loop domains
      util::ShPtrList< LoopDomain> loop_domains;

      // make sure it is defined
      if( !sp_locators.IsDefined())
      {
        BCL_MessageCrt( "No loop domain locators stored with the model");
        return loop_domains;
      }

      // iterate over the locators
      for
      (
        util::ShPtrList< LocatorLoopDomain>::const_iterator
          locator_itr( sp_locators->Begin()), locator_itr_end( sp_locators->End());
        locator_itr != locator_itr_end; ++locator_itr
      )
      {
        // locate the loop domain
        util::ShPtr< LoopDomain> this_domain( ( *locator_itr)->Locate( PROTEIN_MODEL));

        // if collecting not only the unclosed ones
        if( !m_CollectUnclosedOnly)
        {
          // then directly insert it
          loop_domains.PushBack( this_domain);
        }
        else
        {
          // otherwise insert only if the loop is not closed
          if( !LocatorLoopDomain::IsClosed( **locator_itr, *sp_model, m_LoopClosureThreshold))
          {
            loop_domains.PushBack( this_domain);
          }
        }
      }

      return loop_domains;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CollectorLoopDomain::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_CollectUnclosedOnly, ISTREAM);
      io::Serialize::Read( m_LoopClosureThreshold, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CollectorLoopDomain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_CollectUnclosedOnly, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LoopClosureThreshold, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
