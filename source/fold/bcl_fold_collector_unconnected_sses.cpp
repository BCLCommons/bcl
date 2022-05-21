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
#include "fold/bcl_fold_collector_unconnected_sses.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> CollectorUnconnectedSSE::s_Instance
    (
      util::Enumerated< find::CollectorInterface< util::SiPtrList< const assemble::SSE>, assemble::DomainInterface> >::AddInstance
      (
        new CollectorUnconnectedSSE()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorUnconnectedSSE::CollectorUnconnectedSSE() :
      m_Direction(),
      m_TestBondConnected(),
      m_SSTypes(),
      m_IgnoreTermini()
    {

    }

    //! @brief construct from sequence direction
    //! @param DIRECTION the sequence direction is
    //! @param TEST_BOND_CONNECTION
    //! @param SS_TYPES
    //! @param IGNORE_TERMINAL_SSE
    CollectorUnconnectedSSE::CollectorUnconnectedSSE
    (
      const biol::AASequenceFlexibility::SequenceDirection &DIRECTION,
      const bool TEST_BOND_CONNECTION,
      const storage::Set< biol::SSType> &SS_TYPES,
      const bool IGNORE_TERMINAL_SSE
    ) :
      m_Direction( DIRECTION),
      m_TestBondConnected( TEST_BOND_CONNECTION),
      m_SSTypes( SS_TYPES),
      m_IgnoreTermini( IGNORE_TERMINAL_SSE)
    {
    }

    //! @brief clone function
    //! @return pointer to a new CollectorUnconnectedSSE
    CollectorUnconnectedSSE *CollectorUnconnectedSSE::Clone() const
    {
      return new CollectorUnconnectedSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &CollectorUnconnectedSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorUnconnectedSSE::GetAlias() const
    {
      static const std::string s_alias( "CollectorUnconnectedSSE");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorUnconnectedSSE::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Optimization implementation for Monte Carlo Metropolis algorithms.");
      serializer.AddInitializer
      (
        "sequence direction",
        "sequence direction in which to collect the sses",
        io::Serialization::GetAgent( &m_Direction)
      );
      serializer.AddInitializer
      (
        "bond connected",
        "whether to test bond connection",
        io::Serialization::GetAgent( &m_TestBondConnected)
      );
      serializer.AddInitializer
      (
        "sse types",
        "sse types to consider",
        io::Serialization::GetAgent( &m_SSTypes)
      );
      serializer.AddInitializer
      (
        "ignore termini",
        "ignore terminal loops",
        io::Serialization::GetAgent( &m_IgnoreTermini)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief collects unconnected SSEs from the given domain
    //! @param DOMAIN_INTERFACE domain from which to collect unconnected SSEs
    //! @return unconnected SSEs in the given domain
    util::SiPtrList< const assemble::SSE> CollectorUnconnectedSSE::Collect
    (
      const assemble::DomainInterface &DOMAIN_INTERFACE
    ) const
    {
      // get all coil sses
      util::SiPtrVector< const assemble::SSE> sses( DOMAIN_INTERFACE.GetSSEs( m_SSTypes));

      // collect all sses with m_Direction flexibility
      util::SiPtrList< const assemble::SSE> eligible_sses;

      for( util::SiPtrVector< const assemble::SSE>::const_iterator itr( sses.Begin()), itr_end( sses.End()); itr != itr_end; ++itr)
      {
        // get adjacent sses
        const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > adjacent( DOMAIN_INTERFACE.GetAdjacentSSEs( **itr));

        // ignore terminating sses
        if( m_IgnoreTermini)
        {
          const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > neighors( DOMAIN_INTERFACE.GetNeighborSSEs( **itr));

          // skip if termini, either neighbor is not defined
          if( !neighors.First().IsDefined() || !neighors.Second().IsDefined())
          {
            continue;
          }
        }

        // if the direction dependent sse is present and peptide bonded, ignore that coil
        if( m_Direction == biol::AASequenceFlexibility::e_CTerminal)
        {
          // skip if there is an adjacent sse and possibly peptide bond connected
          if
          (
            adjacent.Second().IsDefined()
            && ( !m_TestBondConnected || biol::AABase::AreAminoAcidsPeptideBonded( *( *itr)->GetLastAA(), *adjacent.Second()->GetFirstAA(), true))
          )
          {
            continue;
          }
        }
        else if( m_Direction == biol::AASequenceFlexibility::e_NTerminal)
        {
          // skip if there is an adjacent sse and possibly peptide bond connected
          if
          (
            adjacent.First().IsDefined() &&
            ( !m_TestBondConnected || biol::AABase::AreAminoAcidsPeptideBonded( *adjacent.First()->GetLastAA(), *( *itr)->GetFirstAA(), true))
          )
          {
            continue;
          }
        }
        else
        {
          continue;
        }

        // sse passed all conditions
        eligible_sses.PushBack( *itr);
       }

      return eligible_sses;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
