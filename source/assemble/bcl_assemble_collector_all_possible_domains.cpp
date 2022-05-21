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
#include "assemble/bcl_assemble_collector_all_possible_domains.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_sse_geometry_packing_criteria.h"
#include "assemble/bcl_assemble_sse_geometry_packing_pickers.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_combination.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorAllPossibleDomains::s_Instance
    (
      GetObjectInstances().AddInstance( new CollectorAllPossibleDomains())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the domain size and a packing criteria
    //! @param DOMAIN_SIZE domain size range
    //! @param CONNECTED_DOMAIN bool whether the SSEs in the domain should remain connected
    //! @param PACKING_CRITERIA packing criteria to be used
    CollectorAllPossibleDomains::CollectorAllPossibleDomains
    (
      const math::Range< size_t> DOMAIN_SIZE,
      const bool CONNECTED_DOMAIN,
      const util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > &PACKING_CRITERIA
    ) :
      m_DomainSize( DOMAIN_SIZE),
      m_ForceDomainConnectivity( CONNECTED_DOMAIN),
      m_PackingCriteria( *PACKING_CRITERIA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorAllPossibleDomains
    CollectorAllPossibleDomains *CollectorAllPossibleDomains::Clone() const
    {
      return new CollectorAllPossibleDomains( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorAllPossibleDomains::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return packing criteria function
    //! @return packing criteria function
    const util::ShPtr
    <
      math::FunctionInterfaceSerializable< SSEGeometryPacking, bool>
    > &CollectorAllPossibleDomains::GetDefaultPackingCriteria()
    {
      // construct the vector of criteria vector
      static const util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > s_packing_criteria
      (
        new SSEGeometryPackingCriteriaCombine
        (
          util::ShPtrVector< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >::Create
          (
            util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >
            (
              new SSEGeometryPackingCriteriaDistancePerType()
            ),
            util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> >
            (
              new SSEGeometryPackingCriteriaInteractionWeight
              (
                0.5, math::Comparisons< double>::GetEnums().e_GreaterEqual
              )
            )
          )
        )
      );

      // end
      return s_packing_criteria;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorAllPossibleDomains::GetAlias() const
    {
      static const std::string s_alias( "CollectorAllPossibleDomains");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorAllPossibleDomains::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects all possible domains.");
      serializer.AddInitializer
      (
        "domain size",
        "size of domains randomly chosen",
        io::Serialization::GetAgent( &m_DomainSize)
      );
      serializer.AddInitializer
      (
        "force domain connectivity",
        "whether remaining SSEs should be considered a domain",
        io::Serialization::GetAgent( &m_ForceDomainConnectivity)
      );
      serializer.AddInitializer
      (
        "packing criterion",
        "criterion when the optimization will be terminated",
        io::Serialization::GetAgent( &m_PackingCriteria)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Collect returns all domains in the domain argument
    //! @param SSE_DOMAIN domain from which Domains will be collected
    //! @return returns ShPtrVector of domains found in the given domain
    util::ShPtrVector< Domain> CollectorAllPossibleDomains::Collect( const DomainInterface &SSE_DOMAIN) const
    {
      // get the SSEs
      const util::SiPtrVector< const SSE> sses( SSE_DOMAIN.GetSSEs());

      // if the domain is empty
      if( sses.IsEmpty())
      {
        // return empty vector
        return util::ShPtrVector< Domain>();
      }

      // adjust the potential domain size based on the size of SSE_DOMAIN
      const math::Range< size_t> adjusted_range
      (
        std::min( m_DomainSize.GetMin(), sses.GetSize()),
        std::min( m_DomainSize.GetMax(), sses.GetSize())
      );

      // get the random domain size
      const size_t domain_size( random::GetGlobalRandom().SizeT( adjusted_range));

      // construct a set from the SSEs
      const storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> sse_set( sses.Begin(), sses.End());

      // construct a combination representing the SSEs that are in the domain but not the domain
      const math::Combination< util::SiPtr< const SSE>, SSELessThanNoOverlap> topology_combination
      (
        sse_set,
        sses.GetSize() - domain_size
      );

      // get all possible combinations
      const storage::List< storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> > combinations
      (
        topology_combination.GetAllCombinations()
      );

      // construct a topology from the SSEs
      const Topology complete_topology
      (
        Topology::BuildTopologyGraphFromGeometries
        (
          sses,
          GetSSEGeometryPackingPickers().e_BestInteractionWeight,
          *m_PackingCriteria
        )
      );

      // initialize domain vectors
      util::ShPtrVector< Domain> subdomains;

      // iterate through the list of combinations
      for
      (
        storage::List< storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> >::const_iterator
          combination_itr( combinations.Begin()), combination_itr_end( combinations.End());
        combination_itr != combination_itr_end; ++combination_itr
      )
      {
        // make a hardcopy of the complete graph
        util::ShPtr< Topology::GraphType> sp_domain_graph( complete_topology.GetGraph().HardCopy());
        storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> domain_set( sse_set);

        // iterate through the SSEs in the set
        for
        (
          storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap>::const_iterator
            set_itr( combination_itr->Begin()), set_itr_end( combination_itr->End());
          set_itr != set_itr_end; ++set_itr
        )
        {
          // remove the corresponding vertex
          sp_domain_graph->DeleteVertex
          (
            sp_domain_graph->FindVertex( util::SiPtr< const SSEGeometryInterface>( *set_itr))
          );
          domain_set.Erase( *set_itr);
        }

        // if the graph has no unconnected vertices
        if( sp_domain_graph->IsConnected() && !sp_domain_graph->GetVertices().IsEmpty())
        {
          // if the domain needs to be connected
          if( m_ForceDomainConnectivity)
          {
            // make a hardcopy of the complete graph
            util::ShPtr< Topology::GraphType> sp_protein_graph( complete_topology.GetGraph().HardCopy());

            // iterate through the SSEs in the set
            for
            (
              storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap>::const_iterator
                set_itr( domain_set.Begin()), set_itr_end( domain_set.End());
              set_itr != set_itr_end; ++set_itr
            )
            {
              // remove the corresponding vertex
              sp_protein_graph->DeleteVertex
              (
                sp_protein_graph->FindVertex( util::SiPtr< const SSEGeometryInterface>( *set_itr))
              );
            }

            // move to the next combination if the model is not complete
            if( !sp_protein_graph->IsConnected())
            {
              continue;
            }
          }

          // convert the SiPtr set to a ShPtr set
          storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> sp_sse_set;
          for
          (
            storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap>::const_iterator
              sse_itr( domain_set.Begin()), sse_itr_end( domain_set.End());
            sse_itr != sse_itr_end; ++sse_itr
          )
          {
            sp_sse_set.Insert( util::ShPtr< SSE>( ( *sse_itr)->Clone()));
          }

          // add the domain to the vector
          subdomains.PushBack( util::ShPtr< Domain>( new Domain( sp_sse_set)));
        }
      }

      // end
      return subdomains;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
