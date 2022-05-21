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
#include "chemistry/bcl_chemistry_conformation_comparison_by_substructure.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_symmetry_rmsd.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // anonymous namespace to prevent namespace polution
    namespace
    {
      //! @brief add all instances of this class to the enumerated list
      util::SiPtr< const util::ObjectInterface> AddInstances()
      {
        util::Enumerated< ConformationComparisonInterface>::AddInstance
        (
          new ConformationComparisonBySubstructure( false, graph::CommonSubgraphIsomorphismBase::e_Unconnected)
        );
        util::Enumerated< ConformationComparisonInterface>::AddInstance
        (
          new ConformationComparisonBySubstructure( true, graph::CommonSubgraphIsomorphismBase::e_Unconnected)
        );
        util::Enumerated< ConformationComparisonInterface>::AddInstance
        (
          new ConformationComparisonBySubstructure( false, graph::CommonSubgraphIsomorphismBase::e_Connected)
        );
        util::Enumerated< ConformationComparisonInterface>::AddInstance
        (
          new ConformationComparisonBySubstructure( true, graph::CommonSubgraphIsomorphismBase::e_Connected)
        );
        util::Enumerated< ConformationComparisonInterface>::AddInstance
        (
          new ConformationComparisonBySubstructure( ConformationComparisonBySymmetryRmsd(), graph::CommonSubgraphIsomorphismBase::e_Connected)
        );
        return *util::Enumerated< ConformationComparisonInterface>::AddInstance
        (
          new ConformationComparisonBySubstructure( ConformationComparisonBySymmetryRmsd(), graph::CommonSubgraphIsomorphismBase::e_Unconnected)
        );
      }
    }

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonBySubstructure::s_Instance( AddInstances());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, takes whether to compute tanimoto vs raw values, and whether to enforce a connected solution
    //! @param TANIMOTO whether to compute tanimoto values (vs. raw size of substructure)
    //! @param CONNECTED whether to enforce a connected substructure solution
    ConformationComparisonBySubstructure::ConformationComparisonBySubstructure
    (
      const bool &TANIMOTO,
      const graph::CommonSubgraphIsomorphismBase::SolutionType &TYPE
    ) :
      m_LowerBound( 0.0),
      m_BondComparisonType( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness),
      m_AtomComparisonType( ConformationGraphConverter::e_ElementType),
      m_SolutionType( TYPE),
      m_ComputeTanimoto( TANIMOTO),
      m_CompareH( false),
      m_InteriorWeighted( m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Connected ? false : true),
      m_Exhaustive( false)
    {
    }

    //! @param COMPARISON This metric will be used to compare the resulting substructures rather
    //!        than tanimoto or raw distance
    //! @param TYPE whether to enforce a connected substructure solution
    ConformationComparisonBySubstructure::ConformationComparisonBySubstructure
    (
      const ConformationComparisonInterface &COMPARISON,
      const graph::CommonSubgraphIsomorphismBase::SolutionType &TYPE
    ) :
      m_LowerBound( -std::numeric_limits< double>::max()),
      m_BondComparisonType( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness),
      m_AtomComparisonType( ConformationGraphConverter::e_ElementType),
      m_SolutionType( TYPE),
      m_ComputeTanimoto( true),
      m_CompareH( false),
      m_Metric( COMPARISON),
      m_InteriorWeighted( m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Connected ? false : true),
      m_Exhaustive( false)
    {
    }

    //! virtual copy constructor
    ConformationComparisonBySubstructure *ConformationComparisonBySubstructure::Clone() const
    {
      return new ConformationComparisonBySubstructure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonBySubstructure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonBySubstructure::GetAlias() const
    {
      static const std::string s_tanimoto_connected( "LargestCommonSubstructureTanimoto"),
                               s_tanimoto_unconnected( "LargestCommonDisconnectedSubstructureTanimoto"),
                               s_raw_connected( "LargestCommonSubstructureSize"),
                               s_raw_unconnected( "LargestCommonDisconnectedSubstructureSize"),
                               s_wrapper_connected( "LargestCommonSubstructureMetric"),
                               s_wrapper_unconnected( "LargestCommonDisconnectedSubstructureMetric");

      return !m_Metric.IsDefined()
             ? (
               m_ComputeTanimoto
               ? (
                   m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Connected
                   ? s_tanimoto_connected
                   : s_tanimoto_unconnected
                 )
               : (
                   m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Connected
                   ? s_raw_connected
                   : s_raw_unconnected
                 )
            )
            : (
                   m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Connected
                   ? s_wrapper_connected
                   : s_wrapper_unconnected
              );
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return the RMSD between first and second molecule
    double ConformationComparisonBySubstructure::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      // get graphs for the molecules
      storage::Map< util::SiPtr< const ConformationInterface>, graph::ConstGraph< size_t, size_t> >::iterator
        itr_graph_a( m_Graphs.Find( MOLECULE_A)), itr_graph_b( m_Graphs.Find( MOLECULE_B));

      // check that the graphs were found
      if( itr_graph_a == m_Graphs.End())
      {
        BCL_MessageDbg( "Calling PrepareToCompare inside operator(), this is not thread safe!");
        m_Mutex.Lock();
        Prepare( MOLECULE_A);
        m_Mutex.Unlock();
        itr_graph_a = m_Graphs.Find( MOLECULE_A);
      }
      if( itr_graph_b == m_Graphs.End())
      {
        BCL_MessageDbg( "Calling PrepareToCompare inside operator(), this is not thread safe!");
        m_Mutex.Lock();
        Prepare( MOLECULE_B);
        m_Mutex.Unlock();
        itr_graph_b = m_Graphs.Find( MOLECULE_B);
      }
      BCL_Assert
      (
        itr_graph_a != m_Graphs.End() && itr_graph_b != m_Graphs.End(),
        "PrepareToCompare was not called with all molecules of interest"
      );

      // acquire an isomorphism
      storage::List< graph::CommonSubgraphIsomorphism< size_t, size_t> >::iterator itr_iso( AcquireIsomorphism());
      graph::CommonSubgraphIsomorphism< size_t, size_t> &isomorphism( *itr_iso);

      // set up the isomorphisms with pointers to the graphs
      // copying the graphs is fairly expensive, so we avoid it whenever possible
      isomorphism.SetGraphA
      (
        util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &itr_graph_a->second, false)
      );
      isomorphism.SetGraphB
      (
        util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &itr_graph_b->second, false)
      );

      // determine sum of # of molecules in the isomorphism
      const double sum_molecule_sizes( isomorphism.GetGraphA().GetSize() + isomorphism.GetGraphB().GetSize());

      // determine the minimum # of common atoms needed to satisfy the tanimoto coefficient, if applicable
      const size_t min_number_in_common
      (
        !m_Metric.IsDefined()
        ? (
            m_ComputeTanimoto
            ? m_LowerBound / ( 1 + m_LowerBound) * sum_molecule_sizes
            : m_LowerBound
          )
        : 1
      );

      // get the actual isomorphism, using estimated upper bounds on its size
      isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds(), min_number_in_common);

      // compute the tanimoto coefficient
      const size_t isomorphism_size( isomorphism.GetIsomorphism().GetSize());
      const size_t interior_vertices
      (
        (
          isomorphism.GetSubgraphIsomorphismsOfGraphA().FirstElement().GetIdsOfInteriorVertices().GetSize()
          + isomorphism.GetSubgraphIsomorphismsOfGraphB().FirstElement().GetIdsOfInteriorVertices().GetSize()
        ) / 2
      );

      if( m_Metric.IsDefined())
      {
        FragmentComplete frag_a
        (
          ConformationGraphConverter::CreateAtomsFromGraph
          (
            ConformationGraphConverter::CreateGraphWithAtoms( MOLECULE_A).GetSubgraph
            (
              !m_InteriorWeighted
              ? isomorphism.GetIsomorphism().GetKeysAsVector()
              : isomorphism.GetSubgraphIsomorphismsOfGraphA().FirstElement().GetIdsOfInteriorVertices()
            ),
            false
          ),
          ""
        );
        FragmentComplete frag_b
        (
          ConformationGraphConverter::CreateAtomsFromGraph
          (
            ConformationGraphConverter::CreateGraphWithAtoms( MOLECULE_B).GetSubgraph
            (
              !m_InteriorWeighted
              ? isomorphism.GetIsomorphism().GetMappedValues()
              : isomorphism.GetSubgraphIsomorphismsOfGraphB().FirstElement().GetIdsOfInteriorVertices()
            ),
            false
          ),
          ""
        );

        // release the isomorphsim
        ReleaseIsomorphism( itr_iso);

        m_Metric->Prepare( frag_a);
        m_Metric->Prepare( frag_b);

        // record the old ignore setting
        const bool old_ignore_setting
        (
          ConformationComparisonInterface::GetIgnoreAtomAndBondTypesWhenDeterminingComparability()
        );
        // set it to ignore, necessary since we're now just comparing parts of molecules
        ConformationComparisonInterface::GetIgnoreAtomAndBondTypesWhenDeterminingComparability() = true;
        // compute the metric
        util::Implementation< ConformationComparisonInterface> cci( m_Metric);
        const double value( cci->operator ()( frag_a, frag_b));
        // restore the original setting
        ConformationComparisonInterface::GetIgnoreAtomAndBondTypesWhenDeterminingComparability() = old_ignore_setting;
        return value;
      }

      // simple tanimoto or distance metric

      // release the isomorphsim
      ReleaseIsomorphism( itr_iso);

      const double metric( m_InteriorWeighted ? interior_vertices : isomorphism_size);

      const double comparison_value
      (
        m_ComputeTanimoto
        ? double( metric) / double( sum_molecule_sizes - metric)
        : metric
      );

      return comparison_value < m_LowerBound ? 0.0 : comparison_value;
    }

    //! @brief prepare the class for comparing a conformation
    //! @param MOLECULE the molecule to prepare to compare
    void ConformationComparisonBySubstructure::Prepare( const ConformationInterface &MOLECULE) const
    {
      // create a graph converter
      ConformationGraphConverter graph_maker( m_AtomComparisonType, m_BondComparisonType, !m_CompareH);

      // make the graph and append it
      m_Graphs[ MOLECULE] = graph_maker( MOLECULE);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonBySubstructure::GetSerializer() const
    {
      io::Serializer parameters;
      const std::string measure( m_ComputeTanimoto ? "the tanimoto coefficient of the " : "");
      const std::string connectivity
      (
        m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Connected
        ? "connected"
        : m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Unconnected
          ? "(possibly disconnected)"
          : "set of connected"
      );
      parameters.SetClassDescription
      (
        "Computes " + measure + " largest common " + connectivity + " substructure between two molecules."
      );

      parameters.AddInitializer
      (
        "atom comparison",
        "atom data that is compared to determine whether atoms are equivalent",
        io::Serialization::GetAgent( &m_AtomComparisonType),
        "ElementType"
      );
      parameters.AddInitializer
      (
        "bond comparison",
        "bond data that is compared",
        io::Serialization::GetAgent( &m_BondComparisonType),
        "BondOrderOrAromaticWithRingness"
      );
      parameters.AddInitializer
      (
        "compare h",
        "whether to compare hydrogens too. "
        "This generally slows calculations down, with little effect on the relative results, but may be desired for some applications",
        io::Serialization::GetAgent( &m_CompareH),
        "False"
      );

      if( m_SolutionType != graph::CommonSubgraphIsomorphismBase::e_Connected)
      {
        parameters.AddInitializer
        (
          "exhaustive",
          "If true, uses the exhaustive unconnected isomorphism search algorithm. See Krissinel, Evgeny, Henrick; "
          "\"Common subgraph isomorphism detection by backtracking search\". This is at least O(N^2) slower on chemical "
          "graphs than the default algorithm (iterative finding of the largest common connected substructure) and results "
          "are typically identical, so this should not normally be used except for testing the above algorithm",
          io::Serialization::GetAgent( &m_Exhaustive),
          "False"
        );
      }
      if( m_Metric.IsDefined())
      {
        parameters.AddInitializer
        (
          "metric",
          "Choose how to compare the largest common substructures. "
          "Note that standard RMSD (with superimposition) should not be used unless it is known that the largest substructures "
          "will have no equivalent vertices or performance is more critical than accuracy",
          io::Serialization::GetAgent( &m_Metric),
          "SymmetryRMSD"
        );
        parameters.AddInitializer
        (
          "interior",
          "If true, isomorphismic component compared by the metric will only include the interior vertices in the isomorphism (those which "
          "have all their neighbors in the isomorphism.) "
          "This is often desired for unconnected isomorphisms, for which many vertices may end up in the isomorphism, "
          ", even though most the connected components in the isomorphism are tiny and lack chemical significance. ",
          io::Serialization::GetAgent( &m_InteriorWeighted),
          m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Connected ? "False" : "True"
        );
      }
      else
      {
        parameters.AddInitializer
        (
          "interior",
          "If true, isomorphisms' size is calculated as the number of interior vertices in the isomorphism (those which "
          "have all their neighbors in the isomorphism.) "
          "This is often desired for unconnected isomorphisms, for which many vertices may end up in the isomorphism, "
          ", even though most the connected components in the isomorphism are tiny and lack chemical significance. ",
          io::Serialization::GetAgent( &m_InteriorWeighted),
          m_SolutionType == graph::CommonSubgraphIsomorphismBase::e_Connected ? "False" : "True"
        );
        // the min parameter cannot be used with an arbitrary metric as there is no way to know how many atoms must be
        // matched to satisfy the criterion
        parameters.AddInitializer
        (
          "min",
          "lowest " + std::string( m_ComputeTanimoto ? "tanimoto coefficient of interest" : "number of atoms to consider")
          + "; > 0 -> faster calculation, but less detailed comparison of weakly-similar structures",
          io::Serialization::GetAgentWithMin( &m_LowerBound, 0.0),
          "0.0"
        );
      }
      return parameters;
    }

    //! @brief acquire an isomorphism object
    //! @return iterator to the isomorphism object
    storage::List< graph::CommonSubgraphIsomorphism< size_t, size_t> >::iterator
      ConformationComparisonBySubstructure::AcquireIsomorphism() const
    {
      m_Mutex.Lock();
      // all isomorphisms have been allocated
      if( m_Available.IsEmpty())
      {
        graph::CommonSubgraphIsomorphismBase::SolutionType effective_type( m_SolutionType);
        if( !m_Exhaustive && effective_type != graph::CommonSubgraphIsomorphismBase::e_Connected)
        {
          effective_type = graph::CommonSubgraphIsomorphismBase::e_GreedyUnconnected;
        }
        m_Available.PushBack( graph::CommonSubgraphIsomorphism< size_t, size_t>( effective_type));
      }
      // splice the last hidden vector set from the available list onto the allocated list
      storage::List< graph::CommonSubgraphIsomorphism< size_t, size_t> >::iterator avail_itr( m_Available.Begin());
      m_Allocated.InternalData().splice( m_Allocated.Begin(), m_Available.InternalData(), avail_itr);

      // save the iterator to the now-allocated array set
      avail_itr = m_Allocated.Begin();
      m_Mutex.Unlock();
      return avail_itr;
    }

    //! @brief release a given isomorphism object
    //! @param ITR iterator to the isomorphism object
    void ConformationComparisonBySubstructure::ReleaseIsomorphism
    (
      const storage::List< graph::CommonSubgraphIsomorphism< size_t, size_t> >::iterator &ITR
    ) const
    {
      m_Mutex.Lock();
      // splice the iterator from the allocated back onto available
      m_Available.InternalData().splice( m_Available.Begin(), m_Allocated.InternalData(), ITR);
      m_Mutex.Unlock();
    }

  } // namespace chemistry
} // namespace bcl
