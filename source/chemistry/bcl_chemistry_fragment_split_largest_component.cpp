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
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitLargestComponent::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitLargestComponent())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param MIN_SIZE get the minimum size of largetst component that is desired
    FragmentSplitLargestComponent::FragmentSplitLargestComponent( const size_t MIN_SIZE) :
      m_MinSize( MIN_SIZE)
    {
    }

    //! virtual copy constructor
    FragmentSplitLargestComponent *FragmentSplitLargestComponent::Clone() const
    {
      return new FragmentSplitLargestComponent( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitLargestComponent::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitLargestComponent::GetAlias() const
    {
      static const std::string s_name_largest( "Largest");
      return s_name_largest;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitLargestComponent::GetClassDescription() const
    {
      static const std::string s_desc_largest( "Get largest component of molecule");
      return s_desc_largest;
    }

    //! get the minimum size of a component of interest
    const size_t FragmentSplitLargestComponent::GetMinSize() const
    {
      return m_MinSize;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief returns connected components of a graph that is not connected
    //! @param MOLECULE molecule of interest
    //! @param MOLECULE_GRAPH graph of molecule of interest
    //! @return connected components of a graph that is not connected
    storage::List< storage::Vector< size_t> > FragmentSplitLargestComponent::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {
      if( !graph::Connectivity::IsConnected( MOLECULE_GRAPH)) // !connected -> multiple components
      {
        // get components of the graph
        storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( MOLECULE_GRAPH));
        size_t max_component_size( 0);
        storage::List< storage::Vector< size_t> > largest_components;

        // get largest component if is larger than the minimum size
        for
        (
          storage::List< storage::Vector< size_t> >::const_iterator
            itr( components.Begin()), itr_end( components.End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->GetSize() > max_component_size)
          {
            max_component_size = itr->GetSize();
            largest_components.Reset();
            largest_components.PushBack( *itr);
          }
        }
        if( max_component_size > GetMinSize())
        {
          return largest_components;
        }
      }
      else if( MOLECULE.GetNumberAtoms() > GetMinSize())
      {
        linal::Vector< size_t> molecule_atoms( linal::FillVector< size_t>( MOLECULE.GetNumberAtoms(), 0, 1));
        return storage::List< storage::Vector< size_t> >
        (
          1,
          storage::Vector< size_t>( molecule_atoms.Begin(), molecule_atoms.End())
        );
      }
      return storage::List< storage::Vector< size_t> >();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitLargestComponent::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "splits a macromolecule and returns largest componenet");
      return parameters;
    }

    //! @brief helper function that can be used to convert the output of GetComponentVertices into an ensemble
    //! @param MOLECULE the molecule to split
    //! @param COMPONENTS list of components to create
    //! @param MOLECULE_GRAPH graph of the molecule with atoms
    //! @return a fragment ensemble
    FragmentEnsemble FragmentSplitLargestComponent::ConvertComponentsIntoEnsemble
    (
      const ConformationInterface &MOLECULE,
      const storage::List< storage::Vector< size_t> > &COMPONENTS,
      const ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {
      // ensemble to store components of molecule of interest
      FragmentEnsemble ensemble_of_components;

      // unlike the base classes version, this function preserves molecule properties, since there remains a 1-1
      // relationship between the original and new molecule
      // add all rings to the rings ensemble
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator
          itr( COMPONENTS.Begin()), itr_end( COMPONENTS.End());
        itr != itr_end;
        ++itr
      )
      {
        // if size of component is greater than minimum size then add to ensemble
        if( itr->GetSize() > GetMinSize())
        {
          FragmentComplete new_fragment
          (
            ConformationGraphConverter::CreateAtomsFromGraph( MOLECULE_GRAPH.GetSubgraph( *itr)),
            GetAlias() + " isolated from " + MOLECULE.GetName(),
            MOLECULE.GetStoredProperties().GetMDLProperties()
          );
          ensemble_of_components.PushBack( new_fragment);
        }
      }
      return ensemble_of_components;
    }

  } // namespace chemistry
} // namespace bcl
