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
#include "chemistry/bcl_chemistry_fragment_split_isolate.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitIsolate::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitIsolate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param MIN_SIZE get the minimum size of isolate that is desired
    FragmentSplitIsolate::FragmentSplitIsolate( const size_t MIN_SIZE) :
      m_MinSize( MIN_SIZE)
    {
    }

    //! virtual copy constructor
    FragmentSplitIsolate *FragmentSplitIsolate::Clone() const
    {
      return new FragmentSplitIsolate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitIsolate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitIsolate::GetAlias() const
    {
      static const std::string s_name_isolate( "Isolate");
      return s_name_isolate;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitIsolate::GetClassDescription() const
    {
      static const std::string s_desc_isolate( "Isolate molecules from complexes");
      return s_desc_isolate;
    }

    //! get the minimum size of a component of interest
    const size_t FragmentSplitIsolate::GetMinSize() const
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
    storage::List< storage::Vector< size_t> > FragmentSplitIsolate::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {
      storage::List< storage::Vector< size_t> > indices;
      if( !graph::Connectivity::IsConnected( MOLECULE_GRAPH)) // !connected -> multiple components
      {
        storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( MOLECULE_GRAPH));
        for
        (
          storage::List< storage::Vector< size_t> >::const_iterator itr( components.Begin()), itr_end( components.End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->GetSize() > GetMinSize())
          {
            indices.PushBack( *itr);
          }
        }
        return indices;
      }
      return indices;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitIsolate::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "splits a macromolecule into component parts");
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
