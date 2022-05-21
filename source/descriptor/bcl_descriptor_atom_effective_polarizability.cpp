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
#include "descriptor/bcl_descriptor_atom_effective_polarizability.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "descriptor/bcl_descriptor_atom_polarizability.h"
#include "graph/bcl_graph_connectivity.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new AtomEffectivePolarizability
    AtomEffectivePolarizability *AtomEffectivePolarizability::Clone() const
    {
      return new AtomEffectivePolarizability( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomEffectivePolarizability::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomEffectivePolarizability::GetAlias() const
    {
      static const std::string s_name( "Atom_EffectivePolarizability");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomEffectivePolarizability::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // create the graph, if necessary
      if( m_Graph.GetSize() == size_t( 0))
      {
        chemistry::ConformationGraphConverter graph_maker;
        util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());
        m_Graph = graph_maker( *conformation);
      }

      // find the distances to all other vertices
      util::ShPtr< storage::Vector< size_t> > shptr_distances
      (
        graph::Connectivity::DistancesToOtherVertices( m_Graph, ELEMENT.GetPosition())
      );
      // get a reference to the distances vector
      const storage::Vector< size_t> &distances( *shptr_distances);

      STORAGE( 0) = 0.0;

      // get a new iterator
      Iterator< chemistry::AtomConformationalInterface> itr_other_atom( ELEMENT.Begin());

      // walk over all other atoms
      storage::Vector< size_t>::const_iterator itr_distances( distances.Begin());
      // sum up to polarizability for every atom
      float effective_polarizability( 0.0);
      for( ; itr_other_atom.NotAtEnd(); ++itr_other_atom, ++itr_distances)
      {
        if( ELEMENT != itr_other_atom( 0))
        {
          effective_polarizability += AtomPolarizability::GetPolarizability( *itr_other_atom( 0)) / float( 1 << *itr_distances);
        }
      }
      STORAGE( 0) = effective_polarizability;
    } // Recalculate

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AtomEffectivePolarizability::SetObjectHook()
    {
      m_Graph = graph::ConstGraph< size_t, size_t>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomEffectivePolarizability::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "atomic polarizability smoothed over molecule");
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
