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
#include "descriptor/bcl_descriptor_atom_ring_size.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a bool
    //! @param LARGEST whether to prefer the larger of the rings that the atom is connected to,
    //!        if it is connected to multiple
    AtomRingSize::AtomRingSize ( const bool &LARGEST) :
      m_PreferLargest( LARGEST)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new AtomRingSize
    AtomRingSize *AtomRingSize::Clone() const
    {
      return new AtomRingSize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomRingSize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomRingSize::GetAlias() const
    {
      static const std::string s_max_ring_size( "AtomMaxRingSize"), s_min_ring_size( "AtomMinRingSize");
      return m_PreferLargest ? s_max_ring_size : s_min_ring_size;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomRingSize::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // initialize the storage with the default value for atoms in chains
      const size_t default_value( m_PreferLargest ? s_ChainAtomsMaxRingSize : s_ChainAtomsMinRingSize);
      STORAGE( 0) = default_value;

      // determine whether this atom is even part of a ring
      const size_t number_ring_bonds
      (
        ELEMENT->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
      );
      if( number_ring_bonds == size_t( 0))
      {
        // trivial case, atom is not in a ring
        return;
      }
      if( m_AtomRingSizes.IsEmpty())
      {
        chemistry::ConformationGraphConverter graph_maker
        (
          chemistry::ConformationGraphConverter::e_AtomType,
          chemistry::ConfigurationalBondTypeData::e_IsInRing
        );

        util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());
        // create a graph of the small molecule
        const graph::ConstGraph< size_t, size_t> graph( graph_maker( *conformation));

        // find its rings
        storage::List< graph::Ring> rings( graph::EdgeCoverRingPerception( graph).GetRings());
        m_AtomRingSizes.Resize( graph.GetSize(), default_value);
        for
        (
          storage::List< graph::Ring>::const_iterator itr( rings.Begin()), itr_end( rings.End());
          itr != itr_end;
          ++itr
        )
        {
          // walk through all vertices in the ring, update the atoms ring size according to whether the largest
          // or smallest value is preferred
          const size_t ring_size( itr->GetSize());
          for
          (
            graph::Ring::const_iterator itr_ring( itr->Begin()), itr_ring_end( itr->End());
            itr_ring != itr_ring_end;
            ++itr_ring
          )
          {
            // ring contains this atom, test its size
            if( m_PreferLargest)
            {
              m_AtomRingSizes( *itr_ring) = std::max( m_AtomRingSizes( *itr_ring), ring_size);
            }
            else
            {
              m_AtomRingSizes( *itr_ring) = std::min( m_AtomRingSizes( *itr_ring), ring_size);
            }
          }
        }
      }
      // get the size out of the vector
      STORAGE( 0) = m_AtomRingSizes( ELEMENT.GetPosition());
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AtomRingSize::SetObjectHook()
    {
      m_AtomRingSizes.Reset();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomRingSize::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        std::string( m_PreferLargest ? "Maximum " : "Minimum ") +
        "size of a ring that this atom is part of. For atoms that are not in a ring, returns " +
        util::Format()( size_t( m_PreferLargest ? s_ChainAtomsMaxRingSize : s_ChainAtomsMinRingSize))
        + ". " + std::string
        (
          m_PreferLargest
          ? "The size of the largest ring with no internal rings is returned for atoms that are part of ring systems"
          : ""
        )
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
