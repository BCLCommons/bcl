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
#include "descriptor/bcl_descriptor_molecule_rings.h"

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

    //! @brief constructor from number of steps, and mapped atom property
    MoleculeRings::MoleculeRings
    (
      const chemistry::ConstitutionalBondTypeData::Conjugation &CONJUGATION,
      const bool &MACROCYCLES_ONLY
    )
      :
      m_Conjugation( CONJUGATION),
      m_MacrocyclesOnly( MACROCYCLES_ONLY)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeRings
    MoleculeRings *MoleculeRings::Clone() const
    {
      return new MoleculeRings( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeRings::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name of the property, which is dependent on the conjugation
    storage::Vector< std::string> MoleculeRingNames()
    {
      storage::Vector< std::string> strings( chemistry::ConstitutionalBondTypeData::ConjugationEnum::GetStringVector());
      strings.PushBack( chemistry::ConstitutionalBondTypeData::ConjugationEnum( chemistry::ConstitutionalBondTypeData::s_NumberOfConjugations));
      storage::Vector< std::string> names( 2 * strings.GetSize());
      for
      (
        storage::Vector< std::string>::iterator itr( strings.Begin()), itr_name( names.Begin()), itr_end( strings.End());
        itr != itr_end;
        ++itr, ++itr_name
      )
      {
//        *itr = "N" + *itr + "Rings";
        *( itr_name++) = "N" + *itr + "Rings";
        *itr_name = "N" + *itr + "MacrocyclicRings";
      }
      // When any ring is sought, just call the property NRings
      names( 2 * chemistry::ConstitutionalBondTypeData::e_Any) = "NRings";
      names( 2 * chemistry::ConstitutionalBondTypeData::e_Any + 1) = "NMacrocyclicRings";
      return names;
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeRings::GetAlias() const
    {
      static storage::Vector< std::string> s_properties( MoleculeRingNames());
      return s_properties( 2 * size_t( m_Conjugation) + size_t( m_MacrocyclesOnly));
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeRings::Calculate( linal::VectorReference< float> &STORAGE)
    {
      float number_rings( 0.0);
      chemistry::ConformationGraphConverter graph_maker
      (
        chemistry::ConformationGraphConverter::e_AtomType,
        chemistry::ConfigurationalBondTypeData::e_ConfigurationalBondType
      );

      util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());
      // create a graph of the small molecule
      const graph::ConstGraph< size_t, size_t> graph( graph_maker( *conformation));

      // find its rings
      storage::List< graph::Ring> rings( graph::EdgeCoverRingPerception( graph).GetRings());

      // if conjugation was unspecified, just return the number of rings
      if( m_Conjugation == chemistry::ConstitutionalBondTypeData::e_Any)
      {
        if( !m_MacrocyclesOnly)
        {
          STORAGE = rings.GetSize();
        }
        else
        {
          size_t n_macrocycles( 0);
          for
          (
            storage::List< graph::Ring>::const_iterator itr_ring( rings.Begin()), itr_ring_end( rings.End());
            itr_ring != itr_ring_end;
            ++itr_ring
          )
          {
            if( itr_ring->GetSize() > 8)
            {
              ++n_macrocycles;
            }
          }
          STORAGE = n_macrocycles;
        }
        return;
      }

      // determine how many aromatic, conjugated, and non-conjugated rings there are, all in one loop
      size_t number_aromatic_rings( 0);
      size_t number_conjugated_rings( 0);
      size_t number_non_conjugated_rings( 0);
      size_t number_macrocyclic_aromatic_rings( 0);
      size_t number_macrocyclic_conjugated_rings( 0);
      size_t number_macrocyclic_non_conjugated_rings( 0);

      for
      (
        storage::List< graph::Ring>::iterator itr_ring( rings.Begin()), itr_ring_end( rings.End());
        itr_ring != itr_ring_end;
        ++itr_ring
      )
      {
        // count how many aromatic and non-conjugated bonds there are in the ring
        size_t number_aromatic_bonds( 0);
        size_t number_non_conjugated_bonds( 0);
        for
        (
          graph::Ring::const_iterator itr( ++itr_ring->Begin()), itr_prev( itr_ring->Begin()), itr_end( itr_ring->End());
          itr != itr_end;
          ++itr, ++itr_prev
        )
        {
          // get the bond type index between the atoms indicated by itr and itr_prev
          const size_t bond_type_index( graph.GetEdgeData( *itr, *itr_prev));
          chemistry::ConfigurationalBondType bond_type( bond_type_index);

          // tally up aromatic and non-conjugated bonds (mutually exclusive)
          if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic)
          {
            ++number_aromatic_bonds;
          }
          else if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Nonconjugated)
          {
            ++number_non_conjugated_bonds;
          }
        }

        const size_t bond_type_index( graph.GetEdgeData( itr_ring->FirstElement(), itr_ring->LastElement()));
        chemistry::ConfigurationalBondType bond_type( bond_type_index);
        // tally up aromatic and non-conjugated bonds (mutually exclusive)
        if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic)
        {
          ++number_aromatic_bonds;
        }
        else if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Nonconjugated)
        {
          ++number_non_conjugated_bonds;
        }
        const size_t number_atoms_in_ring( itr_ring->GetSize());
        if( number_non_conjugated_bonds > size_t( 0))
        {
          ++number_non_conjugated_rings;
          if( number_atoms_in_ring > 8) // Macrocycle > 8 atoms
          {
            ++number_macrocyclic_non_conjugated_rings;
          }
        }
        else if( number_aromatic_bonds == number_atoms_in_ring)
        {
          ++number_aromatic_rings;
          if( number_atoms_in_ring > 8) // Macrocycle > 8 atoms
          {
            ++number_macrocyclic_aromatic_rings;
          }
        }
        else
        {
          ++number_conjugated_rings;
          if( number_atoms_in_ring > 8) // Macrocycle > 8 atoms
          {
            ++number_macrocyclic_conjugated_rings;
          }
        }
      }

      // set number_rings to the number of rings with the selected conjugation
      if( m_Conjugation == chemistry::ConstitutionalBondTypeData::e_Aromatic)
      {
        if( !m_MacrocyclesOnly)
        {
          number_rings = number_aromatic_rings;
        }
        else
        {
          number_rings = number_macrocyclic_aromatic_rings;
        }
      }
      else if( m_Conjugation == chemistry::ConstitutionalBondTypeData::e_Conjugated)
      {
        if( !m_MacrocyclesOnly)
        {
          number_rings = number_conjugated_rings;
        }
        else
        {
          number_rings = number_macrocyclic_conjugated_rings;
        }
      }
      else if( m_Conjugation == chemistry::ConstitutionalBondTypeData::e_Nonconjugated)
      {
        if( !m_MacrocyclesOnly)
        {
          number_rings = number_non_conjugated_rings;
        }
        else
        {
          number_rings = number_macrocyclic_non_conjugated_rings;
        }
      }
      STORAGE = number_rings;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeRings::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Number of"
        + ( m_Conjugation != chemistry::ConstitutionalBondTypeData::e_Any
            ? " " + m_Conjugation.GetString()
            : std::string()
          )
        + ( m_MacrocyclesOnly ? " macrocyclic (>8 atoms)" : std::string())
        + " rings in the molecule"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
