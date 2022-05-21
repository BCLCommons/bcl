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
#include "descriptor/bcl_descriptor_molecule_complexity.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeComplexity
    MoleculeComplexity *MoleculeComplexity::Clone() const
    {
      return new MoleculeComplexity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeComplexity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeComplexity::GetAlias() const
    {
      static const std::string s_name( "MoleculeComplexity");
      return s_name;
    }
      
    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MoleculeComplexity::GetNormalSizeOfFeatures() const
    {
      return size_t( 1);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeComplexity::Calculate( linal::VectorReference< float> &STORAGE)
    {

      // Number of stereo centers + 1 (for explanation see reference given in class description)
      size_t n_stereocenters( GetCheminfoProperties().calc_NStereo->SumOverObject( *this->GetCurrentObject())( 0) + 1);

      // Number of atoms
      size_t n_atoms( GetCurrentObject()->GetSize());

      chemistry::ConformationGraphConverter graph_maker
      (
        chemistry::ConformationGraphConverter::e_AtomType,
        chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
      );

      // Current molecule
      util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());

      // create a graph of the small molecule
      const graph::ConstGraph< size_t, size_t> graph( graph_maker( *conformation));

      // find its rings
      storage::List< graph::Ring> rings( graph::EdgeCoverRingPerception( graph).GetRings());

      // vector to mark if a vertex is part of multiple rings
      linal::Vector< size_t> bridge_atoms( graph.GetSize(), size_t( 0));

      // Macrocycles (start at 1)
      size_t n_macrocycles( 1);

      // Step through each ring individually and add up how many rings each vertex is in
      for
      (
        storage::List< graph::Ring>::const_iterator itr_ring( rings.Begin()), itr_ring_end( rings.End());
        itr_ring != itr_ring_end;
        ++itr_ring
      )
      {
        // If there are more than 8 atoms in a ring, it is a macrocycle
        if( itr_ring->GetSize() > 8)
        {
          ++n_macrocycles;
        }

        // Compare all unique pairs of rings
        storage::List< graph::Ring>::const_iterator itr_other_ring( itr_ring);
        ++itr_other_ring;
        for( ; itr_other_ring != itr_ring_end; ++itr_other_ring)
        {

          // Iterate through the second ring, check if the two rings share any atoms
          for
          (
            graph::Ring::const_iterator itr_atom( itr_other_ring->Begin()), itr_atom_end( itr_other_ring->End());
            itr_atom != itr_atom_end;
            ++itr_atom
          )
          {
            if( itr_ring->Contains( *itr_atom))
            {
              bridge_atoms( *itr_atom) = 1;
            }
          }
        }
      }

      // Number of bridging atoms
      size_t n_bridge( bridge_atoms.Sum() + 1);

      // Whether an atom is a spirocyclic atom; if an atom has four ring bonds it pretty much
      // has to be in some form of spirocyclic ring or another
      CheminfoProperty is_spirocyclic_atom
      (
        util::ObjectDataLabel
        (
          "GreaterEqual(lhs=BondTypeCount(property=IsInRing,value=1),rhs=4)"
        )
      );

      // Count the number of spirocyclic atoms
      size_t n_spiro( is_spirocyclic_atom->SumOverObject( *this->GetCurrentObject())( 0) + 1);

      // composite score, see ref
      float score
      (
        log10( n_spiro)
        + log10( n_bridge)
        + log10( n_macrocycles)
        + log10( n_stereocenters)
        + float( std::pow( n_atoms, 1.005) - n_atoms)
      );
      STORAGE = score;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeComplexity::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "calculates the complexity of a molecule (see http://www.jcheminf.com/content/1/1/8)");
      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
