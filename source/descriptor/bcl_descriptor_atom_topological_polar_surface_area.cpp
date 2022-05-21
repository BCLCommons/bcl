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
#include "descriptor/bcl_descriptor_atom_topological_polar_surface_area.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
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
    //! @return pointer to new AtomTopologicalPolarSurface
    AtomTopologicalPolarSurfaceArea *AtomTopologicalPolarSurfaceArea::Clone() const
    {
      return new AtomTopologicalPolarSurfaceArea( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomTopologicalPolarSurfaceArea::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomTopologicalPolarSurfaceArea::GetAlias() const
    {
      static const std::string s_name( "Atom_TopologicalPolarSurfaceArea");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomTopologicalPolarSurfaceArea::Calculate
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

      static const storage::Map< std::string, float> s_TPSAMap( MakeTPSATypeStringToPSAMap());

      // tally up the number of bonds of each type
      size_t number_of_single_bonds( 0);
      size_t number_of_double_bonds( 0);
      size_t number_of_triple_bonds( 0);
      size_t number_of_aromatic_bonds( 0);
      size_t number_of_hydrogen_bonds( ELEMENT->GetNumberofValenceBondsWithOrder( 1));
      bool has_ring_bonds( false);
      for
      (
        storage::Vector< chemistry::BondConformational>::const_iterator
          itr_bonds( ELEMENT->GetBonds().Begin()), itr_bonds_end( ELEMENT->GetBonds().End());
        itr_bonds != itr_bonds_end;
        ++itr_bonds
      )
      {
        // keep track of whether any ringed bonds were seen
        has_ring_bonds = has_ring_bonds || itr_bonds->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_IsInRing) == 1;
        // if the bond is aromatic
        if( itr_bonds->GetBondType()->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic)
        {
          ++number_of_aromatic_bonds;
        }
        else if( itr_bonds->GetBondType()->GetNumberOfElectrons() == 2) // if the bond is nominally a single bond
        {
          // is this bond to a hydrogen?
          if( itr_bonds->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Hydrogen) // bonded to hydrogen
          {
            ++number_of_hydrogen_bonds;
          }
          else // bonded to heavy atom
          {
            ++number_of_single_bonds;
          }
        }
        else if( itr_bonds->GetBondType()->GetNumberOfElectrons() == 4) // if the bond is nominally a float  bond
        {
          ++number_of_double_bonds;
        }
        else if( itr_bonds->GetBondType()->GetNumberOfElectrons() == 6) // if the bond is nominally a triple bond
        {
          ++number_of_triple_bonds;
        }
      }

      // test whether the atom was in a 3-membered ring.  This only needs to be checked when ring bonds were found
      const bool is_in_3_membered_ring
      (
        has_ring_bonds
        && graph::Connectivity::IsInThreeMemberedCycle( m_Graph, ELEMENT.GetPosition())
      );

      // make a string encoding all the information we now have about the type
      const std::string psa_string
      (
        GetTPSATypeString
        (
          ELEMENT->GetElementType()->GetChemicalSymbol(),
          number_of_hydrogen_bonds,
          number_of_single_bonds,
          number_of_double_bonds,
          number_of_triple_bonds,
          number_of_aromatic_bonds,
          is_in_3_membered_ring,
          ELEMENT->GetCharge()
        )
      );

      // look for the string in the map
      storage::Map< std::string, float>::const_iterator map_itr( s_TPSAMap.Find( psa_string));

      // if the string was found in the map, set the corresponding place in the vector to this value
      if( map_itr != s_TPSAMap.End())
      {
        STORAGE( 0) = map_itr->second;
      }
      else // otherwise, set the place in the vector to 0.0
      {
        STORAGE( 0) = 0.0;
      }
    } // Recalculate

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AtomTopologicalPolarSurfaceArea::SetObjectHook()
    {
      m_Graph = graph::ConstGraph< size_t, size_t>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @return a unique string for the given connectivity information
    std::string AtomTopologicalPolarSurfaceArea::GetTPSATypeString
    (
      const std::string &SYMBOL,
      const size_t &NUMBER_HYDROGEN_BONDS,
      const size_t &NUMBER_SINGLE_BONDS_TO_HEAVY_ATOMS,
      const size_t &NUMBER_DOUBLE_BONDS,
      const size_t &NUMBER_TRIPLE_BONDS,
      const size_t &NUMBER_AROMATIC_BONDS,
      const bool &IS_IN_3_MEMBERED_RING,
      const short &CHARGE
    )
    {
      std::ostringstream combiner;
      combiner << SYMBOL << ','
               << NUMBER_HYDROGEN_BONDS << ','
               << NUMBER_SINGLE_BONDS_TO_HEAVY_ATOMS << ','
               << NUMBER_DOUBLE_BONDS << ','
               << NUMBER_TRIPLE_BONDS << ','
               << NUMBER_AROMATIC_BONDS << ','
               << IS_IN_3_MEMBERED_RING << ','
               << CHARGE;
      return combiner.str();
    } // GetTPSATypeString

    //! @brief create the connectivity string to topological polar surface area
    storage::Map< std::string, float>
    AtomTopologicalPolarSurfaceArea::MakeTPSATypeStringToPSAMap()
    {
      static storage::Map< std::string, float> map;
      if( map.IsEmpty())
      {
        // compose the map using data from Ertl, et. al. J. Med. Chem. 2000, 43, 3715
        //                     Sym Hy 1x 2x 3x Ar 3R Chrg     PSA
        map[ GetTPSATypeString( "N", 0, 0, 0, 1, 0, 0, 0)] = 23.79;
        map[ GetTPSATypeString( "N", 0, 1, 1, 0, 0, 0, 0)] = 12.36;
        map[ GetTPSATypeString( "N", 0, 0, 0, 0, 2, 0, 0)] = 12.89;
        map[ GetTPSATypeString( "N", 1, 0, 1, 0, 0, 0, 0)] = 23.85;
        map[ GetTPSATypeString( "N", 0, 3, 0, 0, 0, 1, 0)] = 3.01;
        map[ GetTPSATypeString( "N", 0, 3, 0, 0, 0, 0, 0)] = 3.24;
        map[ GetTPSATypeString( "N", 1, 2, 0, 0, 0, 0, 0)] = 12.03;
        map[ GetTPSATypeString( "N", 1, 2, 0, 0, 0, 1, 0)] = 21.94;
        map[ GetTPSATypeString( "N", 2, 1, 0, 0, 0, 0, 0)] = 26.02;
        map[ GetTPSATypeString( "N", 0, 0, 0, 0, 3, 0, 0)] = 4.41;
        map[ GetTPSATypeString( "N", 0, 1, 0, 0, 2, 0, 0)] = 4.93;
        map[ GetTPSATypeString( "N", 1, 0, 0, 0, 2, 0, 0)] = 15.79;
        map[ GetTPSATypeString( "N", 0, 1, 0, 1, 0, 0, 1)] = 4.36;
        map[ GetTPSATypeString( "N", 0, 2, 1, 0, 0, 0, 1)] = 3.01;
        map[ GetTPSATypeString( "N", 0, 1, 0, 0, 2, 0, 1)] = 3.88;
        map[ GetTPSATypeString( "N", 0, 0, 0, 0, 3, 0, 1)] = 4.1;
        map[ GetTPSATypeString( "N", 1, 1, 1, 0, 0, 0, 1)] = 13.97;
        map[ GetTPSATypeString( "N", 1, 0, 0, 0, 2, 0, 1)] = 14.14;
        map[ GetTPSATypeString( "N", 2, 0, 1, 0, 0, 0, 1)] = 25.59;
        map[ GetTPSATypeString( "N", 1, 3, 0, 0, 0, 0, 1)] = 4.44;
        map[ GetTPSATypeString( "N", 2, 2, 0, 0, 0, 0, 1)] = 16.61;
        map[ GetTPSATypeString( "N", 3, 1, 0, 0, 0, 0, 1)] = 27.64;
        map[ GetTPSATypeString( "N", 0, 0, 1, 1, 0, 0, 0)] = 13.6;
        map[ GetTPSATypeString( "N", 0, 0, 1, 0, 2, 0, 0)] = 8.39;
        map[ GetTPSATypeString( "N", 0, 1, 2, 0, 0, 0, 0)] = 11.68;
        map[ GetTPSATypeString( "O", 0, 1, 0, 0, 0, 0, -1)] = 23.0;
        map[ GetTPSATypeString( "O", 0, 0, 1, 0, 0, 0, 0)] = 17.07;
        map[ GetTPSATypeString( "O", 0, 2, 0, 0, 0, 0, 0)] = 9.23;
        map[ GetTPSATypeString( "O", 0, 2, 0, 0, 0, 1, 0)] = 12.53;
        map[ GetTPSATypeString( "O", 1, 1, 0, 0, 0, 0, 0)] = 20.23;
        map[ GetTPSATypeString( "O", 0, 0, 0, 0, 2, 0, 0)] = 13.14;
        map[ GetTPSATypeString( "P", 0, 1, 1, 0, 0, 0, 0)] = 34.14;
        map[ GetTPSATypeString( "P", 0, 3, 0, 0, 0, 0, 0)] = 13.59;
        map[ GetTPSATypeString( "P", 0, 3, 1, 0, 0, 0, 0)] = 9.81;
        map[ GetTPSATypeString( "P", 1, 2, 1, 0, 0, 0, 0)] = 23.47;
        map[ GetTPSATypeString( "S", 0, 2, 0, 0, 0, 0, 0)] = 25.3;
        map[ GetTPSATypeString( "S", 0, 0, 1, 0, 0, 0, 0)] = 32.09;
        map[ GetTPSATypeString( "S", 0, 2, 1, 0, 0, 0, 0)] = 19.21;
        map[ GetTPSATypeString( "S", 0, 2, 2, 0, 0, 0, 0)] = 8.38;
        map[ GetTPSATypeString( "S", 0, 0, 0, 0, 2, 0, 0)] = 28.24;
        map[ GetTPSATypeString( "S", 0, 0, 1, 0, 2, 0, 0)] = 21.7;
        map[ GetTPSATypeString( "S", 1, 1, 0, 0, 0, 0, 0)] = 38.8;
      }

      return map;
    } // MakeTPSATypeStringToPSAMap

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomTopologicalPolarSurfaceArea::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "see Ertl, et. al. J. Med. Chem. 2000, 43, 3715");
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
