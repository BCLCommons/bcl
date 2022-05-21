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
#include "chemistry/bcl_chemistry_reaction_complete.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "sdf/bcl_sdf_rxn_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ReactionComplete::s_Instance
    (
      GetObjectInstances().AddInstance( new ReactionComplete())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ReactionComplete::ReactionComplete() :
      m_Description(),
      m_NumberReactants( 0),
      m_NumberProducts( 0),
      m_ReactantStructures(),
      m_ProductStructures(),
      m_ValidAtomMapping( false),
      m_ReactiveAtomsReactants(),
      m_ReactiveAtomsProducts()
    {
    }

    //! Complete constructor
    //! @param REACTANTS the reactants
    //! @param PRODUCTS the products
    //! @param REACTIVE_ATOMS_REACTANTS a mapping of reactive atoms to atoms in reactants
    //! @param REACTIVE_ATOMS_PRODUCTS a mapping of reactive atoms to atoms in products
    ReactionComplete::ReactionComplete
    (
      const storage::Vector< ReactionStructure> &REACTANTS,
      const storage::Vector< ReactionStructure> &PRODUCTS,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTIVE_ATOMS_REACTANTS,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTIVE_ATOMS_PRODUCTS,
      const std::string &DESCRIPTION
    ) :
      m_Description( DESCRIPTION),
      m_NumberReactants( REACTANTS.GetSize()),
      m_NumberProducts( PRODUCTS.GetSize()),
      m_ReactantStructures( REACTANTS),
      m_ProductStructures( PRODUCTS),
      m_ValidAtomMapping( true),
      m_ReactiveAtomsReactants(),
      m_ReactiveAtomsProducts()
    {

      if( m_ReactantStructures.IsEmpty())
      {
        BCL_MessageCrt( "Reaction is missing reactants; clearing all data");
        m_ProductStructures.Reset();
        m_ValidAtomMapping = false;
        m_ReactiveAtomsReactants.Reset();
        m_ReactiveAtomsProducts.Reset();
      }

      SetAtomMappings( REACTIVE_ATOMS_REACTANTS, REACTIVE_ATOMS_PRODUCTS);

    }

    void ReactionComplete::SetAtomMappings
    (
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTIVE_ATOMS_REACTANTS,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTIVE_ATOMS_PRODUCTS
    )
    {
      m_ReactiveAtomsReactants = REACTIVE_ATOMS_REACTANTS;
      m_ReactiveAtomsProducts = REACTIVE_ATOMS_PRODUCTS;

      m_ValidAtomMapping = true;
      
      // Make sure that each reactive atom value is associated with an atom in both reactants and products
      size_t n_prod_keys( m_ReactiveAtomsProducts.GetSize());
      size_t n_react_keys( m_ReactiveAtomsReactants.GetSize());
      if( n_react_keys != n_prod_keys)
      {
        BCL_MessageStd( "Different numbers of reactive atoms were specified for reactants and products");
        m_ValidAtomMapping = false;
      }

      for
      ( 
        storage::Map< size_t, storage::Pair< size_t, size_t> >::const_iterator 
          itr_react_map( m_ReactiveAtomsReactants.Begin()), itr_react_map_end( m_ReactiveAtomsReactants.End());
        m_ValidAtomMapping && itr_react_map != itr_react_map_end;
        ++itr_react_map
      )
      {
        if( !m_ReactiveAtomsProducts.Has( itr_react_map->first))
        {
          BCL_MessageStd( "Different reactive atoms were specified for reactants and products");
          m_ValidAtomMapping = false;
          break;
        }
        // TODO add atom index checking
      }
    }

    //! @brief Clone function
    //! @return pointer to new ReactionComplete
    ReactionComplete *ReactionComplete::Clone() const
    {
      return new ReactionComplete( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ReactionComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a description of the reaction
    //! @return the reaction description
    const std::string &ReactionComplete::GetDescription() const
    {
      return m_Description;
    }
      
    //! @brief checks if the reaction has a valid map between reactant and product atoms
    //! @return true if the reaction has a valid (one-to-one) reactive atom mapping
    const bool &ReactionComplete::HasValidAtomMap() const
    {
      return m_ValidAtomMapping;
    }

    //! @brief get the number of atoms involved in the reaction
    size_t ReactionComplete::GetNumberReactiveAtoms() const
    {
      return m_ReactiveAtomsReactants.GetSize();
    }

    //! @brief return the number of atoms
    //! @return the number of atoms
    const size_t &ReactionComplete::GetNumberReactants() const
    {
      return m_NumberReactants;
    }

    //! @brief gets all reactants
    //! @return a vector of reactants
    //const storage::Vector< FragmentComplete> &ReactionComplete::GetReactants() const
    //{
    //  return m_Reactants;
    //}

    //! @brief gets a reactant with an index
    //! @param INDEX the reactant to fetch
    const FragmentComplete &ReactionComplete::GetReactant( const size_t &INDEX) const
    {
      BCL_Assert
      (
        INDEX < m_NumberReactants,
        "Reactant index is outside the range of reactants (" + util::Format()( INDEX)
          + ">" + util::Format()( m_NumberReactants - 1) + ")"
      );

      //return m_Reactants( INDEX);
      return m_ReactantStructures( INDEX).GetFragment();
    }

    //! @brief gets which reactive atoms are present in a given product
    //! @param INDEX the product index of interest
    //! @return a map of reactive atoms (keys) to reactant atom indices (values)
    storage::Map< size_t, size_t> ReactionComplete::GetReactiveAtomsInReactant( const size_t &INDEX) const
    {
      storage::Map< size_t, size_t> atom_mapping;
      BCL_Assert
      (
        INDEX < m_NumberReactants,
        "Requested reactant index is outside available range (" + util::Format()( INDEX)
        + " > " + util::Format()( m_NumberReactants - 1) + ")"
      );

      // Iterate through the reactive atoms map; any atom that is in the desired reactant should be stored
      for
      (
        storage::Map< size_t, storage::Pair< size_t, size_t> >::const_iterator itr_map( m_ReactiveAtomsReactants.Begin()), itr_map_end( m_ReactiveAtomsReactants.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->second.First() == INDEX)
        {
          atom_mapping[ itr_map->first] = itr_map->second.Second();
        }
      }
      return atom_mapping;
    }

    //! @brief return the number of atoms
    //! @return the number of atoms
    const size_t &ReactionComplete::GetNumberProducts() const
    {
      return m_NumberProducts;
    }

    //! @brief gets a product with an index
    //! @param INDEX the product to fetch
    const FragmentComplete &ReactionComplete::GetProduct( const size_t &INDEX) const
    {
      BCL_Assert
      (
        INDEX < m_NumberProducts,
        "Product index is outside the range of products (" + util::Format()( INDEX)
          + ">" + util::Format()( m_NumberProducts - 1) + ")"
      );

      //return m_Products( INDEX);
      return m_ProductStructures( INDEX).GetFragment();
    }

    //! @brief gets which reactive atoms are present in a given product
    //! @param INDEX the product index of interest
    //! @return a map of reactive atoms (keys) to reactant atom indices (values)
    storage::Map< size_t, size_t> ReactionComplete::GetReactiveAtomsInProduct( const size_t &INDEX) const
    {
      storage::Map< size_t, size_t> atom_mapping;
      BCL_Assert
      (
        INDEX < m_NumberProducts,
        "Requested product index is outside available range (" + util::Format()( INDEX)
        + " > " + util::Format()( m_NumberProducts - 1)
      );

      // Iterate through the reactive atoms map; any atom that is in the desired product should be stored
      for
      (
        storage::Map< size_t, storage::Pair< size_t, size_t> >::const_iterator itr_map( m_ReactiveAtomsProducts.Begin()), itr_map_end( m_ReactiveAtomsProducts.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->second.First() == INDEX)
        {
          atom_mapping[ itr_map->first] = itr_map->second.Second();
        }
      }
      return atom_mapping;
    }

    //! @brief reverses the reaction, i.e. makes the products the reactants and vice versa
    void ReactionComplete::Reverse()
    {
      std::swap( m_ReactantStructures, m_ProductStructures);
      std::swap( m_ReactiveAtomsReactants, m_ReactiveAtomsProducts);
      std::swap( m_NumberReactants, m_NumberProducts);
    }

    //! @brief equality operator; reactants/products must have the same configuration and be listed in the same order
    //! @param OTHER another ReactionComplete to compare
    //! @return true if the reactions are equivalent
    bool ReactionComplete::operator ==( const ReactionComplete &OTHER) const
    {
      // TODO: Switch these to check graph equality rather than configuration equality
      for( size_t r( 0); r < m_NumberReactants; ++r)
      {
        ConfigurationSet config_set;
        config_set.Insert( FragmentConfigurationShared( m_ReactantStructures( r).GetFragment()));
        if( config_set.Insert( FragmentConfigurationShared( OTHER.m_ReactantStructures( r).GetFragment())).second)
        {
          return false;
        }
      }
      for( size_t p( 0); p < m_NumberProducts; ++p)
      {
        ConfigurationSet config_set;
        config_set.Insert( FragmentConfigurationShared( m_ProductStructures( p).GetFragment()));
        if( config_set.Insert( FragmentConfigurationShared( OTHER.m_ProductStructures( p).GetFragment())).second)
        {
          return false;
        }
      }
      return m_ReactiveAtomsReactants == OTHER.m_ReactiveAtomsReactants
        && m_ReactiveAtomsProducts == OTHER.m_ReactiveAtomsProducts;
    }
      
    //! @brief write the reaction in RXN format to an output stream
    //! @param OSTREAM the output stream to write to
    //! @return the same output stream that was written to 
    std::ostream &ReactionComplete::WriteRXN( std::ostream &OSTREAM) const
    {
      storage::Vector< storage::Vector< sdf::AtomInfo> > reactant_atom_infos;
      storage::Vector< storage::Vector< sdf::BondInfo> > reactant_bond_infos;

      for( size_t r( 0); r < m_ReactantStructures.GetSize(); ++r)
      {
        reactant_atom_infos.PushBack( m_ReactantStructures( r).GetFragment().GetAtomInfo());
        reactant_bond_infos.PushBack( m_ReactantStructures( r).GetFragment().GetBondInfo());
      }

      storage::Vector< storage::Vector< sdf::AtomInfo> > product_atom_infos;
      storage::Vector< storage::Vector< sdf::BondInfo> > product_bond_infos;

      for( size_t p( 0); p < m_ReactantStructures.GetSize(); ++p)
      {
        product_atom_infos.PushBack( m_ProductStructures( p).GetFragment().GetAtomInfo());
        product_bond_infos.PushBack( m_ProductStructures( p).GetFragment().GetBondInfo());
      }

      return sdf::RXNHandler::WriteToRXN
             (
               OSTREAM,
               GetDescription(),
               reactant_atom_infos,
               reactant_bond_infos,
               product_atom_infos,
               product_bond_infos,
               m_ReactiveAtomsReactants,
               m_ReactiveAtomsProducts
             );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ReactionComplete::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ReactionComplete::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
