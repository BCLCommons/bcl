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

#ifndef BCL_CHEMISTRY_REACTION_COMPLETE_H_
#define BCL_CHEMISTRY_REACTION_COMPLETE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_reaction_structure.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ReactionComplete
    //! @brief Class that contains reaction data
    //! @details Stores multiple reactants and products as FragmentCompletes.  Hypothetically this class stores
    //!          the minimal structures needed to execute a reaction with all preserved atoms mapped (as reactive
    //!          atoms).
    //!
    //! TODO: add the ability to store misc information (e.g. reagents, conditions, etc)
    //!
    //! @see @link example_chemistry_reaction_complete.cpp @endlink
    //! @author geanesar
    //! @date Jan 09, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ReactionComplete :
      public util::ObjectInterface
    {
    public:

    private:

    //////////
    // data //
    //////////

      //! a description of the reaction
      std::string m_Description;

      //! the number of reactants
      size_t m_NumberReactants;

      //! the number of products
      size_t m_NumberProducts;

      //! FragmentCompletes for each reactant
      storage::Vector< ReactionStructure> m_ReactantStructures;

      //! FragmentCompletes for each product
      storage::Vector< ReactionStructure> m_ProductStructures;

      //! whether a valid mapping was given
      bool m_ValidAtomMapping;

      //! a mapping of reactive atoms to their atom identity in the reactants
      //! stored as reactive atom (key) to (reactant index,atom number) (value)
      storage::Map< size_t, storage::Pair< size_t, size_t> > m_ReactiveAtomsReactants;

      //! a mapping of reactive atoms to their atom identity in the products
      //! stored as reactive atom (key) to (product index,atom number) (value)
      storage::Map< size_t, storage::Pair< size_t, size_t> > m_ReactiveAtomsProducts;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ReactionComplete();

      //! Complete constructor
      //! @param REACTANTS the reactants
      //! @param PRODUCTS the products
      //! @param REACTIVE_ATOMS_REACTANTS a mapping of reactive atoms to atoms in reactants
      //! @param REACTIVE_ATOMS_PRODUCTS a mapping of reactive atoms to atoms in products
      ReactionComplete
      (
        const storage::Vector< ReactionStructure> &REACTANTS,
        const storage::Vector< ReactionStructure> &PRODUCTS,
        const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTIVE_ATOMS_REACTANTS,
        const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTIVE_ATOMS_PRODUCTS,
        const std::string &DESCRIPTION = std::string()
      );

      //! @brief Clone function
      //! @return pointer to new ReactionComplete
      ReactionComplete *Clone() const;

    /////////////////
    // data access //
    /////////////////

      void SetAtomMappings
      (
        const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTIVE_ATOMS_REACTANTS,
        const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTIVE_ATOMS_PRODUCTS
      );

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a description of the reaction
      //! @return the reaction description
      const std::string &GetDescription() const;
      
      //! @brief checks if the reaction has a valid, one-to-one, map between reactant and product reactive atoms
      //! @return true if the reaction has a valid reactive atom mapping
      const bool &HasValidAtomMap() const;

      //! @brief get the number of atoms involved in the reaction
      size_t GetNumberReactiveAtoms() const;

      //! @brief return the number of atoms
      //! @return the number of atoms
      const size_t &GetNumberReactants() const;

      //! @brief gets all reactants
      //! @return a vector of reactants
      //const storage::Vector< FragmentComplete> &GetReactants() const;
      const storage::Vector< ReactionStructure> &GetReactantStructures() const
      {
        return m_ReactantStructures;
      }

      storage::Vector< FragmentComplete> GetReactants() const
      {
        storage::Vector< FragmentComplete> mols;
        mols.AllocateMemory( m_ReactantStructures.GetSize());
        
        for( size_t s( 0), end_s( m_ReactantStructures.GetSize()); s < end_s; ++s)
        {
          mols.PushBack( m_ReactantStructures( s).GetFragment());
        }
        return mols;
      }

      //! @brief gets a reactant with an index
      //! @param INDEX the reactant to fetch
      const FragmentComplete &GetReactant( const size_t &INDEX) const;

      //! @brief gets which reactive atoms are present in a given reactant
      //! @param INDEX the product index of interest
      //! @return a map of reactive atoms (keys) to reactant atom indices (values)
      storage::Map< size_t, size_t> GetReactiveAtomsInReactant( const size_t &INDEX) const;

      //! @brief gets which reactive atoms are present in all reactants
      //! @return a map of reactive atoms (keys) to pairs of (reactant index,reactant atom) (values)
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &GetReactiveAtomsAllReactants() const
      {
        return m_ReactiveAtomsReactants;
      }

      //! @brief return the number of atoms
      //! @return the number of atoms
      const size_t &GetNumberProducts() const;

      //! @brief gets all products
      //! @return a vector of products
      //const storage::Vector< FragmentComplete> &GetProducts() const;
      const storage::Vector< ReactionStructure> &GetProductStructures() const
      {
        return m_ProductStructures;
      }

      storage::Vector< FragmentComplete> GetProducts() const
      {
        storage::Vector< FragmentComplete> mols;
        mols.AllocateMemory( m_ProductStructures.GetSize());
        
        for( size_t s( 0), end_s( m_ProductStructures.GetSize()); s < end_s; ++s)
        {
          mols.PushBack( m_ProductStructures( s).GetFragment());
        }
        return mols;
      }

      //! @brief gets a product with an index
      //! @param INDEX the product to fetch
      const FragmentComplete &GetProduct( const size_t &INDEX) const;

      //! @brief gets which reactive atoms are present in a given product
      //! @param INDEX the product index of interest
      //! @return a map of reactive atoms (keys) to product atom indices (values)
      storage::Map< size_t, size_t> GetReactiveAtomsInProduct( const size_t &INDEX) const;

      //! @brief gets which reactive atoms are present in all product
      //! @return a map of reactive atoms (keys) to pairs of (product index,product atom) (values)
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &GetReactiveAtomsAllProducts() const
      {
        return m_ReactiveAtomsProducts;
      }

      //! @brief reverses the reaction, i.e. makes the products the reactants and vice versa
      void Reverse();

    ///////////////
    // operators //
    ///////////////

      //! @brief equality operator
      //! @param OTHER another ReactionComplete to compare
      //! @return true if the reactions are equivalent
      bool operator ==( const ReactionComplete &OTHER) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write the reaction in RXN format to an output stream
      //! @param OSTREAM the output stream to write to
      //! @return the same output stream that was written to 
      std::ostream &WriteRXN( std::ostream &OSTREAM) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ReactionComplete

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_REACTION_COMPLETE_H_
