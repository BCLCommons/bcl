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

//
#ifndef BCL_CHEMISTRY_REACTION_STRUCTURE_H_
#define BCL_CHEMISTRY_REACTION_STRUCTURE_H_

// include the namespace header
#include "bcl_chemistry.h"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_complete.h"
#include "graph/bcl_graph_const_graph.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ReactionStructure
    //! @brief Class that contains molecular configuration data
    //! @details Models stereochemistry, isomeric fragments, chemical adjacency, and aromatic and ring structures
    //!
    //! @see @link example_chemistry_reaction_structure.cpp @endlink
    //! @author geanesar
    //! @date Oct 15, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ReactionStructure :
      public util::ObjectInterface
    {

    public:

      // aromaticity, should be a bitfield
      enum Aromaticity
      {
        e_Aliphatic = 1,
        e_Aromatic = 2,
        e_AliphaticOrAromatic = 3,
        s_NumberAromaticities = 4
      };

    private:

    //////////
    // data //
    //////////

      //! molecular structure
      FragmentComplete m_Fragment;

      //! graph of m_Fragment colored by element type
      graph::ConstGraph< size_t, size_t> m_ElementGraph; //!< graph of the structure

      //! vector specifying allowed aromaticities for each atom
      std::vector< Aromaticity> m_AllowedAromaticity; //!< for queries, allowed aromaticity

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! graph converter used to create graphs for this class (initialized with correct parameters)
      static const ConformationGraphConverter s_Converter;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      ReactionStructure();

      //! @brief constructor with molecular structure and aromaticity specification
      //! @param FRAGMENT to use
      //! @param AROMATICITY aromaticity specification of all atoms
      ReactionStructure
      (
        const FragmentComplete &FRAGMENT,
        const std::vector< Aromaticity> &AROMATICITY = std::vector< Aromaticity>()
      );

      //! @brief Clone function
      //! @return pointer to new ReactionStructure
      ReactionStructure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the internally stored ElementType/BondOrderOrAromatic graph
      const graph::ConstGraph< size_t, size_t> &GetGraph() const;

      //! @brief get the number of atoms in the structure
      const size_t &GetSize() const;

      //! @brief get internal molecular structure data
      const FragmentComplete &GetFragment() const;

      //! @brief sets the aromaticity of a single atom
      //! @param ATOM_NO atom number to change aromaticity of
      //! @param AROMATICITY the aromaticity specification
      void SetAllowedAromaticity( const size_t &ATOM_NO, const Aromaticity &AROMATICITY);

    ///////////////
    // operators //
    ///////////////

      //! @brief get a list of isomorphisms that match this structure
      //! @param CONFORMATION the conformation to search
      //! @return a list of atom-to-atom isomorphisms that match this structure
      storage::Vector< storage::Vector< size_t> > GetMatchingSubstructures( const ConformationInterface &CONFORMATION) const;

      //! @brief determine if this structure contains another
      //! @param OTHER the reaction structure to check for
      //! @return true if this structure contains OTHER (or is equal to)
      bool Contains( const ReactionStructure &OTHER) const;

      //! @brief determine if this structure matches parts of a molecule
      //! @param CONFORMATION the conformation to check
      //! @return true if this structure matches any parts of CONFORMATION
      bool ContainedIn( const ConformationInterface &CONFORMATION) const;

      //! @brief equality operator; graphs are equal sized, isomorphic, and all properties (e.g. aromaticity) are equivalent
      bool operator ==( const ReactionStructure &OTHER) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief helper function to determine if a query aromaticity is allowed according to another aromaticity
      //! @param QUERY the aromaticity to check, will likely be one of e_Aromatic or e_Aliphatic
      //! @param TEST the aromaticity to check against, may be any Aromaticity specification
      //! @return true if QUERY is an allowed aromaticity value accordin to TEST
      //! @details since Aromaticity is implemented as a bitfield, this means that all bits set in QUERY are also set in TEST
      static bool AromaticityMatches( const Aromaticity &QUERY, const Aromaticity &TEST);

    protected:

      //! @brief get subgraph isomorphisms of one graph with graph from this object, including matching aromaticity
      //! @param OTHER_GRAPH the other graph to search
      //! @param OTHER_AROMATICITY aromaticity specification of atoms in OTHER_GRAPH
      //! @return a list of isomorphisms for substructures of OTHER_GRAPH that are equivalent to this one
      storage::Vector< storage::Vector< size_t> > GetMatchingSubstructures
      (
        const graph::ConstGraph< size_t, size_t> &OTHER_GRAPH,
        const std::vector< Aromaticity> &OTHER_AROMATICITY
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ReactionStructure

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_REACTION_STRUCTURE_H_
