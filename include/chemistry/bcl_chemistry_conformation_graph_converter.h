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

#ifndef BCL_CHEMISTRY_CONFORMATION_GRAPH_CONVERTER_H_
#define BCL_CHEMISTRY_CONFORMATION_GRAPH_CONVERTER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_configurational_bond_type_data.h"
#include "bcl_chemistry_conformation_interface.h"
#include "graph/bcl_graph_const_graph.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationGraphConverter
    //! @brief Convert between graphs and small molecules
    //! @details Here, we have a graph with numbered vertices and edges. size_t's represent bond type or element type.
    //! From the graph, construct a small molecule object by constructing its connectivity matrix.
    //!
    //! @see @link example_chemistry_conformation_graph_converter.cpp @endlink
    //! @author mendenjl
    //! @date Jan 24, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationGraphConverter :
      public util::ObjectInterface
    {
    public:

    //////////
    // Enum //
    //////////

      enum AtomComparisonType
      {
        e_Identity,                        //!< 1 for all atoms
        e_ElementType,                     //!< Element type index
        e_AtomType,                        //!< Atom type index
        e_AtomTypeAndChirality,            //!< Chirality index and atom type index packed together
        e_AtomTypeAndComplexRingChirality, //!< Like e_AtomTypeAndChirality, but chirality only packed together for atoms
                                           //!< with at least three ring bonds
        e_AtomTypeAndSymmetry,             //!< Like e_AtomType, but with extra bits to indicate how many other atoms this
                                           //!< atom is symmetric with
        e_AtomTypeAndHasSymmetry,          //!< Like e_AtomTypeAndSymmetry, but only indicates whether it is symmetric with
        e_AtomTypeAndNumberHydrogens,      //!< Atom type and number of hydrogens
        e_AtomTypeAndNumberHydrogensOnRings, //!< Atom type and number of hydrogens on ring atoms
        //! Atom type and number of hydrogens on ring atoms and distinguish multiple hydrogen substituents
        //! This is used to prevent combinatorial explosion when doing graph isomorphism searches with all hydrogens
        e_AtomTypeAndNumberHydrogensOnRingsAndDistinguishHydrogens,
        e_AtomTypeAndDistinguishHydrogens,
        e_CIPPriorityHighToLow,            //!< Priority of the atoms
        e_CouldHaveSubstituents,           //!< Atoms that can be substituted
        s_NumberAtomComparisonTypes
      };

      //! @brief Data as string
      //! @param DATA the data whose name is desired
      //! @return the name as string
      static const std::string &GetAtomComparisonType( const AtomComparisonType &DATA);

      //! DataEnum simplifies the usage of the Data enum of this class
      typedef util::WrapperEnum< AtomComparisonType, &GetAtomComparisonType, s_NumberAtomComparisonTypes>
        AtomComparisonTypeEnum;

      //! Type of atom graph
      typedef graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> t_AtomGraph;

    private:

    //////////
    // data //
    //////////

      AtomComparisonTypeEnum                m_AtomRepresentation; //!< method of representing atoms as size_ts for the graph
      ConfigurationalBondTypeData::DataEnum m_BondRepresentation; //!< method of representing bonds as size_ts in the graph
      bool                                  m_RemoveH;            //!< True if H should be removed from any graphs that are created

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, compares by atom type, chirality, and e_BondOrderAmideWithIsometryOrAromaticWithRingness
      ConformationGraphConverter();

      //! @brief constructor from coloring schemes
      //! @param ATOM_COMPARISON_TYPE means by which to color vertices of a molecule
      //! @param BOND_TYPE_INFO the desired information use for the bonds
      //! @param REMOVE_H true to remove H from any graphs that are created
      ConformationGraphConverter
      (
        const AtomComparisonType &ATOM_COMPARISON_TYPE,
        const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO,
        const bool &REMOVE_H = false
      );

      //! @brief virtual copy constructor
      ConformationGraphConverter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Create a graph of a conformation
      //! @param CONFORMATION a conformation
      //! @return The conformation converted into a graph with the given atom/bond representations
      graph::ConstGraph< size_t, size_t> operator()( const ConformationInterface &CONFORMATION) const;

      //! @brief Create a tree of a conformation, vertex data indicates the number of atoms represented by the point
      //!        1 for atoms in a chain, > 2 for atoms in rings
      //! @param CONFORMATION a conformation
      //! @return The conformation converted into a graph with the given atom/bond representations,
      //!         along with a vector, each index corresponds to atom in the molecule, value in vector is the index in the
      //!         graph corresponding to that point
      storage::Pair< graph::ConstGraph< size_t, size_t>, storage::Vector< size_t> >
        ToTree( const ConformationInterface &CONFORMATION) const;

      //! @brief Create a tree of a conformation, vertex data indicates the number of atoms represented by the point
      //!        1 for atoms in a chain, > 2 for atoms in rings
      //! @param CONFORMATION a conformation
      //! @return The conformation converted into a graph with the given atom/bond representations,
      //!         along with a vector, each index corresponds to atom in the molecule, value in vector is the index in the
      //!         graph corresponding to that point
      size_t CountNonRingVariantIsomorphisms( const ConformationInterface &CONFORMATION) const;

      //! @brief Create an atom vector from a graph with atoms
      //! @param GRAPH a graph created by this converter with CreateGraphWithAtoms
      //! @param RECALCULATE_CONFIGURATION whether to recalculate configuration
      //! @return the generated atom vector
      static AtomVector< AtomComplete> CreateAtomsFromGraph
      (
        const t_AtomGraph &GRAPH,
        const bool &RECALCULATE_CONFIGURATION = true
      );

      //! @brief Given a small molecule, instantiate its graphical representation
      //! @param SMALL_MOL a conformation
      //! @return The resulting graph after the conversion.
      static t_AtomGraph CreateGraphWithAtoms( const ConformationInterface &CONFORMATION);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief get the atom's representation as a size_t using m_AtomRepresentation
      //! @param ATOM atom to retrieve data from
      //! @param TYPE the type of atom comparison that will be used
      //! @return the atom's representation as a size_t using m_AtomRepresentation
      static size_t ConvertAtomData( const AtomConformationalInterface &ATOM, const AtomComparisonType &TYPE);

      //! @brief get the atom's representation as a size_t using m_AtomRepresentation
      //! @param ATOM atom type to retrieve data from
      //! @param TYPE the type of atom comparison that will be used
      //! @return the atom's representation as a size_t using m_AtomRepresentation
      static size_t ConvertAtomTypeData( const AtomType &ATOM, const AtomComparisonType &TYPE);

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

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class ConformationGraphConverter

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_GRAPH_CONVERTER_H_

