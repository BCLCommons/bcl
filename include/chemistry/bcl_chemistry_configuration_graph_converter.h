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

#ifndef BCL_CHEMISTRY_CONFIGURATION_GRAPH_CONVERTER_H_
#define BCL_CHEMISTRY_CONFIGURATION_GRAPH_CONVERTER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configuration_interface.h"
#include "bcl_chemistry_configurational_bond_type_data.h"
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
    //! @class ConfigurationGraphConverter
    //! @brief Convert between graphs and small molecules
    //! @details Here, we have a graph with numbered vertices and edges. size_t's represent bond type or element type.
    //! From the graph, construct a small molecule object by constructing its connectivity matrix.
    //!
    //! @see @link example_chemistry_configuration_graph_converter.cpp @endlink
    //! @author mendenjl
    //! @date Jan 24, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConfigurationGraphConverter :
      public util::ObjectInterface
    {
    public:

    //////////
    // Enum //
    //////////

      enum AtomComparisonType
      {
        e_Identity,                //!< 1 for all atoms
        e_ElementType,             //!< Element type index
        e_AtomType,                //!< Atom type index
        e_AtomTypeAndChirality,    //!< Chirality index and atom type index packed together
        s_NumberAtomComparisonTypes
      };

      //! @brief Data as string
      //! @param DATA the data whose name is desired
      //! @return the name as string
      static const std::string &GetAtomComparisonType( const AtomComparisonType &DATA);

      //! DataEnum simplifies the usage of the Data enum of this class
      typedef util::WrapperEnum< AtomComparisonType, &GetAtomComparisonType, s_NumberAtomComparisonTypes>
        AtomComparisonTypeEnum;

    //////////
    // data //
    //////////

      AtomComparisonTypeEnum                m_AtomRepresentation; //!< method of representing atoms as size_ts for the graph
      ConfigurationalBondTypeData::DataEnum m_BondRepresentation; //!< method of representing bonds as size_ts in the graph

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, compares by atom type, chirality, and e_BondOrderAmideWithIsometryOrAromaticWithRingness
      ConfigurationGraphConverter();

      //! @brief constructor from coloring schemes
      //! @param ATOM_COMPARISON_TYPE means by which to color vertices of a molecule
      //! @param BOND_TYPE_INFO the desired information use for the bonds
      ConfigurationGraphConverter
      (
        const AtomComparisonType &ATOM_COMPARISON_TYPE,
        const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO
      );

      //! @brief virtual copy constructor
      ConfigurationGraphConverter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Create a graph of a configuration
      //! @param CONFIGURATION a configuration
      //! @return The configuration converted into a graph with the given coloring schemes
      graph::ConstGraph< size_t, size_t> operator()( const ConfigurationInterface &CONFIGURATION) const;

    //////////////////////
    // input and output //
    //////////////////////

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

      //! @brief get the atom's representation as a size_t using m_AtomRepresentation
      //! @param ATOM atom to retrieve data from
      //! @return the atom's representation as a size_t using m_AtomRepresentation
      size_t GetAtomData( const AtomConfigurationalInterface &ATOM) const;

    }; // class ConfigurationGraphConverter

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFIGURATION_GRAPH_CONVERTER_H_

