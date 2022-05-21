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

#ifndef BCL_CHEMISTRY_FRAGMENT_GRAPH_MARKER_H_
#define BCL_CHEMISTRY_FRAGMENT_GRAPH_MARKER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configuration_set.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_split_interface.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentGraphMarker
    //! @brief Class that finds conformer cluster centers for all unique conformations of a given configuration
    //!
    //! @see @link example_chemistry_fragment_graph_marker.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Oct 23, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentGraphMarker :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

      //! Set of fragments to color differently
      ConstitutionSet m_SpecialFragments;

      //! Conformation graph converter object
      ConformationGraphConverter m_Converter;

      //! Object that will split a given fragment into different components, each of which should be colored
      //! differently
      util::Implementation< FragmentSplitInterface> m_Splitter;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! @param BOND_TYPE the bond type to be used for graph isomorphism
      //! @param SPLIT_PARAMETER the type of splitting that needs to be done like getting only rings, chain, scaffolds etc
      FragmentGraphMarker
      (
        const ConformationGraphConverter &BOND_TYPE,
        const util::Implementation< FragmentSplitInterface> &SPLIT_PARAMETER
      );

      //! @brief Clone function
      //! @return pointer to new FragmentGraphMarker
      FragmentGraphMarker *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the unique configuration passed through this class
      //! @return unique configuration passed through this class
      const util::ShPtrList< FragmentConstitutionShared> &GetConstitution() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief Create a graph of a conformation
      //! @param CONFORMATION a conformation
      //! @return The conformation converted into a graph with the given atom/bond representations
      graph::ConstGraph< size_t, size_t> operator()( const ConformationInterface &CONFORMATION) const;

      //! @brief create a vector of graphs of conformations that are passed as an ensemble
      //! @param ENSEMBLE ensemble of conformations whose graphs are desired
      //! @return a vector of graphs of a conformations that are passed as an ensemble
      storage::Vector< graph::ConstGraph< size_t, size_t> > operator()( const FragmentEnsemble &ENSEMBLE) const;

      //! @brief Create a graph of a conformation and store configuration in a configurationset
      //! @param CONFORMATION a conformation
      //! @return The conformation converted into a graph with the given atom/bond representations
      graph::ConstGraph< size_t, size_t> Insert( const ConformationInterface &ENSEMBLE);

      //! @brief create a vector of graphs of conformations that are passed as an ensemble and storing the ensemble in a configuratin set
      //! @param ENSEMBLE ensemble of conformations whose graphs are desired
      //! @return create a vector of graphs of conformations that are passed as an ensemble
      storage::Vector< graph::ConstGraph< size_t, size_t> > Insert( const FragmentEnsemble &ENSEMBLE);

    /////////////////////
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

    }; // class FragmentGraphMarker

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_GRAPH_MARKER_H_
