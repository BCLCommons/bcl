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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_FRAGMENTS_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_FRAGMENTS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_split_interface.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonByFragments
    //! @brief This class compares molecules according to the number of matching fragments between molecules
    //! @details Comparisons are done using molecular graphs of fragments 
    //!
    //! @see @link example_chemistry_conformation_comparison_by_fragments.cpp @endlink
    //! @author aguilaji, geanesar
    //! @date Dec 19, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonByFragments :
      public ConformationComparisonInterface
    {

    //////////
    // data //
    //////////

      //! single instances of this class
      static const util::SiPtr< const util::ObjectInterface> s_RawInstance;
      static const util::SiPtr< const util::ObjectInterface> s_TanimotoInstance;

      //! graphs of fragments for each molecule
      mutable storage::Map
      < 
        util::SiPtr< const ConformationInterface>, 
        storage::Vector< graph::ConstGraph< size_t, size_t> >
      > m_FragmentGraphs;
      
      //! how to compare atoms of graphs
      ConformationGraphConverter::AtomComparisonTypeEnum m_AtomComparison;
      
      //! how to compare bonds of graphs
      ConfigurationalBondTypeData::DataEnum m_BondComparison;

      //! the graph maker
      ConformationGraphConverter m_GraphConverter;

      //! Flag to compute tanimoto coefficient (vs. raw # of fragments in common)
      bool m_ComputeTanimoto;

      //! A vector of split interfaces to break a molecule into fragments
      storage::Vector< util::Implementation< FragmentSplitInterface> > m_Splits;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor, takes whether to compute tanimoto vs raw values, and whether to enforce a connected solution
      //! @param TANIMOTO whether to compute tanimoto values (vs. raw size of substructure)
      ConformationComparisonByFragments( const bool &TANIMOTO);

      //! @brief copy constructor
      //! @return a pointer to a copy of this class
      ConformationComparisonByFragments *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief gets the splits that this class will use
      //! @return a vector of split implementations
      const storage::Vector< util::Implementation< FragmentSplitInterface> > &GetSplits() const;

      //! @brief adds a split to the class
      //! @param SPLIT the spliting method to add
      void AddSplit( const FragmentSplitInterface &SPLIT);

    ////////////////
    // operations //
    ////////////////

      //! @brief find the RMSD between two conformations
      //! @param MOLECULE_B - the fragment to align against
      //! @param MOLECULE_B - the fragment being aligned
      //! @return the RMSD between MOLECULE_A and MOLECULE_B
      double operator()
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief prepare the class for comparing a conformation
      //! @param MOLECULE the molecule to prepare to compare
      void Prepare( const ConformationInterface &MOLECULE) const;

      //! @brief makes graphs of fragments of a molecule
      //! @param MOLECULE the molecule to fragment
      //! @return a vector of unique graphs that represent fragments
      storage::Vector< graph::ConstGraph< size_t, size_t> > MakeFragmentGraphs
      (
        const ConformationInterface &MOLECULE
      ) const;

      //! @brief determine if a graph is contained in a vector of graphs, i.e. same size and an isomorphism exists
      //! @param GRAPH the graph to check for
      //! @param GRAPH_VECTOR the vector of graphs to check
      //! @return the index of the first graph that matches, or GRAPH_VECTOR.GetSize() if not found
      size_t GraphVectorFind
      ( 
        const graph::ConstGraph< size_t, size_t> &GRAPH, 
        const storage::Vector< graph::ConstGraph< size_t, size_t> > &GRAPH_VECTOR
      ) const;

    protected:
      
      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_FRAGMENTS_H_
