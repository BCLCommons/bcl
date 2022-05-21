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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_SUBSTRUCTURE_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_SUBSTRUCTURE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configurational_bond_type_data.h"
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonBySubstructure
    //! @brief This class compares the largest common substructure of molecules
    //!
    //! @see @link example_chemistry_conformation_comparison_by_substructure.cpp @endlink
    //! @author mendenjl
    //! @date Jun 22, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonBySubstructure :
      public ConformationComparisonInterface
    {

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! Graphs of all molecules
      mutable storage::Map< util::SiPtr< const ConformationInterface>, graph::ConstGraph< size_t, size_t> > m_Graphs;

      //! Lower bound (for tanimoto coefficient or # atoms for raw comparison
      double m_LowerBound;

      //! Scheme used for comparing whether two bonds are equivalent
      ConfigurationalBondTypeData::DataEnum m_BondComparisonType;

      //! Scheme used for comparing whether two atoms are equivalent
      ConformationGraphConverter::AtomComparisonTypeEnum m_AtomComparisonType;

      //! Flag to require a connected substructure
      graph::CommonSubgraphIsomorphismBase::SolutionType m_SolutionType;

      //! Flag to compute tanimoto coefficient (vs. raw # of atoms in common)
      bool m_ComputeTanimoto;

      //! whether to consider hydrogens.  This generally just slows things down, but may be important for some applications
      bool m_CompareH;

      //! Metric for comparison; if left undefined, computes tanimoto or distance
      util::Implementation< ConformationComparisonInterface> m_Metric;

      //! Whether to compare by number of interior vertices vs. number of vertices.
      //! The difference is small for connected subgraphs, but significant for unconnected subgraphs
      bool m_InteriorWeighted;

      //! Whether to use the exhaustive unconnected vertex algorithm.  This is roughly a multiple of O(N^2) slower on
      //! chemical graphs, which makes it unacceptably slow for most uses...but maybe someone wants it.
      bool m_Exhaustive;

      //! Mutex for accessing m_Allocated and m_Available
      mutable sched::Mutex m_Mutex;

      //! Pool of isomorphism-computing objects, stored because they are very large and excessive memory reallocations
      //! can be avoided by maintaining such a pool, without making it a data member directly, which would make this
      //! class's operator() not thread safe
      mutable storage::List< graph::CommonSubgraphIsomorphism< size_t, size_t> > m_Allocated;
      mutable storage::List< graph::CommonSubgraphIsomorphism< size_t, size_t> > m_Available;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor, takes whether to compute tanimoto vs raw values, and whether to enforce a connected solution
      //! @param TANIMOTO whether to compute tanimoto values (vs. raw size of substructure)
      //! @param TYPE whether to enforce a connected substructure solution
      ConformationComparisonBySubstructure
      (
        const bool &TANIMOTO,
        const graph::CommonSubgraphIsomorphismBase::SolutionType &TYPE
      );

      //! @param COMPARISON This metric will be used to compare the resulting substructures rather
      //!        than tanimoto or raw distance
      //! @param TYPE whether to enforce a connected substructure solution
      ConformationComparisonBySubstructure
      (
        const ConformationComparisonInterface &COMPARISON,
        const graph::CommonSubgraphIsomorphismBase::SolutionType &TYPE
      );

      //! virtual copy constructor
      ConformationComparisonBySubstructure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

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

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief acquire an isomorphism object
      //! @return iterator to the isomorphism object
      storage::List< graph::CommonSubgraphIsomorphism< size_t, size_t> >::iterator AcquireIsomorphism() const;

      //! @brief release a given isomorphism object
      //! @param ITR iterator to the isomorphism object
      void ReleaseIsomorphism
      (
        const storage::List< graph::CommonSubgraphIsomorphism< size_t, size_t> >::iterator &ITR
      ) const;

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_SUBSTRUCTURE_H_
