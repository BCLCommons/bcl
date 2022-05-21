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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_SYMMETRY_RMSD_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_SYMMETRY_RMSD_H_

// include the namespace header
#include "bcl_chemistry.h"
#include "signal/bcl_signal_slots.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configurational_bond_type_data.h"
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "quality/bcl_quality_rmsd_preprocessor.h"
#include "storage/bcl_storage_hash_map.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonBySymmetryRmsd
    //! @brief This class is designed to be used for determining and comparing 3D structures for molecules by taking into.
    //! account automorphisms.
    //! @see @link example_chemistry_conformation_comparison_by_symmetry_rmsd.cpp @endlink
    //! @author kothiwsk
    //! @date Mar 11, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonBySymmetryRmsd :
      public ConformationComparisonInterface,
      public signal::Slots
    {
      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_RealSpaceInstance;

      //! bool, whether to perform super-imposition
      bool m_Superimpose;

      //! bool, to consider H
      bool m_ConsiderH;

      //! atom comparison type
      ConformationGraphConverter::AtomComparisonTypeEnum m_AtomComparison;

      //! bond comparison type
      ConfigurationalBondTypeData::DataEnum m_BondComparison;

      mutable storage::Vector< storage::Vector< size_t> > m_Isomorphisms;

      mutable storage::Vector< AtomType> m_PreviousAtomTypes;
      mutable storage::Vector< graph::UndirectedEdge< size_t> > m_PreviousEdgeTypes;
      mutable std::string m_PreviousMolA;

      //! Isomorphism limit. For highly symmetric molecules, limiting the number of isomorphisms will substantially speed up
      //! the calculation, usually with only minimal loss of precision
      size_t m_NumberIsomorphismLimit;

      //! hash of molecules to preprocessed RMSD stuff
      mutable storage::HashMap< size_t, storage::Vector< quality::RMSDPreprocessor> > m_Cache;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, accepts bool that indicates whether to superimpose (default) or not
      ConformationComparisonBySymmetryRmsd( const bool &SUPERIMPOSE = true, const size_t &ISO_LIMIT = 1000);

      //! virtual copy constructor
      ConformationComparisonBySymmetryRmsd *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief returns the number of isomorphisms
      //! @return the number of isomorphisms
      size_t GetNumberIsomorphisms() const;

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

      //! @brief prepare the cache entry for a given molecule
      const storage::Vector< quality::RMSDPreprocessor> &GetCachedInfo( const ConformationInterface &MOLECULE) const;

      //! @brief remove results for this object from the cache
      //! @param ARGUMENT argument that should be removed
      void RemoveResultFromCache( const ConformationInterface &ARGUMENT);

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_BY_SYMMETRY_RMSD_H_
