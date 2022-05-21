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

#ifndef BCL_CHEMISTRY_ROTAMER_LIBRARY_INTERFACE_H_
#define BCL_CHEMISTRY_ROTAMER_LIBRARY_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_const_graph.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RotamerLibraryInterface
    //! @brief is a handler for storing and retrieving rotamer libraries from a file system.
    //!
    //! @details
    //!
    //! @see @link example_chemistry_rotamer_library_interface.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Jul 01, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RotamerLibraryInterface :
      public util::SerializableInterface
    {

    public:

      //! typename for counts at a given set of bond angles
      typedef storage::List< storage::Pair< linal::Matrix< double>, math::RunningAverage< linal::Vector< double> > > >
        t_BondAnglesWithCounts;

      typedef storage::Triplet
      <
        ConformationGraphConverter::AtomComparisonTypeEnum,
        size_t,
        storage::Vector< storage::Pair< size_t, size_t> >
      > t_BondAngleMapKey;

      // Key <- atom-type <> Vector< bond type, atom type>
      // Value <- Matrix with unit-vector coordinates of all atoms after the first.
      // The first atom is always moved to 1 0 0, second atom is moved such that it is 0 in
      typedef storage::Map< t_BondAngleMapKey, t_BondAnglesWithCounts> t_BondAngleMap;

      //! @return pointer to new RotamerLibraryInterface
      virtual RotamerLibraryInterface *Clone() const = 0;

      //! @return returns default rotamer library to be used
      static const std::string &GetDefault();

    /////////////////
    // data access //
    /////////////////

      //! @brief detect whether the rotamer library is defined
      virtual bool IsDefined() const = 0;

      //! @brief get the directed-graph whose vertices are constitution ids and child nodes are substructures of parent node
      virtual graph::ConstGraph< size_t, size_t> RetrieveRotamerSubstructureTree() const = 0;

      //! @brief get the directed-graph whose vertices are constitution ids and child nodes are substructures of parent node
      virtual graph::ConstGraph< size_t, size_t> RetrieveRotamerRequirementsGraph() const = 0;

      //! @brief retrive constitution associated with id.
      //! @param IDS
      virtual storage::Vector< graph::ConstGraph< size_t, size_t> > RetrieveConstitution( const storage::Vector< size_t> &IDS) const = 0;

      //! @brief retrieve configurations associated with the given constitution id.
      //! @param CONSTITUTION_ID
      virtual storage::Vector< storage::Set< size_t> > RetrieveConfigurationMapping() const = 0;

      //! @brief get constitutions that are root in the substructure tree i.e. they are not contained in any other fragments
      virtual storage::Vector< graph::ConstGraph< size_t, size_t> > GetRootConstitutions() const = 0;

      //! @brief get the bond angle map
      virtual t_BondAngleMap GetBondAngleMap() const = 0;

      virtual FragmentEnsemble RetrieveAssociatedConfigurations( const storage::Set< size_t> &IDS) const = 0;

      //! @brief create the substructure tree of unique constitutions for the  given rotamer library
      //! @param CONFORMATIONS
      virtual void Create( FragmentEnsemble &CONFORMATIONS) const;

      //! @brief get constitutions that are root in the substructure tree i.e. they are not contained in any other fragments
      virtual graph::ConstGraph< size_t, size_t> CreateConstitutionGraphFromString
      (
        const storage::Vector< std::string> &GRAPH_VERTEX_EDGE_DATA,
        const storage::Vector< std::string> &MOLECULE_INFO
      ) const;

    private:

      //! @brief helper function that creates the actual substructure tree
      //! @param ROTAMER_SUBSTRUCTURE_TREE
      //! @param CONSTITUTIONS
      //! @param CONFORMATION_TO_CONSTITUTION_MAPPING
      //! @param CONFORMATIONS
      virtual void CreateImpl
      (
        const size_t &NUMBER_OF_ROOT_NODES,
        const storage::Vector< storage::Set< size_t> > &ROTAMER_SUBSTRUCTURE_TREE,
        const util::ShPtrVector< FragmentConstitutionShared> &CONSTITUTIONS,
        const storage::Vector< storage::Set< size_t> > &CONSTITUTION_TO_CONFORMATION_MAPPING,
        const FragmentEnsemble &CONFORMATIONS,
        const storage::Vector< size_t> &PARENTS = storage::Vector< size_t>()
      ) const = 0;

    }; // class RotamerLibraryInterface

  } // namespace chemistry

} // namespace bcl

#endif // BCL_CHEMISTRY_ROTAMER_LIBRARY_INTERFACE_H_
