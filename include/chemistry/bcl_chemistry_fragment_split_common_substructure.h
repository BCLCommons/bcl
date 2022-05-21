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

#ifndef BCL_CHEMISTRY_FRAGMENT_SPLIT_COMMON_SUBSTRUCTURE_H_
#define BCL_CHEMISTRY_FRAGMENT_SPLIT_COMMON_SUBSTRUCTURE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_split_interface.h"
#include "util/bcl_util_enumerated.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentSplitCommonSubstructure
    //! @brief Splits molecules into the largest common substructure they possess relative to molecules of an input file
    //!
    //! @see @link example_chemistry_fragment_split_largest_common_substructure.cpp @endlink
    //! @author geanesar, brownbp1
    //! @date Mar 14, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentSplitCommonSubstructure :
      public FragmentSplitInterface
    {

    private:

    //////////
    // data //
    //////////

        //! the conformation graph converter to use
        ConformationGraphConverter m_Converter;

        //! the filename to read molecules from
        std::string m_File;

        //! the type of atom comparison to perform
        ConformationGraphConverter::AtomComparisonTypeEnum m_AtomComparison;

        //! the type of bond comparison to perform
        ConfigurationalBondTypeData::DataEnum m_BondComparison;

        //! obtain the complement to the largest common substructure instead
        bool m_Complement;

        //! whether molecules have been read or not
        mutable bool m_FileWasRead;

        //! the molecules that are read from the file, used to compare incoming molecules
        mutable FragmentEnsemble m_Molecules;

        //! graphs of the molecules that are read from file
        mutable storage::Vector< graph::ConstGraph< size_t, size_t> > m_MoleculeGraphs;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      FragmentSplitCommonSubstructure *Clone() const;

      //! @brief constructor
      //! @param GET_RINGS true if rings are desired, false if chains are desired
      //! @param MIN_SIZE get the minimum size of ring or chain that is desired
      FragmentSplitCommonSubstructure
      (
        const std::string &FILENAME = "",
        const ConformationGraphConverter::AtomComparisonType ATOM_COMPARISON = ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data BOND_COMPARISON = ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness,
        const bool &COMPLEMENT = false
      );

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief Get a description for what this class does (used when writing help)
      //! @return a description for what this class does (used when writing help)
      const std::string &GetClassDescription() const;

      //! @brief gets the minimum size of fragments
      //! @return the minimum size of fragments
      const size_t GetMinSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns an ensemble of fragments of a molecule
      //! @param CONFORMATION molecule of interest
      //! @return an ensemble of common substructures relative to those in a file
      FragmentEnsemble operator()( const ConformationInterface &CONFORMATION) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief reads in molecules from a given file if it is necessary
      void ReadFile() const;

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

#endif // BCL_CHEMISTRY_FRAGMENT_SPLIT_COMMON_SUBSTRUCTURE_H_
