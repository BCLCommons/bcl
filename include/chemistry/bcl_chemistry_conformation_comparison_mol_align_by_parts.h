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

/*
 *bcl_chemistry_conformation_comparison_mol_align_by_parts.h
 *
 *  Created on: Jun 26, 2015
 *      Author: brownbp1
 */

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_MOL_ALIGN_BY_PARTS_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_MOL_ALIGN_BY_PARTS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "bcl_chemistry_conformation_comparison_psi_field.h"
#include "bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_sample_conformations.h"
#include "descriptor/bcl_descriptor_dataset_builder.h"
#include "storage/bcl_storage_triplet.h"
namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonMolAlignByParts
    //! @brief This class performs iterative alignment of molecules with restraints placed on conformer generation and
    //! atom positions; these restraints are iteratively updated throughout the alignment procedure
    //!
    //! @see @link example_chemistry_conformation_comparison_mol_align_by_parts.cpp @endlink
    //! @author brownbp1
    //! @date Aug 15, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonMolAlignByParts :
      public ConformationComparisonPsiField
    {

    protected:

    //////////
    // data //
    //////////

      //! mobile atom indices
      std::string m_MobileAtoms;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      ConformationComparisonMolAlignByParts();

      //! virtual copy constructor
      ConformationComparisonMolAlignByParts *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief returns the atom indices that can be sampled during alignment
      //! @return the mobile atom indices
      const std::string &GetMobileAtoms() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief make conformers of input molecules, randomly pair conformers
      //! @brief align via multi-tiered MCM and find best property RMSD between input molecules
      //! @param MOLECULE_A - first molecule from which conformational library will be generated
      //! @param MOLECULE_B - second molecule from which conformational library will be generated
      //! @return best conformer pairings and property rmsd value
      double operator()
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B
      ) const;

      //! @brief align two small molecule objects and find the property RMSD
      //! @param CONFORMERS_MOL_A - conformers generated from input molecule A
      //! @param CONFORMERS_MOL_B - conformers generated from input molecule B
      storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > MolAlignByParts
      (
        const FragmentEnsemble &CONFORMERS_MOL_A,
        const FragmentEnsemble &CONFORMERS_MOL_B
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set the number of conformer pairs
      void SetSampleByParts( const std::string &MOBILE_ATOM_INDICES);

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_MOL_ALIGN_BY_PARTS_H_
