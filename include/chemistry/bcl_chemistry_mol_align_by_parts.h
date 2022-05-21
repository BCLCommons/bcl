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

#ifndef BCL_CHEMISTRY_MOL_ALIGN_BY_PARTS_H_
#define BCL_CHEMISTRY_MOL_ALIGN_BY_PARTS_H_

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
    //! @class MolAlignByParts
    //! @brief This class performs iterative alignment of molecules with restraints placed on conformer generation and
    //! atom positions; these restraints are iteratively updated throughout the alignment procedure
    //!
    //! @see @link example_chemistry_mol_align_by_parts.cpp @endlink
    //! @author brownbp1
    //! @date Aug 16, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MolAlignByParts :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

//    private:

      //! mobile atom indices
      storage::Vector< size_t> m_MobileAtoms;

      //! fixed (immobilized) atom indices
      storage::Vector< size_t> m_FixedAtoms;

      //! masked (non-scored) atom indices
      storage::Vector< size_t> m_TargetMaskedAtoms;

      //! visible (scored) atom indices
      storage::Vector< size_t> m_TargetVisibleAtoms;

      //! alignment trials per scaffold
      size_t m_TrialsPerScaffold;

      //! molalign object if alignment is to allow movement of more than internal DOFs
      ConformationComparisonPsiField m_MolAlign;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      MolAlignByParts();

      //! virtual copy constructor
      MolAlignByParts *Clone() const;

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
      const storage::Vector< size_t> &GetMobileAtoms() const;

      //! @brief returns the atom indices that cannot be sampled during alignment
      //! @return the fixed atom indices
      const storage::Vector< size_t> &GetFixedAtoms() const;

      //! @brief returns the atom indices that are hidden from scoring during alignment
      //! @return the masked atom indices
      const storage::Vector< size_t> &GetMaskedAtoms() const;

      //! @brief returns the atom indices that can scored during alignment
      //! @return the visible atom indices
      const storage::Vector< size_t> &GetVisibleAtoms() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief align two small molecule objects and find the property RMSD
      //! @param CONFORMERS_MOL_A - conformers generated from input molecule A
      //! @param CONFORMERS_MOL_B - conformers generated from input molecule B
      storage::Vector< storage::Triplet< FragmentComplete, FragmentEnsemble, double> > Align
      (
        const FragmentComplete &TARGET_MOL,
        const FragmentEnsemble &SCAFFOLD_ENS
      );

      //! @brief initialize the sample by parts atoms from member data on current molecule
      void InitializeSampleByParts( FragmentComplete &MOL);

      //! @brief initialize masking from member data
      //! @param MOL - target molecule
      //! @param MASK_FIXED_ATOMS - only allow scoring on mobile atoms; mask the fixed atoms
      void InitializeAtomMasks
      (
        const FragmentComplete &MOL,
        const bool &MASK_FIXED_ATOMS
      );

      //! @brief set atom properties and remove size zero properties
      //! @param PROPERTIES descriptor label to use
      descriptor::Combine< AtomConformationalInterface, float> PrepareProperties( const util::ObjectDataLabel &PROPERTIES) const;

      //! @brief set the atom indices that can be sampled during alignment
      void SetMobileAtoms( const storage::Vector< size_t> &MOBILE_ATOM_INDICES);

      //! @brief set the atom indices that are fixed during alignment
      void SetFixedAtoms( const storage::Vector< size_t> &FIXED_ATOM_INDICES);

      //! @brief set the atom indices that are hidden from scoring during alignment
      void SetMaskedAtoms( const storage::Vector< size_t> &MASKED_ATOM_INDICES);

      //! @brief set the atom indices that can be scored during alignment
      void SetVisibleAtoms( const storage::Vector< size_t> &VISIBLE_ATOM_INDICES);

      //! @brief combine a bunch of fragments without destroying internal connectivity
      FragmentComplete CreateMutantFragment( const FragmentEnsemble &ENSEMBLE) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief write a BCL vector in a horizontal string format for easier reading
      std::string PrintClean( const storage::Vector< size_t> &INPUT) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MOL_ALIGN_BY_PARTS_H_
