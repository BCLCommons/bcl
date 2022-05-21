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

#ifndef BCL_DESCRIPTOR_MOL_ALIGN_PHARM_SCORE_H_
#define BCL_DESCRIPTOR_MOL_ALIGN_PHARM_SCORE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "bcl_descriptor_molecule_3da_smooth_sign_code.h"
#include "bcl_descriptor_molecule_3da_smooth_sign_occlusion_code.h"
#include "bcl_descriptor_window.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "util/bcl_util_implementation.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MolAlignPharmScore
    //! @brief code object for 3DARealSpace
    //! @details Calculates the 3DASign of the contacting residues between a protein and a ligand
    //!  such that the 3DA distance bins are partitioned across contact distance.
    //!
    //! @see @link example_descriptor_molecule_3DA_pair_convolution.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Sep 20, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MolAlignPharmScore :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      //! Scaffold molecule filename
      std::string m_ScaffoldFilename;

      //! Alignment type (rigid vs. flexible)
      std::string m_AlignmentType;

      //! Scaffold molecule for alignment
      chemistry::FragmentComplete m_Scaffold;

      //! Molecular alignment engine
      chemistry::ConformationComparisonPsiFlexField m_Aligner;

      //! Sample conformations for flexible alignment
      chemistry::SampleConformations m_SampleConformations;

      //! Voxel grid for the scaffold molecule
      chemistry::VoxelGridAtom m_VoxelGridScaffold;

      //! Maximum atom distance for determining mutually matching atoms and neighbors
      float m_MaxAtomDistance;

      //! Properties to compute as part of the alignment
      Combine< chemistry::AtomConformationalInterface, float> m_Properties;

      //! Weights for each atom property in the overall correlation
      linal::Vector< float> m_PropertyWeights;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MolAlignPharmScore();

      //! @brief virtual copy constructor
      MolAlignPharmScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 3;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

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

    }; // class MolAlignPharmScore

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOL_ALIGN_PHARM_SCORE_H_
