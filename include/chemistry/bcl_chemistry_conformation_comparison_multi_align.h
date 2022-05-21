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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_MULTI_ALIGN_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_MULTI_ALIGN_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_conformation_comparison_multi_align.h"
#include "bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "bcl_chemistry_conformation_comparison_psi_field.h"
#include "bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_sample_conformations.h"
#include "descriptor/bcl_descriptor_prediction.h"
#include "storage/bcl_storage_triplet.h"
namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonMultiAlign
    //! @brief This class is designed to be an expanded version of ConformationComparisonPsiField to
    //!        include pharmacophore mapping
    //!
    //! @see @link example_chemistry_conformation_comparison_multi_align.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date December 31, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonMultiAlign :
      public ConformationComparisonPsiFlexField
      {
    private:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! name of output file
      std::string m_FileTagMulti;

      //!
      std::string m_MovieLigFitFilename;

      //! additional instance of alignment class with re-calculated hypothesis molecule properties
      ConformationComparisonPsiFlexField m_InstanceTwo;

      //! QSAR model with which to select first alignment conformers
      descriptor::Prediction< AtomConformationalInterface> m_Predictor;

      //! binding pocket for the target ligands
      util::SiPtr< FragmentComplete> m_BindingPocket;

      //!
      linal::Vector3D m_StartPos;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConformationComparisonMultiAlign();

      //! virtual copy constrctor
      ConformationComparisonMultiAlign *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      double operator()
      (
        const FragmentEnsemble &MOLECULE_LIST,
        const FragmentEnsemble &POCKET_LIST
      ) const;

    //////////////////////
    // helper functions //
    /////////////////////

      linal::Vector< double> GetAxisAlignedBoundingBox( FragmentComplete &TARGET) const;

      storage::Triplet< FragmentComplete, FragmentComplete, double> AlignFullyFlex( const FragmentEnsemble &ENS_A, const FragmentEnsemble &ENS_B) const;

      FragmentComplete BondAlignInfEns
      (
        FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const storage::Vector< storage::Triplet< size_t, size_t, float> > &INDICES
          ) const;

      //Used to rapidly align molecules in conformational co-space
      storage::Pair< storage::Vector< storage::Triplet< size_t, size_t, float> >, storage::Pair< FragmentComplete, FragmentComplete> >
      AlignMolecularEnsemble
      (
        const FragmentEnsemble &ENSEMBLE_A,
        const FragmentEnsemble &ENSEMBLE_B
      ) const;

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

      }; // ConformationComparisonMultiAlign
  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_MULTI_ALIGN_H_
