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

#ifndef BCL_CHEMISTRY_LIGAND_POCKET_FIT_SCORE_H_
#define BCL_CHEMISTRY_LIGAND_POCKET_FIT_SCORE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "bcl_chemistry_conformation_comparison_property_rmsd_x.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "linal/bcl_linal_vector_3d.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_function_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_metropolis.h"
#include "mc/bcl_mc_temperature_default.h"
#include "mc/bcl_mc_temperature_interface.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_function.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_result_threshold.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LigandPocketFitScore
    //! @brief This class is designed to be used for superimposing and subsequently comparing 3D structures for
    //!        molecules based on euclidean and property distance using rigid body transformations.
    //!
    //! @see @link example_chemistry_ligand_pocket_fit_score.cpp @endlink
    //! @author brownbp1
    //! @date Jan 06, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LigandPocketFitScore :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>
    {

    private:

      //! per-residue b-factor of the binding pocket
      mutable storage::Vector< double> m_BFactor;

      //! weight anchoring poses to reference coordinate / binding pocket centroid
      double m_CentroidWeight;

      //! molecule fragment representation of binding pocket
      FragmentComplete m_Pocket;

      //! reference coordinate given as alternative to m_Pocket center
      linal::Vector3D m_StartPosition;

      //! Atom properties to include in the distance
      mutable storage::Vector< descriptor::CheminfoProperty> m_Properties;

      //! Number of atoms to compare in correlation score
      size_t m_NumberOfAtomsToAlign;

      //! Weights for each atom property in the overall correlation
      storage::Vector< double> m_PropertyWeights;

      //! Correlation lengths: distance over which the property influence should decay to half its value
      storage::Vector< double> m_PropertyCorrelationLengths;

    protected:
      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LigandPocketFitScore();

      //! @brief pocket constructor
      explicit LigandPocketFitScore
      (
        const FragmentComplete &POCKET
      );

      //! @brief constructor specify alternative reference coordinate
      //! when pocket center is not desired (i.e. when the full protein protein is given)
      LigandPocketFitScore
      (
        const FragmentComplete &POCKET,
        const linal::Vector3D &START_POSITION
      );

      //! virtual copy constructor
      LigandPocketFitScore *Clone() const;

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

      //! @brief align two small molecule objects and find the property RMSD
      //! @param MOLECULE - molecule to be fit into pocket
      //! @return the molecule in its new pose relative to the static pocket
      double operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief align two small molecule objects and find the property RMSD
      //! @param MOLECULE - molecule to be fit into pocket
      //! @param POCKET - pocket into which MOLECULE is being geometrically fit
      //! @return the molecule in its new pose relative to the static pocket
      double operator()
      (
        const FragmentComplete &MOLECULE,
        const FragmentComplete &POCKET
      ) const;

      ///////////////////
      // Helper Functions
      ///////////////////

      double CollisionScore
      (
        const FragmentComplete &MOL,
        const FragmentComplete &POCK,
        const storage::Vector< double> &BFACTOR = storage::Vector< double>(),
        const bool &STARTING_POSE_BIAS = false
      ) const;

      double PropertyCorrelationScore
      (
        const FragmentComplete &MOL,
        const FragmentComplete &POCK
      ) const;

      //! @brief prepare the class for comparing a conformation
      //! @param MOLECULE the molecule to prepare to compare
      void Prepare( const ConformationInterface &MOLECULE) const;

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

#endif // BCL_CHEMISTRY_LIGAND_POCKET_FIT_SCORE_H_

