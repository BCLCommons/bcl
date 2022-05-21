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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_PSI_FLEX_FIELD_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_PSI_FLEX_FIELD_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "bcl_chemistry_conformation_comparison_psi_field.h"
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
    //! @class ConformationComparisonPsiFlexField
    //! @brief This class is designed to be an expanded version of ConformationComparisonPsiField to
    //!        include comparisons between libraries of molecule conformers.
    //!
    //! @see @link example_chemistry_conformation_comparison_psi_flex_field.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Jul 21, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonPsiFlexField :
      public ConformationComparisonPsiField
    {
    protected:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! number of sampled conformer pairs
      size_t m_PairNumber;

      //! number of best final pairs desired
      size_t m_NumberBestPairsFromEachTraj;

      //! number of MCM 2 iterations
      size_t m_FilterIterations;

      //! max number of unimproved MCM 2 iterations
      size_t m_FilterLimit;

      //! number of MCM 3 iterations
      size_t m_RefinementIterations;

      //! max number of unimproved MCM 3 iterations
      size_t m_RefinementLimit;

      //! number of trajectories to run for best fraction
      size_t m_NumberFlexibleTrajectories;

      //! fraction of best conformer pairs to keep on first filter
      float m_FractionFilterInitial;

      //! fraction of best conformer pairs to keep with each subsequent iteration
      float m_FractionFilterIterations;

      //! align molecule b using only the input conformer
      bool m_RigidMoleculeB;

      //! conformation sampler to generate conformations
      SampleConformations m_SampleConformations;

      //! properties for conformer pair selection
      descriptor::Combine< AtomConformationalInterface, float> m_ConfPairProperties;

      //! assigns cached properties to conformer ensembles
      mutable descriptor::DatasetBuilder< AtomConformationalInterface> m_DatasetBuilder;

    public:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class TripletComparison
      //! @brief Small class that compares the third element of two triplets
      //!
      //! @remarks example unnecessary
      //! @author brownbp1, mendenjl
      //! @date Jul 21, 2015
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class TripletComparison
      {
      public:
        bool operator()( const storage::Triplet< FragmentComplete, FragmentComplete, double> &FIRST, const storage::Triplet< FragmentComplete, FragmentComplete, double> &SECOND)
        {
          return FIRST.Third() < SECOND.Third();
        }
      };

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      ConformationComparisonPsiFlexField();

      //! @brief full constructor
      ConformationComparisonPsiFlexField
      (
        const ConformationComparisonPsiField &PSIFIELD,
        const size_t &N_PAIRS,
        const size_t &N_BEST_PAIRS_EACH_TRAJ,
        const size_t &FILTER_ITERATIONS,
        const size_t &FILTER_LIMIT,
        const size_t &REFINEMENT_ITERATIONS,
        const size_t &REFINMENT_LIMIT,
        const size_t &N_FLEXIBLE_TRAJ,
        const float &FRACTION_FILTER_INITIAL,
        const float &FRACTION_FILTER_ITERATIONS,
        const bool &RIGID_MOL_B,
        const SampleConformations &SAMPLE_CONFS
      );

      //! virtual copy constructor
      ConformationComparisonPsiFlexField *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief return the sample conformations object used for flexible alignment
      //! @return the sample conformers object
      const SampleConformations GetSampleConfs() const;

      //! @brief return the number of conformer pairs
      //! @return the number of conformer pairs
      const size_t GetNumberConformerPairs() const;

      //! @brief return the number of iterations during iterative filter steps
      //! @return the number of filter iterations
      const size_t GetNumberFilterIterations() const;

      //! @brief return the number of max unimproved iterations allowed in the MCM during filter steps
      //! @return the number of max unimproved filter iterations
      const size_t GetNumberFilterMaxUnimproved() const;

      //! @brief return the number of refinement iterations post-filtering
      //! @return the number of refinement iterations
      const size_t GetNumberRefinementIterations() const;

      //! @brief return the number of max unimproved iterations allowed during the refinement MCM
      //! @return the number of max unimproved refinement iterations
      const size_t GetNumberRefinementMaxUnimproved() const;

      //! @brief return the fraction of conformer pairs filtered out after first alignment round
      //! @return the fraction filtered initially
      const float GetFractionFilterInitial() const;

      //! @brief return the fraction of conformer pairs filtered out at each successive alignment iteration
      //! @return the fraction filtered iteratively
      const float GetFractionFilterIterative() const;

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
      //! @param PAIRS - number of starting molecule A/B conformer pairs
      //! @param ITERATIONS - number of MCM tier 1 iterations
      //! @param MAX_UNIMPROVED - number of allowed unimproved iterations
      //! @param FILTER_ITERATIONS - number of MCM tier 2 iterations
      //! @param FILTER_LIMIT - number of allowed unimproved MCM tier 2 iterations
      //! @param REFINEMENT_ITERATIONS - number of MCM tier 3 iterations
      //! @param REFINEMENT_LIMIT - number of allowed unimproved MCM tier 3 iterations
      //! @param FRACTION_FILTER_INITIAL -fraction of conformer pairs filtered out after MCM 1
      //! @param FRACTION_FILTER_ITERATIONS - fraction of conformer pairs filtered out after each round of MCM 2
      //! @return the new molecules A and B, and the RMSDX between the two aligned molecules
      storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > FieldOptimizeOrientationFlex
      (
        const FragmentEnsemble &CONFORMERS_MOL_A,
        const FragmentEnsemble &CONFORMERS_MOL_B,
        const size_t &PAIRS,
        const size_t &ITERATIONS,
        const size_t &MAX_UNIMPROVED,
        const size_t &FILTER_ITERATIONS,
        const size_t &FILTER_LIMIT,
        const size_t &REFINEMENT_ITERATIONS,
        const size_t &REFINEMENT_LIMIT,
        const float &FRACTION_FILTER_INITIAL,
        const float &FRACTION_FILTER_ITERATIONS,
        const bool &CENTER_ACONFS_ON_BCONFS = true,
        const FragmentEnsemble &POCKETS = FragmentEnsemble()
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief align two small molecule objects and find the property RMSD
      //! @param MOL - molecule for which to generate conformers
      //! @return conformer ensemble of the input molecule
      FragmentEnsemble MakeConformers( const FragmentComplete &MOL) const;

      static linal::Vector< double> PickConformerPairs( const FragmentComplete &MOL, const FragmentEnsemble &ENS);

      //! @brief set the sample conformations object
      void SetSampleConformations( const SampleConformations &SAMPLE_CONFS);

      //! @brief set the number of conformer pairs
      void SetNumberConformerPairs( const size_t N_CONF_PAIRS);

      //! @brief set the number of iterations during iterative filter steps
      void SetNumberFilterIterations( const size_t N_FILTER_ITERATIONS);

      //! @brief set the number of max unimproved iterations allowed in the MCM during filter steps
      void SetNumberFilterMaxUnimproved( const size_t N_FILTER_MAX_UNIMPROVED);

      //! @brief set the number of refinement iterations post-filtering
      void SetNumberRefinementIterations( const size_t N_REFINE_ITERATIONS);

      //! @brief set the number of max unimproved iterations allowed during the refinement MCM
      void SetNumberRefinementMaxUnimproved( const size_t N_REFINE_MAX_UNIMPROVED);

      //! @brief set the fraction of conformer pairs filtered out after first alignment round
      void SetFractionFilterInitial( const float FILTER_FRACTION_INITIAL);

      //! @brief set the fraction of conformer pairs filtered out at each successive alignment iteration
      void SetFractionFilterIterative( const float FILTER_FRACTION_ITERATIVE);

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

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_PSI_FLEX_FIELD_H_
