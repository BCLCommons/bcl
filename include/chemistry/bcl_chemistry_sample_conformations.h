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

#ifndef BCL_CHEMISTRY_SAMPLE_CONFORMATIONS_H_
#define BCL_CHEMISTRY_SAMPLE_CONFORMATIONS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_descriptor_to_score_adaptor.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_mutate_clash_resolver.h"
#include "bcl_chemistry_search_fragment_library_from_tree.h"
#include "bcl_chemistry_small_molecule_fragment_mapping.h"
#include "descriptor/bcl_descriptor_base.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SampleConformations
    //! @brief Used for attaching a grow fragment to a base fragment
    //!
    //! @see @link example_chemistry_sample_conformations.cpp @endlink
    //! @author kothiwsk, mendenjl, brownbp1
    //! @date Oct 14, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SampleConformations :
      public util::SerializableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! object that searches for fragments of the molecule of interest from the rotamer library
      util::ShPtr< SearchFragmentLibraryFromTree>                 m_SearchFragmentLibrary;

      //! method for comparing conformations
      util::Implementation< ConformationComparisonInterface>      m_ConformationComparison;

      //! tolerance for conformation comparer
      double                                                      m_Tolerance;

      //! number of conformations desired
      size_t                                                      m_NumberConformations;

      //! number of iterations for searching conformers of a molecule
      size_t                                                      m_MaxIterations;

      //! change chirality
      bool                                                        m_ChangeChirality;

      //! relative probability of just making random dihedral swaps
      double                                                      m_RandomDihedralChangeProbability;

      //! boolean to specify if 3D structure needs to be generated from scratch
      bool                                                        m_Generate3D;

      //! whether to cluster the returned conformations so as to nominally optimize coverage of conformational space.
      bool                                                        m_Cluster;

      //! Number of conformations to choose entirely based on diversity from existing structures
      size_t                                                      m_ClusterNBasedOnDiversity;

      //! Optional filename where all scores of all fragments can be written out
      std::string                                                 m_OutputIndividualScoresFilename;

      //! true to allow sampling of dihedral angles
      bool                                                        m_SampleDihedralAngles;

      //! true to allow sampling of ring conformations
      bool                                                        m_SampleRingConformations;

      //! true to allow sampling of bond angles and lengths
      bool                                                        m_SampleBondAnglesAndLengths;

      //! clash resolver
      double                                                      m_MaxNumberClashResolutionCycles;

      //! Absolute maximum clash score
      double                                                      m_MaxClashScore;

      //! How to handle SampleByParts indices set in SDF MDL property line
      enum SampleByPartsOptions
      {
        e_IgnoreAll,      //! Do not use SampleByparts at all
        e_IgnoreOnlyMDL,  //! Ignore MDL; use only command-line set options
        e_UseOnlyMDL,     //! Use only MDL; ignore command-line options
        e_UseAll          //! Combine command-line settings with MDL specification
      };
      SampleByPartsOptions m_SampleByPartsHandler;

      //! MDL property string for sample by parts indices
      std::string                                                 m_MDL;

      //! MDL property line sample by parts;
      mutable linal::Vector< size_t>                              m_MDLIndicesToSample;

      //! Explicitly set indices to include in sample by parts
      mutable storage::Vector< size_t>                            m_IndicesToSample;
      std::string                                                 m_IndicesToSampleString;

      //! Fragments whose atom indices to include in sample by parts
      mutable storage::Vector< size_t>                            m_FragmentsToSample;
      std::string                                                 m_ReferenceFragments;

      //! Combined SampleByParts indices
      mutable linal::Vector< size_t>                              m_SampleByPartsIndices;

      //! mutex for file reading
      mutable sched::Mutex                                        m_Mutex;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SampleConformations();

      //! @brief constructor from provided arguments
      //! @param ROTAMER_LIBRARY rotamer library which contains fragment conformations used
      //! @param CONFORMATION_COMPARER method to be used for comparing conformations
      //! @param TOLERANCE tolerance to be used for differentiating conforme
      //! @param NUMBER_CONFORMATIONS number of conformations desired
      //! @param MAX_ITERATIONS maximum number of iterations to perform
      //! @param CHANGE_CHIRALITY if type of chirality should be changed at chiral centers
      //! @param RANDOM_DIHEDRAL_CHANGE_CHANCE probability of swapping a dihedral at random
      //! @param GENERATE_3D whether or not 3D structure should be generated from scratch
      //! @param CLASH_TOLERANCE tolerance for vdw overlaps of atoms
      //! @param CLUSTER perform leader-follower clustering during sampling
      //! @param NUMBER_CLASH_RESOLUTION_CYCLES cycles of clash resolver
      //! @oaram SCORES_FILENAME optional file where all scores of all fragments can be written
      SampleConformations
      (
        const RotamerLibraryInterface &ROTAMER_LIBRARY,
        const std::string &CONFORMATION_COMPARER,
        const double &TOLERANCE,
        const size_t NUMBER_CONFORMATIONS,
        const size_t &MAX_ITERATIONS,
        const bool CHANGE_CHIRALITY,
        const double &RANDOM_DIHEDRAL_CHANGE_CHANCE = double( 0.0),
        bool GENERATE_3D = false,
        double CLASH_TOLERANCE = 0.1,
        bool CLUSTER = true,
        double NUMBER_CLASH_RESOLUTION_CYCLES = 0.1,
        std::string SCORES_FILENAME = ""
      );

      //! @brief clone constructor
      SampleConformations *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetAlias() const;

      //! fragment rotamer library
      const util::Implementation< RotamerLibraryInterface> &GetRotamerLibrary() const;

      //! method for comparing conformations
      const std::string &GetConfComparer() const;

      //! tolerance for conformation comparer
      const double &GetTolerance() const;

      //! number of conformations desired
      const size_t &GetMaxNumConfs() const;

      //! number of iterations for searching conformers of a molecule
      const size_t &GetMaxIterations() const;

      //! change chirality
      const bool &GetIfChangeChirality() const;

      //! relative probability of just making random dihedral swaps
      const double &GetRandDihedChangeProb() const;

      //! boolean to specify if 3D structure needs to be generated from scratch
      const bool &GetIfGenerate3D() const;

      //! whether to cluster the returned conformations so as to nominally optimize coverage of conformational space.
      const bool &GetIfCluster() const;

      //! Number of conformations to choose entirely based on diversity from existing structures
      const size_t &GetClusterNDiversity() const;

      //! clash resolver
      const double &GetClashResolution() const;

      //! Absolute maximum clash score
      const double &GetMaxClashScore() const;

    ////////////////
    // operations //
    ////////////////

      //! method for comparing conformations
      void SetConfComparer( const std::string &COMPARER);

      //! tolerance for conformation comparer
      void SetTolerance( const double TOLERANCE);

      //! number of conformations desired
      void SetMaxNumConfs( const size_t N_MAX_CONFS);

      //! number of iterations for searching conformers of a molecule
      void SetMaxIterations( const size_t N_MAX_ITERATIONS);

      //! change chirality
      void SetIfChangeChirality( const bool CHANGE_CHIRALITY);

      //! relative probability of just making random dihedral swaps
      void SetRandDihedChangeProb( const double RAND_DIHED_PROB);

      //! boolean to specify if 3D structure needs to be generated from scratch
      void SetIfGenerate3D( const bool GENERATE_3D);

      //! whether to cluster the returned conformations so as to nominally optimize coverage of conformational space.
      void SetIfCluster( const bool CLUSTER);

      //! Number of conformations to choose entirely based on diversity from existing structures
      void SetClusterNDiversity( const size_t N_FORCED_DIVERSE_CLUSTERS);

      //! clash resolver
      void SetClashResolution( const double CLASH_RESOLUTION);

      //! Absolute maximum clash score
      void SetMaxClashScore( const double CLASH_SCORE);

      //! @brief change the output score filename
      void SetOutputScoreFilename( const std::string &FILENAME);

      //! @brief set the sampling types to limit or expand the types of conformational sampling will be added
      //! @param SAMPLE_DIHEDRALS whether dihedral angles can be substantially changed (substantially outside the 30 degree bins)
      //! @param SAMPLE_RINGS whether to sample ring conformations
      //! @param SAMPLE_BOND_ANGLES whether to sample bond angles (and lengths)
      //! @param SAMPLE_CHIRALITY whether to sample chirality
      void SetSamplingPreferences
      (
        const bool &SAMPLE_DIHEDRALS,
        const bool &SAMPLE_RINGS,
        const bool &SAMPLE_BOND_ANGLES,
        const bool &SAMPLE_CHIRALITY
      );

      //! @brief set sample by parts indices from MDL property line
      void SetSampleByPartsMDL
      (
        const FragmentComplete &MOLECULE,
        const std::string &MDL = "SampleByParts"
      ) const;

      //! @brief set sample by parts indices
      void SetSampleByPartsIndices( const storage::Vector< size_t> &INDICES_TO_SAMPLE) const;

      //! @brief set sample by parts indices
      void SetSampleByPartsIndices( const std::string &INDICES_TO_SAMPLE) const;

      //! @brief set the sample by parts atoms based on fragments
      void SetSampleByPartsFragments
      (
        const FragmentComplete &FRAGMENT,
        const FragmentEnsemble &FRAGMENTS_TO_SAMPLE,
        const bool COMPLEMENT
      ) const;

      //! @brief set the sample by parts atoms based on fragments
      void SetSampleByPartsFragments
      (
        const FragmentComplete &FRAGMENT,
        const std::string &FRAGMENTS_TO_SAMPLE,
        const bool COMPLEMENT
      ) const;

      //! @brief sets the final sample by parts atoms used by operator() from member data
      void PoolSampleByPartsIndices() const;

      //! @brief set the filter criteria for determining what to keep from the final ensemble
      //! @param ENSEMBLE conformers to be filtered
      //! @param STARTING_MOL molecule (conformer) against which the ensemble will be filtered
      //! @param CONFORMATION_COMPARER the measure of similarity to use
      //! @param COMPARISON are the conformers less, greater, etc. than TOLERANCE according
      //! to CONFORMATION_COMPARER to starting conformer
      //! @param TOLERANCE comparison threshold for conformer relation
      //! @param ATOMS_TO_COMPARE atom subset to compare
      FragmentEnsemble FilterConformerEnsemble
      (
        const FragmentEnsemble &ENSEMBLE,
        const FragmentComplete &STARTING_MOL,
        const std::string &CONFORMATION_COMPARER,
        const math::Comparisons< float>::Comparison &COMPARISON,
        const float &TOLERANCE,
        const storage::Vector< size_t> &ATOMS_TO_COMPARE
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief returns an ensemble of generated conformers
      //! @param MOLECULE molecule of interest whose conformations need to be sampled
      //! @return an ensemble of generated conformers
      storage::Pair< FragmentEnsemble, FragmentEnsemble> operator()( const FragmentComplete &MOLECULE) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief applies mutates to generate conformers for molecule
      //! @param STARTING_ARGUMENT the argument that will be the first for sampling
      //! @param CONFORMATION_SET conformation set to store results that are produced on application of mutates
      //! @param MUTATE mutates that need to be applied to the molecule of interest
      //! @param SCORE_FUNCTION object that scores a the mutate result
      //! @param MAX_ITERATIONS maximum number of iterations during monte carlo sampling
      //! @param TEMPERATURE temperature to make decision of accepting a move
      //! @param MAX_RESULTS desired number of new mutate results that are accepted/improved
      void Sampling
      (
        const FragmentComplete &STARTING_ARGUMENT,
        FragmentEnsemble &CONFORMATION_SET,
        const util::ShPtr< math::MutateInterface< FragmentComplete> > &MUTATE,
        const FragmentProbabilityScore &SCORE,
        const size_t &MAX_ITERATIONS,
        const double &CLASH_TOLERANCE,
        bool &SHARED_CHOOSE_ONLY_BEST,
        AtomClashScore &SHARED_CLASH_SCORE
      ) const;

      //! @brief determine if the coordinates of atoms of a molecule are defined
      //! @param MOLECULE molecule to check for
      //! @result true if coordinates are defined, false otherwise
      bool HasUndefinedCoordinates( const ConformationInterface &MOLECULE) const;

      //! @brief determine if the coordinates of atoms of a molecule are defined
      //! @param MOLECULE molecule whose rings are desired
      //! @param RINGS small fragment isomorphisms of rings that need to be checked for wholeness
      //! @result true if coordinates are defined, false otherwise
      bool CheckForWholeRings
      (
        const FragmentComplete &MOLECULE,
        const util::ShPtr< FragmentEnsemble> &FRAGMENTS
      ) const;

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

    protected:

      //! @brief get the relative probability that a conformation was the closest-to-native
      //!        estimated from the best conf score and its conf score
      //! @param BEST_CONF_SCORE best conf-score observed for this molecule
      //! @param CONF_SCORE actual conf score for the molecule
      static double GetRelativeProbabilityOfConformerBeingClosestToNative( const double &BEST_CONF_SCORE, const double &CONF_SCORE);

    }; // class SampleConformations

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_SAMPLE_CONFORMATIONS_H_
