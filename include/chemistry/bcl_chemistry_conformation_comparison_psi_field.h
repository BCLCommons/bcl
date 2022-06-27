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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_PSI_FIELD_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_PSI_FIELD_H_

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
#include "math/bcl_math_cubic_spline_damped.h"
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
    //! @class ConformationComparisonPsiField
    //! @brief This class is designed to be used for superimposing and subsequently comparing 3D structures for
    //!        molecules based on euclidean and property distance using rigid body transformations.
    //!
    //! @see @link example_chemistry_conformation_comparison_psi_field.cpp @endlink
    //! @author brownbp1, geanesar, mendenjl
    //! @date Jun 13, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonPsiField :
      public ConformationComparisonPropertyFieldCorrelation
    {
    protected:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! perform align to scaffold in addition to approximator mutates
      bool m_AlignToScaffold;

      //! perform an initial random rotation of molecule A prior to alignment
      bool m_InitialRotation;

      //! number of MCM 1 iterations
      size_t m_IterationNumber;

      //! max number of unimproved MCM 1 iterations
      size_t m_LimitNumber;

      //! the number of trajectories for rigid alignment
      size_t m_NumberRigidTrajectories;

      //! upper distance range for selection of maximum atom distance
      float m_AtomDistanceUpperLimit;

      //! lower distance range for selection of maximum atom distance
      float m_AtomDistanceLowerLimit;

      //! option whether or not to output molecule A
      std::string m_OutputA;

      //! option whether or not to output molecule B
      std::string m_OutputB;

      //! number of unique alignments to output in order of best to worst
      size_t m_NumberOutputSDF;

      //! tolerance for real space symmetry rmsd to compare poses in angstroms
      double m_PoseTolerance;

      //! threshold below which the alignment is accepted as a valid pose
      double m_PoseScoreThreshold;

      //! alignment move relative probabilities
      double m_FlipProb;
      double m_BigRotProb;
      double m_BondSwapProb;
      double m_BondAlignProb;
      double m_BondAlign3Prob;
      double m_BondAlignInfProb;
      double m_SmallRotProb;
      double m_SmallTranslateProb;
      double m_ConformerSwapProb;

      //! Save snapshots from MCM trajectory as .sdf file
      std::string m_Movie;

      //! exclude atoms from molecule A during scoring
      storage::Vector< size_t> m_ExclusionIndicesA;
      std::string m_ExclusionAtomsA;
      mutable storage::Vector< size_t> m_KeepIndicesA;

      //! exclude atoms from molecule B during scoring
      storage::Vector< size_t> m_ExclusionIndicesB;
      std::string m_ExclusionAtomsB;
      mutable storage::Vector< size_t> m_KeepIndicesB;

      // TODO: not fully implemented
      //! Map of statistical pair potential
      storage::Map< storage::Pair< AtomType, AtomType>, math::CubicSplineDamped> m_PairPotentials;

      // TODO: not fully implemented
      //! Filename for statistical pair potential
      std::string m_PotentialFilename;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConformationComparisonPsiField();

      //! @brief full constructor
      ConformationComparisonPsiField
      (
        const bool ALIGN_TO_SCAFFOLD,
        const size_t ITERATIONS,
        const size_t LIMIT,
        const size_t N_RIGID_TRAJ,
        const float ATOM_DIST_UPPER,
        const float ATOM_DIST_LOWER,
        const size_t N_OUTPUTS,
        const double POSE_TOLERANCE,
        const double POSE_SCORE_THRESHOLD
      );

      //! @brief full constructor, default probabilities
      ConformationComparisonPsiField
      (
        const bool ALIGN_TO_SCAFFOLD,
        const size_t ITERATIONS,
        const size_t LIMIT,
        const size_t N_RIGID_TRAJ,
        const float ATOM_DIST_UPPER,
        const float ATOM_DIST_LOWER,
        const size_t N_OUTPUTS,
        const double POSE_TOLERANCE,
        const double POSE_SCORE_THRESHOLD,
        const storage::Vector< size_t> &EXCLUSION_INDICES_A,
        const storage::Vector< size_t> &EXCLUSION_INDICES_B
      );

      //! @brief full constructor
      ConformationComparisonPsiField
      (
        const ConformationComparisonPropertyFieldCorrelation &COMPARER,
        const bool ALIGN_TO_SCAFFOLD,
        const size_t ITERATIONS,
        const size_t LIMIT,
        const size_t N_RIGID_TRAJ,
        const float ATOM_DIST_UPPER,
        const float ATOM_DIST_LOWER,
        const size_t N_OUTPUTS,
        const double POSE_TOLERANCE,
        const double POSE_SCORE_THRESHOLD,
        const float FLIP_PROB,
        const float BIG_ROT_PROB,
        const float BOND_SWAP_PROB,
        const float BOND_ALIGN_PROB,
        const float BOND_ALIGN3_PROB,
        const float BOND_ALIGN_INF_PROB,
        const float SMALL_ROT_PROB,
        const float SMALL_TRANS_PROB,
        const float CONFORMER_SWAP_PROB,
        const storage::Vector< size_t> &EXCLUSION_INDICES_A,
        const storage::Vector< size_t> &EXCLUSION_INDICES_B
      );

      //! virtual copy constructor
      ConformationComparisonPsiField *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief whether to perform substructure-based alignments as part of sampling
      //! @return true if performing substructure-based alignments, else false
      const bool GetAlignToScaffold() const;

      //! @brief whether to perform an initial random rotation of molecule A
      //! @return true if initial rotation is set, else false
      const bool GetInitialRotation() const;

      //! @brief return the number of iterations during the initial alignment
      //! @return the number of iterations
      const size_t GetNumberIterations() const;

      //! @brief return the number of max unimproved iterations allowed in MCM
      //! @return the number of max unimproved iterations
      const size_t GetNumberMaxUnimproved() const;

      //! @brief return the number of rigid alignment trajectories to perform
      //! @return the number of trajectories
      const size_t GetNumberRigidTraj() const;

      //! @brief return the upper limit on the max atom distance range
      //! @return the upper limit on max atom distance range
      const float GetAtomDistUpper() const;

      //! @brief return the lower limit on the max atom distance range
      //! @return the lower limit on max atom distance range
      const float GetAtomDistLower() const;

      //! @brief returns the prefix for the aligned molecule A output SDF file
      //! @return the prefix of the output file for molecule A alignment
      const std::string &GetOutputAPrefix() const;

      //! @brief returns the prefix for the aligned molecule B output SDF file
      //! @return the prefix of the output file for molecule B alignment
      const std::string &GetOutputBPrefix() const;

      //! @brief return the number of maximum unique output alignments
      //! @return the number of outputs
      const size_t GetNumberMaxOutputs() const;

      //! @brief return the tolerance in Angstroms differentiating unique poses
      //! @return the tolerance
      const double GetPoseTolerance() const;

      //! @brief return the threshold in score units accepting a pose
      //! @return the score threshold
      const double GetPoseScoreThreshold() const;

      //! @brief returns the name of the alignment movie output file
      //! @return the movie name
      const std::string &GetMovieName() const;

      //! @brief returns the exclusion indices for molecule A
      //! @return a vector containing indices to be excluded from alignment in molecule A
      const storage::Vector< size_t> &GetExclusionIndicesA() const;

      //! @brief returns the indices to be kept for molecule A
      //! @return a vector containing indices to be included in alignment in molecule A
      const storage::Vector< size_t> &GetKeepIndicesA() const;

      //! @brief returns the exclusion indices for molecule B
      //! @return a vector containing indices to be excluded from alignment in molecule B
      const storage::Vector< size_t> &GetExclusionIndicesB() const;

      //! @brief returns the indices to be kept for molecule B
      //! @return a vector containing indices to be included in alignment in molecule B
      const storage::Vector< size_t> &GetKeepIndicesB() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief align two small molecule objects and find the property RMSD
      //! @param MOLECULE_A - first molecule being aligned
      //! @param MOLECULE_B - second molecule being aligned
      //! @return the RMSDX between first and second molecule
      double operator()
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B
      ) const;

      //! @brief align two small molecule objects and find the property RMSD
      //! @param MOLECULE_A - first molecule being aligned
      //! @param MOLECULE_B - second molecule being aligned
      //! @param ITERATIONS - number of MCM tier 1 iterations
      //! @param MAX_UNIMPROVED - number of allowed unimproved MCM tier 1 iterations
      //! @param RECENTER_FIRST - recenter molecules before MCM sampling
      //! @param MAX_ATOM_DIST - max atom distance determining atom pairs during scoring
      //! @param MOLECULE_A_CONFORMERS = molecule A conformer ensemble
      //! @param POCKET_BOX = bounding box for pocket
      //! @return the new molecule A and RMSDX between the two aligned molecules
      storage::Pair< FragmentComplete, double> FieldOptimizeOrientation
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B,
        const size_t &ITERATIONS,
        const size_t &MAX_UNIMPROVED,
        const bool &RECENTER_FIRST = true,
        const float &MAX_ATOM_DIST = 1.00,
        const FragmentEnsemble &MOLECULE_A_CONFORMERS = FragmentEnsemble(),
        const storage::Vector< size_t> &EXCLUSION_ATOMS_A = storage::Vector< size_t>(),
        const storage::Vector< size_t> &EXCLUSION_ATOMS_B = storage::Vector< size_t>(),
        const FragmentEnsemble &POCKETS = FragmentEnsemble()
      ) const;

      //! @brief align two small molecule objects and find the property RMSD
      //! @param MOLECULE_A - first molecule being aligned
      //! @param MOLECULE_B - second molecule being aligned
      //! @param ITERATIONS - number of MCM tier 1 iterations
      //! @param MAX_UNIMPROVED - number of allowed unimproved MCM tier 1 iterations
      //! @param RECENTER_FIRST - recenter molecules before MCM sampling
      //! @param MAX_ATOM_DIST - max atom distance determining atom pairs during scoring
      //! @param MOLECULE_A_CONFORMERS = molecule A conformer ensemble
      //! @param POCKET_BOX = bounding box for pocket
      //! @return the new molecule A and RMSDX between the two aligned molecules
      storage::Vector< storage::Pair< FragmentComplete, double> > PseudoOperator
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B,
        const size_t &ITERATIONS,
        const size_t &MAX_UNIMPROVED,
        const bool &RECENTER_FIRST = true,
        const float &MAX_ATOM_DIST = 1.00,
        const FragmentEnsemble &MOLECULE_A_CONFORMERS = FragmentEnsemble(),
        const storage::Vector< size_t> &EXCLUSION_ATOMS_A = storage::Vector< size_t>(),
        const storage::Vector< size_t> &EXCLUSION_ATOMS_B = storage::Vector< size_t>(),
        const FragmentEnsemble &POCKETS = FragmentEnsemble()
      ) const;

      ///////////////////
      // Helper Functions
      ///////////////////

      //! @brief whether to perform substructure-based alignments as part of sampling
      void SetAlignToScaffold( const bool ALIGN_TO_SCAFFOLD);

      //! @brief whether to perform an initial random rotation of molecule A
      void SetInitialRotation( const bool INITIAL_ROTATION);

      //! @brief set the number of iterations during the initial alignment
      void SetNumberIterations( const size_t N_ITERATIONS);

      //! @brief set the number of max unimproved iterations allowed in MCM
      void SetNumberMaxUnimproved( const size_t N_MAX_UNIMPROVED);

      //! @brief set the number of rigid alignment trajectories to perform
      void SetNumberRigidTraj( const size_t N_RIGID_TRAJ);

      //! @brief set the upper limit on the max atom distance range
      void SetAtomDistUpper( const float ATOM_DIST_UPPER);

      //! @brief set the lower limit on the max atom distance range
      void SetAtomDistLower( const float ATOM_DIST_LOWER);

      //! @brief set the prefix for the aligned molecule A output SDF file
      void SetOutputAPrefix( const std::string &OUTPUT_A);

      //! @brief set the prefix for the aligned molecule B output SDF file
      void SetOutputBPrefix( const std::string &OUTPUT_B);

      //! @brief set the number of maximum unique output alignments
      void SetNumberMaxOutputs( const size_t N_MAX_OUTPUTS);

      //! @brief set the tolerance in Angstroms differentiating unique poses
      void SetPoseTolerance( const double POSE_TOLERANCE);

      //! @brief set the threshold in score units accepting a pose
      void SetPoseScoreThreshold( const double SCORE_THRESHOLD);

      //! @brief set the relative probability of performing a Flip move during alignment
      void SetFlipProb( const double FLIP_PROB);

      //! @brief set the relative probability of performing a BigRotation move during alignment
      void SetBigRotProb( const double BIG_ROT_PROB);

      //! @brief set the relative probability of performing a BondSwap move during alignment
      void SetBondSwapProb( const double BOND_SWAP_PROB);

      //! @brief set the relative probability of performing a BondSlign move during alignment
      void SetBondAlignProb( const double BOND_ALIGN_PROB);

      //! @brief set the relative probability of performing a BondAlign3 move during alignment
      void SetBondAlign3Prob( const double BOND_ALIGN3_PROB);

      //! @brief set the relative probability of performing a BondAlignInf move during alignment
      void SetBondAlignInfProb( const double BOND_ALIGN_INF_PROB);

      //! @brief set the relative probability of performing a SmallRotation move during alignment
      void SetSmallRotProb( const double SMALL_ROT_PROB);

      //! @brief set the relative probability of performing a SmallTranslation move during alignment
      void SetSmallTransProb( const double SMALL_TRANS_PROB);

      //! @brief set the relative probability of performing a ConformerSwap move during alignment
      void SetConfSwapProb( const double CONF_SWAP_PROB);

      //! @brief set the exclusion indices for molecule A
      void SetExclusionIndicesA( const storage::Vector< size_t> &EXCLUSION_INDICES_A);

      //! @brief set the indices to be kept for molecule A
      void SetKeepIndicesA( const storage::Vector< size_t> &KEEP_INDICES_A);

      //! @brief set the exclusion indices for molecule B
      void SetExclusionIndicesB( const storage::Vector< size_t> &EXCLUSION_INDICES_B);

      //! @brief set the indices to be kept for molecule B
      void SetKeepIndicesB( const storage::Vector< size_t> &KEEP_INDICES_B);

      //! @brief set the name of the alignment movie output file
      void SetMovieName( const std::string &MOVIE_NAME);

      //! @brief perform deterministic substructure-based alignment
      //! @param MOLECULE_A - the molecule being superimposed
      //! @param MOLECULE_B - the scaffold molecule
      //! @param ATOM_TYPE - atom coloring used for substructure alignment
      //! @param BOND_TYPE - bond coloring used for substructure alignment
      //! @return true if successful
      bool ArbitraryScaffoldAlignment
      (
        FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
        const ConfigurationalBondTypeData::Data &BOND_TYPE,
        const storage::Vector< size_t> &TARGET_MOL_INDICES = storage::Vector< size_t>(),
        const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES = storage::Vector< size_t>()
      ) const;

      //! @brief obtain atoms that can be perturbed during alignment
      //! @param MOLECULE_A - first input molecule
      //! @param MOLECULE_B - second input molecule
      //! @param EXCLUSION_INDICES_A - indices of atoms in molecule A to hide during alignment
      //! @param EXCLUSION_INDICES_B - indices of atoms in molecule B to hide during alignment
      //! @return void
      static void GetNonMaskedAtoms
      (
        const FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const storage::Vector< size_t> &EXCLUSION_INDICES_A,
        const storage::Vector< size_t> &EXCLUSION_INDICES_B,
        storage::Vector< size_t> &KEEP_INDICES_A,
        storage::Vector< size_t> &KEEP_INDICES_B
      );

      //! @brief superimpose molecules at all mutually matching atom pairs
      //! @param MOLECULE_A - first input molecule
      //! @param MOLECULE_B - second input molecule
      //! @return true if successful
      static bool BondAlignInf
      (
        FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const storage::Vector< size_t> &KEEP_INDICES_A,
        const storage::Vector< size_t> &KEEP_INDICES_B
      );

      //! @brief scoring function definition for MCM engine
      //! @param MOLECULE_A - first input molecule
      //! @param MOLECULE_B - second input molecule
      //! @param ATOM_DISTANCE - max atom distance determining atom pairs for scoring
      //! @return property RMSD score to MCM engine
      double GetBaseField
      (
        const FragmentComplete &MOL_A,
        const FragmentComplete &MOL_B,
        const float &ATOM_DISTANCE = util::GetUndefined< float>(),
        const storage::Vector< size_t> &EXCLUDE_ATOMS_A = storage::Vector< size_t>(),
        const storage::Vector< size_t> &EXCLUDE_ATOMS_B = storage::Vector< size_t>(),
        const FragmentEnsemble &POCKETS = FragmentEnsemble()
      ) const;

      //! @brief scoring atom pair complementarity between molecule and pocket neighbor atoms
      //! @param MOL_ATOM - first input atom
      //! @param POCKET_ATOM - second input atom
      //! @param DISTANCE - distance between atom pairs
      //! @return atom type complementarity score
      double CheckAtomTypeComplementarity
      (
        const util::SiPtr<const AtomConformationalInterface> &MOL_ATOM,
        const util::SiPtr<const AtomConformationalInterface> &POCKET_ATOM,
        const double &DISTANCE
      ) const;

      //! @brief scoring atom pair complementarity between molecule and pocket neighbor atoms
      //! @param MOL_ATOM - first input atom
      //! @param POCKET_ATOM - second input atom
      //! @param DISTANCE - distance between atom pairs
      //! @return atom type complementarity score
      double CheckPairPotential
      (
        const util::SiPtr<const AtomConformationalInterface> &MOL_ATOM,
        const util::SiPtr<const AtomConformationalInterface> &POCKET_ATOM,
        const double &DISTANCE
      ) const;

      //! @brief find graph isomorphism between the scaffold and small molecule graphs
      //! @param MOLECULE - the molecule being superimposed
      //! @param SCAFFOLD_GRAPH - the scaffold graph
      storage::Map< size_t, size_t> FindIsomorphism
      (
        const FragmentComplete &MOLECULE,
        graph::ConstGraph< size_t, size_t> &SCAFFOLD_GRAPH
      ) const;

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ScoreByField
      //! @brief This helper class is designed to be used to score two molecules by property Mahalanobis distance
      //!        using a ConformationComparisonPropertyFieldCorrelation object
      //!
      //! @author brownbp1, geanesar, mendenjl
      //! @date Jun 23, 2015
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class ScoreByField :
        public math::FunctionInterfaceSerializable< FragmentComplete, double>
      {
      private:
        // Molecule_B constant throughout
        FragmentComplete m_MoleculeB;

        // a pointer to the parent ConformationComparisonPsiField that will instantiate this class
        util::SiPtr< const ConformationComparisonPsiField> m_ParentComparer;

        //! simple point to vector of properties
        util::SiPtr< const storage::Vector< descriptor::CheminfoProperty> > m_Properties;

        //! simple pointer to tracker
        util::SiPtr< const opti::Tracker< FragmentComplete, double> > m_ScoreTracker;

        //! atom indices to exclude from scoring / atom pair matching
        storage::Vector< size_t> m_ExclusionIndicesA;
        storage::Vector< size_t> m_ExclusionIndicesB;

        //! distance cutoff for matched atom pairs
        float m_MaxAtomDistance;

        //! if a pocket is defined then use its voxel grid to compute collisions
        FragmentEnsemble m_Pockets;

      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief constructor
        ScoreByField
        (
          const FragmentComplete &MOLECULE_B,
          const util::SiPtr< const ConformationComparisonPsiField> &PARENT,
          const float &MAX_ATOM_DISTANCE,
          const storage::Vector< size_t> &EXCLUSION_ATOMS_A,
          const storage::Vector< size_t> &EXCLUSION_ATOMS_B,
          const FragmentEnsemble &POCKETS
        );

        //! virtual copy constructor
        ScoreByField *Clone() const;

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! the class name as const ref std::string
        const std::string &GetClassIdentifier() const;

        void SetTracker( const opti::Tracker< FragmentComplete, double> &TRACKER)
        {
          m_ScoreTracker = TRACKER;
        }

        // input
        std::istream &Read( std::istream &ISTREAM);

        // output
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      ////////////////
      // operations //
      ////////////////

        //! @brief score two molecules by property rmsd
        //! @return a score for the approximator
        double operator()( const FragmentComplete &MOLECULE_A) const;
      };

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class TransformationMutation
      //! @brief This helper class is designed to make one of two superimposed molecules undergo a random
      //!        move, e.g. rotation, translation, bondalign, etc.
      //!
      //! @author brownbp1, geanesar, mendenjl
      //! @date Jun 23, 2015
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class TransformationMutation :
        public math::MutateInterface< FragmentComplete>
      {
      private:
          util::SiPtr< const FragmentEnsemble> m_MolEnsA;
          util::SiPtr< const FragmentComplete> m_MolBPtr;
          util::SiPtr< const ConformationComparisonPsiField> m_Parent;
          util::SiPtr< const ConformationComparisonPsiField::ScoreByField> m_ScoreKnower;
          std::string m_MovieFilename;
          util::SiPtr< const opti::Tracker< FragmentComplete, double> > m_Tracker;
          mutable std::string m_Scheme;
      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

         //! @default constructor
         TransformationMutation
         (
           const FragmentEnsemble &MOLECULE_A_CONFORMERS,
           const FragmentComplete &MOL_B,
           const util::SiPtr< const ConformationComparisonPsiField> &PARENT,
           const util::SiPtr< const ConformationComparisonPsiField::ScoreByField> &SCORE_KNOWER,
           const std::string &MOVIE_FILENAME
         );

         //! virtual copy constructor
         TransformationMutation *Clone() const;

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! the class name as const ref std::string
        const std::string &GetClassIdentifier() const;
        //! @brief gets the scheme for this mutate
        //! @return the scheme for this mutate
        virtual const std::string &GetScheme() const
        {
          return m_Scheme;
        }

        //! @brief sets tracker for approximator mutator
        //! @param TRACKER - tracker with argument FragmentComplete and result double
        //! return void
        void SetTracker( const opti::Tracker< FragmentComplete, double> &TRACKER)
        {
          m_Tracker = util::SiPtr< const opti::Tracker< FragmentComplete, double> >( TRACKER);
        }

        // input
        std::istream &Read( std::istream &ISTREAM);

        // output
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      ////////////////
      // operations //
      ////////////////

        //! @brief makes a random rotation or translation to molecule A
        //! @return the mutation result for the approximator
        math::MutateResult< FragmentComplete> operator()( const FragmentComplete &MOLECULE_A) const;
      };

    private:

      //! @brief selects the method with which to link the fragments
      //! @return a string indicating the method to use
      std::string ChooseAlignmentMove( const bool &CONFS) const;

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

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_PSI_FIELD_H_
