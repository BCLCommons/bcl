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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_real_space_rmsd.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_symmetry_rmsd.h"
#include "chemistry/bcl_chemistry_conformation_comparison_multi_align.h"
#include "chemistry/bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_fragment_align_to_scaffold.h"
#include "chemistry/bcl_chemistry_ligand_pocket_fit_score.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "quality/bcl_quality_rmsd.h"
#include "score/bcl_score.h"
#include "storage/bcl_storage_pair.h"

using bcl::io::File;

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonPsiField::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonPsiField())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ConformationComparisonPsiField::ConformationComparisonPsiField() :
      ConformationComparisonPropertyFieldCorrelation(),
      m_AlignToScaffold( false),
      m_InitialRotation( false),
      m_IterationNumber( 400),
      m_LimitNumber( 160),
      m_NumberRigidTrajectories( 5),
      m_AtomDistanceUpperLimit( 1.15),
      m_AtomDistanceLowerLimit( 0.70),
      m_NumberOutputSDF( 1),
      m_PoseTolerance( 0.125),
      m_PoseScoreThreshold( 2.0),
      m_FlipProb( 0.06),
      m_BigRotProb( 0.06),
      m_BondSwapProb( 0.06),
      m_BondAlignProb( 0.30),
      m_BondAlign3Prob( 0.20),
      m_BondAlignInfProb( 0.10),
      m_SmallRotProb( 0.12),
      m_SmallTranslateProb( 0.12),
      m_ConformerSwapProb( 0.06),
      m_ExclusionIndicesA( storage::Vector< size_t>()),
      m_ExclusionIndicesB( storage::Vector< size_t>()),
      m_PairPotentials( storage::Map< storage::Pair< AtomType, AtomType>, math::CubicSplineDamped>())
    {
    }

    //! @brief full constructor
    ConformationComparisonPsiField::ConformationComparisonPsiField
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
    ) :
      ConformationComparisonPropertyFieldCorrelation(),
      m_AlignToScaffold( ALIGN_TO_SCAFFOLD),
      m_InitialRotation( false),
      m_IterationNumber( ITERATIONS),
      m_LimitNumber( LIMIT),
      m_NumberRigidTrajectories( N_RIGID_TRAJ),
      m_AtomDistanceUpperLimit( ATOM_DIST_UPPER),
      m_AtomDistanceLowerLimit( ATOM_DIST_LOWER),
      m_NumberOutputSDF( N_OUTPUTS),
      m_PoseTolerance( POSE_TOLERANCE),
      m_PoseScoreThreshold( POSE_SCORE_THRESHOLD),
      m_FlipProb( 0.06),
      m_BigRotProb( 0.06),
      m_BondSwapProb( 0.06),
      m_BondAlignProb( 0.30),
      m_BondAlign3Prob( 0.20),
      m_BondAlignInfProb( 0.10),
      m_SmallRotProb( 0.12),
      m_SmallTranslateProb( 0.12),
      m_ConformerSwapProb( 0.06),
      m_ExclusionIndicesA( storage::Vector< size_t>()),
      m_ExclusionIndicesB( storage::Vector< size_t>()),
      m_PairPotentials( storage::Map< storage::Pair< AtomType, AtomType>, math::CubicSplineDamped>())
    {
    }

    //! @brief full constructor, default probabilities
    ConformationComparisonPsiField::ConformationComparisonPsiField
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
    ) :
      ConformationComparisonPropertyFieldCorrelation(),
      m_AlignToScaffold( ALIGN_TO_SCAFFOLD),
      m_InitialRotation( false),
      m_IterationNumber( ITERATIONS),
      m_LimitNumber( LIMIT),
      m_NumberRigidTrajectories( N_RIGID_TRAJ),
      m_AtomDistanceUpperLimit( ATOM_DIST_UPPER),
      m_AtomDistanceLowerLimit( ATOM_DIST_LOWER),
      m_NumberOutputSDF( N_OUTPUTS),
      m_PoseTolerance( POSE_TOLERANCE),
      m_PoseScoreThreshold( POSE_SCORE_THRESHOLD),
      m_FlipProb( 0.06),
      m_BigRotProb( 0.06),
      m_BondSwapProb( 0.06),
      m_BondAlignProb( 0.30),
      m_BondAlign3Prob( 0.20),
      m_BondAlignInfProb( 0.10),
      m_SmallRotProb( 0.12),
      m_SmallTranslateProb( 0.12),
      m_ConformerSwapProb( 0.06),
      m_ExclusionIndicesA( EXCLUSION_INDICES_A),
      m_ExclusionIndicesB( EXCLUSION_INDICES_B),
      m_PairPotentials( storage::Map< storage::Pair< AtomType, AtomType>, math::CubicSplineDamped>())
    {
    }

    //! @brief full constructor
    ConformationComparisonPsiField::ConformationComparisonPsiField
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
    ) :
      ConformationComparisonPropertyFieldCorrelation( COMPARER),
      m_AlignToScaffold( ALIGN_TO_SCAFFOLD),
      m_InitialRotation( false),
      m_IterationNumber( ITERATIONS),
      m_LimitNumber( LIMIT),
      m_NumberRigidTrajectories( N_RIGID_TRAJ),
      m_AtomDistanceUpperLimit( ATOM_DIST_UPPER),
      m_AtomDistanceLowerLimit( ATOM_DIST_LOWER),
      m_NumberOutputSDF( N_OUTPUTS),
      m_PoseTolerance( POSE_TOLERANCE),
      m_PoseScoreThreshold( POSE_SCORE_THRESHOLD),
      m_FlipProb( FLIP_PROB),
      m_BigRotProb( BIG_ROT_PROB),
      m_BondSwapProb( BOND_SWAP_PROB),
      m_BondAlignProb( BOND_ALIGN_PROB),
      m_BondAlign3Prob( BOND_ALIGN3_PROB),
      m_BondAlignInfProb( BOND_ALIGN_INF_PROB),
      m_SmallRotProb( SMALL_ROT_PROB),
      m_SmallTranslateProb( SMALL_TRANS_PROB),
      m_ConformerSwapProb( CONFORMER_SWAP_PROB),
      m_ExclusionIndicesA( EXCLUSION_INDICES_A),
      m_ExclusionIndicesB( EXCLUSION_INDICES_B),
      m_PairPotentials( storage::Map< storage::Pair< AtomType, AtomType>, math::CubicSplineDamped>())
    {
    }

    //! virtual copy constructor
    ConformationComparisonPsiField *ConformationComparisonPsiField::Clone() const
    {
      return new ConformationComparisonPsiField( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConformationComparisonPsiField::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ConformationComparisonPsiField::GetAlias() const
    {
      static std::string s_name( "PsiField");
      return s_name;
    }

    //! @brief whether to perform substructure-based alignments as part of sampling
    //! @return true if performing substructure-based alignments, else false
    const bool ConformationComparisonPsiField::GetAlignToScaffold() const
    {
      return m_AlignToScaffold;
    }

    //! @brief whether to perform an initial random rotation of molecule A
    //! @return true if initial rotation is set, else false
    const bool ConformationComparisonPsiField::GetInitialRotation() const
    {
      return m_InitialRotation;
    }

    //! @brief return the number of iterations during the initial alignment
    //! @return the number of iterations
    const size_t ConformationComparisonPsiField::GetNumberIterations() const
    {
      return m_IterationNumber;
    }

    //! @brief return the number of max unimproved iterations allowed in MCM
    //! @return the number of max unimproved iterations
    const size_t ConformationComparisonPsiField::GetNumberMaxUnimproved() const
    {
      return m_LimitNumber;
    }

    //! @brief return the number of rigid alignment trajectories to perform
    //! @return the number of trajectories
    const size_t ConformationComparisonPsiField::GetNumberRigidTraj() const
    {
      return m_NumberRigidTrajectories;
    }

    //! @brief return the upper limit on the max atom distance range
    //! @return the upper limit on max atom distance range
    const float ConformationComparisonPsiField::GetAtomDistUpper() const
    {
      return m_AtomDistanceUpperLimit;
    }

    //! @brief return the lower limit on the max atom distance range
    //! @return the lower limit on max atom distance range
    const float ConformationComparisonPsiField::GetAtomDistLower() const
    {
      return m_AtomDistanceLowerLimit;
    }

    //! @brief returns the prefix for the aligned molecule A output SDF file
    //! @return the prefix of the output file for molecule A alignment
    const std::string &ConformationComparisonPsiField::GetOutputAPrefix() const
    {
      return m_OutputA;
    }

    //! @brief returns the prefix for the aligned molecule B output SDF file
    //! @return the prefix of the output file for molecule B alignment
    const std::string &ConformationComparisonPsiField::GetOutputBPrefix() const
    {
      return m_OutputB;
    }

    //! @brief return the number of maximum unique output alignments
    //! @return the number of outputs
    const size_t ConformationComparisonPsiField::GetNumberMaxOutputs() const
    {
      return m_NumberOutputSDF;
    }

    //! @brief return the tolerance in Angstroms differentiating unique poses
    //! @return the tolerance
    const double ConformationComparisonPsiField::GetPoseTolerance() const
    {
      return m_PoseTolerance;
    }

    //! @brief return the threshold in score units accepting a pose
    //! @return the score threshold
    const double ConformationComparisonPsiField::GetPoseScoreThreshold() const
    {
      return m_PoseScoreThreshold;
    }

    //! @brief returns the name of the alignment movie output file
    //! @return the movie name
    const std::string &ConformationComparisonPsiField::GetMovieName() const
    {
      return m_Movie;
    }

    //! @brief returns the exclusion indices for molecule A
    //! @return a vector containing indices to be excluded from alignment in molecule A
    const storage::Vector< size_t> &ConformationComparisonPsiField::GetExclusionIndicesA() const
    {
      return m_ExclusionIndicesA;
    }

    //! @brief returns the indices to be kept for molecule A
    //! @return a vector containing indices to be included in alignment in molecule A
    const storage::Vector< size_t> &ConformationComparisonPsiField::GetKeepIndicesA() const
    {
      return m_KeepIndicesA;
    }

    //! @brief returns the exclusion indices for molecule B
    //! @return a vector containing indices to be excluded from alignment in molecule B
    const storage::Vector< size_t> &ConformationComparisonPsiField::GetExclusionIndicesB() const
    {
      return m_ExclusionIndicesB;
    }

    //! @brief returns the indices to be kept for molecule B
    //! @return a vector containing indices to be included in alignment in molecule B
    const storage::Vector< size_t> &ConformationComparisonPsiField::GetKeepIndicesB() const
    {
      return m_KeepIndicesB;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the property RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_B - second molecule being aligned
    //! @return the RMSD between first and second molecule

    double ConformationComparisonPsiField::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      // we need to get the keep_indices from exclusion_indices for the scaffold-based alignment
      storage::Vector< size_t> keep_indices_a, keep_indices_b;
      GetNonMaskedAtoms( MOLECULE_A, MOLECULE_B, m_ExclusionIndicesA, m_ExclusionIndicesB, keep_indices_a, keep_indices_b);
      m_KeepIndicesA = keep_indices_a;
      m_KeepIndicesB = keep_indices_b;

      // manage max atom distance
      float dmax_increment( ( m_AtomDistanceUpperLimit - m_AtomDistanceLowerLimit) / ( m_NumberRigidTrajectories - 1));
      storage::Vector< storage::Pair< FragmentComplete, double> > traj_tracker;
      storage::Vector< double> distances;

      //perform m_NumberRigidTrajectories number of independent rigid alignments
      for( size_t traj( 0); traj < m_NumberRigidTrajectories; ++traj)
      {
        // cache properties MOLECULE_A
        FragmentComplete mol_a( MOLECULE_A);
        mol_a.ShareCache( MOLECULE_A);

        //center mol_a on MOLECULE_B
        mol_a.Translate( MOLECULE_B.GetCenter() - mol_a.GetCenter());

        //randomly rotate off of MOLECULE_B
        if( m_InitialRotation)
        {
          math::RotationMatrix3D rot;
          rot.SetRand();
          linal::Vector3D mol_a_centered( mol_a.GetCenter());
          mol_a.Translate( -mol_a_centered);
          mol_a.Rotate( rot);
          mol_a.Translate( mol_a_centered);
        }

        //perform alignments
        storage::Pair< FragmentComplete, double> aligned;

        float max_atom_distance( m_AtomDistanceLowerLimit + ( dmax_increment * traj));
        distances.PushBack( max_atom_distance);
        aligned = FieldOptimizeOrientation( mol_a, MOLECULE_B, m_IterationNumber, m_LimitNumber, false, max_atom_distance, FragmentEnsemble(),
          m_ExclusionIndicesA, m_ExclusionIndicesB);

        //keep track of independent alignment trajectories
        storage::Pair< FragmentComplete, double> best_mol( aligned.First(), aligned.Second());
        traj_tracker.PushBack( best_mol);
      }

      if( m_AlignToScaffold)
      {
        // cache properties MOLECULE_A
        FragmentComplete mol_a( MOLECULE_A);
        mol_a.ShareCache( MOLECULE_A);
        storage::Vector< storage::Pair< ConformationGraphConverter::AtomComparisonTypeEnum, ConfigurationalBondTypeData::DataEnum> > ats_options;

        static storage::List< util::Stopwatch> stopwatches;
        if( stopwatches.IsEmpty())
        {
          for( size_t e_atom( 0); e_atom < ConformationGraphConverter::s_NumberAtomComparisonTypes; ++e_atom)
          {
            for( size_t e_bond( 1); e_bond < ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness; ++e_bond)
            {
              stopwatches.PushBack
              (
                util::Stopwatch
                (
                  ConformationGraphConverter::GetAtomComparisonType( ConformationGraphConverter::AtomComparisonType( e_atom))
                  + "_" + ConfigurationalBondTypeData::GetDataName( static_cast< ConfigurationalBondTypeData::Data>( e_bond))
                  + " Timing",
                  util::Time( 10, 0),
                  util::Message::e_Verbose,
                  true,
                  false
                )
              );
            }
          }
        }
        auto itr_stop( stopwatches.Begin());
        for( size_t e_atom( 0); e_atom < ConformationGraphConverter::s_NumberAtomComparisonTypes; ++e_atom)
        {
          for( size_t e_bond( 1); e_bond < ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness; ++e_bond, ++itr_stop)
          {
            itr_stop->Start();
            ArbitraryScaffoldAlignment
            (
              mol_a,
              MOLECULE_B,
              static_cast< ConformationGraphConverter::AtomComparisonType>( e_atom),
              static_cast< ConfigurationalBondTypeData::Data>( e_bond),
              m_KeepIndicesA,
              m_KeepIndicesB
            );

            float score_a
            (
              this->PropertyDistanceScorer
              (
                mol_a,
                MOLECULE_B,
                m_AtomDistanceUpperLimit,
                m_ExclusionIndicesA,
                m_ExclusionIndicesB
              ).First()
            );
            traj_tracker.PushBack();
            traj_tracker.LastElement().First() = mol_a;
            traj_tracker.LastElement().Second() = score_a;
            itr_stop->Stop();
          }
        }
      }

      // Formerly just evaluated the lower, middle, and upper bounds on the max atom distances;
      // now averages across all sampled distances
      storage::Vector< storage::Pair< double, size_t> > sorter( traj_tracker.GetSize());
      size_t tracker_index( 0);
      for
      (
        storage::Vector< storage::Pair< FragmentComplete, double >>::iterator itr( traj_tracker.Begin()), itr_end( traj_tracker.End());
        itr < itr_end;
        ++itr, ++tracker_index
      )
      {
        math::RunningAverage< double> sum( 0.0);
        for( size_t i( 0); i < distances.GetSize(); ++i)
        {
          sum += PropertyDistanceScorer( itr->First(), MOLECULE_B, distances( i), m_ExclusionIndicesA, m_ExclusionIndicesB).First();
        }
        double ave( sum.GetAverage());
        sorter( tracker_index).First() = ave;
        itr->First().StoreProperty( "PropertyFieldDistance", util::Format()( ave));
        sorter( tracker_index).Second() = tracker_index;
      }

      // Sort molecules by best
      sorter.Sort( std::less< storage::Pair< double, size_t> >());
      storage::Vector< size_t> best_order( traj_tracker.GetSize());
      for( size_t i( 0); i < traj_tracker.GetSize(); ++i)
      {
        best_order( i) = sorter( i).Second();
      }
      traj_tracker.Reorder( best_order);

      // prune similar molecules
      ConformationComparisonBySymmetryRmsd compare_symrmsd( false, 100); // consider real space difference
      storage::Vector< storage::Pair< FragmentComplete, double> > pruned_pair_scores;

      // add first molecule to new molecules collection
      pruned_pair_scores.PushBack( traj_tracker( 0));

      // check other alignment poses; if they deviate by m_PoseTolerance angstroms then they are different and save
      storage::Vector< size_t> score_indices;
      for( size_t i( 0); i < traj_tracker.GetSize(); ++i)
      {
        // loop over the saved poses
        bool same( false);
        for( size_t j( 0); j < pruned_pair_scores.GetSize(); ++j)
        {
          if( compare_symrmsd( traj_tracker( i).First(), pruned_pair_scores( j).First()) <= m_PoseTolerance)
          {
            // outer loop contains something same as already in new collection
            same = true;
            break;
          }
        }
        if( !same && sorter( i).First() < m_PoseScoreThreshold)
        {
          pruned_pair_scores.PushBack
          (
            storage::Pair< FragmentComplete, double>(
              std::make_pair
              (
                traj_tracker( i).First(),
                sorter( i).First()
              )
            )
          );
          score_indices.PushBack( i);
        }
      }

      // write results
      static sched::Mutex s_mutex;
      s_mutex.Lock();
      for( size_t best_pairs_index( 0); best_pairs_index < m_NumberOutputSDF; ++best_pairs_index)
      {
        if( best_pairs_index < pruned_pair_scores.GetSize())
        {
          if( !m_OutputA.empty())
          {
            io::OFStream out;
            io::File::MustOpenOFStream( out, m_OutputA + "_" + std::to_string( best_pairs_index) + ".sdf", std::ios::app);
            pruned_pair_scores( best_pairs_index).First().StoreProperty( "PropertyFieldDistance", util::Format()( pruned_pair_scores( best_pairs_index).Second()));
            pruned_pair_scores( best_pairs_index).First().WriteMDL( out);
            io::File::CloseClearFStream( out);
          }
          if( !m_OutputB.empty())
          {
            io::OFStream out;
            io::File::MustOpenOFStream( out, m_OutputB + "_" + std::to_string( best_pairs_index) + ".sdf", std::ios::app);
            FragmentComplete mol_b( MOLECULE_B);
            mol_b.WriteMDL( out);
            io::File::CloseClearFStream( out);
          }
        }
        s_mutex.Unlock();
      }
      return traj_tracker( 0).Second();
    }

    storage::Pair< FragmentComplete, double> ConformationComparisonPsiField::FieldOptimizeOrientation
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B,
      const size_t &ITERATIONS,
      const size_t &MAX_UNIMPROVED,
      const bool &RECENTER_FIRST,
      const float &MAX_ATOM_DIST,
      const FragmentEnsemble &MOLECULE_A_CONFORMERS,
      const storage::Vector< size_t> &EXCLUSION_ATOMS_A,
      const storage::Vector< size_t> &EXCLUSION_ATOMS_B,
      const FragmentEnsemble &POCKETS
    ) const
    {

      // cache properties MOLECULE_B
      FragmentComplete mol_b( MOLECULE_B);
      mol_b.ShareCache( MOLECULE_B);

      // cache properties MOLECULE_A
      FragmentComplete mol_a( MOLECULE_A);
      mol_a.ShareCache( MOLECULE_A);

      // score object for use by the approximator
      util::SiPtr< const ConformationComparisonPsiField> compare_ptr( this);
      ScoreByField score( mol_b, compare_ptr, MAX_ATOM_DIST, EXCLUSION_ATOMS_A, EXCLUSION_ATOMS_B, POCKETS);

      // mutation object for use by the approximator
      TransformationMutation mutator( MOLECULE_A_CONFORMERS, mol_b, compare_ptr, score, m_Movie);

      // approximator termination criteria
      opti::CriterionCombine< FragmentComplete, double> criterion_combine;

      // terminate if max iterations reached
      opti::CriterionNumberIterations< FragmentComplete, double> maximum_number_iterations( ITERATIONS);
      criterion_combine.InsertCriteria( maximum_number_iterations);

      // end early if we meet our criteria, which in this case is a perfect alignment
      opti::CriterionResultThreshold< FragmentComplete, double> criterion_threshold( 0.0);
      criterion_combine.InsertCriteria( criterion_threshold);

      // terminate if we fail to improve after specified iterations
      opti::CriterionUnimproved< FragmentComplete, double> unimproved_threshold( MAX_UNIMPROVED);
      criterion_combine.InsertCriteria( unimproved_threshold);

      // Temperature control object, used in the metropolis criteria (constant temperature throughout, no exponential decay)
      util::ShPtr< mc::TemperatureInterface> sp_temperature( new mc::TemperatureAccepted( 0.5, 0.01, ITERATIONS, 1.0, 10));
      mc::Metropolis< double> metropolis( sp_temperature, true);

      // setting up the approximator for module 1
      mc::Approximator< FragmentComplete, double> approximator
      (
        score,
        mutator,
        metropolis,
        criterion_combine,
        mol_a
      );

      util::ShPtr< TransformationMutation> sp_mutation( approximator.GetMutator());
      sp_mutation->SetTracker( approximator.GetTracker());
      util::ShPtr< ScoreByField> sp_score( approximator.GetObjective());

      // run the approximator
      approximator.Approximate();

      // get the score of the best molecule
      return storage::Pair< FragmentComplete, double>
          (
              approximator.GetTracker().GetBest()->First(),
              approximator.GetTracker().GetBest()->Second()
          );
    }

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
    storage::Vector< storage::Pair< FragmentComplete, double> > ConformationComparisonPsiField::PseudoOperator
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B,
      const size_t &ITERATIONS,
      const size_t &MAX_UNIMPROVED,
      const bool &RECENTER_FIRST,
      const float &MAX_ATOM_DIST,
      const FragmentEnsemble &MOLECULE_A_CONFORMERS,
      const storage::Vector< size_t> &EXCLUSION_ATOMS_A,
      const storage::Vector< size_t> &EXCLUSION_ATOMS_B,
      const FragmentEnsemble &POCKETS
    ) const
    {
      // we need to get the keep_indices from exclusion_indices for the scaffold-based alignment
      storage::Vector< size_t> keep_indices_a, keep_indices_b;
      GetNonMaskedAtoms( MOLECULE_A, MOLECULE_B, m_ExclusionIndicesA, m_ExclusionIndicesB, keep_indices_a, keep_indices_b);
      m_KeepIndicesA = keep_indices_a;
      m_KeepIndicesB = keep_indices_b;

      // manage max atom distance
      float dmax_increment( ( m_AtomDistanceUpperLimit - m_AtomDistanceLowerLimit) / ( m_NumberRigidTrajectories - 1));
      storage::Vector< storage::Pair< FragmentComplete, double> > traj_tracker;
      storage::Vector< double> distances;

      //perform m_NumberRigidTrajectories number of independent rigid alignments
      for( size_t traj( 0); traj < m_NumberRigidTrajectories; ++traj)
      {
        // cache properties MOLECULE_A
        FragmentComplete mol_a( MOLECULE_A);
        mol_a.ShareCache( MOLECULE_A);

        //center mol_a on MOLECULE_B
        mol_a.Translate( MOLECULE_B.GetCenter() - mol_a.GetCenter());

        //randomly rotate off of MOLECULE_B
        if( m_InitialRotation)
        {
          math::RotationMatrix3D rot;
          rot.SetRand();
          linal::Vector3D mol_a_centered( mol_a.GetCenter());
          mol_a.Translate( -mol_a_centered);
          mol_a.Rotate( rot);
          mol_a.Translate( mol_a_centered);
        }

        //perform alignments
        storage::Pair< FragmentComplete, double> aligned;

        float max_atom_distance( m_AtomDistanceLowerLimit + ( dmax_increment * traj));
        distances.PushBack( max_atom_distance);
        aligned = FieldOptimizeOrientation( mol_a, MOLECULE_B, m_IterationNumber, m_LimitNumber, false, max_atom_distance, FragmentEnsemble(),
          m_ExclusionIndicesA, m_ExclusionIndicesB);

        //keep track of independent alignment trajectories
        storage::Pair< FragmentComplete, double> best_mol( aligned.First(), aligned.Second());
        traj_tracker.PushBack( best_mol);
      }

      if( m_AlignToScaffold)
      {
        // cache properties MOLECULE_A
        FragmentComplete mol_a( MOLECULE_A);
        mol_a.ShareCache( MOLECULE_A);
        storage::Vector< storage::Pair< ConformationGraphConverter::AtomComparisonTypeEnum, ConfigurationalBondTypeData::DataEnum>> ats_options;

        static storage::List< util::Stopwatch> stopwatches;
        if( stopwatches.IsEmpty())
        {
          for( size_t e_atom( 0); e_atom < ConformationGraphConverter::s_NumberAtomComparisonTypes; ++e_atom)
          {
            for( size_t e_bond( 1); e_bond < ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness; ++e_bond)
            {
              stopwatches.PushBack
              (
                util::Stopwatch
                (
                  ConformationGraphConverter::GetAtomComparisonType( ConformationGraphConverter::AtomComparisonType( e_atom))
                  + "_" + ConfigurationalBondTypeData::GetDataName( static_cast< ConfigurationalBondTypeData::Data>( e_bond))
                  + " Timing",
                  util::Time( 10, 0),
                  util::Message::e_Verbose,
                  true,
                  false
                )
              );
            }
          }
        }
        auto itr_stop( stopwatches.Begin());
        for( size_t e_atom( 0); e_atom < ConformationGraphConverter::s_NumberAtomComparisonTypes; ++e_atom)
        {
          for( size_t e_bond( 1); e_bond < ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness; ++e_bond, ++itr_stop)
          {
            itr_stop->Start();
            ArbitraryScaffoldAlignment
            (
              mol_a,
              MOLECULE_B,
              static_cast< ConformationGraphConverter::AtomComparisonType>( e_atom),
              static_cast< ConfigurationalBondTypeData::Data>( e_bond),
              m_KeepIndicesA,
              m_KeepIndicesB
            );

            float score_a
            (
              this->PropertyDistanceScorer
              (
                mol_a,
                MOLECULE_B,
                m_AtomDistanceUpperLimit,
                m_ExclusionIndicesA,
                m_ExclusionIndicesB
              ).First()
            );
            traj_tracker.PushBack();
            traj_tracker.LastElement().First() = mol_a;
            traj_tracker.LastElement().Second() = score_a;
            itr_stop->Stop();
          }
        }
      }

      // Formerly just evaluated the lower, middle, and upper bounds on the max atom distances;
      // now averages across all sampled distances
      storage::Vector< storage::Pair< double, size_t>> sorter( traj_tracker.GetSize());
      size_t tracker_index( 0);
      for
      (
          storage::Vector< storage::Pair< FragmentComplete, double >>::iterator itr( traj_tracker.Begin()), itr_end( traj_tracker.End());
          itr < itr_end;
          ++itr, ++tracker_index
      )
      {
        math::RunningAverage< double> sum( 0.0);
        for( size_t i( 0); i < distances.GetSize(); ++i)
        {
          sum += PropertyDistanceScorer( itr->First(), MOLECULE_B, distances( i), m_ExclusionIndicesA, m_ExclusionIndicesB).First();
        }
        double ave( sum.GetAverage());
        sorter( tracker_index).First() = ave;
        itr->First().StoreProperty( "PropertyFieldDistance", util::Format()( ave));
        sorter( tracker_index).Second() = tracker_index;
      }

      // Sort molecules by best
      sorter.Sort( std::less< storage::Pair< double, size_t>>());
      storage::Vector< size_t> best_order( traj_tracker.GetSize());
      for( size_t i( 0); i < traj_tracker.GetSize(); ++i)
      {
        best_order( i) = sorter( i).Second();
      }
      traj_tracker.Reorder( best_order);

      // prune similar molecules
      ConformationComparisonBySymmetryRmsd compare_symrmsd( false, 100); // consider real space difference
      storage::Vector< storage::Pair< FragmentComplete, double> > pruned_pair_scores;

      // add first molecule to new molecules collection
      pruned_pair_scores.PushBack( traj_tracker( 0));

      // check other alignment poses; if they deviate by m_PoseTolerance angstroms then they are different and save
      storage::Vector< size_t> score_indices;
      for( size_t i( 0); i < traj_tracker.GetSize(); ++i)
      {
        // loop over the saved poses
        bool same( false);
        for( size_t j( 0); j < pruned_pair_scores.GetSize(); ++j)
        {
          if( compare_symrmsd( traj_tracker( i).First(), pruned_pair_scores( j).First()) <= m_PoseTolerance)
          {
            // outer loop contains something same as already in new collection
            same = true;
            break;
          }
        }
        if( !same && sorter( i).First() < m_PoseScoreThreshold)
        {
          pruned_pair_scores.PushBack
          (
            storage::Pair< FragmentComplete, double>(
              std::make_pair
              (
                traj_tracker( i).First(),
                sorter( i).First()
              )
            )
          );
          score_indices.PushBack( i);
        }
      }

      // output the best alignment pairs and their rmsdx values
      return pruned_pair_scores;
    }

    ///////////////////
    // Helper Functions
    ///////////////////

    //! @brief whether to perform substructure-based alignments as part of sampling
    void ConformationComparisonPsiField::SetAlignToScaffold( const bool ALIGN_TO_SCAFFOLD)
    {
      m_AlignToScaffold = ALIGN_TO_SCAFFOLD;
    }

    //! @brief whether to perform an initial random rotation of molecule A
    void ConformationComparisonPsiField::SetInitialRotation( const bool INITIAL_ROTATION)
    {
      m_InitialRotation = INITIAL_ROTATION;
    }

    //! @brief set the number of iterations during the initial alignment
    void ConformationComparisonPsiField::SetNumberIterations( const size_t N_ITERATIONS)
    {
      m_IterationNumber = N_ITERATIONS;
    }

    //! @brief set the number of max unimproved iterations allowed in MCM
    void ConformationComparisonPsiField::SetNumberMaxUnimproved( const size_t N_MAX_UNIMPROVED)
    {
      m_LimitNumber = N_MAX_UNIMPROVED;
    }

    //! @brief return the number of rigid alignment trajectories to perform
    void ConformationComparisonPsiField::SetNumberRigidTraj( const size_t N_RIGID_TRAJ)
    {
      m_NumberRigidTrajectories = N_RIGID_TRAJ;
    }

    //! @brief return the upper limit on the max atom distance range
    void ConformationComparisonPsiField::SetAtomDistUpper( const float ATOM_DIST_UPPER)
    {
      m_AtomDistanceUpperLimit = ATOM_DIST_UPPER;
    }

    //! @brief return the lower limit on the max atom distance range
    void ConformationComparisonPsiField::SetAtomDistLower( const float ATOM_DIST_LOWER)
    {
      m_AtomDistanceLowerLimit = ATOM_DIST_LOWER;
    }

    //! @brief returns the prefix for the aligned molecule A output SDF file
    void ConformationComparisonPsiField::SetOutputAPrefix( const std::string &OUTPUT_A)
    {
      m_OutputA = OUTPUT_A;
    }

    //! @brief returns the prefix for the aligned molecule B output SDF file
    void ConformationComparisonPsiField::SetOutputBPrefix( const std::string &OUTPUT_B)
    {
      m_OutputB = OUTPUT_B;
    }

    //! @brief return the number of maximum unique output alignments
    void ConformationComparisonPsiField::SetNumberMaxOutputs( const size_t N_MAX_OUTPUTS)
    {
      m_NumberOutputSDF = N_MAX_OUTPUTS;
    }

    //! @brief return the tolerance in Angstroms differentiating unique poses
    void ConformationComparisonPsiField::SetPoseTolerance( const double POSE_TOLERANCE)
    {
      m_PoseTolerance = POSE_TOLERANCE;
    }

    //! @brief return the threshold in score units accepting a pose
    void ConformationComparisonPsiField::SetPoseScoreThreshold( const double SCORE_THRESHOLD)
    {
      m_PoseScoreThreshold = SCORE_THRESHOLD;
    }

    //! @brief set the relative probability of performing a Flip move during alignment
    void ConformationComparisonPsiField::SetFlipProb( const double FLIP_PROB)
    {
      m_FlipProb = FLIP_PROB;
    }

    //! @brief set the relative probability of performing a BigRotation move during alignment
    void ConformationComparisonPsiField::SetBigRotProb( const double BIG_ROT_PROB)
    {
      m_BigRotProb = BIG_ROT_PROB;
    }

    //! @brief set the relative probability of performing a BondSwap move during alignment
    void ConformationComparisonPsiField::SetBondSwapProb( const double BOND_SWAP_PROB)
    {
      m_BondSwapProb = BOND_SWAP_PROB;
    }

    //! @brief set the relative probability of performing a BondSlign move during alignment
    void ConformationComparisonPsiField::SetBondAlignProb( const double BOND_ALIGN_PROB)
    {
      m_BondAlignProb = BOND_ALIGN_PROB;
    }

    //! @brief set the relative probability of performing a BondAlign3 move during alignment
    void ConformationComparisonPsiField::SetBondAlign3Prob( const double BOND_ALIGN3_PROB)
    {
      m_BondAlign3Prob = BOND_ALIGN3_PROB;
    }

    //! @brief set the relative probability of performing a BondAlignInf move during alignment
    void ConformationComparisonPsiField::SetBondAlignInfProb( const double BOND_ALIGN_INF_PROB)
    {
      m_BondAlignInfProb = BOND_ALIGN_INF_PROB;
    }

    //! @brief set the relative probability of performing a SmallRotation move during alignment
    void ConformationComparisonPsiField::SetSmallRotProb( const double SMALL_ROT_PROB)
    {
      m_SmallRotProb = SMALL_ROT_PROB;
    }

    //! @brief set the relative probability of performing a SmallTranslation move during alignment
    void ConformationComparisonPsiField::SetSmallTransProb( const double SMALL_TRANS_PROB)
    {
      m_SmallTranslateProb = SMALL_TRANS_PROB;
    }

    //! @brief set the relative probability of performing a ConformerSwap move during alignment
    void ConformationComparisonPsiField::SetConfSwapProb( const double CONF_SWAP_PROB)
    {
      m_ConformerSwapProb = CONF_SWAP_PROB;
    }

    //! @brief returns the exclusion indices for molecule A
    void ConformationComparisonPsiField::SetExclusionIndicesA( const storage::Vector< size_t> &EXCLUSION_INDICES_A)
    {
      m_ExclusionIndicesA = EXCLUSION_INDICES_A;
    }

    //! @brief returns the indices to be kept for molecule A
    void ConformationComparisonPsiField::SetKeepIndicesA( const storage::Vector< size_t> &KEEP_INDICES_A)
    {
      m_KeepIndicesA = KEEP_INDICES_A;
    }

    //! @brief returns the exclusion indices for molecule B
    void ConformationComparisonPsiField::SetExclusionIndicesB( const storage::Vector< size_t> &EXCLUSION_INDICES_B)
    {
      m_ExclusionIndicesB = EXCLUSION_INDICES_B;
    }

    //! @brief returns the indices to be kept for molecule B
    void ConformationComparisonPsiField::SetKeepIndicesB( const storage::Vector< size_t> &KEEP_INDICES_B)
    {
      m_KeepIndicesB = KEEP_INDICES_B;
    }

    //! @brief returns the name of the alignment movie output file
    void ConformationComparisonPsiField::SetMovieName( const std::string &MOVIE_NAME)
    {
      m_Movie = MOVIE_NAME;
    }

    bool ConformationComparisonPsiField::ArbitraryScaffoldAlignment
    (
      FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE,
      const storage::Vector< size_t> &TARGET_MOL_INDICES,
      const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES
    ) const
    {
      FragmentAlignToScaffold align_to_scaffold( ATOM_TYPE, BOND_TYPE, size_t( 3));
      return align_to_scaffold.AlignToScaffold( MOLECULE_A, MOLECULE_B, TARGET_MOL_INDICES, SCAFFOLD_MOL_INDICES);
    }

    //! @brief obtain atoms that can be perturbed during alignment
    //! @param MOLECULE_A - first input molecule
    //! @param MOLECULE_B - second input molecule
    //! @param EXCLUSION_INDICES_A - indices of atoms in molecule A to hide during alignment
    //! @param EXCLUSION_INDICES_B - indices of atoms in molecule B to hide during alignment
    //! @return void
    void ConformationComparisonPsiField::GetNonMaskedAtoms
    (
      const FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const storage::Vector< size_t> &EXCLUSION_INDICES_A,
      const storage::Vector< size_t> &EXCLUSION_INDICES_B,
      storage::Vector< size_t> &KEEP_INDICES_A,
      storage::Vector< size_t> &KEEP_INDICES_B
    )
    {
      for( size_t a_i( 0), a_sz( MOLECULE_A.GetSize()); a_i < a_sz; ++a_i)
      {
        // do not include the atoms we said we wanted to exclude
        bool exclude( false);
        for( size_t ea_i( 0), ea_sz( EXCLUSION_INDICES_A.GetSize()); ea_i < ea_sz; ++ea_i)
        {
          if( a_i == EXCLUSION_INDICES_A( ea_i))
          {
            exclude = true;
            break;
          }
        }
        // if we did not find it in the exclusion indices, add it
        if( !exclude)
        {
          KEEP_INDICES_A.PushBack( a_i);
        }
      }
      for( size_t b_i( 0), b_sz( MOLECULE_B.GetSize()); b_i < b_sz; ++b_i)
      {
        // do not include the atoms we said we wanted to exclude
        bool exclude( false);
        for( size_t eb_i( 0), eb_sz( EXCLUSION_INDICES_B.GetSize()); eb_i < eb_sz; ++eb_i)
        {
          if( b_i == EXCLUSION_INDICES_B( eb_i))
          {
            exclude = true;
            break;
          }
        }
        // if we did not find it in the exclusion indices, add it
        if( !exclude)
        {
          KEEP_INDICES_B.PushBack( b_i);
        }
      }
    }

    bool ConformationComparisonPsiField::BondAlignInf
    (
      FragmentComplete &MOLECULE_A,
      const FragmentComplete &MOLECULE_B,
      const storage::Vector< size_t> &KEEP_INDICES_A,
      const storage::Vector< size_t> &KEEP_INDICES_B
    )
    {
//      storage::Vector< size_t> atom_indices_aligned_a, atom_indices_aligned_b;
//      ConformationComparisonPsiField::GetAlignedAtoms( MOLECULE_A, MOLECULE_B, atom_indices_aligned_a, atom_indices_aligned_b);
      auto atom_indices_aligned( ConformationComparisonPsiField::GetAlignedAtoms( MOLECULE_A, MOLECULE_B, KEEP_INDICES_A, KEEP_INDICES_B));

      // get the aligned atoms
      const size_t n_atoms_a( MOLECULE_A.GetNumberAtoms());
      const size_t n_atoms_b( MOLECULE_B.GetNumberAtoms());

      iterate::Generic< const AtomConformationalInterface> r_itr_a( MOLECULE_A.GetAtomsIterator());
      storage::Vector< storage::Pair< double, size_t> > mol_a_best_match( n_atoms_a, storage::Pair< double, size_t>( util::GetUndefined< double>(), n_atoms_b));
      storage::Vector< storage::Pair< double, size_t> > mol_b_best_match( n_atoms_b, storage::Pair< double, size_t>( util::GetUndefined< double>(), n_atoms_a));

//      // for each atom in MOLECULE_A find the nearest atom in MOLECULE_B
//      for( size_t mol_a_pos( 0); mol_a_pos < n_atoms_a; ++r_itr_a, ++mol_a_pos)
//      {
//        // make sure this is an allowed molecule a index
//        bool keep( false);
//        for( size_t keep_a_i( 0), keep_a_sz( KEEP_INDICES_A.GetSize()); keep_a_i < keep_a_sz; ++keep_a_i)
//        {
//          if( mol_a_pos == KEEP_INDICES_A( keep_a_i))
//          {
//            keep = true;
//            break;
//          }
//        }
//        if(!keep)
//        {
//          continue;
//        }
//
//        size_t mol_b_pos( 0);
//        const linal::Vector3D &a_pos( r_itr_a->GetPosition());
//        for
//        (
//          iterate::Generic< const AtomConformationalInterface> r_itr_b( MOLECULE_B.GetAtomsIterator());
//          mol_b_pos < n_atoms_b;
//          ++r_itr_b, ++mol_b_pos
//        )
//        {
//          // make sure this is an allowed molecule b index
//          bool keep( false);
//          for( size_t keep_b_i( 0), keep_b_sz( KEEP_INDICES_B.GetSize()); keep_b_i < keep_b_sz; ++keep_b_i)
//          {
//            if( mol_b_pos == KEEP_INDICES_B( keep_b_i))
//            {
//              keep = true;
//              break;
//            }
//          }
//          if(!keep)
//          {
//            continue;
//          }
//
//          const float distance( linal::SquareDistance( r_itr_b->GetPosition(), a_pos));
//          if( mol_a_best_match( mol_a_pos).Second() == n_atoms_b || distance < mol_a_best_match( mol_a_pos).First())
//          {
//            mol_a_best_match( mol_a_pos).Second() = mol_b_pos;
//            mol_a_best_match( mol_a_pos).First() = distance;
//          }
//          if( mol_b_best_match( mol_b_pos).Second() == n_atoms_a || distance < mol_b_best_match( mol_b_pos).First())
//          {
//            mol_b_best_match( mol_b_pos).Second() = mol_a_pos;
//            mol_b_best_match( mol_b_pos).First() = distance;
//          }
//        }
//      }
//      atom_indices_aligned_a.Reset();
//      atom_indices_aligned_b.Reset();
//
//      // save the indices of the mutually closest atom pairs from MOLECULE_A and MOLECULE_B
//      for( size_t mol_a_pos( 0); mol_a_pos < n_atoms_a; ++mol_a_pos)
//      {
//        if( mol_a_best_match( mol_a_pos).Second() < n_atoms_b && mol_b_best_match( mol_a_best_match( mol_a_pos).Second()).Second() == mol_a_pos)
//        {
//          atom_indices_aligned_a.PushBack( mol_a_pos);
//          atom_indices_aligned_b.PushBack( mol_a_best_match( mol_a_pos).Second());
//        }
//      }

//      if( atom_indices_aligned_a.GetSize() < 3)
      if( atom_indices_aligned.First().GetSize() < 3)
      {
        return false;
      }

      util::SiPtrVector< const linal::Vector3D> matched_atom_pos_a, matched_atom_pos_b;
//      for( size_t aligned_pos( 0), n_aligned( atom_indices_aligned_a.GetSize()); aligned_pos < n_aligned; ++aligned_pos)
//      {
//        matched_atom_pos_a.PushBack( MOLECULE_A.GetAtomVector()( atom_indices_aligned_a( aligned_pos)).GetPosition());
//        matched_atom_pos_b.PushBack( MOLECULE_B.GetAtomVector()( atom_indices_aligned_b( aligned_pos)).GetPosition());
//      }
      for( size_t aligned_pos( 0), n_aligned( atom_indices_aligned.First().GetSize()); aligned_pos < n_aligned; ++aligned_pos)
      {
        matched_atom_pos_a.PushBack( MOLECULE_A.GetAtomVector()( atom_indices_aligned.First()( aligned_pos)).GetPosition());
        matched_atom_pos_b.PushBack( MOLECULE_B.GetAtomVector()( atom_indices_aligned.Second()( aligned_pos)).GetPosition());
      }

      quality::RMSD calculator;
      math::TransformationMatrix3D transformed_coords( calculator.CalculateSuperimposition( matched_atom_pos_a, matched_atom_pos_b));
      MOLECULE_A.Transform( transformed_coords);
      return true;
    }

    // scoring function definition for ScoreByField class to create Score for approximator
    double ConformationComparisonPsiField::GetBaseField
    (
      const FragmentComplete &MOL_A,
      const FragmentComplete &MOL_B,
      const float &ATOM_DISTANCE,
      const storage::Vector< size_t> &EXCLUDE_ATOMS_A,
      const storage::Vector< size_t> &EXCLUDE_ATOMS_B,
      const FragmentEnsemble &POCKETS
    ) const
    {
      if( !POCKETS.IsEmpty() && util::IsDefined( ATOM_DISTANCE))
      {
        //Compute molecule alignment score
        double molecule_alignment_score( ConformationComparisonPropertyFieldCorrelation::PropertyDistanceScorer( MOL_A, MOL_B, ATOM_DISTANCE, EXCLUDE_ATOMS_A, EXCLUDE_ATOMS_B).First());

        //Build voxel grids for the pocket and the molecule
        VoxelGridAtom voxel_grid_pocket( 4.0), voxel_grid_mol( 4.0);
        voxel_grid_mol.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( MOL_A.GetAtomsIterator(), MOL_A.GetAtomsIterator().End()));
        voxel_grid_pocket.SetObjects( util::SiPtrVector< const AtomConformationalInterface>( POCKETS.Begin()->GetAtomsIterator(), POCKETS.Begin()->GetAtomsIterator().End()));
        auto neighbors( voxel_grid_mol.GetNeighborsIn( voxel_grid_pocket, 4.0));

        //Check for atom overlap between all neighbor pairs
        for( auto itr_neighbors( neighbors.Begin()), itr_neighbors_end( neighbors.End()); itr_neighbors != itr_neighbors_end; ++itr_neighbors)
        {
          //Get distance between atoms
          const double distance( itr_neighbors->Third());

          //Compute csd_vdw radii
          double mol_csd_vdw_radius( itr_neighbors->First()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));
          double pocket_csd_vdw_radius( itr_neighbors->Second()->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));

          //Compute collision score
          double collision_score( mol_csd_vdw_radius + pocket_csd_vdw_radius - distance);

          //Check atom pair complementarity
//          double mol_size_norm( -molecule_alignment_score / ( ( MOL_A.GetNumberAtoms() - MOL_A.GetNumberHydrogens() + MOL_B.GetNumberAtoms() - MOL_B.GetNumberHydrogens()) / 2.0));
          // double best_complementarity_score( mol_size_norm * 0.30), best_pairpotential_score( mol_size_norm * 0.20);
//          molecule_alignment_score += std::max(best_complementarity_score,this->CheckAtomTypeComplementarity(itr_neighbors->First(),itr_neighbors->Second(),distance));
//          molecule_alignment_score += std::max(best_pairpotential_score,this->CheckPairPotential(itr_neighbors->First(),itr_neighbors->Second(),-collision_score));

          //Penalize alignments with collisions
          if( collision_score < 1.40)
          {
            molecule_alignment_score += 0.0;
          }
          else
          {
            molecule_alignment_score += 5.0;
          }
        }

        static LigandPocketFitScore mol_pocket_complementarity_scorer;
        molecule_alignment_score += mol_pocket_complementarity_scorer.PropertyCorrelationScore( MOL_A, *POCKETS.Begin());
        return molecule_alignment_score;
      }
      else if( util::IsDefined( ATOM_DISTANCE))
      {
        return ConformationComparisonPropertyFieldCorrelation::PropertyDistanceScorer( MOL_A, MOL_B, ATOM_DISTANCE, EXCLUDE_ATOMS_A, EXCLUDE_ATOMS_B).First();
      }
      else
      {
        return ConformationComparisonPropertyFieldCorrelation::operator()( MOL_A, MOL_B);
      }
    }

    double ConformationComparisonPsiField::CheckAtomTypeComplementarity
    (
      const util::SiPtr<const AtomConformationalInterface> &MOL_ATOM,
      const util::SiPtr<const AtomConformationalInterface> &POCKET_ATOM,
      const double &DISTANCE
    ) const
    {
      // Determine if atom in ligand is hydrogen bond donor or acceptor
      double bonus(0.0);
      bool h_donor(false), h_acceptor(false);
      if( MOL_ATOM->GetElementType()->GetAtomicNumber() == size_t( 7) || MOL_ATOM->GetElementType()->GetAtomicNumber() == size_t( 8))
      {
        if( MOL_ATOM->GetNumberCovalentlyBoundHydrogens() > size_t( 0))
        {
          h_donor = true;
        }
        else if( MOL_ATOM->GetNumberElectronsInValenceBonds() < 2 * MOL_ATOM->GetNumberValenceBonds())
        {
          h_donor = true;
        }
        else
        {
          h_acceptor = true;
        }
      }

      // Determine if atom in binding pocket is hydrogen bond donor or acceptor
      if( POCKET_ATOM->GetElementType()->GetAtomicNumber() == size_t( 7) || POCKET_ATOM->GetElementType()->GetAtomicNumber() == size_t( 8))
      {
        if( POCKET_ATOM->GetNumberCovalentlyBoundHydrogens() > size_t( 0))
        {
          h_donor = true;
        }
        else if( POCKET_ATOM->GetNumberElectronsInValenceBonds() < 2 * POCKET_ATOM->GetNumberValenceBonds())
        {
          h_donor = true;
        }
        else
        {
          h_acceptor = true;
        }
      }

      // the maximum bonus per complementary atom pair is 0.75
      double polar_max_bonus( 0.75); // nonpolar_max_bonus( 0.10);
      if( h_acceptor && h_donor)
      {
        if( DISTANCE > double( 3.00) && DISTANCE < double( 4.00))
        {
          bonus += polar_max_bonus;
        }
        else if( DISTANCE > 4.00)
        {
          bonus += 0.0;
        }
        else if( DISTANCE > 2.00)
        {
          bonus += polar_max_bonus * ( 3.00 - DISTANCE);
        }
        else
        {
          bonus += 0.0;
        }
      }

      return -bonus;
    }

    double ConformationComparisonPsiField::CheckPairPotential
    (
      const util::SiPtr<const AtomConformationalInterface> &MOL_ATOM,
      const util::SiPtr<const AtomConformationalInterface> &POCKET_ATOM,
      const double &DISTANCE
    ) const
    {
      auto itr( m_PairPotentials.Find( storage::Pair< AtomType, AtomType>( MOL_ATOM->GetAtomType(), POCKET_ATOM->GetAtomType())));
      if( itr == m_PairPotentials.End())
      {
        return 0;
      }

      return itr->second( DISTANCE);
    }

    // ScoreByField class

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    ConformationComparisonPsiField::ScoreByField::ScoreByField
    (
      const FragmentComplete &MOLECULE_B,
      const util::SiPtr< const ConformationComparisonPsiField> &PARENT,
      const float &MAX_ATOM_DISTANCE,
      const storage::Vector< size_t> &EXCLUSION_ATOMS_A,
      const storage::Vector< size_t> &EXCLUSION_ATOMS_B,
      const FragmentEnsemble &POCKETS
    ) :
      m_MoleculeB( MOLECULE_B),
      m_ParentComparer( PARENT),
      m_MaxAtomDistance( MAX_ATOM_DISTANCE),
      m_ExclusionIndicesA( EXCLUSION_ATOMS_A),
      m_ExclusionIndicesB( EXCLUSION_ATOMS_B),
      m_Pockets( POCKETS)
    {
    }

    //! virtual copy constructor
    ConformationComparisonPsiField::ScoreByField *ConformationComparisonPsiField::ScoreByField::Clone() const
    {
      return new ScoreByField( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonPsiField::ScoreByField::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    // input
    std::istream &ConformationComparisonPsiField::ScoreByField::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    // output
    std::ostream &ConformationComparisonPsiField::ScoreByField::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief score two molecules by property rmsd
    //! @return a score for the approximator

    double ConformationComparisonPsiField::ScoreByField::operator()( const FragmentComplete &MOLECULE_A) const
    {
      return m_ParentComparer->GetBaseField( MOLECULE_A, m_MoleculeB, m_MaxAtomDistance, m_ExclusionIndicesA, m_ExclusionIndicesB, m_Pockets);
    }

    //TransformationMutation class

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

     //! @default constructor
     ConformationComparisonPsiField::TransformationMutation::TransformationMutation
     (
       const FragmentEnsemble &MOLECULE_A_CONFORMERS,
       const FragmentComplete &MOL_B,
       const util::SiPtr< const ConformationComparisonPsiField> &PARENT,
       const util::SiPtr< const ConformationComparisonPsiField::ScoreByField> &SCORE_KNOWER,
       const std::string &MOVIE_FILENAME
     ) :
         m_MolEnsA( MOLECULE_A_CONFORMERS.GetMolecules().IsEmpty() ?
             util::SiPtr<const FragmentEnsemble>() :
             util::SiPtr<const FragmentEnsemble>( MOLECULE_A_CONFORMERS)),
         m_MolBPtr( MOL_B),
         m_Parent( PARENT),
         m_ScoreKnower( SCORE_KNOWER),
         m_MovieFilename( MOVIE_FILENAME)
     {
     }

     //! virtual copy constructor
     ConformationComparisonPsiField::TransformationMutation *ConformationComparisonPsiField::TransformationMutation::Clone() const
     {
       return new TransformationMutation( *this);
     }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonPsiField::TransformationMutation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    // input
    std::istream &ConformationComparisonPsiField::TransformationMutation::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    // output
    std::ostream &ConformationComparisonPsiField::TransformationMutation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief makes a random rotation or translation to molecule A
    //! @return the mutation result for the approximator
    math::MutateResult< FragmentComplete> ConformationComparisonPsiField::TransformationMutation::operator()( const FragmentComplete &MOLECULE_A) const
    {
      // no point in this if there are no keep indices
      BCL_Assert
      (
        m_Parent->m_KeepIndicesA.GetSize(),
        "There are no atoms in molecule A for alignment!"
      );
      BCL_Assert
      (
        m_Parent->m_KeepIndicesB.GetSize(),
        "There are no atoms in molecule B for alignment!"
      );

      FragmentComplete molecule_a( MOLECULE_A);
      static const storage::Vector< float> flip_angles( storage::Vector< float>::Create( 90.0, 120.0, 180.0));
      std::string chosen_scheme( m_Parent->ChooseAlignmentMove( m_MolEnsA.IsDefined()));

      if( chosen_scheme == "Flip")
      {
        linal::Vector3D axis;
        axis.SetRandomTranslation( 1.0);
        size_t rand_flip_angle( random::GetGlobalRandom().Random< size_t>( 0, 2));
        math::RotationMatrix3D flip( axis, flip_angles( rand_flip_angle));
        linal::Vector3D molecule_a_centered( molecule_a.GetCenter());
        molecule_a.Translate( -molecule_a_centered);
        molecule_a.Rotate( flip);
        molecule_a.Translate( molecule_a_centered);
        m_Scheme = "Flip ";
      }

      // rotate molecule randomly with magnitude between 0 and 180 degrees
      else if( chosen_scheme == "BigRot")
      {
        math::RotationMatrix3D rot;
        rot.SetRand();
        linal::Vector3D molecule_a_centered( molecule_a.GetCenter());
        molecule_a.Translate( -molecule_a_centered);
        molecule_a.Rotate( rot);
        molecule_a.Translate( molecule_a_centered);
        m_Scheme = "BigRot ";
      }
      // bondswap
      else if( chosen_scheme == "BondSwap")
      {
        const storage::Vector< sdf::BondInfo> ori_bonds( molecule_a.GetBondInfo());
        if( ori_bonds.GetSize() > size_t( 1))
        {
          size_t tries( 0);
          for( ; tries < 100; ++tries)
          {
            const size_t bond_id_a( random::GetGlobalRandom().Random( ori_bonds.GetSize() - 1));
            const size_t bond_id_b( random::GetGlobalRandom().Random( ori_bonds.GetSize() - 2));
            const sdf::BondInfo &bond_a( ori_bonds( bond_id_a));
            const sdf::BondInfo &bond_b( ori_bonds( bond_id_b >= bond_id_a ? bond_id_b + 1 : bond_id_b));

            // make sure that the atoms in the chosen bonds are in keepindices
            if
            (
                m_Parent->m_KeepIndicesA.Find( bond_a.GetAtomIndexLow()) < m_Parent->m_KeepIndicesA.GetSize() &&
                m_Parent->m_KeepIndicesA.Find( bond_a.GetAtomIndexHigh()) < m_Parent->m_KeepIndicesA.GetSize() &&
                m_Parent->m_KeepIndicesB.Find( bond_b.GetAtomIndexLow()) < m_Parent->m_KeepIndicesB.GetSize() &&
                m_Parent->m_KeepIndicesB.Find( bond_b.GetAtomIndexHigh()) < m_Parent->m_KeepIndicesB.GetSize()
            )
            {
              coord::LineSegment3D line_a(molecule_a.GetAtomVector()(bond_a.GetAtomIndexLow()).GetPosition(), molecule_a.GetAtomVector()(bond_a.GetAtomIndexHigh()).GetPosition());
              coord::LineSegment3D line_b(molecule_a.GetAtomVector()(bond_b.GetAtomIndexLow()).GetPosition(), molecule_a.GetAtomVector()(bond_b.GetAtomIndexHigh()).GetPosition());

              if( random::GetGlobalRandom().Boolean())
              {
                line_b = coord::LineSegment3D( molecule_a.GetAtomVector()( bond_b.GetAtomIndexHigh()).GetPosition(), molecule_a.GetAtomVector()( bond_b.GetAtomIndexLow()).GetPosition());
              }
              math::TransformationMatrix3D bondswap_matrix( line_a, line_b);
              molecule_a.Transform( bondswap_matrix);

              m_Scheme = "BondSwap ";
            }
          }
        }
        else
        {
          return this->operator ()( MOLECULE_A);
        }
      }
      // bondalign
      else if( chosen_scheme == "BondAlign")
      {
        float minimum_distance( math::GetHighestBoundedValue< float>());
        size_t mol_b_atom_index_tracker( util::GetUndefinedSize_t());

        // get a random allowed atom from molecule A
        auto mol_a_atom_indices( m_Parent->m_KeepIndicesA);
        mol_a_atom_indices.Shuffle();
        size_t atom_a_rand_index( mol_a_atom_indices( 0));
        const AtomComplete &atom_a( molecule_a.GetAtomVector()( atom_a_rand_index));

        // get the closest atom in molecule B to the chosen atom in molecule A that is allowed by keepindices
        for( size_t mol_b_atom_index( 0); mol_b_atom_index < m_MolBPtr->GetSize(); ++mol_b_atom_index)
        {
          // make sure this is an allowed molecule b index
          bool keep( false);
          for( size_t keep_b_i( 0), keep_b_sz( m_Parent->m_KeepIndicesB.GetSize()); keep_b_i < keep_b_sz; ++keep_b_i)
          {
            if( mol_b_atom_index == m_Parent->m_KeepIndicesB( keep_b_i))
            {
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }

          //compute the distance between atom from molecule_a and all atoms in molecule_b
          float current_distance = linal::Distance( atom_a.GetPosition(), m_MolBPtr->GetAtomVector()( mol_b_atom_index).GetPosition());

          //identify closest atom in molecule b and keep its index
          if( current_distance < minimum_distance)
          {
            minimum_distance = current_distance;
            mol_b_atom_index_tracker = mol_b_atom_index;
          }
        }

        // if there is no such atom in molecule B then return input
        if( !util::IsDefined( mol_b_atom_index_tracker))
        {
          return this->operator ()( MOLECULE_A);
        }

        // get atom b
        const AtomComplete &atom_b( m_MolBPtr->GetAtomVector()( mol_b_atom_index_tracker));

        // get the bond indices for each of the chosen atoms and randomize them
        storage::Vector< size_t> atom_a_bonds( atom_a.GetBonds().GetSize()), atom_b_bonds( atom_b.GetBonds().GetSize());
        for( size_t i( 0), sz( atom_a_bonds.GetSize()); i < sz; ++i)
        {
          atom_a_bonds( i) = i;
        }
        for( size_t i( 0), sz( atom_b_bonds.GetSize()); i < sz; ++i)
        {
          atom_b_bonds( i) = i;
        }
        atom_a_bonds.Shuffle();
        atom_b_bonds.Shuffle();

        // pick a bonded atom from atom a that is contained within keepindices
        size_t bonded_a_index( util::GetUndefinedSize_t());
        for
        (
            auto atom_a_itr( atom_a_bonds.Begin()), atom_a_itr_end( atom_a_bonds.End());
            atom_a_itr != atom_a_itr_end;
            ++atom_a_itr
        )
        {
          // make sure this is an allowed molecule a index
          bool keep( false);
          for( size_t keep_a_i( 0), keep_a_sz( m_Parent->m_KeepIndicesA.GetSize()); keep_a_i < keep_a_sz; ++keep_a_i)
          {
            size_t i( molecule_a.GetAtomVector().GetAtomIndex( atom_a.GetBonds()( *atom_a_itr).GetTargetAtom()));
            if( i == m_Parent->m_KeepIndicesA( keep_a_i))
            {
              bonded_a_index = i;
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }
        }
        if( !util::IsDefined( bonded_a_index))
        {
          return this->operator ()( MOLECULE_A);
        }

        // pick a bonded atom from atom b that is contained within keepindices
        size_t bonded_b_index( util::GetUndefinedSize_t());
        for
        (
            auto atom_b_itr( atom_b_bonds.Begin()), atom_b_itr_end( atom_b_bonds.End());
            atom_b_itr != atom_b_itr_end;
            ++atom_b_itr
        )
        {
          // make sure this is an allowed molecule a index
          bool keep( false);
          for( size_t keep_b_i( 0), keep_b_sz( m_Parent->m_KeepIndicesB.GetSize()); keep_b_i < keep_b_sz; ++keep_b_i)
          {
            size_t i( m_MolBPtr->GetAtomVector().GetAtomIndex( atom_b.GetBonds()( *atom_b_itr).GetTargetAtom()));
            if( i == m_Parent->m_KeepIndicesB( keep_b_i))
            {
              bonded_b_index = i;
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }
        }
        if( !util::IsDefined( bonded_b_index))
        {
          return this->operator ()( MOLECULE_A);
        }

        // get the bonded atoms
        const AtomConformationalInterface &bonded_a( molecule_a.GetAtomVector()( bonded_a_index));
        const AtomConformationalInterface &bonded_b( m_MolBPtr->GetAtomVector()( bonded_b_index));

        // create a line segment for each bond
        coord::LineSegment3D line_a (atom_a.GetPosition(), bonded_a.GetPosition());
        coord::LineSegment3D line_b (atom_b.GetPosition(), bonded_b.GetPosition());

        // align the molecule
        math::TransformationMatrix3D bondalign_matrix( line_b, line_a);
        molecule_a.Transform( bondalign_matrix);
        m_Scheme = "BondAlign ";
      }
      // bondalign3
      else if( chosen_scheme == "BondAlign3")
      {
        // useless if user only has diatomic molecules or less than 3 unmasked atoms
        if( molecule_a.GetNumberAtoms() < 3)
        {
          return this->operator ()( MOLECULE_A);
        }

        // initialize data
        bool correct_bond_number = false;
        float minimum_distance( math::GetHighestBoundedValue< float>());
        size_t mol_b_atom_index_tracker( util::GetUndefined< size_t>());
        size_t atom_a_rand_index( util::GetUndefined< size_t>());

        //if the chosen atom does not make at least 2 bonds then choose another atom
        size_t tries( 0);
        while( correct_bond_number == false && tries < 100)
        {
          // pick a random atom from molecule a that is allowed by keepindices
          auto mol_a_atom_indices( m_Parent->m_KeepIndicesA);
          mol_a_atom_indices.Shuffle();
          atom_a_rand_index = mol_a_atom_indices( 0);
          const AtomComplete &atom_a_1( molecule_a.GetAtomVector()( atom_a_rand_index));
          if( atom_a_1.GetAtomType()->GetNumberBonds() < 2)
          {
            // BCL_MessageStd("Atom makes < 2 bonds, choosing another atom");
          }
          else
          {
            correct_bond_number = true;
          }
          ++tries;
        }
        // require that this is defined or else quit
        if( !util::IsDefined( atom_a_rand_index))
        {
          return this->operator ()( MOLECULE_A);
        }

        // now we have atom a
        const AtomComplete &atom_a_1( molecule_a.GetAtomVector()( atom_a_rand_index));

        // find the closest atom in molecule b from atom a which makes at least 2 bonds
        size_t chosen_bond_a1_index, chosen_bond_a2_index, chosen_bond_b1_index, chosen_bond_b2_index;
        for( size_t mol_b_atom_index( 0); mol_b_atom_index < m_MolBPtr->GetSize(); ++mol_b_atom_index)
        {
          // make sure this is an allowed molecule b index
          bool keep( false);
          for( size_t keep_b_i( 0), keep_b_sz( m_Parent->m_KeepIndicesB.GetSize()); keep_b_i < keep_b_sz; ++keep_b_i)
          {
            if( mol_b_atom_index == m_Parent->m_KeepIndicesB( keep_b_i))
            {
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }

          //compute the distance between atoms from molecule_a and all atoms in molecule_b
          float current_distance = linal::Distance( atom_a_1.GetPosition(), m_MolBPtr->GetAtomVector()( mol_b_atom_index).GetPosition());

          //identify closest atom in molecule b which also makes at least 2 bonds
          if( current_distance < minimum_distance)
          {
            const AtomComplete &atom_b_1( m_MolBPtr->GetAtomVector()( mol_b_atom_index));
            if( atom_b_1.GetAtomType()->GetNumberBonds() >= 2)
            {
              minimum_distance = current_distance;
              mol_b_atom_index_tracker = mol_b_atom_index;
            }
            else
            {
              continue;
            }
          }
        }
        // require that this is defined or else quit
        if( !util::IsDefined( mol_b_atom_index_tracker))
        {
          return this->operator ()( MOLECULE_A);
        }

        //Nearest neighbor in molecule_b to original random atom from molecule_a
        const AtomComplete &atom_b_1( m_MolBPtr->GetAtomVector()( mol_b_atom_index_tracker));

        // get the bond indices for each of the chosen atoms and randomize them
        storage::Vector< size_t> atom_a_bonds( atom_a_1.GetBonds().GetSize()), atom_b_bonds( atom_b_1.GetBonds().GetSize());
        for( size_t i( 0), sz( atom_a_bonds.GetSize()); i < sz; ++i)
        {
          atom_a_bonds( i) = i;
        }
        for( size_t i( 0), sz( atom_b_bonds.GetSize()); i < sz; ++i)
        {
          atom_b_bonds( i) = i;
        }
        atom_a_bonds.Shuffle();
        atom_b_bonds.Shuffle();

        // pick a bonded atom from atom a that is contained within keepindices
        size_t bonded_a_index( util::GetUndefinedSize_t());
        for
        (
            auto atom_a_itr( atom_a_bonds.Begin()), atom_a_itr_end( atom_a_bonds.End());
            atom_a_itr != atom_a_itr_end;
            ++atom_a_itr
        )
        {
          // make sure this is an allowed molecule a index
          bool keep( false);
          for( size_t keep_a_i( 0), keep_a_sz( m_Parent->m_KeepIndicesA.GetSize()); keep_a_i < keep_a_sz; ++keep_a_i)
          {
            size_t i( molecule_a.GetAtomVector().GetAtomIndex( atom_a_1.GetBonds()( *atom_a_itr).GetTargetAtom()));
            if( i == m_Parent->m_KeepIndicesA( keep_a_i))
            {
              bonded_a_index = i;
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }
        }
        if( !util::IsDefined( bonded_a_index))
        {
          return this->operator ()( MOLECULE_A);
        }

        // pick a bonded atom from atom b that is contained within keepindices
        size_t bonded_b_index( util::GetUndefinedSize_t());
        for
        (
            auto atom_b_itr( atom_b_bonds.Begin()), atom_b_itr_end( atom_b_bonds.End());
            atom_b_itr != atom_b_itr_end;
            ++atom_b_itr
        )
        {
          // make sure this is an allowed molecule a index
          bool keep( false);
          for( size_t keep_b_i( 0), keep_b_sz( m_Parent->m_KeepIndicesB.GetSize()); keep_b_i < keep_b_sz; ++keep_b_i)
          {
            size_t i( m_MolBPtr->GetAtomVector().GetAtomIndex( atom_b_1.GetBonds()( *atom_b_itr).GetTargetAtom()));
            if( i == m_Parent->m_KeepIndicesB( keep_b_i))
            {
              bonded_b_index = i;
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }
        }
        if( !util::IsDefined( bonded_b_index))
        {
          return this->operator ()( MOLECULE_A);
        }

        // get the bonded atoms
        const AtomConformationalInterface &atom_a_2( molecule_a.GetAtomVector()( bonded_a_index));
        const AtomConformationalInterface &atom_b_2( m_MolBPtr->GetAtomVector()( bonded_b_index));

        // now get the second bonded atom for each molecule
        size_t bonded_a_index_two( util::GetUndefinedSize_t());
        for
        (
            auto atom_a_itr( atom_a_bonds.Begin()), atom_a_itr_end( atom_a_bonds.End());
            atom_a_itr != atom_a_itr_end;
            ++atom_a_itr
        )
        {
          // make sure this is an allowed molecule a index
          bool keep( false);
          for( size_t keep_a_i( 0), keep_a_sz( m_Parent->m_KeepIndicesA.GetSize()); keep_a_i < keep_a_sz; ++keep_a_i)
          {
            size_t i( molecule_a.GetAtomVector().GetAtomIndex( atom_a_1.GetBonds()( *atom_a_itr).GetTargetAtom()));
            if( i == m_Parent->m_KeepIndicesA( keep_a_i) && i != bonded_a_index)
            {
              bonded_a_index_two = i;
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }
        }
        if( !util::IsDefined( bonded_a_index_two))
        {
          return this->operator ()( MOLECULE_A);
        }

        // pick a bonded atom from atom b that is contained within keepindices
        size_t bonded_b_index_two( util::GetUndefinedSize_t());
        for
        (
            auto atom_b_itr( atom_b_bonds.Begin()), atom_b_itr_end( atom_b_bonds.End());
            atom_b_itr != atom_b_itr_end;
            ++atom_b_itr
        )
        {
          // make sure this is an allowed molecule a index
          bool keep( false);
          for( size_t keep_b_i( 0), keep_b_sz( m_Parent->m_KeepIndicesB.GetSize()); keep_b_i < keep_b_sz; ++keep_b_i)
          {
            size_t i( m_MolBPtr->GetAtomVector().GetAtomIndex( atom_b_1.GetBonds()( *atom_b_itr).GetTargetAtom()));
            if( i == m_Parent->m_KeepIndicesB( keep_b_i) && i != bonded_b_index)
            {
              bonded_b_index_two = i;
              keep = true;
              break;
            }
          }
          if( !keep)
          {
            continue;
          }
        }
        if( !util::IsDefined( bonded_b_index_two))
        {
          return this->operator ()( MOLECULE_A);
        }

        // get the bonded atoms
        const AtomConformationalInterface &atom_a_3( molecule_a.GetAtomVector()( bonded_a_index_two));
        const AtomConformationalInterface &atom_b_3( m_MolBPtr->GetAtomVector()( bonded_b_index_two));

        util::SiPtrVector< const linal::Vector3D> mol_a_atoms, mol_b_atoms;
        mol_a_atoms.PushBack( util::ToSiPtr( atom_a_1.GetPosition()));
        mol_a_atoms.PushBack( util::ToSiPtr( atom_a_2.GetPosition()));
        mol_a_atoms.PushBack( util::ToSiPtr( atom_a_3.GetPosition()));
        const bool swap_pos
        (
          linal::Distance( atom_a_2.GetPosition(), atom_b_2.GetPosition()) + linal::Distance( atom_a_3.GetPosition(), atom_b_3.GetPosition())
          > linal::Distance( atom_a_2.GetPosition(), atom_b_3.GetPosition()) + linal::Distance( atom_a_3.GetPosition(), atom_b_2.GetPosition())
        );
        mol_b_atoms.PushBack( util::ToSiPtr( atom_b_1.GetPosition()));
        mol_b_atoms.PushBack( util::ToSiPtr( swap_pos ? atom_b_3.GetPosition() : atom_b_2.GetPosition()));
        mol_b_atoms.PushBack( util::ToSiPtr( swap_pos ? atom_b_2.GetPosition() : atom_b_3.GetPosition()));

        quality::RMSD calculator;
        math::TransformationMatrix3D transformed_coords( calculator.CalculateSuperimposition( mol_a_atoms, mol_b_atoms));
        molecule_a.Transform( transformed_coords);

        m_Scheme = "BondAlign3 ";
      }

      else if( chosen_scheme == "BondAlignInf")
      {
        // superimposes matched atom pairs
        if( !ConformationComparisonPsiField::BondAlignInf( molecule_a, *m_MolBPtr, m_Parent->m_KeepIndicesA, m_Parent->m_KeepIndicesB))
        {
          return this->operator ()( MOLECULE_A);
        }
        m_Scheme = "BondAlignInf ";
      }

      // rotate molecule randomly with magnitude between 0 and 5 degrees
      else if( chosen_scheme == "SmallRot")
      {
        math::RotationMatrix3D rot;
        rot.SetRand( 0.08);
        linal::Vector3D molecule_a_centered( molecule_a.GetCenter());
        molecule_a.Translate( -molecule_a_centered);
        molecule_a.Rotate( rot);
        molecule_a.Translate( molecule_a_centered);
        m_Scheme = "SmallRot ";
      }

      else if( chosen_scheme == "SmallTranslate")
      {
        // translate molecule a random distance between 0 and 1 angstroms
        const double mag( random::GetGlobalRandom().Random( 1.0));
        linal::Vector3D trans( mag, 0.0, 0.0);
        molecule_a.Translate( trans.Rotate( math::RotationMatrix3D().SetRand()));
        m_Scheme = "SmallTranslate " + util::Format()( mag);
      }

      //whole conformational swap
      else if( chosen_scheme == "ConformerSwap")
      {
        FragmentEnsemble::const_iterator mol_itr( m_MolEnsA->Begin());
        size_t pos( random::GetGlobalRandom().Random< double>( m_MolEnsA->GetSize()));
        std::advance( mol_itr, pos);
        util::SiPtrVector< const linal::Vector3D> scaffold_coords( molecule_a.GetAtomCoordinates());
        util::SiPtrVector< const linal::Vector3D> molecule_coords( mol_itr->GetAtomCoordinates());
        math::TransformationMatrix3D transform( quality::RMSD::SuperimposeCoordinates( scaffold_coords, molecule_coords));
        molecule_a = *mol_itr;
        molecule_a.Transform( transform);
        m_Scheme = "ConformerSwap ";
      }

      //And now output the transformed molecule with the rest of the mutate info
      util::ShPtr< FragmentComplete> mutation( new FragmentComplete( molecule_a));
      mutation->ShareCache( MOLECULE_A);
      math::MutateResult< FragmentComplete> result( mutation, *this);

      // Outputs each transformation to an .sdf file
      // Primarily used for debug purposes
      if( !m_MovieFilename.empty())
      {
        io::OFStream out;
        io::File::MustOpenOFStream( out, m_MovieFilename + ".sdf.gz", std::ios::app);
        mutation->WriteMDL( out);
      }

      return result;
    }

    //! @brief selects the method with which to link the fragments
    //! @return a string indicating the method to use
    std::string ConformationComparisonPsiField::ChooseAlignmentMove( const bool &CONFS) const
    {

        // determine what type of linkage to perform
        float sum_probs( 0.0);

        // different probabilities depending on whether or not we are extending the fragment internally
        CONFS ?
        sum_probs = m_FlipProb + m_BigRotProb + m_BondSwapProb + m_BondAlignProb + m_BondAlign3Prob + m_BondAlignInfProb + m_SmallRotProb + m_SmallTranslateProb + m_ConformerSwapProb :
        sum_probs = m_FlipProb + m_BigRotProb + m_BondSwapProb + m_BondAlignProb + m_BondAlign3Prob + m_BondAlignInfProb + m_SmallRotProb + m_SmallTranslateProb;

        // if everything is zero then return an empty string
        if( !sum_probs)
        {
          return std::string();
        }

        // options
        storage::Vector< std::string> options;
        options.PushBack( "Flip");
        options.PushBack( "BigRot");
        options.PushBack( "BondSwap");
        options.PushBack( "BondAlign");
        options.PushBack( "BondAlign3");
        options.PushBack( "BondAlignInf");
        options.PushBack( "SmallRot");
        options.PushBack( "SmallTranslate");
        if( CONFS)
        {
          options.PushBack( "ConformerSwap");
        }

        // probs
        storage::Vector< float> probs;
        probs.PushBack( m_FlipProb);
        probs.PushBack( m_BigRotProb);
        probs.PushBack( m_BondSwapProb);
        probs.PushBack( m_BondAlignProb);
        probs.PushBack( m_BondAlign3Prob);
        probs.PushBack( m_BondAlignInfProb);
        probs.PushBack( m_SmallRotProb);
        probs.PushBack( m_SmallTranslateProb);
        if( CONFS)
        {
          probs.PushBack( m_ConformerSwapProb);
        }

        // obtain a link method
        float rand( random::GetGlobalRandom().Random< float>( 0, sum_probs));
        float cumulative_weighted_sum( 0.0);
        size_t i( 0);
        for( ; i < options.GetSize(); ++i)
        {
          cumulative_weighted_sum += probs( i);
          if( rand < cumulative_weighted_sum)
          {
            break;
          }
        }

        // return the link method string
        return options( i);
    }

     //! @brief return parameters for member data that are set up from the labels
     //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonPsiField::GetSerializer() const
    {
      io::Serializer parameters( ConformationComparisonPropertyFieldCorrelation::GetSerializer());
      parameters.SetClassDescription( "Aligns two molecules and finds best property RMSD based on ConformationComparisonPropertyFieldCorrelation scoring.");

      parameters.AddInitializer
      (
        "iterations",
        "number_of_mc_iterations",
        io::Serialization::GetAgent( &m_IterationNumber),
        "400"
      );

      parameters.AddInitializer
      (
        "max_unimproved_steps",
        "limit number unimproved mc iterations",
        io::Serialization::GetAgent( &m_LimitNumber),
        "160"
      );

      parameters.AddInitializer
      (
        "number_rigid_trajectories",
        "number of mc trajectories for rigid (single conformer) alignment; "
        "each will be initialized with a uniform max atom distance between the lower and upper limits",
        io::Serialization::GetAgent( &m_NumberRigidTrajectories),
        "5"
      );

      parameters.AddInitializer
      (
        "number_outputs",
        "number of best pairs to output as SDF files",
        io::Serialization::GetAgent( &m_NumberOutputSDF),
        "1"
      );

      parameters.AddInitializer
      (
        "align_to_scaffold",
        "include AlignToScaffold as a move in addition to mutates",
        io::Serialization::GetAgent( &m_AlignToScaffold),
        "false"
      );

      parameters.AddInitializer
      (
        "initial_rand_rotation",
        "perform an initial random rotation within the range 0 - 180 degrees; "
        "may be useful to assist with benchmarking",
        io::Serialization::GetAgent( &m_InitialRotation),
        "false"
      );

      parameters.AddOptionalInitializer
      (
        "output_aligned_mol_a",
        "if true output molecule(s) A pose with the specified name",
        io::Serialization::GetAgent( &m_OutputA)
      );

      parameters.AddOptionalInitializer
      (
        "output_aligned_mol_b",
        "if true output molecule(s) B pose with the specified name",
        io::Serialization::GetAgent( &m_OutputB)
      );

      parameters.AddInitializer
      (
        "exclusion_indices_a",
        "atom indices (0-indexed) to be ignored during alignment scoring in molecule A",
        io::Serialization::GetAgent( &m_ExclusionAtomsA),
        ""
      );

      parameters.AddInitializer
      (
        "exclusion_indices_b",
        "atom indices (0-indexed) to be ignored during alignment scoring in molecule B",
        io::Serialization::GetAgent( &m_ExclusionAtomsB),
        ""
      );

      parameters.AddInitializer
      (
        "pose_tolerance",
        "real space symmetry rmsd tolerance to discriminate unique poses for output",
        io::Serialization::GetAgent( &m_PoseTolerance),
        "0.125"
      );

      parameters.AddInitializer
      (
        "pose_score_threshold",
        "threshold below which the alignment is accepted as a valid pose",
        io::Serialization::GetAgent( &m_PoseScoreThreshold),
        "2.0"
      );

      parameters.AddInitializer
      (
        "flip_prob",
        "relative probability of performing a Flip alignment move",
        io::Serialization::GetAgent( &m_FlipProb),
        "0.06"
      );

      parameters.AddInitializer
      (
        "big_rot_prob",
        "relative probability of performing a BigRotation alignment move",
        io::Serialization::GetAgent( &m_BigRotProb),
        "0.06"
      );

      parameters.AddInitializer
      (
        "bond_swap_prob",
        "relative probability of performing a BondSwap alignment move",
        io::Serialization::GetAgent( &m_BondSwapProb),
        "0.06"
      );

      parameters.AddInitializer
      (
        "bond_align_prob",
        "relative probability of performing a BondAlign alignment move",
        io::Serialization::GetAgent( &m_BondAlignProb),
        "0.30"
      );

      parameters.AddInitializer
      (
        "bond_align3_prob",
        "relative probability of performing a BondAlign3 alignment move",
        io::Serialization::GetAgent( &m_BondAlign3Prob),
        "0.20"
      );

      parameters.AddInitializer
      (
        "bond_align_inf_prob",
        "relative probability of performing a BondAlignInf alignment move",
        io::Serialization::GetAgent( &m_BondAlignInfProb),
        "0.10"
      );

      parameters.AddInitializer
      (
        "small_rot_prob",
        "relative probability of performing a SmallRotation alignment move",
        io::Serialization::GetAgent( &m_SmallRotProb),
        "0.12"
      );

      parameters.AddInitializer
      (
        "small_trans_prob",
        "relative probability of performing a SmallTranslation alignment move",
        io::Serialization::GetAgent( &m_SmallTranslateProb),
        "0.12"
      );

      parameters.AddInitializer
      (
        "conf_swap_prob",
        "relative probability of performing a ConformerSwap alignment move",
        io::Serialization::GetAgent( &m_ConformerSwapProb),
        "0.06"
      );

      parameters.AddOptionalInitializer
      (
        "movie",
        "string - save snapshots from MCM trajectory - WARNING - large file sizes easily accumulate",
        io::Serialization::GetAgent( &m_Movie)
      );

      parameters.AddInitializer
      (
        "potential_map_filename",
        "string - filename of pair potential map",
        io::Serialization::GetAgent( &m_PotentialFilename),
        "new_potentials.txt"
      );

      return parameters;
    }
    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ConformationComparisonPsiField::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      for( size_t best_pairs_index( 0); 1; ++best_pairs_index)
      {
        // remove/overwrite files with identical names
        const std::string filename_A( m_OutputA + "_" + std::to_string( best_pairs_index) + ".sdf");
        io::DirectoryEntry entry_A( filename_A);
        if( entry_A.DoesExist())
        {
          entry_A.Remove();
        }
        else
        {
          break;
        }
        const std::string filename_B( m_OutputB + "_" + std::to_string( best_pairs_index) + ".sdf");
        io::DirectoryEntry entry_B( filename_B);
        if( entry_B.DoesExist())
        {
          entry_B.Remove();
        }
        else
        {
          break;
        }

      }

      // read exclusion atoms
      if( m_ExclusionAtomsA.size())
      {
        m_ExclusionIndicesA = util::SplitStringToNumerical< size_t>( m_ExclusionAtomsA);
      }
      if( m_ExclusionAtomsB.size())
      {
        m_ExclusionIndicesB = util::SplitStringToNumerical< size_t>( m_ExclusionAtomsB);
      }

      // initialize pair potentials
      storage::Map< storage::Pair< AtomType, AtomType>, linal::Vector< double> > mapped_atom_pair_potentials_vec;
      io::IFStream file;
      io::File::MustOpenIFStream( file, score::Score::AddHistogramPath( m_PotentialFilename));
      io::Serialize::Read( mapped_atom_pair_potentials_vec, file);
      io::File::CloseClearFStream( file);
      size_t itr_index( 0);
      for
      (
          auto mapped_itr( mapped_atom_pair_potentials_vec.Begin()), mapped_itr_end( mapped_atom_pair_potentials_vec.End());
          mapped_itr != mapped_itr_end;
          ++mapped_itr, ++itr_index
      )
      {
        m_PairPotentials[ mapped_itr->first] = math::CubicSplineDamped().Train( 0.0, 0.25, mapped_itr->second, 0, 0);
      }

      return ConformationComparisonPropertyFieldCorrelation::ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }

  } // namespace chemistry
} // namespace bcl
