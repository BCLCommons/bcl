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
#include "chemistry/bcl_chemistry_sample_conformations.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_vdw_score.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_rmsd.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_symmetry_rmsd.h"
#include "chemistry/bcl_chemistry_fragment_align_to_scaffold.h"
#include "chemistry/bcl_chemistry_fragment_probability_score.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_mutate_chirality.h"
#include "chemistry/bcl_chemistry_mutate_fragment.h"
#include "chemistry/bcl_chemistry_mutate_multi_fragment.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_mutate_combine.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_sum_function.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

// Uncomment the next line to perform profiling of the tree search
#define BCL_PROFILE_SampleConformations
#ifdef BCL_PROFILE_SampleConformations
#include "util/bcl_util_stopwatch.h"
#endif

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SampleConformations::s_Instance
    (
      GetObjectInstances().AddInstance( new SampleConformations())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SampleConformations::SampleConformations() :
      m_SearchFragmentLibrary( new SearchFragmentLibraryFromTree()),
      m_ConformationComparison( ConformationComparisonBySymmetryRmsd()),
      m_Tolerance( 0.25),
      m_NumberConformations( 200),
      m_MaxIterations( 200),
      m_ChangeChirality( false),
      m_RandomDihedralChangeProbability( 0.01),
      m_Generate3D( true),
      m_Cluster( false),
      m_SampleDihedralAngles( true),
      m_SampleRingConformations( true),
      m_SampleBondAnglesAndLengths( true),
      m_MaxNumberClashResolutionCycles( 2),
      m_MaxClashScore( 0.1),
      m_MDL( "SampleByParts")
    {
    }

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
    SampleConformations::SampleConformations
    (
      const RotamerLibraryInterface &ROTAMER_LIBRARY,
      const std::string &CONFORMATION_COMPARER,
      const double &TOLERANCE,
      const size_t NUMBER_CONFORMATIONS,
      const size_t &MAX_ITERATIONS,
      const bool CHANGE_CHIRALITY,
      const double &RANDOM_DIHEDRAL_CHANGE_CHANCE,
      bool GENERATE_3D,
      double CLASH_TOLERANCE,
      bool CLUSTER,
      double NUMBER_CLASH_RESOLUTION_CYCLES,
      std::string SCORES_FILENAME
    ) :
      m_SearchFragmentLibrary( new SearchFragmentLibraryFromTree( ROTAMER_LIBRARY)),
      m_ConformationComparison( CONFORMATION_COMPARER),
      m_Tolerance( TOLERANCE),
      m_NumberConformations( NUMBER_CONFORMATIONS),
      m_MaxIterations( MAX_ITERATIONS),
      m_ChangeChirality( CHANGE_CHIRALITY),
      m_RandomDihedralChangeProbability( RANDOM_DIHEDRAL_CHANGE_CHANCE),
      m_Generate3D( GENERATE_3D),
      m_Cluster( CLUSTER),
      m_OutputIndividualScoresFilename( SCORES_FILENAME),
      m_SampleDihedralAngles( true),
      m_SampleRingConformations( true),
      m_SampleBondAnglesAndLengths( true),
      m_MaxNumberClashResolutionCycles( NUMBER_CLASH_RESOLUTION_CYCLES),
      m_MaxClashScore( CLASH_TOLERANCE),
      m_MDL( "SampleByParts")
    {
      if( m_Cluster && m_MaxIterations < 2 * m_NumberConformations)
      {
        BCL_MessageStd( "Disabled conformation clustering as it is futile without sampling at least twice the number of molecules desired");
        m_Cluster = false;
      }
      if( !m_Tolerance && !m_Cluster)
      {
        m_ConformationComparison.Reset();
      }
      
    }

    //! @brief clone constructor
    SampleConformations *SampleConformations::Clone() const
    {
      return new SampleConformations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SampleConformations::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SampleConformations::GetAlias() const
    {
      static const std::string s_name_isolate( "SampleConformations");
      return s_name_isolate;
    }

    //! object that searches for fragments of the molecule of interest from the rotamer library
    const util::Implementation< RotamerLibraryInterface> &SampleConformations::GetRotamerLibrary() const
    {
      return m_SearchFragmentLibrary->GetRotamerLibrary();
    }

    //! method for comparing conformations
    const std::string &SampleConformations::GetConfComparer() const
    {
      return m_ConformationComparison->GetAlias();
    }

    //! tolerance for conformation comparer
    const double &SampleConformations::GetTolerance() const
    {
      return m_Tolerance;
    }

    //! number of conformations desired
    const size_t &SampleConformations::GetMaxNumConfs() const
    {
      return m_NumberConformations;
    }

    //! number of iterations for searching conformers of a molecule
    const size_t &SampleConformations::GetMaxIterations() const
    {
      return m_MaxIterations;
    }

    //! change chirality
    const bool &SampleConformations::GetIfChangeChirality() const
    {
      return m_ChangeChirality;
    }

    //! relative probability of just making random dihedral swaps
    const double &SampleConformations::GetRandDihedChangeProb() const
    {
      return m_RandomDihedralChangeProbability;
    }

    //! boolean to specify if 3D structure needs to be generated from scratch
    const bool &SampleConformations::GetIfGenerate3D() const
    {
      return m_Generate3D;
    }

    //! whether to cluster the returned conformations so as to nominally optimize coverage of conformational space.
    const bool &SampleConformations::GetIfCluster() const
    {
      return m_Cluster;
    }

    //! Number of conformations to choose entirely based on diversity from existing structures
    const size_t &SampleConformations::GetClusterNDiversity() const
    {
      return m_ClusterNBasedOnDiversity;
    }

    //! clash resolver
    const double &SampleConformations::GetClashResolution() const
    {
      return m_MaxNumberClashResolutionCycles;
    }

    //! Absolute maximum clash score
    const double &SampleConformations::GetMaxClashScore() const
    {
      return m_MaxClashScore;
    }

  ////////////////
  // operations //
  ////////////////

    //! method for comparing conformations
    void SampleConformations::SetConfComparer( const std::string &COMPARER)
    {
      m_ConformationComparison = COMPARER;
    }

    //! tolerance for conformation comparer
    void SampleConformations::SetTolerance( const double TOLERANCE)
    {
      m_Tolerance = TOLERANCE;
    }

    //! number of conformations desired
    void SampleConformations::SetMaxNumConfs( const size_t N_MAX_CONFS)
    {
      m_NumberConformations = N_MAX_CONFS;
    }

    //! number of iterations for searching conformers of a molecule
    void SampleConformations::SetMaxIterations( const size_t N_MAX_ITERATIONS)
    {
      m_MaxIterations = N_MAX_ITERATIONS;
    }

    //! change chirality
    void SampleConformations::SetIfChangeChirality( const bool CHANGE_CHIRALITY)
    {
      m_ChangeChirality = CHANGE_CHIRALITY;
    }

    //! relative probability of just making random dihedral swaps
    void SampleConformations::SetRandDihedChangeProb( const double RAND_DIHED_PROB)
    {
      m_RandomDihedralChangeProbability = RAND_DIHED_PROB;
    }

    //! boolean to specify if 3D structure needs to be generated from scratch
    void SampleConformations::SetIfGenerate3D( const bool GENERATE_3D)
    {
      m_Generate3D = GENERATE_3D;
    }

    //! whether to cluster the returned conformations so as to nominally optimize coverage of conformational space.
    void SampleConformations::SetIfCluster( const bool CLUSTER)
    {
      m_Cluster = CLUSTER;
    }

    //! Number of conformations to choose entirely based on diversity from existing structures
    void SampleConformations::SetClusterNDiversity( const size_t N_FORCED_DIVERSE_CLUSTERS)
    {
      m_ClusterNBasedOnDiversity = N_FORCED_DIVERSE_CLUSTERS;
    }

    //! clash resolver
    void SampleConformations::SetClashResolution( const double CLASH_RESOLUTION)
    {
      m_MaxNumberClashResolutionCycles = CLASH_RESOLUTION;
    }

    //! Absolute maximum clash score
    void SampleConformations::SetMaxClashScore( const double CLASH_SCORE)
    {
      m_MaxClashScore = CLASH_SCORE;
    }

    //! @brief change the output score filename
    void SampleConformations::SetOutputScoreFilename( const std::string &FILENAME)
    {
      m_OutputIndividualScoresFilename = FILENAME;
    }

    //! @brief set the sampling types to limit or expand the types of conformational sampling will be added
    //! @param SAMPLE_DIHEDRALS whether dihedral angles can be substantially changed (substantially outside the 30 degree bins)
    //! @param SAMPLE_RINGS whether to sample ring conformations
    //! @param SAMPLE_BOND_ANGLES whether to sample bond angles (and lengths)
    //! @param SAMPLE_CHIRALITY whether to sample chirality
    void SampleConformations::SetSamplingPreferences
    (
      const bool &SAMPLE_DIHEDRALS,
      const bool &SAMPLE_RINGS,
      const bool &SAMPLE_BOND_ANGLES,
      const bool &SAMPLE_CHIRALITY
    )
    {
      m_SampleDihedralAngles = SAMPLE_DIHEDRALS;
      m_SampleRingConformations = SAMPLE_RINGS;
      m_SampleBondAnglesAndLengths = SAMPLE_BOND_ANGLES;
      m_ChangeChirality = SAMPLE_CHIRALITY;
    }

    //! @brief set sample by parts indices from MDL property line
    void SampleConformations::SetSampleByPartsMDL( const FragmentComplete &MOLECULE, const std::string &MDL) const
    {
      m_MDLIndicesToSample = linal::Vector< size_t>( MOLECULE.GetStoredProperties().GetMDLPropertyAsVector( MDL));
    }

    //! @brief set sample by parts indices
    void SampleConformations::SetSampleByPartsIndices( const storage::Vector< size_t> &INDICES_TO_SAMPLE) const
    {
      m_IndicesToSample = INDICES_TO_SAMPLE;
    }

    //! @brief set sample by parts indices
    void SampleConformations::SetSampleByPartsIndices( const std::string &INDICES_TO_SAMPLE) const
    {
      if( INDICES_TO_SAMPLE.size())
      {
        m_IndicesToSample.Reset();
        m_IndicesToSample = util::SplitStringToNumerical< size_t>( INDICES_TO_SAMPLE);
      }
    }

    //! @brief set the sample by parts atoms based on fragments
    void SampleConformations::SetSampleByPartsFragments
    (
      const FragmentComplete &FRAGMENT,
      const FragmentEnsemble &FRAGMENTS_TO_SAMPLE,
      const bool COMPLEMENT
    ) const
    {
      m_FragmentsToSample = FragmentTrackMutableAtoms::SetMutableFragments( FRAGMENT, FRAGMENTS_TO_SAMPLE, COMPLEMENT);
    }

    //! @brief set the sample by parts atoms based on fragments
    void SampleConformations::SetSampleByPartsFragments
    (
      const FragmentComplete &FRAGMENT,
      const std::string &FRAGMENTS_TO_SAMPLE,
      const bool COMPLEMENT
    ) const
    {
      if( FRAGMENTS_TO_SAMPLE.size())
      {
        m_Mutex.Lock();
        io::IFStream input;
        io::File::MustOpenIFStream( input, FRAGMENTS_TO_SAMPLE);
        FragmentEnsemble reference_mol;
        reference_mol.ReadMoreFromMdl( input);
        io::File::CloseClearFStream( input);
        m_FragmentsToSample = FragmentTrackMutableAtoms::SetMutableFragments( FRAGMENT, reference_mol, COMPLEMENT);
        m_Mutex.Unlock();
      }
    }

    //! @brief sets the final sample by parts atoms used by operator() from member data
    void SampleConformations::PoolSampleByPartsIndices() const
    {
      // get the MDL property values and store unique
      storage::Set< size_t> sbp_indices_uniq;
      if( m_MDLIndicesToSample.GetSize())
      {
        sbp_indices_uniq.InsertElements( m_MDLIndicesToSample.Begin(), m_MDLIndicesToSample.End());
      }

      // add the explicitly set indices
      if( m_IndicesToSample.GetSize())
      {
        sbp_indices_uniq.InsertElements( m_IndicesToSample.Begin(), m_IndicesToSample.End());
      }

      // add the fragment indices
      if( m_FragmentsToSample.GetSize())
      {
        sbp_indices_uniq.InsertElements( m_FragmentsToSample.Begin(), m_FragmentsToSample.End());
      }

      // set final
      m_SampleByPartsIndices = linal::Vector< size_t>( sbp_indices_uniq.Begin(), sbp_indices_uniq.End());
    }

    //! @brief set the filter criteria for determining what to keep from the final ensemble
    //! @param ENSEMBLE conformers to be filtered
    //! @param STARTING_MOL molecule (conformer) against which the ensemble will be filtered
    //! @param CONFORMATION_COMPARER the measure of similarity to use
    //! @param COMPARISON are the conformers less, greater, etc. than TOLERANCE according
    //! to CONFORMATION_COMPARER to starting conformer
    //! @param TOLERANCE comparison threshold for conformer relation
    //! @param ATOMS_TO_COMPARE atom subset to compare
    FragmentEnsemble SampleConformations::FilterConformerEnsemble
    (
      const FragmentEnsemble &ENSEMBLE,
      const FragmentComplete &STARTING_MOL,
      const std::string &CONFORMATION_COMPARER,
      const math::Comparisons< float>::Comparison &COMPARISON,
      const float &TOLERANCE,
      const storage::Vector< size_t> &ATOMS_TO_COMPARE
    )
    {
      // we will return a different ensemble than the one we started with
      FragmentEnsemble ensemble( ENSEMBLE);
      FragmentComplete starting_mol( STARTING_MOL);

      // initialize a comparison object with the provided string
      util::Implementation< ConformationComparisonInterface> conf_comparer;
      if( CONFORMATION_COMPARER.size())
      {
        conf_comparer = util::Implementation< ConformationComparisonInterface>( CONFORMATION_COMPARER);
      }
      else if( m_ConformationComparison.IsDefined())
      {
        conf_comparer = util::Implementation< ConformationComparisonInterface>( m_ConformationComparison.GetAlias());
      }
      else
      {
        BCL_MessageStd( "No comparison type specified for FilterConformerEnsemble; defaulting to RMSD");
        conf_comparer = util::Implementation< ConformationComparisonInterface>( ConformationComparisonByRmsd());
      }

      // to do this as expected we really need to align the conformers fully or by
      // the atoms being compared
      FragmentAlignToScaffold ats
      (
        ConformationGraphConverter::AtomComparisonType::e_ElementType,
        ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness,
        size_t( 2)
      );

      // If we are only comparing a subset of atoms then remove undesired atoms from atom vectors
      if( ATOMS_TO_COMPARE.GetSize())
      {
        // reference molecule
        AtomVector< AtomComplete> reference_atom_v( starting_mol.GetAtomVector());
        reference_atom_v.Reorder( ATOMS_TO_COMPARE);
        starting_mol = FragmentComplete( reference_atom_v, starting_mol.GetName());

        // all conformers in ensemble
        for
        (
            auto conf_itr( ensemble.Begin()), conf_itr_end( ensemble.End());
            conf_itr != conf_itr_end;
            ++conf_itr
        )
        {
          AtomVector< AtomComplete> atom_v( conf_itr->GetAtomVector());
          atom_v.Reorder( ATOMS_TO_COMPARE);
          *conf_itr = FragmentComplete( atom_v, conf_itr->GetName());
        }
      }

      // compare each conformer to reference
      FragmentEnsemble filtered_ensemble;
      storage::Vector< FragmentComplete> ens( ENSEMBLE.Begin(), ENSEMBLE.End());
      size_t conf_index( 0);
      for
      (
          auto conf_itr( ensemble.Begin()), conf_itr_end( ensemble.End());
          conf_itr != conf_itr_end;
          ++conf_itr, ++conf_index
      )
      {
        // align all atoms
        ats.AlignToScaffold( *conf_itr, starting_mol);

        float comparison_value( ( *conf_comparer)( *conf_itr, starting_mol));
        if( ( **COMPARISON)( comparison_value, TOLERANCE))
        {
          // true
          ats.AlignToScaffold( ens( conf_index), starting_mol);
          ens( conf_index).GetStoredPropertiesNonConst().SetMDLProperty( "FilterEnsembleConformerComparer", conf_comparer.GetAlias());
          ens( conf_index).GetStoredPropertiesNonConst().SetMDLProperty( "FilterEnsembleComparerType", COMPARISON.GetName());
          ens( conf_index).GetStoredPropertiesNonConst().SetMDLProperty( "FilterEnsembleComparerValue", util::Format()( comparison_value));
          filtered_ensemble.PushBack( ens( conf_index));
        }
      }
      return filtered_ensemble;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief returns an ensemble of generated conformers
    //! @param MOLECULE molecule of interest whose conformations need to be sampled
    //! @return an ensemble of generated conformers
    storage::Pair< FragmentEnsemble, FragmentEnsemble> SampleConformations::operator()( const FragmentComplete &MOLECULE) const
    {
      if( MOLECULE.GetNumberValences())
      {
        BCL_MessageCrt
        (
          "Molecule contains unsatisfied valences; ConfScore will be very unreliable and conformers will suffer "
          "greatly in accuracy"
        );
      }
      storage::Set< size_t> parts_to_sample;

      // if 2D structure is given generate a 3D conformation for the molecule of interest
      FragmentComplete generated_mol( MOLECULE);

      // Insert a conformation score of 1 (highest possible = worst score), which will be used if a conformation cannot
      // be generated; otherwise it will be overwritten by the actual score
      generated_mol.StoreProperty
      (
        "ConfScore",
        storage::Vector< double>( 1, double( 1.0))
      );
      generated_mol.StoreProperty
      (
        "ConfScorePBest",
        storage::Vector< double>( 1, double( 1.0))
      );
      generated_mol.GetStoredPropertiesNonConst().RemoveConformationalDescriptors();

      // test for disconnected molecule
      if( !graph::Connectivity::IsConnected( ConformationGraphConverter()( MOLECULE)))
      {
        BCL_MessageCrt( "Cannot sample conformations for a disconnected molecule, returning the input conformation");
        FragmentEnsemble ensemble;
        ensemble.PushBack( generated_mol);
        return storage::Pair< FragmentEnsemble, FragmentEnsemble>( ensemble, FragmentEnsemble());
      }

      // collect sample by parts indices
      SetSampleByPartsMDL( MOLECULE, m_MDL);
      SetSampleByPartsIndices( m_IndicesToSampleString);
      SetSampleByPartsFragments( MOLECULE, m_ReferenceFragments, false);
      PoolSampleByPartsIndices();

      // rotamer bin signature of a particular rotamer
      const bool are_sampling_by_parts( m_SampleByPartsIndices.GetSize());
      for
      (
        linal::Vector< size_t>::const_iterator itr_indices( m_SampleByPartsIndices.Begin()), itr_indices_end( m_SampleByPartsIndices.End());
        itr_indices != itr_indices_end;
        ++itr_indices
      )
      {
        size_t index( *itr_indices - size_t( 0));
        parts_to_sample.Insert( index);
      }
      MutateChirality chirality_changer( generated_mol);
      if( m_Generate3D && !are_sampling_by_parts)
      {
        generated_mol.IdealizeGeometry();
        AtomVector< AtomComplete> atoms( generated_mol.GetAtomVector());
        chirality_changer.ApplyChirality( atoms, MOLECULE.GetAtomVector(), true);
        chirality_changer.ApplyDoubleBondIsometry( atoms, MOLECULE.GetAtomVector(), true);
        generated_mol = FragmentComplete( atoms, MOLECULE.GetName());
        // Insert a conformation score of 1 (highest possible = worst score), which will be used if a conformation cannot
        // be generated; otherwise it will be overwritten by the actual score
        generated_mol.StoreProperty
        (
          "ConfScore",
          storage::Vector< double>( 1, double( 1.0))
        );
        generated_mol.StoreProperty
        (
          "ConfScorePBest",
          storage::Vector< double>( 1, double( 1.0))
        );
        generated_mol.GetStoredPropertiesNonConst().RemoveConformationalDescriptors();
      }
      else if( m_SampleBondAnglesAndLengths)
      {
        generated_mol.StandardizeBondLengths();
      }

      #ifdef BCL_PROFILE_SampleConformations
      static util::Stopwatch s_find_frags( "FindFrags", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_find_frags.Start();
      #endif

      // search fragments for molecule of interest from the rotamer library
      util::ShPtrVector< SmallMoleculeFragmentIsomorphism> fragments_iso_lib
      (
        m_SampleDihedralAngles || m_SampleRingConformations || m_Generate3D
        ? m_SearchFragmentLibrary->FindFragmentsOfMolecule( generated_mol)
        : util::ShPtrVector< SmallMoleculeFragmentIsomorphism>()
      );

      util::ShPtrVector< RotamerDihedralBondData> pick_iso_saved;

      #ifdef BCL_PROFILE_SampleConformations
      s_find_frags.Stop();
      #endif

      // make a vector to store fragments that have been found for the molecule of interest. Used for creating
      util::ShPtr< FragmentEnsemble> tmp_fragments( new FragmentEnsemble);

      // store fragments in a temporary vector
      for
      (
          util::ShPtrVector< SmallMoleculeFragmentIsomorphism>::const_iterator
          itr_isomorph( fragments_iso_lib.Begin()), itr_isomorph_end( fragments_iso_lib.End());
          itr_isomorph != itr_isomorph_end;
          ++itr_isomorph
      )
      {
        tmp_fragments->PushBack( ( *itr_isomorph)->GetFragment());
      }

      generated_mol.UpdateH();

      // initialize shared variables for the mutates, one to indicate whether the mutate should only select the optimal/most probable
      // selection (so as to rapidly locate the optima) and another to hold the clash scoring object
      util::ShPtr< bool> mutate_to_best_sp( new bool( true));
      util::ShPtr< AtomClashScore> clash_score_sp( new AtomClashScore( true));
      util::ShPtr< AtomClashScore> clash_score_sp_no_h( new AtomClashScore( false));

      #ifdef BCL_PROFILE_SampleConformations
      static util::Stopwatch s_map_frags( "MapFragmentIsomorphisms", util::Message::e_Standard, true, false);
      s_map_frags.Start();
      #endif

      // create a graph where atoms are colored by atom type or chirality and edges by whether edge is a bond ring or not
      ConformationGraphConverter graph_maker
      (
        ConformationGraphConverter::e_AtomType,
        ConfigurationalBondTypeData::e_IsInRing
      );

      // create a graph for the molecule of interest
      const graph::ConstGraph< size_t, size_t> molecule_graph( graph_maker( generated_mol));

      // get mapping of fragments to the molecule of interest
      util::ShPtrVector< RotamerDihedralBondData> updated_bond_mapping
      (
        SmallMoleculeFragmentMapping().MapFragmentIsomorphisms( generated_mol, fragments_iso_lib, parts_to_sample)
      );
      util::ShPtrVector< BondAngleAssignment> bond_angle_assignments
      (
        m_SearchFragmentLibrary->GetBondAngleAssignments( generated_mol, !m_ChangeChirality)
      );

      auto dihedrals( PriorityDihedralAngles().GetDihedralEdges( generated_mol));
      storage::Vector< size_t> dihedral_counts( dihedrals.GetSize(), size_t( 0));
      size_t degrees_freedom_chain( 0), degrees_freedom_ring( 0);
      {
        // because we are interested in sampling conformations, remove dihedral bond data that correspond to rigid rings,
        // which do not matter for sampling or scoring
        auto itr_place( updated_bond_mapping.Begin());
        for
        (
          auto itr_next( updated_bond_mapping.Begin()), itr_end( updated_bond_mapping.End());
          itr_next != itr_end;
          ++itr_next
        )
        {
          if( !m_SampleDihedralAngles && ( *itr_next)->GetFragment().GetNumberDihedralChainBonds())
          {
            continue;
          }
          bool was_added_ring( false);
          if( m_Generate3D && !( *itr_next)->GetFragment().GetNumberDihedralChainBonds())
          {
            bool did_something( false);
            MutateFragment cur_fragment
            (
              *itr_next,
              molecule_graph,
              generated_mol,
              chirality_changer,
              false,
              false,
              false,
              false
            );
            for( size_t tries( 0), mx_tries( 10); tries < mx_tries && !did_something; ++tries)
            {
              math::MutateResult< FragmentComplete> mutated( cur_fragment( generated_mol));
              if( mutated.GetArgument().IsDefined())
              {
                generated_mol = *mutated.GetArgument();
                did_something = true;
              }
            }
            // try again without constraints on chirality, which may be messed up due to bad bond angles (should
            // just do the bond angles part first in this case, but oh well)
            if( !did_something)
            {
              cur_fragment = MutateFragment
              (
                *itr_next,
                molecule_graph,
                generated_mol,
                chirality_changer,
                true,
                false,
                false,
                false
              );
              for( size_t tries( 0), mx_tries( 10); tries < mx_tries && !did_something; ++tries)
              {
                math::MutateResult< FragmentComplete> mutated( cur_fragment( generated_mol));
                if( mutated.GetArgument().IsDefined())
                {
                  generated_mol = *mutated.GetArgument();
                  did_something = true;
                }
              }
              if( did_something)
              {
                BCL_MessageStd
                (
                  "Adding ring fragment required modifying chirality. This will usually be fixed in a later step"
                );
              }
            }
            if( !did_something)
            {
              BCL_MessageVrb( "Failed to add ring fragment: ");
              continue;
            }
            else
            {
              was_added_ring = true;
            }
          }
          for
          (
            auto itr_cb( ( *itr_next)->GetCenterBondIsomorphism().Begin()),
                 itr_cb_end( ( *itr_next)->GetCenterBondIsomorphism().End());
            itr_cb != itr_cb_end;
            ++itr_cb
          )
          {
            ++dihedral_counts( *itr_cb);
          }
          if( !m_SampleRingConformations && ( *itr_next)->GetFragment().GetNumberRingBonds())
          {
            continue;
          }
          // also skip rigid rings that are purely planar; they can aid in neither scoring nor sampling and so it's a large
          // waste of sampling to swap them in and out
          if
          (
            !( *itr_next)->GetFragment().GetNumberDihedralChainBonds()
            && !( *itr_next)->ContainsRingRotamers()
            && m_MaxIterations
          )
          {
            if( was_added_ring)
            {
              continue;
            }
            BCL_MessageDbg( "Constant Fragment");
            MutateFragment cur_fragment
            (
              *itr_next,
              molecule_graph,
              generated_mol,
              chirality_changer,
              false,
              true,
              true
            );
            cur_fragment.GetChooseBestOnlySharedPointer() = mutate_to_best_sp;
            math::MutateResult< FragmentComplete> mutated( cur_fragment( generated_mol));
            if( mutated.GetArgument().IsDefined())
            {
              double clash_before( ( *clash_score_sp_no_h)( generated_mol));
              double clash_after( ( *clash_score_sp_no_h)( *mutated.GetArgument()));
              if( clash_after <= clash_before + 0.0001)
              {
                generated_mol = *mutated.GetArgument();
                BCL_MessageDbg( "Constant Fragment Accepted");
                continue;
              }
              BCL_MessageDbg( "Constant Fragment Rejected " + util::Format()( clash_before) + " " + util::Format()( clash_after));
            }
            else
            {
              BCL_MessageDbg( "Constant Fragment Rejected; could not create");
            }
          }
          if( ( *itr_next)->ContainsRings() && !( *itr_next)->GetFragment().GetNumberDihedralChainBonds())
          {
            ++degrees_freedom_ring;
          }
          if( itr_place != itr_next)
          {
            *itr_place = *itr_next;
          }
          ++itr_place;
        }
        const size_t new_rot_bond_data_size( size_t( std::distance( updated_bond_mapping.Begin(), itr_place)));
        if( new_rot_bond_data_size != updated_bond_mapping.GetSize())
        {
          updated_bond_mapping.Resize( new_rot_bond_data_size);
        }
      }

      storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > chain_uncovered;
      storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > ring_uncovered;

      {
        size_t dpos( 0);
        for
        (
          auto itr_cnts( dihedral_counts.Begin()), itr_cnts_end( dihedral_counts.End());
          itr_cnts != itr_cnts_end;
          ++itr_cnts, ++dpos
        )
        {
          if
          (
            !parts_to_sample.IsEmpty()
            &&
            (
              !parts_to_sample.Contains( dihedrals( dpos).GetIndexLow())
              ||
              !parts_to_sample.Contains( dihedrals( dpos).GetIndexHigh())
            )
          )
          {
            continue;
          }

          if
          (
            MOLECULE.GetAtomVector()( dihedrals( dpos).GetIndexLow()).GetAtomType()->GetFormsOnlyLinearBonds()
            || MOLECULE.GetAtomVector()( dihedrals( dpos).GetIndexHigh()).GetAtomType()->GetFormsOnlyLinearBonds()
          )
          {
            continue;
          }
          if( !dihedrals( dpos).GetEdgeData()->IsBondInRing())
          {
            ++degrees_freedom_chain;
          }
          if( !*itr_cnts)
          {
            if( !dihedrals( dpos).GetEdgeData()->IsBondInRing())
            {
              chain_uncovered.PushBack( dihedrals( dpos));
            }
            else
            {
              ring_uncovered.PushBack( dihedrals( dpos));
            }
          }
        }
      }
      #ifdef BCL_PROFILE_SampleConformations
      s_map_frags.Stop();
      #endif
      if( m_Generate3D && !ring_uncovered.IsEmpty() && !are_sampling_by_parts)
      {
        BCL_MessageCrt( "No conformations available for one or more rings.");
        FragmentEnsemble ensemble;
        return storage::Pair< FragmentEnsemble, FragmentEnsemble>( ensemble, *tmp_fragments);
      }
      if( m_Generate3D)
      {
        generated_mol.StandardizeBondLengths();
        generated_mol.UpdateH();
      }

      #ifdef BCL_PROFILE_SampleConformations
      s_map_frags.Stop();
      #endif

      // sum function to combine scores
      FragmentProbabilityScore fps( updated_bond_mapping, m_ChangeChirality);

      // intialize an object probability distribution to store the default mutates and weights
      util::ShPtrVector< MutateDihedralsInterface> dihedral_mutates;
      util::ShPtrVector< MutateDihedralsInterface> bnd_mutates;
      util::ShPtrVector< MutateDihedralsInterface> dihedral_wiggler_mutates;

      for
      (
        util::ShPtrVector< RotamerDihedralBondData>::const_iterator
          itr_ref( updated_bond_mapping.Begin()), itr_ref_end( updated_bond_mapping.End());
        itr_ref != itr_ref_end;
        ++itr_ref
      )
      {
        if
        (
          !( *itr_ref)->GetFragment().GetNumberDihedralChainBonds()
          && ( *itr_ref)->GetFragment().GetNumberChainBonds()
          && m_MaxIterations
        )
        {
          continue;
        }
        dihedral_mutates.PushBack
        (
          util::ShPtr< MutateDihedralsInterface>
          (
            new MutateFragment
            (
              *itr_ref,
              molecule_graph,
              generated_mol,
              chirality_changer,
              m_ChangeChirality,
              false,
              true
            )
          )
        );
        dihedral_mutates.LastElement()->GetChooseBestOnlySharedPointer() = mutate_to_best_sp;
      }
      if( chain_uncovered.GetSize())
      {
        BCL_MessageStd
        (
          "# of dihedral angles not covered by existing knowledge-based potential: "
          + util::Format()( chain_uncovered.GetSize())
        );
      }
      for
      (
        auto itr_unc( chain_uncovered.Begin()), itr_unc_end( chain_uncovered.End());
        itr_unc != itr_unc_end;
        ++itr_unc
      )
      {
        if( m_SampleDihedralAngles || m_MaxIterations > size_t( 1))
        {
          dihedral_mutates.PushBack
          (
            util::ShPtr< MutateDihedralsInterface>
            (
              new MutateDihedralBond
              (
                MOLECULE.GetAtomVector()( itr_unc->GetIndexLow()),
                MOLECULE.GetAtomVector()( itr_unc->GetIndexHigh()),
                molecule_graph,
                MOLECULE,
                itr_unc->GetEdgeData()->GetNumberOfElectrons() > size_t( 2)
                || itr_unc->GetEdgeData()->GetConjugation() == ConstitutionalBondTypeData::e_Amide,
                !m_SampleDihedralAngles
              )
            )
          );
          dihedral_mutates.LastElement()->GetChooseBestOnlySharedPointer() = mutate_to_best_sp;
        }
      }

      auto all_edges( PriorityDihedralAngles::GetDihedralEdges( MOLECULE));
      storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > rotatable_edges;
      rotatable_edges.AllocateMemory( all_edges.GetSize());
      for( auto itr_edges( all_edges.Begin()), itr_edges_end( all_edges.End()); itr_edges != itr_edges_end; ++itr_edges)
      {
        if( itr_edges->GetEdgeData()->IsBondInRing() || itr_edges->GetEdgeData()->GetNumberOfElectrons() > size_t( 4))
        {
          continue;
        }
        if
        (
          MOLECULE.GetAtomVector()( itr_edges->GetIndexLow()).GetAtomType()->GetFormsOnlyLinearBonds()
          || MOLECULE.GetAtomVector()( itr_edges->GetIndexHigh()).GetAtomType()->GetFormsOnlyLinearBonds()
        )
        {
          continue;
        }
        if( itr_edges->GetEdgeData()->GetNumberOfElectrons() == size_t( 4))
        {
          if( m_ChangeChirality)
          {
            rotatable_edges.PushBack( *itr_edges);
          }
        }
        else
        {
          rotatable_edges.PushBack( *itr_edges);
        }
      }
      for( auto itr_edges( rotatable_edges.Begin()), itr_edges_end( rotatable_edges.End()); itr_edges != itr_edges_end; ++itr_edges)
      {
        if
        (
          !parts_to_sample.IsEmpty()
          && ( !parts_to_sample.Contains( itr_edges->GetIndexLow()) || !parts_to_sample.Contains( itr_edges->GetIndexHigh()))
        )
        {
          continue;
        }
        if( m_RandomDihedralChangeProbability && m_SampleDihedralAngles)
        {
          bnd_mutates.PushBack
          (
            util::ShPtr< MutateDihedralsInterface>
            (
              new MutateDihedralBond
              (
                MOLECULE.GetAtomVector()( itr_edges->GetIndexLow()),
                MOLECULE.GetAtomVector()( itr_edges->GetIndexHigh()),
                molecule_graph,
                MOLECULE,
                itr_edges->GetEdgeData()->GetNumberOfElectrons() > size_t( 2)
                || itr_edges->GetEdgeData()->GetConjugation() == ConstitutionalBondTypeData::e_Amide,
                !m_SampleDihedralAngles
              )
            )
          );
          bnd_mutates.LastElement()->GetChooseBestOnlySharedPointer() = mutate_to_best_sp;
        }
      }

      util::ShPtrList< MutateDihedralsInterface> mutate_bond_angles, mutate_bond_angles_refinement;

      size_t degrees_freedom_bond_angles( 0), constant_bond_angles( 0);
      for( auto itr_bas( bond_angle_assignments.Begin()), itr_bas_end( bond_angle_assignments.End()); itr_bas != itr_bas_end; ++itr_bas)
      {
        if( !parts_to_sample.IsEmpty() && !parts_to_sample.Contains( ( *itr_bas)->GetCentralAtomIndex()))
        {
          continue;
        }
        // check for constant bond angle
//        if( m_SampleBondAnglesAndLengths && ( *itr_bas)->GetBondAngleCounts().GetSize() == size_t( 1))
//        {
//          // fix the bond angle if needs fixing
//          if( ( *itr_bas)->GetFreeEnergy( generated_mol) <= 0.0)
//          {
//            BCL_MessageStd( "Existing bond angle was only that has ever been seen for " + util::Format()( ( *itr_bas)->GetCentralAtomIndex()));
//            ++constant_bond_angles;
//            continue;
//          }
//          MutateBondAngles changeit
//          (
//            *itr_bas,
//            generated_mol,
//            chirality_changer,
//            m_ChangeChirality,
//            false
//          );
//          changeit.GetChooseBestOnlySharedPointer() = mutate_to_best_sp;
//          auto mutres( changeit( generated_mol));
//          if( mutres.GetArgument().IsDefined())
//          {
//            generated_mol = *mutres.GetArgument();
//            ++constant_bond_angles;
//            continue;
//          }
//        }
        if( m_Generate3D)
        {
          MutateBondAngles changeit
          (
            *itr_bas,
            generated_mol,
            chirality_changer,
            m_ChangeChirality,
            false
          );
          changeit.GetChooseBestOnlySharedPointer() = mutate_to_best_sp;
          auto mutres( changeit( generated_mol));
          if( mutres.GetArgument().IsDefined())
          {
            generated_mol = *mutres.GetArgument();
          }
        }
        if( m_SampleBondAnglesAndLengths)
        {
          ++degrees_freedom_bond_angles;
        }
        mutate_bond_angles.PushBack
        (
          util::ShPtr< MutateBondAngles>
          (
            new MutateBondAngles
            (
              *itr_bas,
              generated_mol,
              chirality_changer,
              m_ChangeChirality,
              false
            )
          )
        );
        mutate_bond_angles.LastElement()->GetChooseBestOnlySharedPointer() = mutate_to_best_sp;
        mutate_bond_angles_refinement.PushBack
        (
          util::ShPtr< MutateBondAngles>
          (
            new MutateBondAngles
            (
              *itr_bas,
              generated_mol,
              chirality_changer,
              m_ChangeChirality,
              false,
              true
            )
          )
        );
      }

      // create mutate to change all dihedral angles
      MutateMultiFragment full_dihedral_mutate
      (
        dihedral_mutates,
        bnd_mutates,
        m_RandomDihedralChangeProbability,
        all_edges.GetSize()
      );

      util::ShPtr< MutateClashResolver> mcr( new MutateClashResolver( m_MaxNumberClashResolutionCycles));
      mcr->Setup
      (
        MOLECULE,
        util::ShPtrVector< MutateBondAngles>( mutate_bond_angles_refinement.Begin(), mutate_bond_angles_refinement.End()),
        clash_score_sp
      );

      // create mutate for all bond angles
      util::ShPtrList< math::MutateInterface< FragmentComplete> > mutate_bond_angles_math( mutate_bond_angles.Begin(), mutate_bond_angles.End());
      math::MutateCombine< FragmentComplete> bond_angle_mutate
      (
        m_SampleBondAnglesAndLengths
        ? mutate_bond_angles_math
        : util::ShPtrList< math::MutateInterface< FragmentComplete> >(),
        true,
        "bond angle mutation"
      );

      // create mutate interface object from object probability distribution
      const util::ShPtr< math::MutateCombine< FragmentComplete> > mutate
      (
        m_ChangeChirality
        ? new math::MutateCombine< FragmentComplete>
          (
            chirality_changer,
            bond_angle_mutate,
            full_dihedral_mutate,
            *mcr,
            true,
            "Ring-Chain-Angles combined"
          )
        : new math::MutateCombine< FragmentComplete>
          (
            bond_angle_mutate,
            full_dihedral_mutate,
            *mcr,
            true,
            "Ring-Chain-Angles combined"
          )
      );
      BCL_MessageStd
      (
        "Found " + util::Format()( degrees_freedom_chain)
        + " rotatable bonds, " + util::Format()( degrees_freedom_ring)
        + " flexible rings, " + util::Format()( degrees_freedom_bond_angles) + " changeable bond angles; "
        + ( constant_bond_angles ? " and " + util::Format()( constant_bond_angles) + " inflexible bond angles" : std::string())
      );

      generated_mol.GetStoredPropertiesNonConst().Merge( MOLECULE.GetStoredProperties());
      generated_mol.GetStoredPropertiesNonConst().RemoveConformationalDescriptors();

      // start with an initial minimization with no temperature, just to minimize the structure
      FragmentEnsemble sorted_ensemble;

      // monte carlo sampling to apply different mutates to the molecule of interest for conformation sampling
      const size_t n_iterations( m_MaxIterations);
      Sampling
      (
        generated_mol,
        sorted_ensemble,
        mutate,
        fps,
        n_iterations,
        m_MaxClashScore,
        *mutate_to_best_sp,
        *clash_score_sp
      );
      sorted_ensemble.Shuffle();
      sorted_ensemble.Sort( storage::Vector< std::string>::Create( "ConfScore", "ConfOneFourInteractionScore"));

      if( sorted_ensemble.IsEmpty())
      {
        BCL_MessageCrt( "No conformations were sampled, outputting the input conformation");
        if( !m_Generate3D)
        {
          generated_mol.StoreProperty
          (
            "ConfScore",
            storage::Vector< double>( 1, double( 1.0))
          );
          generated_mol.StoreProperty
          (
            "ConfScorePBest",
            storage::Vector< double>( 1, double( 1.0))
          );
          generated_mol.StoreProperty
          (
            "ConfTolerance",
            storage::Vector< double>( 1, m_Tolerance)
          );
          sorted_ensemble.PushBack( generated_mol);
        }
        else
        {
          // if 2D structure is given generate a 3D conformation for the molecule of interest
          FragmentComplete generated_mol_b( MOLECULE);

          // Insert a conformation score of 1 (highest possible = worst score), which will be used if a conformation cannot
          // be generated; otherwise it will be overwritten by the actual score
          generated_mol_b.StoreProperty
          (
            "ConfScore",
            storage::Vector< double>( 1, double( 1.0))
          );
          generated_mol_b.StoreProperty
          (
            "ConfScorePBest",
            storage::Vector< double>( 1, double( 1.0))
          );
          generated_mol_b.StoreProperty
          (
            "ConfTolerance",
            storage::Vector< double>( 1, m_Tolerance)
          );
          generated_mol_b.GetStoredPropertiesNonConst().RemoveConformationalDescriptors();
          sorted_ensemble.PushBack( generated_mol_b);
        }
      }

      // create an ensemble for storing different conformations of a molecule
      FragmentEnsemble top_scoring_ensemble;

      // scorer can be used to store rotamer information, would could be used to train a machine learner to identify better
      // ligand poses
      //FragmentProbabilityScore scorer( updated_bond_mapping);
      // create a conformation set to store generated conformations
      const double conf_score_best( sorted_ensemble.GetMolecules().FirstElement().GetMDLPropertyAsVector( "ConfScore")( 0));
      storage::Vector< double> score_differences;
      score_differences.AllocateMemory( m_NumberConformations);
      bool do_simple_filtering( !m_Cluster || !m_ConformationComparison.IsDefined());
      double effective_tolerance( m_Tolerance);
      size_t n_preclustering_removed( 0);
      if( do_simple_filtering)
      {
        #ifdef BCL_PROFILE_SampleConformations
        static util::Stopwatch s_cset( "ConformationComparison", util::Time( 1, 0), util::Message::e_Standard, true, false);
        s_cset.Start();
        #endif
        size_t number_conformations( 0);
        util::Implementation< ConformationComparisonInterface> precomparer;
        // symmetry rmsd can be very expensive to compute and if the user is giving a broad tolerance it is faster
        // to test with RMSD first to discard conformers that readily align to other molecules
        if( m_ConformationComparison.IsDefined())
        {
          if( m_ConformationComparison->GetAlias() == "SymmetryRMSD")
          {
            precomparer = util::Implementation< ConformationComparisonInterface>( "RMSD");
          }
          else if( m_ConformationComparison->GetAlias() == "SymmetryRealSpaceRMSD")
          {
            precomparer = util::Implementation< ConformationComparisonInterface>( "RealSpaceRMSD");
          }
          if( precomparer.IsDefined())
          {
            m_ConformationComparison->operator ()( generated_mol, generated_mol);
            util::SiPtr< const ConformationComparisonBySymmetryRmsd> symmetryrmsd( &*m_ConformationComparison);
            if( symmetryrmsd->GetNumberIsomorphisms() <= size_t( 2))
            {
              precomparer.Reset();
            }
          }
        }
        else
        {
          size_t n_ensemble( sorted_ensemble.GetSize());
          if( n_ensemble > m_NumberConformations)
          {
            sorted_ensemble.Shuffle();
            auto itr_s( sorted_ensemble.Begin());
            std::advance( itr_s, m_NumberConformations);
            sorted_ensemble.GetMolecules().Remove( itr_s, sorted_ensemble.End());
            sorted_ensemble.Sort( storage::Vector< std::string>::Create( "ConfScore", "ConfOneFourInteractionScore"));
          }
        }
        for
        (
          auto itr( sorted_ensemble.Begin()), itr_end( sorted_ensemble.End());
          itr != itr_end;
        )
        {
          bool conformation_found( false);
          double cscore( itr->GetMDLPropertyAsVector( "ConfScore")( 0));
          // check with pre-comparer
          if( precomparer.IsDefined())
          {
            for
            (
              auto itr_accepted( top_scoring_ensemble.GetMolecules().ReverseBegin()),
                   itr_accepted_end( top_scoring_ensemble.GetMolecules().ReverseEnd());
              itr_accepted != itr_accepted_end;
              ++itr_accepted
            )
            {
              if( m_Cluster && cscore != itr_accepted->GetMDLPropertyAsVector( "ConfScore")( 0))
              {
                break;
              }
              if( ( *precomparer)( *itr_accepted, *itr) <= effective_tolerance)
              {
                conformation_found = true;
                itr_accepted->GetStoredPropertiesNonConst().SetMDLProperty
                (
                  "BCLConfClusterSize",
                  util::Format()( itr_accepted->GetStoredProperties().GetMDLPropertyAsVector( "BCLConfClusterSize")( 0) + size_t( 1))
                );
                ++n_preclustering_removed;
                break;
              }
            }
          }
          if( !conformation_found && m_ConformationComparison.IsDefined())
          {
            // We try going from the end of the list first, because the list is sorted by ConfScore, and molecules with
            // similar conf scores are more likely to be similar
            for
            (
              auto itr_accepted( top_scoring_ensemble.GetMolecules().ReverseBegin()),
                   itr_accepted_end( top_scoring_ensemble.GetMolecules().ReverseEnd());
              itr_accepted != itr_accepted_end;
              ++itr_accepted
            )
            {
              if( m_Cluster && cscore != itr_accepted->GetMDLPropertyAsVector( "ConfScore")( 0))
              {
                break;
              }
              if( ( *m_ConformationComparison)( *itr_accepted, *itr) <= effective_tolerance)
              {
                itr_accepted->GetStoredPropertiesNonConst().SetMDLProperty
                (
                  "BCLConfClusterSize",
                  util::Format()( itr_accepted->GetStoredProperties().GetMDLPropertyAsVector( "BCLConfClusterSize")( 0) + size_t( 1))
                );
                conformation_found = true;
                ++n_preclustering_removed;
                break;
              }
            }
          }
          if( !conformation_found)
          {
            ++number_conformations;
            itr->SetName( MOLECULE.GetName());
            itr->GetStoredPropertiesNonConst().Merge( generated_mol.GetStoredProperties());
            //  scorer.AddRotamerMiscPropertyInfo( **itr, 200);
            auto itr_prev( itr);
            ++itr;
            top_scoring_ensemble.GetMolecules().InternalData().splice
            (
              top_scoring_ensemble.GetMolecules().InternalData().end(),
              sorted_ensemble.GetMolecules().InternalData(),
              itr_prev
            );
            top_scoring_ensemble.GetMolecules().LastElement().GetStoredPropertiesNonConst().SetMDLProperty
            (
              "BCLConfClusterSize",
              std::string( "1")
            );
            double conf_score( top_scoring_ensemble.GetMolecules().LastElement().GetMDLPropertyAsVector( "ConfScore")( 0));
            score_differences.PushBack( GetRelativeProbabilityOfConformerBeingClosestToNative( conf_score_best, conf_score));
            if( ( !m_Cluster || !m_ConformationComparison.IsDefined()) && score_differences.GetSize() == m_NumberConformations)
            {
              break;
            }
          }
          else
          {
            itr->GetChangeSignal().Emit( *itr);
            ++itr;
          }
        }
        double score_sum( math::Statistics::Sum( score_differences.Begin(), score_differences.End(), 0.0));
        auto itr_score( score_differences.Begin());
        for
        (
          auto itr( top_scoring_ensemble.Begin()), itr_end( top_scoring_ensemble.End());
          itr != itr_end;
          ++itr, ++itr_score
        )
        {
          itr->StoreProperty( "ConfScorePBest", linal::Vector< float>( size_t( 1), float( *itr_score / score_sum)));
        }
        #ifdef BCL_PROFILE_SampleConformations
        s_cset.Stop();
        #endif
      }
      if( m_Cluster)
      {
        if( do_simple_filtering)
        {
          sorted_ensemble.GetMolecules().InternalData().swap( top_scoring_ensemble.GetMolecules().InternalData());
          top_scoring_ensemble.GetMolecules().Reset();
        }
        #ifdef BCL_PROFILE_SampleConformations
        static util::Stopwatch s_clust( "Clustering", util::Time( 1, 0), util::Message::e_Standard, true, false);
        s_clust.Start();
        #endif
        size_t n_conformations_initial( sorted_ensemble.GetSize());
        if( n_conformations_initial > m_NumberConformations)
        {
          // pre-cluster setup. First, for every molecule that has the same score as any other molecule, filter by SymmetryRMSD
          storage::Vector< size_t> is_covered( n_conformations_initial, size_t( 0));
          linal::Vector< double> tmp_vals( n_conformations_initial, double( 0.0));
          linal::Vector< double> tmp_vals_srt( n_conformations_initial, double( 0.0));
          storage::Vector< storage::Pair< double, double> > tmp_vals_srt2( n_conformations_initial);
          linal::Vector< double> confps( n_conformations_initial, double( 0.0));

          const double start_end_ratiob( 0.9425 * math::Pow( double( n_conformations_initial), -0.231446));
          const double end_valueb( 2.0 * start_end_ratiob / ( double( n_conformations_initial) * ( 1.0 + start_end_ratiob)));
          const double start_valueb( end_valueb / start_end_ratiob);
          const double slopeb( ( end_valueb - start_valueb) / ( n_conformations_initial - 1.0));
          {
            size_t i( 0);
            for( auto itr_mol( sorted_ensemble.Begin()), itr_mol_end( sorted_ensemble.End()); itr_mol != itr_mol_end; ++itr_mol, ++i)
            {
              //double conf_score( itr_mol->GetMDLPropertyAsVector( "ConfScore")( 0));
              confps( i) = start_valueb + i * slopeb;
              auto cluster_sizes( itr_mol->GetStoredProperties().GetMDLPropertyAsVector( "BCLConfClusterSize"));
              if( cluster_sizes.GetSize())
              {
                confps( i) *= cluster_sizes( 0);
              }
            }
            confps.SetToSum( 1.0);
          }

          storage::Vector< double> tolerances;
          storage::Vector< double> confscores;
          storage::Vector< double> confscorepsum;
          linal::Vector< double> min_rmsd_any( n_conformations_initial, double( 10000.0));
          tolerances.AllocateMemory( m_NumberConformations);
          confscores.AllocateMemory( m_NumberConformations);
          auto itr_covered( is_covered.Begin());
          size_t n_outer( 0), n_saved( 0), n_uncovered( n_conformations_initial);
          size_t target_n_saved( m_NumberConformations);
          auto itr( sorted_ensemble.Begin()), itr_end( sorted_ensemble.End());
          double p_uncovered( double( target_n_saved) / double( m_NumberConformations));
          for( ; itr != itr_end && n_saved < target_n_saved; ++itr, ++itr_covered, ++n_outer)
          {
            if( *itr_covered)
            {
              continue;
            }
            tmp_vals.Shrink( n_uncovered);
            tmp_vals_srt.Shrink( n_uncovered);
            tmp_vals_srt2.Resize( n_uncovered);
            size_t n_b( 0), n_c( 0);
            // We try going from the end of the list first, because the list is sorted by ConfScore, and molecules with
            // similar conf scores are more likely to be similar
            for
            (
              auto itr_b( sorted_ensemble.Begin()),
                   itr_b_end( sorted_ensemble.End());
              itr_b != itr_b_end;
              ++itr_b, ++n_b
            )
            {
              if( !is_covered( n_b))
              {
                tmp_vals_srt2( n_c).First() = tmp_vals( n_c) = ( *m_ConformationComparison)( *itr_b, *itr);
                tmp_vals_srt2( n_c).Second() = confps( n_b);
                min_rmsd_any( n_b) = std::min( min_rmsd_any( n_b), tmp_vals( n_c));
                ++n_c;
              }
            }
            std::stable_sort( tmp_vals_srt2.Begin(), tmp_vals_srt2.End());
            double p_desired( p_uncovered / double( target_n_saved - n_saved));
            size_t p_pos( 0);
            double new_tolerance( 0.0);
            while( p_desired > 0.0 && p_pos < tmp_vals_srt2.GetSize())
            {
              p_desired -= tmp_vals_srt2( p_pos).Second();
              new_tolerance = tmp_vals_srt2( p_pos).First();
              ++p_pos;
            }
            new_tolerance = std::max( new_tolerance, m_Tolerance);
            tolerances.PushBack( new_tolerance);
            double sumconfps( 0.0);
            ++n_saved;
            n_c = 0;
            for( size_t i_b( 0), i_c( 0), n_uncovered_prev( n_uncovered); i_b < n_uncovered_prev; ++i_b, ++i_c)
            {
              while( is_covered( i_c))
              {
                ++i_c;
              }
              if( tmp_vals( i_b) <= new_tolerance)
              {
                --n_uncovered;
                is_covered( i_c) += 1;
                sumconfps += confps( i_c);
              }
            }
            p_uncovered -= sumconfps;
            confscorepsum.PushBack( sumconfps);
            top_scoring_ensemble.PushBack( *itr);
            top_scoring_ensemble.GetMolecules().LastElement().SetName( MOLECULE.GetName());
            top_scoring_ensemble.GetMolecules().LastElement().GetStoredPropertiesNonConst().Merge( generated_mol.GetStoredProperties());
            confscores.PushBack( util::ConvertStringToNumericalValue< double>( top_scoring_ensemble.GetMolecules().LastElement().GetMDLProperty( "ConfScore")));
          }

          bool did_diversity
          (
            n_uncovered
            || ( n_saved < m_NumberConformations && tolerances.Find( m_Tolerance) < tolerances.GetSize())
          );
          if( did_diversity)
          {
            if( tolerances.Find( m_Tolerance) < tolerances.GetSize())
            {
              itr = sorted_ensemble.Begin();
              min_rmsd_any = 10000.0;
              n_uncovered = n_conformations_initial;
              itr_covered = is_covered.Begin();
              is_covered.SetAllElements( 0);
              n_outer = size_t( 0);
              for
              (
                auto itr_top( top_scoring_ensemble.Begin()), itr_top_end( top_scoring_ensemble.End());
                itr_top != itr_top_end;
                ++itr_top
              )
              {
                size_t n_inner( 0);
                for
                (
                  auto itr_so( sorted_ensemble.Begin()), itr_so_end( sorted_ensemble.End());
                  itr_so != itr_so_end;
                  ++itr_so, ++n_inner
                )
                {
                  if( !min_rmsd_any( n_inner))
                  {
                    continue;
                  }
                  min_rmsd_any( n_inner)
                    = std::min( min_rmsd_any( n_inner), ( *m_ConformationComparison)( *itr_top, *itr_so));
                  if( !min_rmsd_any( n_inner))
                  {
                    is_covered( n_inner) = 1;
                    --n_uncovered;
                    break;
                  }
                }
              }
            }
            for( ; n_saved < m_NumberConformations; ++n_saved)
            {
              // find the max rmsd of all uncovered conformations
              size_t best_index( n_outer);
              auto itr_best( itr);
              bool did_something( true);
              {
                double max_rmsd( 0.0);
                auto itr_nxt( itr);
                auto itr_cvr( itr_covered);
                for( size_t local_index( n_outer); local_index < n_conformations_initial; ++local_index, ++itr_nxt, ++itr_cvr)
                {
                  if( *itr_cvr)
                  {
                    continue;
                  }
                  if( min_rmsd_any( local_index) > max_rmsd)
                  {
                    max_rmsd = min_rmsd_any( local_index);
                    itr_best = itr_nxt;
                    best_index = local_index;
                  }
                }
                if( m_Tolerance && max_rmsd > m_Tolerance)
                {
                  did_something = false;
                  break;
                }
                tolerances.PushBack( max_rmsd);
                top_scoring_ensemble.PushBack( *itr_best);
                top_scoring_ensemble.GetMolecules().LastElement().SetName( MOLECULE.GetName());
                top_scoring_ensemble.GetMolecules().LastElement().GetStoredPropertiesNonConst().Merge( generated_mol.GetStoredProperties());
                confscores.PushBack( util::ConvertStringToNumericalValue< double>( top_scoring_ensemble.GetMolecules().LastElement().GetMDLProperty( "ConfScore")));
                is_covered( best_index) = 1;
              }
              if( n_saved + 1 < m_NumberConformations && did_something)
              {
                auto itr_nxt( itr);
                auto itr_cvr( itr_covered);
                for( size_t local_index( n_outer); local_index < n_conformations_initial; ++local_index, ++itr_nxt, ++itr_cvr)
                {
                  if( *itr_cvr)
                  {
                    continue;
                  }
                  min_rmsd_any( local_index) = std::min( ( *m_ConformationComparison)( *itr_best, *itr_nxt), min_rmsd_any( local_index));
                }
              }
            }
            confscores.Sort( std::less< double>());
            top_scoring_ensemble.Sort( storage::Vector< std::string>::Create( "ConfScore", "ConfOneFourInteractionScore", "BCLConfClusterSize"));
          }
          size_t actual_n_conf( confscores.GetSize());
          const double start_end_ratio( 0.9425 * math::Pow( double( actual_n_conf), -0.231446));
          const double end_value( 2.0 * start_end_ratio / ( double( actual_n_conf) * ( 1.0 + start_end_ratio)));
          const double start_value( end_value / start_end_ratio);
          const double slope( ( end_value - start_value) / ( actual_n_conf - 1.0));
          size_t nth( 0);
          for
          (
            auto itr( top_scoring_ensemble.Begin()), itr_end( top_scoring_ensemble.End());
            itr != itr_end;
            ++itr, ++nth
          )
          {
            itr->StoreProperty( "ConfScorePBest", linal::Vector< float>( size_t( 1), float( start_value + nth * slope)));
            itr->StoreProperty( "ConfTolerance", linal::Vector< float>( size_t( 1), float( tolerances( nth))));
          }
        }
        else
        {
          top_scoring_ensemble.GetMolecules().InternalData().swap( sorted_ensemble.GetMolecules().InternalData());
          for
          (
            auto itr( top_scoring_ensemble.GetMolecules().Begin()), itr_end( top_scoring_ensemble.GetMolecules().End());
            itr != itr_end;
            ++itr
          )
          {
            itr->SetName( MOLECULE.GetName());
            itr->GetStoredPropertiesNonConst().Merge( generated_mol.GetStoredProperties());
            double conf_score( itr->GetMDLPropertyAsVector( "ConfScore")( 0));
            score_differences.PushBack( GetRelativeProbabilityOfConformerBeingClosestToNative( conf_score_best, conf_score));
          }
          double score_sum( math::Statistics::Sum( score_differences.Begin(), score_differences.End(), 0.0));
          auto itr_score( score_differences.Begin());
          for
          (
            auto itr( top_scoring_ensemble.Begin()), itr_end( top_scoring_ensemble.End());
            itr != itr_end;
            ++itr, ++itr_score
          )
          {
            itr->StoreProperty( "ConfScorePBest", linal::Vector< float>( size_t( 1), float( *itr_score / score_sum)));
            itr->StoreProperty( "ConfTolerance", linal::Vector< float>( size_t( 1), float( m_Tolerance)));
          }
        }
        #ifdef BCL_PROFILE_SampleConformations
        s_clust.Stop();
        #endif
      }
      FragmentSplitRigid rigid( 3);
      ConformationGraphConverter::t_AtomGraph mol_graph( ConformationGraphConverter::CreateGraphWithAtoms( generated_mol));
      auto rigid_components( rigid.GetComponentVertices( generated_mol, mol_graph));
      storage::Vector< size_t> largest_rigid_fragment;
      for( auto itr_rigid( rigid_components.Begin()), itr_rigid_end( rigid_components.End()); itr_rigid != itr_rigid_end; ++itr_rigid)
      {
        if( itr_rigid->GetSize() > largest_rigid_fragment.GetSize())
        {
          largest_rigid_fragment = *itr_rigid;
        }
      }
      if( ( !m_SampleDihedralAngles && !m_SampleRingConformations) || largest_rigid_fragment.GetSize() < size_t( 3))
      {
        largest_rigid_fragment = storage::CreateIndexVector( generated_mol.GetSize());
      }
      if( largest_rigid_fragment.GetSize() >= size_t( 3))
      {
        bool mol_was_okay( !HasUndefinedCoordinates( MOLECULE));
        util::SiPtrVector< const linal::Vector3D> coords( ( mol_was_okay ? MOLECULE : generated_mol).GetAtomCoordinates());
        coords.Reorder( largest_rigid_fragment);
        for( auto itr( top_scoring_ensemble.Begin()), itr_end( top_scoring_ensemble.End()); itr != itr_end; ++itr)
        {
          auto rotamer_coords( itr->GetAtomCoordinates());
          rotamer_coords.Reorder( largest_rigid_fragment);
          auto transform( quality::RMSD::SuperimposeCoordinates( coords, rotamer_coords));
          itr->Transform( transform);
        }
      }
      if( !m_OutputIndividualScoresFilename.empty())
      {
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputIndividualScoresFilename);
        for( auto itr( top_scoring_ensemble.Begin()), itr_end( top_scoring_ensemble.End()); itr != itr_end; ++itr)
        {
          auto ens( fps.GetScoreComponents( *itr));
          ens.WriteMDL( output);
        }
        io::File::CloseClearFStream( output);
      }
      return storage::Pair< FragmentEnsemble, FragmentEnsemble>( top_scoring_ensemble, *tmp_fragments);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief applies mutates in monte carlo fasion to generate conformers for molecule
    //! @param STARTING_ARGUMENT the argument that will be the first for sampling
    //! @param CONFORMATION_SET conformation set to store results that are produced on application of mutates
    //! @param MUTATE mutates that need to be applied to the molecule of interest
    //! @param SCORE_FUNCTION object that scores a the mutate result
    //! @param MAX_ITERATIONS maximum number of iterations during monte carlo sampling
    //! @param TEMPERATURE temperature to make decision of accepting a move
    //! @param MAX_RESULTS desired number of new mutate results that are accepted/improved
    void SampleConformations::Sampling
    (
      const FragmentComplete &STARTING_ARGUMENT,
      FragmentEnsemble &CONFORMATION_SET,
      const util::ShPtr< math::MutateInterface< FragmentComplete> > &MUTATE,
      const FragmentProbabilityScore &SCORE,
      const size_t &MAX_ITERATIONS,
      const double &CLASH_TOLERANCE,
      bool &SHARED_CHOOSE_ONLY_BEST,
      AtomClashScore &SHARED_CLASH_SCORE
    ) const
    {
      #ifdef BCL_PROFILE_SampleConformations
      static util::Stopwatch s_cset( "Sampling Total", util::Time( 1000, 0), util::Message::e_Standard, true, false);
      s_cset.Start();
      #endif
      // opti terminate object
      opti::CriterionCombine< FragmentComplete, double> sp_terminate;
      math::RunningAverage< double> aveclashscore;

      if( MAX_ITERATIONS == size_t( 0))
      {
        // if only scores are desired
        CONFORMATION_SET.PushBack( STARTING_ARGUMENT);
        CONFORMATION_SET.GetMolecules().Last()->StoreProperty( "ConfScore", storage::Vector< double>( 1, SCORE( STARTING_ARGUMENT)));
        CONFORMATION_SET.GetMolecules().Last()->StoreProperty( "ConfOneFourInteractionScore", storage::Vector< double>( 1, SCORE.Get14InteractionScore( STARTING_ARGUMENT)));
        CONFORMATION_SET.GetMolecules().Last()->StoreProperty( "ConfClashScore", storage::Vector< double>( 1, SHARED_CLASH_SCORE( STARTING_ARGUMENT)));
      }
      else
      {
        size_t n_clashed( 0);

        size_t n_clashed_in_a_row( 0);
        FragmentComplete clashed_but_best;
        util::SiPtr< const FragmentComplete> to_mutate( STARTING_ARGUMENT);
        for( size_t iteration( 0); iteration < MAX_ITERATIONS && n_clashed_in_a_row < 100; ++iteration)
        {
          if( iteration < size_t( 2) && n_clashed_in_a_row < 2)
          {
            SHARED_CHOOSE_ONLY_BEST = true;
          }
          if( !CONFORMATION_SET.IsEmpty())
          {
            to_mutate =
              &*random::GetGlobalRandom().Iterator( CONFORMATION_SET.Begin(), CONFORMATION_SET.End(), CONFORMATION_SET.GetSize());
          }
          math::MutateResult< FragmentComplete> result( ( *MUTATE)( *to_mutate));
          SHARED_CHOOSE_ONLY_BEST = false;
          if( !result.GetArgument().IsDefined())
          {
            ++n_clashed_in_a_row;
            continue;
          }
          const FragmentComplete &molecule_conf_gen( *result.GetArgument());

          double main_clash_score( SHARED_CLASH_SCORE( molecule_conf_gen));
          if( main_clash_score > CLASH_TOLERANCE)
          {
            --iteration;
            ++n_clashed;
            ++n_clashed_in_a_row;
            continue;
          }

          CONFORMATION_SET.PushBack( molecule_conf_gen);
          CONFORMATION_SET.GetMolecules().Last()->StoreProperty( "ConfScore", storage::Vector< double>( 1, SCORE( molecule_conf_gen) + 100.0 * main_clash_score));
          CONFORMATION_SET.GetMolecules().Last()->StoreProperty( "ConfOneFourInteractionScore", storage::Vector< double>( 1, SCORE.Get14InteractionScore( molecule_conf_gen)));
          CONFORMATION_SET.GetMolecules().Last()->StoreProperty( "ConfClashScore", storage::Vector< double>( 1, main_clash_score));
          n_clashed_in_a_row = 0;
        }
        if( n_clashed > MAX_ITERATIONS / 2)
        {
          BCL_MessageStd( " # rejected due to severe clashes: " + util::Format()( n_clashed));
        }
      }
      #ifdef BCL_PROFILE_SampleConformations
      s_cset.Stop();
      #endif
    }

    //! @brief determine if the coordinates of atoms of a molecule are defined
    //! @param MOLECULE molecule to check for
    //! @result true if coordinates are defined, false otherwise
    bool SampleConformations::HasUndefinedCoordinates( const ConformationInterface &MOLECULE) const
    {
      bool has_undefined_coords( false);
      size_t zero_coord_atoms( 0);
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( MOLECULE.GetAtomsIterator().Begin());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        linal::Vector3D atom_coords( itr_atoms->GetPosition());
        if( math::EqualWithinTolerance( atom_coords.X(), 0) && math::EqualWithinTolerance( atom_coords.Y(), 0))
        {
          if( math::EqualWithinTolerance( atom_coords.Z(), 0))
          {
            ++zero_coord_atoms;
          }
        }
        if( zero_coord_atoms > 1)
        {
          return true;
        }

        // check that coordinates are defined
        if( !itr_atoms->GetPosition().IsDefined())
        {
          // undefined coordinates; geometry is bad
          BCL_MessageVrb
          (
            "Bad Geometry! Atom with undefined position (type: " + itr_atoms->GetAtomType().GetName() + ")"
          );
          has_undefined_coords = true;
          break;
        }

        // for the 3D check, skip atoms with a defined z coordinate
        if( !math::EqualWithinTolerance( itr_atoms->GetPosition().Z(), 0))
        {
          continue;
        }

        const storage::Vector< BondConformational> &connected_atoms( itr_atoms->GetBonds());

        // atoms with 0-2 bonds could have all neighboring atoms with a 0 z-coordinate
        if( connected_atoms.GetSize() < 3)
        {
          continue;
        }

        // atoms with 3 bonds could have all neighboring atoms with a 0 z-coordinate if the geometry is not SP3
        if
        (
          connected_atoms.GetSize() == 3
          && itr_atoms->GetAtomType()->GetHybridOrbitalType() != GetHybridOrbitalTypes().e_SP3
        )
        {
          continue;
        }

        bool had_3d_neighbor( false);
        // check that at least one neighbor has a non-zero z-coordinate
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()), itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          if( !math::EqualWithinTolerance( itr_connected->GetTargetAtom().GetPosition().Z(), 0))
          {
            had_3d_neighbor = true;
            break;
          }
        }
        if( !had_3d_neighbor)
        {
          const std::string &name( itr_atoms->GetAtomType().GetName());
          if( connected_atoms.GetSize() == 3)
          {
            BCL_MessageVrb
            (
              "Bad Geometry! SP3 atom with (type: " + name + "); all neighbors have 0 z-coordinate"
            );
          }
          else
          {
            BCL_MessageVrb
            (
              "Bad Geometry! Atom with 4 bonds (type: " + name + "); all neighbors have 0 z-coordinate"
            );
          }
          has_undefined_coords = true;
        }
      }
      return has_undefined_coords;
    }

    //! @brief determine if the coordinates of atoms of a molecule are defined
    //! @param MOLECULE molecule whose rings are desired
    //! @param MOLECULE_RINGS rings contained in molecule
    //! @param RINGS small fragment isomorphisms of rings that need to be checked for wholeness
    //! @result true if coordinates are defined, false otherwise
    bool SampleConformations::CheckForWholeRings
    (
      const FragmentComplete &MOLECULE,
      const util::ShPtr< FragmentEnsemble> &FRAGMENTS
    ) const
    {
      ConformationGraphConverter graph_maker
      (
        ConformationGraphConverter::e_AtomType,
        ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
      );

      graph::SubgraphIsomorphism< size_t, size_t> csi_substructure;

      storage::Vector< storage::List< storage::Pair< FragmentComplete, graph::ConstGraph< size_t, size_t> > > > rings_in_fragments( MOLECULE.GetSize() + 1);
      ConstitutionSet unique_rings;
      // store fragments in a temporary vector
      for
      (
        storage::List< FragmentComplete>::const_iterator
          itr_fragments( FRAGMENTS->Begin()), itr_fragments_end( FRAGMENTS->End());
        itr_fragments != itr_fragments_end;
        ++itr_fragments
      )
      {
        FragmentEnsemble fragment_rings( FragmentSplitRings( true, 3)( *itr_fragments));
        for
        (
          storage::List< FragmentComplete>::const_iterator
            itr_frings( fragment_rings.Begin()), itr_frings_end( fragment_rings.End());
          itr_frings != itr_frings_end;
          ++itr_frings
        )
        {
          if( unique_rings.Insert( FragmentConstitutionShared( *itr_frings)).second)
          {
            rings_in_fragments( itr_frings->GetNumberAtoms()).PushBack
            (
              storage::Pair< FragmentComplete, graph::ConstGraph< size_t, size_t> >( *itr_frings, graph_maker( *itr_frings))
            );
          }
        }
      }

      // Search for rings within the molecule
      FragmentEnsemble molecule_rings( FragmentSplitRings( true, 3)( MOLECULE));
      storage::Vector< storage::List< storage::Pair< FragmentComplete, graph::ConstGraph< size_t, size_t> > > > rings_in_molecule( MOLECULE.GetSize() + 1);
      ConstitutionSet unique_molecule_rings;

      // store fragments in a temporary vector
      for
      (
        storage::List< FragmentComplete>::const_iterator
        itr_molrings( molecule_rings.Begin()), itr_molrings_end( molecule_rings.End());
          itr_molrings != itr_molrings_end;
        ++itr_molrings
      )
      {
        if( unique_molecule_rings.Insert( FragmentConstitutionShared( *itr_molrings)).second)
        {
          rings_in_molecule( itr_molrings->GetNumberAtoms()).PushBack
          (
            storage::Pair< FragmentComplete, graph::ConstGraph< size_t, size_t> >( *itr_molrings, graph_maker( *itr_molrings))
          );
        }
      }

      bool all_rings_found( true);
      for
      (
        storage::Vector< storage::List< storage::Pair< FragmentComplete, graph::ConstGraph< size_t, size_t> > > >::const_reverse_iterator
          itr_molrings( rings_in_molecule.ReverseBegin()), itr_molrings_end( rings_in_molecule.ReverseEnd());
          itr_molrings != itr_molrings_end;
        ++itr_molrings
      )
      {
        for
        (
          storage::List< storage::Pair< FragmentComplete, graph::ConstGraph< size_t, size_t> > >::const_iterator
            itr_molring( itr_molrings->Begin()), itr_molring_end( itr_molrings->End());
          itr_molring != itr_molring_end;
          ++itr_molring
        )
        {
          size_t mol_ring_size( itr_molring->First().GetNumberAtoms());
          csi_substructure.SetGraph( itr_molring->Second());
          bool ring_found( false);
          for
          (
            storage::Vector< storage::List< storage::Pair< FragmentComplete, graph::ConstGraph< size_t, size_t> > > >::const_reverse_iterator
              itr_fragrings( rings_in_fragments.ReverseBegin()), itr_fragrings_end( rings_in_fragments.ReverseEnd());
            itr_fragrings != itr_fragrings_end;
            ++itr_fragrings
          )
          {
            for
            (
              storage::List< storage::Pair< FragmentComplete, graph::ConstGraph< size_t, size_t> > >::const_iterator
                itr_fragring( itr_fragrings->Begin()), itr_fragring_end( itr_fragrings->End());
              itr_fragring != itr_fragring_end;
              ++itr_fragring
            )
            {
              size_t frag_ring_size( itr_fragring->First().GetNumberAtoms());
              if( mol_ring_size != frag_ring_size)
              {
                continue;
              }
              csi_substructure.SetSubgraph( itr_fragring->Second());
              if( csi_substructure.FindIsomorphism())
              {
                ring_found = true;
                break;
              }
            }
            if( ring_found)
            {
              break;
            }
          }
          if( !ring_found)
          {
            return false;
          }
        }
      }
      return all_rings_found;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SampleConformations::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "sample conformations of molecules.\n"
        "To keep part of the molecule rigid, set an MDL property SampleByParts with the MDL (0-indexed) atom indices that can "
        "be moved, e.g.\n"
        "> <SampleByParts>\n"
        "0 4 6 8 10\n"
        "This would indicate that only dihedrals containing exclusively of atoms 1, 5, 7, 9, and 11 in the "
        "SD file can be rotated. "
      );

      parameters.AddInitializer
      (
        "rotamer_library",
        "flag for rotamer library to use. Crystallographic open database is used by default",
        io::Serialization::GetAgent( &m_SearchFragmentLibrary->GetRotamerLibrary()),
        "File"
      );

      parameters.AddInitializer
      (
        "conformation_comparer",
        "method to compare conformers with",
        io::Serialization::GetAgent( &m_ConformationComparison),
        "SymmetryRMSD"
      );

      parameters.AddInitializer
      (
        "generate_3D",
        "true to completely discard the input conformation and rebuild a 3D conformer from scratch. ",
        io::Serialization::GetAgent( &m_Generate3D),
        "false"
      );

      parameters.AddInitializer
      (
        "tolerance",
        "amount of tolerance allowed between two conformers",
        io::Serialization::GetAgent( &m_Tolerance),
        "0.25"
      );
      parameters.AddInitializer
      (
        "max_cr_iterations",
        "max number of clash resolution iterations performed (as fraction of # dihedrals & bond angles)."
        "More iterations will produced less clashed molecules but may "
        "take much longer. Often it is faster to just make more molecules and discard the excessively "
        "clashed ones if they can't be resolved in a small fraction of moves (e.g. 10% of the total dihedrals + bond angles",
        io::Serialization::GetAgent( &m_MaxNumberClashResolutionCycles),
        "0.1"
      );
      parameters.AddInitializer
      (
        "max_avg_clash",
        "max average clash allowed for conformations =sum(vdw_overlap for all atoms in molecule) / atoms in molecule."
        " larger values will yield potentially more clashed conformers for some molecules",
        io::Serialization::GetAgent( &m_MaxClashScore),
        "0.1"
      );
      parameters.AddInitializer
      (
        "max_iterations",
        "number of iterations",
        io::Serialization::GetAgent( &m_MaxIterations),
        "1000"
      );
      parameters.AddInitializer
      (
        "change_chirality",
        "change chirality during conformation sampling",
        io::Serialization::GetAgent( &m_ChangeChirality),
        "false"
      );
      parameters.AddInitializer
      (
        "max_conformations",
        "maximum number conformations to generate per molecule",
        io::Serialization::GetAgent( &m_NumberConformations),
        "200"
      );
      parameters.AddInitializer
      (
        "cluster",
        "Set to true to select up to max_conformations so as to maximize coverage of conformational space by performing "
        "greedy clustering on all returned conformations. This option should only be turned off when generating boltzmann-like "
        "ensembles of conformers typically. Otherwise, better speed and rmsd recovery will result just by decreasing # of "
        "iterations",
        io::Serialization::GetAgent( &m_Cluster),
        "True"
      );
      parameters.AddInitializer
      (
        "relative_random_dihedral_change_weight",
        "Weight for random dihededral changes relative to fragment based changes. This is desireable for exhaustive sampling "
        "or when parts of a fragment are very poorly represented in the fragment library",
        io::Serialization::GetAgent( &m_RandomDihedralChangeProbability),
        "0.0"
      );

      parameters.AddInitializer
      (
        "sample_by_parts_mdl",
        "MDL property from which to extract sample by parts atom indices from an SDF",
        io::Serialization::GetAgent( &m_MDL),
        "SampleByParts"
      );

      parameters.AddInitializer
      (
        "sample_by_parts_indices",
        "sample dihedral angles containing these atom indices (0-indexed) during conformer generation",
        io::Serialization::GetAgent( &m_IndicesToSampleString),
        ""
      );

      parameters.AddInitializer
      (
        "sample_by_parts_fragments",
        "allow atoms that match these fragments in substructure comparisons to be "
        "sampled",
        io::Serialization::GetAgent( &m_ReferenceFragments),
        ""
      );

//      parameters.AddInitializer
//      (
//        "sample_by_parts_mode",
//        "how to handle atom indices specified by 'sample_by_parts' routines; \n"
//        "IgnoreAll - ignore all settings, perform conformer sampling on all dihedrals; \n"
//        "IgnoreOnlyMDL - ignore all MDL properties in SDF and only use command-line settings; \n"
//        "UseOnlyMDL - ignore command-line settings; use only MDL; \n"
//        "UseAll - combine atom indices in MDL with additional indices specified on command-line",
//        io::Serialization::GetAgent( m_SampleByPartsHandler),
//        SampleByPartsOptions::e_UseAll
//      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool SampleConformations::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }
      if
      (
        !m_SearchFragmentLibrary->GetRotamerLibrary().IsDefined()
        || !m_SearchFragmentLibrary->GetRotamerLibrary()->IsDefined()
      )
      {
        ERR_STREAM << "Rotamer library was undefined";
        return false;
      }
      if( m_Cluster && m_MaxIterations < 2 * m_NumberConformations)
      {
        m_Cluster = false;
        BCL_MessageStd( "Clustering disabled as it is futile when making fewer than twice the desired number of molecules");
      }
      if( !m_Tolerance && !m_Cluster)
      {
        m_ConformationComparison.Reset();
      }
      return true;
    }

    //! @brief get the relative probability that a conformation was the closest-to-native
    //!        estimated from the best conf score and its conf score
    //! @param BEST_CONF_SCORE best conf-score observed for this molecule
    //! @param CONF_SCORE actual conf score for the molecule
    double SampleConformations::GetRelativeProbabilityOfConformerBeingClosestToNative( const double &BEST_CONF_SCORE, const double &CONF_SCORE)
    {
      return math::Pow( CONF_SCORE - BEST_CONF_SCORE + 0.002, -0.608);
    }
  } // namespace chemistry
} // namespace bcl
