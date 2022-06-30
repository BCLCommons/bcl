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

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_align_to_scaffold.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_stochastic_pose_optimizer.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_score_function_generic.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_combine.h"
#include "descriptor/bcl_descriptor_molecule_similarity.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_const_function.h"
#include "math/bcl_math_mutate_combine.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "math/bcl_math_mutate_repeat.h"
#include "math/bcl_math_running_min_max.h"
#include "math/bcl_math_template_instantiations.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "mc/bcl_mc_temperature_default.h"
#include "mc/bcl_mc_temperature_interface.h"
#include "model/bcl_model_retrieve_interface.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_function.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_skipped_steps.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "opti/bcl_opti_tracker.h"
#include "random/bcl_random_uniform_distribution.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf_mdl_handler.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_call_stack.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeFit
    //! @brief Application for fitting molecules against other molecules (alignment), protein binding-pockets (dock),
    //! and/or pharmacophore maps (pharm map). Approaches can be used individually or in custom combinations.
    //!
    //! @author brownbp1
    //! @date 05/18/2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeFit :
      public InterfaceRelease
    {

    private:

      // ThreadManager needs access to private nested classes
      friend class ThreadManager;
      friend class Worker;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ThreadManager
      //! @brief manages threads during molecule fitting (one input molecule per thread)
      //!
      //! @author brownbp1
      //! @date May 18, 2020
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class ThreadManager :
        public util::ObjectInterface
      {

      private:

      ///////////
      // Data //
      ///////////

          // Basic data for general use
          const size_t                                        m_Threads;                    // Number of threads
          size_t                                              m_CurrentMolIndex;            // Index of the molecule to fit
          const size_t                                        m_NumberMolsToFit;            // Number of molecules read in to be fit
          size_t                                              m_NumberMolsFit;              // Number of molecules that have been fit
          storage::Vector< chemistry::FragmentComplete>       m_Molecules;                  // The molecules which have been built and are ready for output
          io::OFStream                                        m_OutputStream;               // Output file to write molecules to
          sched::Mutex                                        m_Mutex;                      // Lock for updating Workers
          opti::Tracker< chemistry::FragmentComplete, double> m_OptiGoal;                   // Optimization criterion
          std::string                                         m_Routine;                    // The fit routine to perform

          // Data for conformer generation
          chemistry::SampleConformations                      m_SampleConfs;                // Base small molecule conformer generator
          chemistry::SampleConformations                      m_SampleConfsLocal;           // Local small molecule conformer generator
          storage::Vector< size_t>                            m_LocalSamplingPreferences;   // Sampling preferences for local conformer generator

          // Data for small molecule alignment
          chemistry::FragmentAlignToScaffold                  m_AlignScaffolds;              // Align molecules to most similar scaffold according to substructure
          util::Implementation
          <
          chemistry::ConformationComparisonInterface
          >                                                   m_Comparer;                    // Comparison to use for substructure-based alignment
          chemistry::ConformationComparisonPsiField           m_RigidAlign;                  // Align molecules to scaffolds based on property similarity (rigid)
          chemistry::ConformationComparisonPsiFlexField       m_MolAlign;                    // Align molecules to scaffolds based on property similarity (flexible)
          const size_t                                        m_AlignmentSolutions;          // Alignment solutions per scaffold
          const bool                                          m_RefineAlignment;             // Refine alignment with local conformer alignments
          const bool                                          m_SkipInitialAlignmentFlag;        // Refine starting pose with local conformer alignments

          // Data for structure-based pose refinement
          chemistry::FragmentStochasticPoseOptimizer          m_PoseOptimizer;               // Rotate/translate/confswap ligand to improve binding pose interactions
          descriptor::CheminfoProperty                        m_PoseDependentPropertyScorer; // Property from which to build the objective function
          const std::string                                   m_PoseDependentMDLProperty;    // Pose-dependent scoring receptor MDL property
          const std::string                                   m_ReceptorFilename;            // Pose-dependent scoring receptor filename
          const size_t                                        m_PoseSamplingIterations;      // Iterations for sampling the protein-ligand pose
          const size_t                                        m_PoseSamplingCycles;          // Number of MCM cycles to perform during pose optimization
          const size_t                                        m_RefinementIterations;        // Refinement iterations for pose sampling
          const float                                         m_MetropolisTemperature;       // Temperature of the Metropolis criterion

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //!
          //! @class Worker
          //! @brief runs the threads
          //!
          //! @author brownbp1
          //! @date May 18, 2020
          //!
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          struct Worker
          {

            // Basic data for general use
            std::string                                                   m_WorkerRoutine;             // The fit routine to perform
            util::ShPtr< storage::Vector< chemistry::FragmentComplete>>   m_InputFragments;            // Input molecules
            util::ShPtr< storage::Vector< chemistry::FragmentComplete>>   m_ScaffoldFragments;         // Scaffold fragments against which to perform alignments
            size_t                                                        m_CurrentWorkerMolIndex;     // index of the molecule to fit
            util::SiPtr< ThreadManager>                                   m_ThreadManager;             // Pointer to the parent thread manager; needed to update Worker

            // Data for small molecule alignment
            chemistry::FragmentAlignToScaffold                            m_WorkerAlignScaffolds;      // Align molecules to scaffold(s)
            chemistry::ConformationComparisonPsiField                     m_WorkerRigidAlign;          // Align molecules to scaffolds based on property similarity (rigid)
            chemistry::ConformationComparisonPsiFlexField                 m_WorkerMolAlign;            // Align molecules to scaffolds based on property similarity (flexible)
            util::Implementation
            <
            chemistry::ConformationComparisonInterface
            >                                                             m_WorkerComparer;            // Comparison to use for substructure-based alignment

            // Data for conformer generation
            chemistry::SampleConformations                                m_WorkerSampleConfs;         // Generate small molecule conformational ensemble
            chemistry::SampleConformations                                m_WorkerSampleConfsLocal;    // Generate local small molecule conformational ensemble

            // Data for structure-based pose refinement
            chemistry::FragmentStochasticPoseOptimizer                    m_WorkerPoseOptimizer;       // Rotate/translate/confswap ligand to improve binding pose interactions

        //////////////
        // operator //
        //////////////

            //! @brief Run the thread to perform the chosen fit operations
            void RunThread()
            {
              // run threads
              do
              {
                // alignment
                if( m_WorkerRoutine == "0")
                {
                  SmallMoleculeAlignment( ( *m_InputFragments)( m_CurrentWorkerMolIndex), *m_ScaffoldFragments);
                }
                else if( m_WorkerRoutine == "1")
                {
                  MolAlignThenLocalDock( ( *m_InputFragments)( m_CurrentWorkerMolIndex), *m_ScaffoldFragments);
                }
                else if( m_WorkerRoutine == "2")
                {
                  if( !m_ThreadManager->m_SkipInitialAlignmentFlag)
                  {
                    MCSEnsembleAlignment
                    (
                      ( *m_InputFragments)( m_CurrentWorkerMolIndex),
                      *m_ScaffoldFragments,
                      m_WorkerComparer,
                      false
                    );
                  }
                  if( m_ThreadManager->m_RefineAlignment)
                  {
                    MCSEnsembleAlignment
                    (
                      ( *m_InputFragments)( m_CurrentWorkerMolIndex),
                      *m_ScaffoldFragments,
                      m_WorkerComparer,
                      true
                    );
                  }
                }
                else if( m_WorkerRoutine == "3")
                {
                  PoseSensitiveMCSAlignment( ( *m_InputFragments)( m_CurrentWorkerMolIndex), *m_ScaffoldFragments);
                }
                else if( m_WorkerRoutine == "4")
                {
                  // pose optimization
                  ProteinLigandDocking( ( *m_InputFragments)( m_CurrentWorkerMolIndex), true);
                }
                else
                {
                  BCL_Assert( false, "No method selected. Exiting...");
                }

                // add fit molecules to molecule vector and increment index
                m_ThreadManager->m_Mutex.Lock();
                if( m_CurrentWorkerMolIndex < m_ThreadManager->GetNumberMoleculesToFit())
                {
                  m_ThreadManager->AddMolecule( (*m_InputFragments)( m_CurrentWorkerMolIndex), m_CurrentWorkerMolIndex);
                  m_ThreadManager->UpdateWorker( *this);
                  m_ThreadManager->m_Mutex.Unlock();
                }
                else
                {
                  // gotta unlock this since it won't happen in UpdateWorker
                  m_ThreadManager->m_Mutex.Unlock();
                }
              } while( m_ThreadManager->GetNumberMoleculesFit() < m_ThreadManager->GetNumberMoleculesToFit() - m_ThreadManager->GetNumberThreads() + 2);
            } // RunThread()

        ////////////////
        // operations //
        ////////////////

            //! @brief Generate small molecule conformational ensemble
            //! @param MOLECULE molecule of interest
            //! @param LOCAL restrict to local conformer sampling
            chemistry::FragmentEnsemble GenerateConformers
            (
              const chemistry::FragmentComplete &MOLECULE,
              const bool &LOCAL
            )
            {
              // generate conformers
              chemistry::FragmentEnsemble ens;
              if( LOCAL)
              {
                ens = m_WorkerSampleConfsLocal( MOLECULE).First();
              }
              else
              {
                ens = m_WorkerSampleConfs( MOLECULE).First();
              }

              if( ens.GetSize())
              {
                return ens;
              }
              else
              {
                return chemistry::FragmentEnsemble( storage::List< chemistry::FragmentComplete>( size_t( 1), MOLECULE));
              }
            }

            //! @brief Perform small molecule alignments
            //! @param MOLECULE molecule to be aligned
            //! @param SCAFFOLDS molecules against with MOLECULE will be aligned
            void SmallMoleculeMCSAlignment
            (
              chemistry::FragmentComplete &MOLECULE,
              const storage::Vector< chemistry::FragmentComplete> &SCAFFOLDS
            )
            {
              // compute the similarity of our current molecule to each of the scaffold molecules
              chemistry::FragmentEnsemble scaffolds( storage::List< chemistry::FragmentComplete>( SCAFFOLDS.Begin(), SCAFFOLDS.End()));
              static descriptor::MoleculeSimilarity similarity( "LargestCommonDisconnectedSubstructureTanimoto", scaffolds);
              linal::Vector< float> similarities( similarity.SumOverObject( MOLECULE));

              // find the most similar scaffold
              float most_similar( similarities.Max());
              size_t most_similar_index( 0);
              for( auto itr( similarities.Begin()), itr_end( similarities.End()); itr != itr_end; ++itr, ++most_similar_index)
              {
                if( *itr == most_similar)
                {
                  break;
                }
              }

              // align molecules
              m_WorkerAlignScaffolds.AlignToScaffold( MOLECULE, SCAFFOLDS( most_similar_index));
            }

            void SmallMoleculeAlignment
            (
              chemistry::FragmentComplete &MOLECULE,
              const storage::Vector< chemistry::FragmentComplete> &SCAFFOLDS
            )
            {
              // will need to collect all of the final poses
              chemistry::FragmentEnsemble final_poses;

              // Perform property-bases small molecule alignment
              storage::Vector< // indexes scaffold
              storage::Vector< // indexes alignment solution
                storage::Triplet< chemistry::FragmentComplete, chemistry::FragmentComplete, double // single alignment result
                > > > aligned_mols( PropertyBasedAlignment( MOLECULE, SCAFFOLDS));

              // Find the best alignment
              float min_rmsdx( math::GetHighestBoundedValue< float>());
              for
              (
                  size_t scaffold_index( 0), n_scaffolds( aligned_mols.GetSize());
                  scaffold_index < n_scaffolds;
                  ++scaffold_index
              )
              {
                // go over each alignment solution per scaffold
                for
                (
                    size_t alignment_index( 0), n_alignments( aligned_mols( scaffold_index).GetSize());
                    alignment_index < n_alignments;
                    ++alignment_index
                )
                {
                  // make a reference to our mol to shorthand this
                  chemistry::FragmentComplete &mol( aligned_mols( scaffold_index)( alignment_index).First());
                  float mol_rmsdx( mol.GetStoredProperties().GetMDLPropertyAsVector( "RMSDX")( 0));

                  // save the best alignment
                  if( mol_rmsdx < min_rmsdx)
                  {
                    min_rmsdx = mol_rmsdx;
                    MOLECULE = mol;
                  }
                } // end current alignment solution
              } // end alignment solutions for current scaffold
            }

            //! @brief Perform small molecule alignments
            //! @param MOLECULE molecule to be aligned
            //! @param SCAFFOLDS molecules against with MOLECULE will be aligned
            void MCSEnsembleAlignment
            (
              chemistry::FragmentComplete &MOLECULE,
              const storage::Vector< chemistry::FragmentComplete> &SCAFFOLDS,
              const util::Implementation< chemistry::ConformationComparisonInterface> COMPARER,
              const bool LOCAL
            )
            {
              // generate conformers
              storage::Vector< storage::Pair< chemistry::FragmentComplete, float> > best_confs( SCAFFOLDS.GetSize());
              chemistry::FragmentEnsemble main_ensemble( GenerateConformers( MOLECULE, LOCAL));

              // align molecules to each scaffold
              size_t scaffold_index( 0);
              for
              (
                  auto scaffold_itr( SCAFFOLDS.Begin()), scaffold_itr_end( SCAFFOLDS.End());
                  scaffold_itr != scaffold_itr_end;
                  ++scaffold_itr, ++scaffold_index
              )
              {
                // this is stupid; I should re-design this so that I don't have to copy the ensemble
                // just unroll the 2d vector into 1d
                chemistry::FragmentEnsemble ensemble( main_ensemble);

                // perform alignment
                storage::Vector< storage::Pair< bool, float> > alignment_score
                (
                  m_WorkerAlignScaffolds.AlignEnsembleToScaffold( ensemble, *scaffold_itr, COMPARER)
                );

                // get the best conformer from this set of alignments
                ensemble.Sort( COMPARER->GetAlias());
                size_t best_conf_index( ensemble.GetMolecules().FirstElement().GetStoredPropertiesNonConst().GetMDLPropertyAsVector( "ConformerIndex")( 0));

                // track the best conformer for this scaffold and its MolAlign score
                best_confs( scaffold_index) = std::make_pair
                (
                  ensemble.GetMolecules().FirstElement(),
                  alignment_score( best_conf_index).First()
                );
              }

              // get the best scoring conformer alignment
              float best_alignment( best_confs( 0).Second());
              float best_alignment_i( 0);
              for
              (
                  size_t aln_i( 1), aln_i_sz( SCAFFOLDS.GetSize());
                  aln_i < aln_i_sz;
                  ++aln_i
              )
              {
                if( best_confs( aln_i).Second() < best_alignment)
                {
                  best_alignment = best_confs( aln_i).Second();
                  best_alignment_i = aln_i;
                }
              }
              MOLECULE = best_confs( best_alignment_i).First();
            }

            //! @brief Perform small molecule alignments
            //! @param MOLECULE molecule to be aligned
            //! @param SCAFFOLDS molecules against with MOLECULE will be aligned
            void PoseSensitiveMCSAlignment
            (
              chemistry::FragmentComplete &MOLECULE,
              const storage::Vector< chemistry::FragmentComplete> &SCAFFOLDS
            )
            {
              // prepare molecule
              MOLECULE.GetStoredPropertiesNonConst().SetMDLProperty
              (
                m_ThreadManager->m_PoseDependentMDLProperty,
                m_ThreadManager->m_ReceptorFilename
              );

              // generate conformers
              chemistry::FragmentEnsemble ensemble( GenerateConformers( MOLECULE, false));

              // align molecules to each scaffold
              float best_score( m_ThreadManager->m_PoseDependentPropertyScorer->SumOverObject( MOLECULE)( 0));
              for
              (
                  auto scaffold_itr( SCAFFOLDS.Begin()), scaffold_itr_end( SCAFFOLDS.End());
                  scaffold_itr != scaffold_itr_end;
                  ++scaffold_itr
              )
              {
                // go over each conformer of current molecule
                for
                (
                    auto conf_itr( ensemble.Begin()), conf_itr_end( ensemble.End());
                    conf_itr != conf_itr_end;
                    ++conf_itr
                )
                {
                  // perform alignment
                  storage::Pair< bool, float> alignment_score
                  (
                    m_WorkerAlignScaffolds.PoseSensitiveAlignToScaffold
                    (
                      *conf_itr,
                      *scaffold_itr,
                      m_ThreadManager->m_PoseDependentPropertyScorer,
                      m_ThreadManager->m_PoseDependentMDLProperty,
                      m_ThreadManager->m_ReceptorFilename
                    )
                  );
                  // if improved score, take it
                  if( alignment_score.First())
                  {
                    if
                    (
                        ( m_ThreadManager->m_OptiGoal.GetImprovementType() == opti::e_SmallerIsBetter && alignment_score.Second() < best_score) ||
                        ( m_ThreadManager->m_OptiGoal.GetImprovementType() == opti::e_LargerIsBetter && alignment_score.Second() > best_score)
                    )
                    {
                      MOLECULE = *conf_itr;
                    }
                  }
                }
              }
            }

            storage::Vector< // indexes scaffold
            storage::Vector< // indexes alignment solution
              storage::Triplet< chemistry::FragmentComplete, chemistry::FragmentComplete, double> // single alignment result
              >
            > PropertyBasedAlignment
            (
              chemistry::FragmentComplete &MOLECULE,
              const storage::Vector< chemistry::FragmentComplete> &SCAFFOLDS
            )
            {
              // output messages for users
              BCL_MessageStd( "Preparing small molecule ligand for flexible property-based alignment...");

              // generate conformers of the molecule
              BCL_MessageStd( "Generating small molecule ligand conformers...");
              chemistry::FragmentEnsemble ensemble( GenerateConformers( MOLECULE, false));

              // loop over each scaffold
              size_t n_scaffolds( SCAFFOLDS.GetSize());
              storage::Vector< storage::Vector< storage::Triplet< chemistry::FragmentComplete, chemistry::FragmentComplete, double> > > all_alignments( n_scaffolds);
              for
              (
                  size_t scaffold_index( 0);
                  scaffold_index < n_scaffolds;
                  ++scaffold_index
              )
              {
                // cast scaffold to FragmentEnsemble
                chemistry::FragmentEnsemble scaffold( storage::List< chemistry::FragmentComplete>( size_t( 1), SCAFFOLDS( scaffold_index)));

                // for each molecule-scaffold pair, we generate multiple solutions
                BCL_MessageStd( "Aligning molecule to scaffold: " + util::Format()( scaffold_index) + " (0-indexed)");
                storage::Vector< storage::Triplet< chemistry::FragmentComplete, chemistry::FragmentComplete, double> > aligned_mols
                (
                  m_WorkerMolAlign.FieldOptimizeOrientationFlex
                  (
                    ensemble,
                    scaffold,
                    m_ThreadManager->m_MolAlign.GetNumberConformerPairs(),
                    m_ThreadManager->m_MolAlign.GetNumberIterations(),
                    m_ThreadManager->m_MolAlign.GetNumberMaxUnimproved(),
                    m_ThreadManager->m_MolAlign.GetNumberFilterIterations(),
                    m_ThreadManager->m_MolAlign.GetNumberFilterMaxUnimproved(),
                    m_ThreadManager->m_MolAlign.GetNumberRefinementIterations(),
                    m_ThreadManager->m_MolAlign.GetNumberRefinementMaxUnimproved(),
                    m_ThreadManager->m_MolAlign.GetFractionFilterInitial(),
                    m_ThreadManager->m_MolAlign.GetFractionFilterIterative()
                  )
                );
                // add our alignment solutions to collection
                all_alignments( scaffold_index) = aligned_mols;
                BCL_MessageStd( "Completed alignment to scaffold: " + util::Format()( scaffold_index) + " (0-indexed)");
                BCL_MessageStd( "Maximum number of alignment solutions requested: " + util::Format()( m_ThreadManager->m_AlignmentSolutions));
                BCL_MessageStd( "Number of alignment solutions returned: " + util::Format()( aligned_mols.GetSize()));
              }
              return all_alignments;
            }

            //! @brief Perform protein-ligand docking
            //! @param MOLECULE molecule to be docked
            //! @param LOCAL whether to restrict to local refinement
            void ProteinLigandDocking
            (
              chemistry::FragmentComplete &MOLECULE,
              bool LOCAL
            )
            {
              // prepare molecule
              BCL_MessageStd( "Preparing small molecule ligand for docking...");
              MOLECULE.GetStoredPropertiesNonConst().SetMDLProperty
              (
                m_ThreadManager->m_PoseDependentMDLProperty,
                m_ThreadManager->m_ReceptorFilename
              );

              // optimize pose
              BCL_MessageStd( "Beginning Monte Carlo - Metropolis docking cycles!");
              BCL_MessageStd( "Optimization improvement type: " + util::Format()( m_ThreadManager->m_OptiGoal.GetImprovementType().GetString()));
              util::ShPtr< chemistry::FragmentComplete> opti_mol;
              if( LOCAL)
              {
                // generate local conformers
                BCL_MessageStd( "Generating small molecule ligand conformers...");
                chemistry::FragmentEnsemble ensemble( GenerateConformers( MOLECULE, true));

                // run local refinement
                opti_mol = m_WorkerPoseOptimizer.RunMCMDocker( ensemble, true, false);
              }
              else
              {
                // generate conformers
                BCL_MessageStd( "Generating small molecule ligand conformers...");
                chemistry::FragmentEnsemble ensemble( GenerateConformers( MOLECULE, false));

                // run docking
                opti_mol = m_WorkerPoseOptimizer.RunMCMDocker( ensemble);
              }
              if( opti_mol.IsDefined() && opti_mol->GetSize())
              {
                BCL_MessageStd( "Saving best pose from docking...");
                MOLECULE = *opti_mol;
              }
            }

            //! @param MOLECULE molecule to be aligned and then docked
            //! @param SCAFFOLDS molecules against which MOLECULE will be aligned
            void MCSAlignThenDock
            (
              chemistry::FragmentComplete &MOLECULE,
              const storage::Vector< chemistry::FragmentComplete> &SCAFFOLDS
            )
            {
              SmallMoleculeMCSAlignment( MOLECULE, SCAFFOLDS);
              ProteinLigandDocking( MOLECULE, true);
            }

            //! @param MOLECULE molecule to be aligned then docked
            //! @param SCAFFOLDS molecules against which MOLECULE will be aligned
            void MolAlignThenLocalDock
            (
              chemistry::FragmentComplete &MOLECULE,
              const storage::Vector< chemistry::FragmentComplete> &SCAFFOLDS
            )
            {
              // NOTE - recommend low cycles high refinement iterations

              // will need to collect all of the final poses
              chemistry::FragmentEnsemble final_poses;

              // Perform property-bases small molecule alignment
              storage::Vector< // indexes scaffold
              storage::Vector< // indexes alignment solution
                storage::Triplet< chemistry::FragmentComplete, chemistry::FragmentComplete, double // single alignment result
                > > > aligned_mols( PropertyBasedAlignment( MOLECULE, SCAFFOLDS));

              // Use each one of the alignment results for a refine docking run
              for
              (
                  size_t scaffold_index( 0), n_scaffolds( aligned_mols.GetSize());
                  scaffold_index < n_scaffolds;
                  ++scaffold_index
              )
              {
                // go over each alignment solution per scaffold
                for
                (
                    size_t alignment_index( 0), n_alignments( aligned_mols( scaffold_index).GetSize());
                    alignment_index < n_alignments;
                    ++alignment_index
                )
                {
                  // make a reference to our mol to shorthand this
                  chemistry::FragmentComplete &mol( aligned_mols( scaffold_index)( alignment_index).First());

                  // prepare molecule by setting up pose-dependent properties
                  BCL_MessageStd( "Preparing small molecule ligand for docking...");
                  mol.GetStoredPropertiesNonConst().SetMDLProperty
                  (
                    m_ThreadManager->m_PoseDependentMDLProperty,
                    m_ThreadManager->m_ReceptorFilename
                  );

                  // generate local conformers for the alignment solution
                  BCL_MessageStd( "Generating local small molecule ligand conformational ensemble...");
                  chemistry::FragmentEnsemble ensemble( GenerateConformers( mol, true));

                  // optimize pose
                  BCL_MessageStd( "Beginning Monte Carlo - Metropolis local docking refinement cycles!");
                  BCL_MessageStd( "Optimization improvement type: " + util::Format()( m_ThreadManager->m_OptiGoal.GetImprovementType().GetString()));
                  util::ShPtr< chemistry::FragmentComplete> opti_mol( m_WorkerPoseOptimizer.RunMCMDocker( ensemble, true, false));
                  if( opti_mol.IsDefined() && opti_mol->GetSize())
                  {
                    mol = *opti_mol;
                  }

                  // save pose
                  final_poses.PushBack( mol);
                } // end current alignment solution
              } // end alignment solutions for current scaffold

              // sort all of the pose-refined alignment solutions
              final_poses.Sort( "FLD_Score");

              // select best pose to save
              if( m_ThreadManager->m_OptiGoal.GetImprovementType() == opti::e_LargerIsBetter)
              {
                MOLECULE = final_poses.GetMolecules().LastElement();
              }
              else
              {
                MOLECULE = final_poses.GetMolecules().FirstElement();
              }
            } // end MolAlignThenLocalDock

          }; // struct Worker

      ////////////////
      // operations //
      ////////////////

          //! @brief Updates worker threads
          //! @param WORKER worker thread whose index is incremented
          void UpdateWorker( ThreadManager::Worker &WORKER)
          {
            // increase worker molecule index
            WORKER.m_CurrentWorkerMolIndex = m_CurrentMolIndex;

            // increase the molecule fit count
            ++m_CurrentMolIndex;
            ++m_NumberMolsFit;
          } // UpdateWorker

      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

          //! brief constructor
          ThreadManager
          (
            const util::ShPtr< storage::Vector< chemistry::FragmentComplete> >       &INPUT_FRAGMENTS,             // Input molecules
            const std::string                                                       &OUTPUT_FILENAME,             // Output SDF filename
            const size_t                                                             NUMBER_THREADS,              // Number of threads (from scheduler)
            const util::ShPtr< storage::Vector< chemistry::FragmentComplete> >       &SCAFFOLD_FRAGMENTS,          // Scaffold fragments against which to align inputs
            const descriptor::CheminfoProperty                                      &PROPERTY_SCORER,             // Pose-dependent scoring function
            const std::string                                                       &POSE_DEPENDENT_MDL_PROPERTY, // Pose-dependent scoring MDL property receptor
            const std::string                                                       &RECEPTOR_FILENAME,           // Pose-dependent scoring receptor filename
            const std::string                                                       &ROUTINE,                     // Molecule fit routine
            const size_t                                                             POSE_SAMPLING_ITERATIONS,    // Iterations for pose sampling
            const size_t                                                             POSE_CYCLES,                 // Number of MCM cycles to perform
            const size_t                                                             REFINEMENT_ITERATIONS,       // Refinement iterations
            const float                                                              METROPOLIS_TEMP,             // Temperature of the Metropolis criterion
            const chemistry::FragmentAlignToScaffold                                &MCS_ALIGN,                   // MCS alignment object
            const chemistry::ConformationComparisonPsiField                          &RIGID_ALIGN,                  // Rigid MolAlign
            const chemistry::ConformationComparisonPsiFlexField                      &FLEX_ALIGN,                  // Flex MolAlign
            const size_t                                                             ALIGNMENT_SOLUTIONS,         // Alignment solutions per scaffold
            const opti::Tracker< chemistry::FragmentComplete, double>               &OPTI_GOAL,                   // Optimization criterion
            const chemistry::SampleConformations                                    &SAMPLE_CONFS,                // Sample conformers object
            const storage::Vector< size_t>                                          &LOCAL_SAMPLE_PREFS,          // Local sampling preferences
            const util::Implementation< chemistry::ConformationComparisonInterface> &COMPARER,
            const bool                                                               REFINE_ALIGNMENT,            // Refine an initial alignment with local conformer alignments
            const bool                                                               SKIP_INITIAL_ALIGNMENT       // Perform alignment refinement on input pose only
          ) :
            m_Threads( std::min( NUMBER_THREADS, INPUT_FRAGMENTS->GetSize())),
            m_CurrentMolIndex( NUMBER_THREADS - 1),
            m_NumberMolsToFit( INPUT_FRAGMENTS->GetSize()),
            m_NumberMolsFit( 0),
            m_Molecules( INPUT_FRAGMENTS->GetSize()),
            m_PoseDependentPropertyScorer( PROPERTY_SCORER),
            m_PoseDependentMDLProperty( POSE_DEPENDENT_MDL_PROPERTY),
            m_ReceptorFilename( RECEPTOR_FILENAME),
            m_Routine( ROUTINE),
            m_PoseSamplingIterations( POSE_SAMPLING_ITERATIONS),
            m_PoseSamplingCycles( POSE_CYCLES),
            m_RefinementIterations( REFINEMENT_ITERATIONS),
            m_MetropolisTemperature( METROPOLIS_TEMP),
            m_AlignScaffolds( MCS_ALIGN),
            m_RigidAlign( RIGID_ALIGN),
            m_MolAlign( FLEX_ALIGN),
            m_AlignmentSolutions( ALIGNMENT_SOLUTIONS),
            m_OptiGoal( OPTI_GOAL),
            m_SampleConfs( SAMPLE_CONFS),
            m_LocalSamplingPreferences( LOCAL_SAMPLE_PREFS),
            m_Comparer( COMPARER),
            m_RefineAlignment( REFINE_ALIGNMENT),
            m_SkipInitialAlignmentFlag( SKIP_INITIAL_ALIGNMENT)
          {
            // prepare output filestream
            io::File::MustOpenOFStream( m_OutputStream, OUTPUT_FILENAME);

            // prepare operations to pass to worker
            SetupConformerGenerator(); // conformer generation
            SetupMolAlign();           // property-based small molecule alignment
            SetupPoseOptimizer();      // protein-ligand dock

            // initialize one worker per thread
            std::vector< Worker> workers( m_Threads);
            size_t current_mol_index( 0);
            for
            (
                std::vector< Worker>::iterator itr( workers.begin()), end( workers.end());
                itr != end;
                ++itr, ++current_mol_index
            )
            {
              // create worker object
              Worker &worker_ref( *itr);
              worker_ref.m_WorkerRoutine                    = m_Routine;                     // Set the routine for the worker
              worker_ref.m_InputFragments                   = INPUT_FRAGMENTS;               // Input molecules
              worker_ref.m_ScaffoldFragments                = SCAFFOLD_FRAGMENTS;            // Scaffold molecules for alignment
              worker_ref.m_CurrentWorkerMolIndex            = current_mol_index;             // Current molecule index
              worker_ref.m_ThreadManager                    = this;                          // Copy of pointer to thread manager
              worker_ref.m_WorkerSampleConfs                = m_SampleConfs;                 // Conformer generator
              worker_ref.m_WorkerSampleConfsLocal           = m_SampleConfsLocal;            // Local conformer generator
              worker_ref.m_WorkerAlignScaffolds             = m_AlignScaffolds;              // AlignToScaffold
              worker_ref.m_WorkerRigidAlign                 = m_RigidAlign;                  // MolAlign (rigid)
              worker_ref.m_WorkerMolAlign                   = m_MolAlign;                    // MolAlign (flexible)
              worker_ref.m_WorkerPoseOptimizer              = m_PoseOptimizer;               // Docking
              worker_ref.m_WorkerComparer                   = m_Comparer;                    // Alignment score
            }

            // Allocate space for jobs
            util::ShPtrVector< sched::JobInterface> jobs;
            jobs.AllocateMemory( m_Threads);

            // assign group id for job priority (all jobs same priority)
            const size_t group_id( 0);

            // add threads to jobs
            for( size_t proc_number( 0); proc_number < m_Threads; ++proc_number)
            {
              Worker &worker_ref( workers[ proc_number]);
              jobs.PushBack
              (
                util::ShPtr< sched::JobInterface>
                (
                  new sched::ThunkJob< Worker, void>
                  (
                    group_id,
                    worker_ref,
                    &Worker::RunThread,
                    sched::JobInterface::e_READY,
                    NULL
                  )
                )
              );
              // Submit the jobs to the scheduler
              sched::GetScheduler().RunJob( jobs.LastElement());
            }

            // join all the jobs
            for( size_t proc_number( 0); proc_number < m_Threads; ++proc_number)
            {
              sched::GetScheduler().Join( jobs( proc_number));
            }

            // Write output
            for( size_t mol_index( 0); mol_index < m_Molecules.GetSize(); ++mol_index)
            {
              m_Molecules( mol_index).WriteMDL( m_OutputStream);
            }

            // Close output
            io::File::CloseClearFStream( m_OutputStream);
          }; // ThreadManager()

          //! @brief clone function
          ThreadManager *Clone() const
          {
            BCL_Exit( "ThreadManager cannot be cloned.", -1);
            return NULL;
          }

          //! @brief Get class identifier string
          const std::string &GetClassIdentifier() const
          {
            return GetStaticClassName( *this);
          }

      /////////////////
      // data access //
      /////////////////

          // Return the current molecule index
          size_t GetCurrentMolIndex()
          {
            return m_CurrentMolIndex;
          }

          // Return the number of molecules that have been fit
          size_t GetNumberMoleculesFit()
          {
            return m_NumberMolsFit;
          }

          // Return the number of molecules to fit
          size_t GetNumberMoleculesToFit()
          {
            return m_NumberMolsToFit;
          }

          // Return the number of threads being used
          size_t GetNumberThreads()
          {
            return m_Threads;
          }

          // Return pocket filename
          std::string GetReceptorFilename()
          {
            return m_ReceptorFilename;
          }

          // Return molecules
          storage::Vector< chemistry::FragmentComplete> &GetMolecules()
          {
            return m_Molecules;
          }

      //////////////////////
      // helper functions //
      //////////////////////

          // Add fit molecules to output vector
          void AddMolecule
          (
            const chemistry::FragmentComplete &MOLECULE,
            const size_t &MOL_INDEX
          )
          {
            m_Molecules( MOL_INDEX) = MOLECULE;
          }

          // Setup the small molecule conformer generator
          void SetupConformerGenerator()
          {
            // make a local sampler from the base object
            m_SampleConfsLocal = m_SampleConfs;
            m_SampleConfsLocal.SetSamplingPreferences
            (
              bool( m_LocalSamplingPreferences( 0)),
              bool( m_LocalSamplingPreferences( 1)),
              bool( m_LocalSamplingPreferences( 2)),
              bool( m_LocalSamplingPreferences( 3))
            );
          }

          // Setup MolAlign object
          void SetupMolAlign()
          {
            // setup the property-based alignment object
            if( 1)
            {
              m_MolAlign.SetSampleConformations( m_SampleConfs);
            }
            else
            {
              m_MolAlign.SetSampleConformations( m_SampleConfsLocal);
            }
          }

          // Setup protein-ligand pose refinement object
          void SetupPoseOptimizer()
          {
            // setup pose optimization object
            m_PoseOptimizer = chemistry::FragmentStochasticPoseOptimizer
            (
              m_PoseDependentPropertyScorer,
              storage::Vector< float>(),
              m_PoseDependentMDLProperty,
              m_ReceptorFilename,
              m_PoseSamplingIterations,
              m_PoseSamplingCycles,
              m_RefinementIterations,
              m_MetropolisTemperature,
              m_OptiGoal,
              float( 5.0)
            );
          }

      //////////////////
      // input/output //
      /////////////////

      protected:

          std::istream &Read( std::istream &INSTREAM)
          {
            return INSTREAM;
          }

          std::ostream &Write( std::ostream &OUTSTREAM, const size_t INDENT) const
          {
            return OUTSTREAM;
          }

      }; // class ThreadManager

    //////////
    // data //
    //////////

      //! flag for input molecules to be fit
      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

      //! flag for defining output filename,
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! flag to control input base fragment
      util::ShPtr< command::FlagInterface> m_ScaffoldFragmentsFlag;

      //! flag to set the resolution of atom-type comparisons for substructure-based alignment
      util::ShPtr< command::FlagInterface> m_AtomComparisonFlag;

      //! flag to set resolution of bond-type comparisons for substructure-based alignment
      util::ShPtr< command::FlagInterface> m_BondComparisonFlag;

      //! flag to specify the sample conformations (BCL::Conf) object settings
      util::ShPtr< command::FlagInterface> m_SampleConfsFlag;

      //! flag to specify degrees of freedom for local conformer sampling
      util::ShPtr< command::FlagInterface> m_LocalSamplingPrefsFlag;

      //! flag to refine an alignment with local conformer alignments
      util::ShPtr< command::FlagInterface> m_RefineAlignmentFlag;

      //! flag to skip the initial alignment and just refine the input pose with local conformer alignments
      util::ShPtr< command::FlagInterface> m_SkipInitialAlignmentFlag;

      //! flag to specify the conformation comparison type
      util::ShPtr< command::FlagInterface> m_ConformerComparerFlag;

      //! flag for the pose-dependent score function
      util::ShPtr< command::FlagInterface> m_PoseDependentScoringFunctionFlag;

      //! flag for the number of iterations to be used in pose-dependent optimization
      util::ShPtr< command::FlagInterface> m_PoseSamplingIterationsFlag;

      //! flag for the number of MCM cycles to perform during pose-dependent optimization
      util::ShPtr< command::FlagInterface> m_PoseSamplingCyclesFlag;

      //! flag for the number of MCM cycles to perform during pose-dependent refinement
      util::ShPtr< command::FlagInterface> m_RefinementIterationsFlag;

      //! flag for the Metroplis criterion temperature
      util::ShPtr< command::FlagInterface> m_MetropolisTemperatureFlag;

      //! flag to remove SampleByParts prior to generating conformers
      util::ShPtr< command::FlagInterface> m_RemoveSampleByPartsFlag;

      //! flag to reverse directionality of optimization (from smaller to larger)
      util::ShPtr< command::FlagInterface> m_LargerScoreBetter;

      //! flag for the number of alignment solutions to carry forward into pose refinement
      util::ShPtr< command::FlagInterface> m_AlignmentSolutionsPerScaffoldFlag;

      //! flag to specify what type of fit routine to perform
      util::ShPtr< command::FlagInterface> m_RoutineFlag;

      //! flag to specify the MolAlign object settings
      util::ShPtr< command::FlagInterface> m_MolAlignRigidFlag;

      //! flag to specify the MolAlign object settings
      util::ShPtr< command::FlagInterface> m_MolAlignFlexFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeFit();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MoleculeFit *Clone() const
      {
        return new MoleculeFit( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Fits other molecules against one another utilizing several routines.";
      }

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const
      {
        static std::string s_read_me =
          "MoleculeFit fit is an in-development parallelizable application that implements a general framework to perform "
          "small molecule alignment, protein-ligand docking, and pharmacophore fitting. Currently, several routines are "
          "available:"
            "0 - MolAlign - performs property-based small molecule flexible alignment of the input molecules to the scaffold molecules and returns the pose of the single best alignment \n"
            "1 - MolAlign with structure-based refinement - performs MolAlign routine followed by local refinement of protein-ligand interactions with PLC-DNN and MCM-guided perturbations \n"
            "2 - Maximum common substructure alignment - enumerates MCS alignments to scaffolds with a conformational ensemble; chooses the best pose by MolAlign score \n"
            "3 - Pose-sensitive maximum common substructure alignment - enumerates MCS alignments to scaffolds with a conformational ensemble; chooses the best pose by PLC-DNN score \n"
            "4 - Protein-ligand refinement - performs refinement of protein-ligand complex using MCM-guided perturbations and PLC-DNN for pose-dependent scoring \n"
            "5 - Protein-ligand docking - docks a flexible small molecule to a rigid protein receptor using MCM-guided perturbations and PLC-DNN for pose-dependent scoring";
        return s_read_me;
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

      ////////////////////
      // common options //
      ////////////////////

        // add command line options to add/remove hydrogens
        sdf::AddMoleculeIOPrefFlags( *sp_cmd);

        // add AtomMdlLine to molecule
        sdf::MdlHandler::AddAtomMdlLineFlag( *sp_cmd);

      /////////////////////
      // special options //
      /////////////////////

        sp_cmd->AddFlag( m_InputFilenamesFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_ScaffoldFragmentsFlag);
        sp_cmd->AddFlag( m_AtomComparisonFlag);
        sp_cmd->AddFlag( m_BondComparisonFlag);
        sp_cmd->AddFlag( m_SampleConfsFlag);
        sp_cmd->AddFlag( m_LocalSamplingPrefsFlag);
        sp_cmd->AddFlag( m_RefineAlignmentFlag);
        sp_cmd->AddFlag( m_SkipInitialAlignmentFlag);
        sp_cmd->AddFlag( m_ConformerComparerFlag);
        sp_cmd->AddFlag( m_MolAlignRigidFlag);
        sp_cmd->AddFlag( m_MolAlignFlexFlag);
        sp_cmd->AddFlag( m_PoseDependentScoringFunctionFlag);
        sp_cmd->AddFlag( m_PoseSamplingIterationsFlag);
        sp_cmd->AddFlag( m_PoseSamplingCyclesFlag);
        sp_cmd->AddFlag( m_RefinementIterationsFlag);
        sp_cmd->AddFlag( m_MetropolisTemperatureFlag);
        sp_cmd->AddFlag( m_RemoveSampleByPartsFlag);
        sp_cmd->AddFlag( m_LargerScoreBetter);
        sp_cmd->AddFlag( m_AlignmentSolutionsPerScaffoldFlag);
        sp_cmd->AddFlag( m_RoutineFlag);

//
//        m_MolAlignIterationsFlag
//        m_MolAlignMaxUnimprovedFlag
//        m_MolAlignFilterIterationsFlag
//        m_MolAlignFilterMaxUnimprovedFlag
//        m_MolAlignRefinementIterationsFlag
//        m_MolAlignRefinementMaxUnimprovedFlag
//        m_MolAlignRigidTrajFlag
//        m_MolAlignFlexTrajFlag
//        m_MolAlignMaxAtomDistanceUpperFlag
//        m_MolAlignMaxAtomDistanceLowerFlag
//        m_MolAlignPoseToleranceFlag
//        m_MolAlignPoseThresholdFlag
//        m_MolAlignConfPairsFlag
//        m_MolAlignBestPairsFlag
//        m_MolAlignFractionFilterInitialFlag
//      m_MolAlignFractionFilterIterativeFlag
//        m_MolAlignRigidScaffoldFlag

      ///////////////////
      // default flags //
      ///////////////////

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const
      {

      //////////////////////////
      // read molecule inputs //
      //////////////////////////

        // open input stream for molecules
        io::IFStream input;

        // read in input molecules
        BCL_Assert( m_InputFilenamesFlag->GetFlag(), "Must provide input files!");
        storage::Vector< std::string> filenames( m_InputFilenamesFlag->GetStringList());
        util::ShPtr< chemistry::FragmentEnsemble> sp_inputs( new chemistry::FragmentEnsemble);

        for
        (
          storage::Vector< std::string>::const_iterator
            itr( filenames.Begin()), itr_end( filenames.End());
          itr != itr_end;
          ++itr
        )
        {
          // read in grow fragments ensemble
          io::File::MustOpenIFStream( input, *itr);
          sp_inputs->ReadMoreFromMdl( input, sdf::e_Maintain);
          io::File::CloseClearFStream( input);
        }
        storage::Vector< chemistry::FragmentComplete> inputs_v;
        for
        (
            auto itr( sp_inputs->Begin()), itr_end( sp_inputs->End());
            itr != itr_end;
            ++itr
        )
        {
          // remove SampleByParts MDL property if desired
          if( m_RemoveSampleByPartsFlag->GetFlag())
          {
            itr->GetStoredPropertiesNonConst().RemoveProperty( std::string( "SampleByParts"));
          }
          inputs_v.PushBack( *itr);
        }
        util::ShPtr< storage::Vector< chemistry::FragmentComplete >> sp_inputs_v( new storage::Vector< chemistry::FragmentComplete>( inputs_v));

        // check for scaffold fragment in alignment-based methods
        if
        (
            m_RoutineFlag->GetFirstParameter()->GetValue() == "0" ||
            m_RoutineFlag->GetFirstParameter()->GetValue() == "1" ||
            m_RoutineFlag->GetFirstParameter()->GetValue() == "2"
        )
        {
          BCL_Assert
          (
            m_ScaffoldFragmentsFlag->GetFlag(),
            "Requested an alignment method but did not provide a reference molecule for alignment! Exiting..."
          );
        }
        // read in scaffold molecules
        util::ShPtr< storage::Vector< chemistry::FragmentComplete >> sp_scaffold_pool_v;
        if( m_ScaffoldFragmentsFlag->GetFlag())
        {
          filenames.Reset();
          filenames = m_ScaffoldFragmentsFlag->GetStringList();
          util::ShPtr< chemistry::FragmentEnsemble> sp_scaffold_pool( new chemistry::FragmentEnsemble);

          for
          (
              storage::Vector< std::string>::const_iterator
              itr( filenames.Begin()), itr_end( filenames.End());
              itr != itr_end;
              ++itr
          )
          {
            // read in grow fragments ensemble
            io::File::MustOpenIFStream( input, *itr);
            sp_scaffold_pool->ReadMoreFromMdl( input, sdf::e_Maintain);
            io::File::CloseClearFStream( input);
          }
          storage::Vector< chemistry::FragmentComplete> scaffold_pool_v;
          for
          (
              auto itr( sp_scaffold_pool->Begin()), itr_end( sp_scaffold_pool->End());
              itr != itr_end;
              ++itr
          )
          {
            // remove SampleByParts MDL property if desired
            if( m_RemoveSampleByPartsFlag->GetFlag())
            {
              itr->GetStoredPropertiesNonConst().RemoveProperty( std::string( "SampleByParts"));
            }
            scaffold_pool_v.PushBack( *itr);
          }
          sp_scaffold_pool_v = util::ShPtr< storage::Vector< chemistry::FragmentComplete >>( new storage::Vector< chemistry::FragmentComplete>( scaffold_pool_v));
        }

        // close input stream for molecules
        io::File::CloseClearFStream( input);

      /////////////////////////
      // parse the arguments //
      /////////////////////////

        // try to read cheminfo property scorer
        descriptor::CheminfoProperty property_scorer;
        std::string mdl_property, receptor_filename;
        if( m_PoseDependentScoringFunctionFlag->GetFlag())
        {
          property_scorer = m_PoseDependentScoringFunctionFlag->GetParameterList()( 0)->GetValue();
          mdl_property = m_PoseDependentScoringFunctionFlag->GetParameterList()( 1)->GetValue();
          receptor_filename = m_PoseDependentScoringFunctionFlag->GetParameterList()( 2)->GetValue();
        }
        else
        {
          // flat energy surface
          property_scorer = "Constant(0.0)";
          mdl_property = "";
          receptor_filename = "";
        }

        // set substructure comparison method
        chemistry::ConformationGraphConverter::AtomComparisonTypeEnum atom_data;
        chemistry::ConfigurationalBondTypeData::DataEnum bond_data;
        if( m_AtomComparisonFlag->GetFlag())
        {
          atom_data = chemistry::ConformationGraphConverter::AtomComparisonTypeEnum
              (
                m_AtomComparisonFlag->GetFirstParameter()->GetValue()
              );
        }
        else
        {
          atom_data = chemistry::ConformationGraphConverter::e_ElementType;
        }
        if( m_BondComparisonFlag->GetFlag())
        {
          bond_data = chemistry::ConfigurationalBondTypeData::DataEnum
              (
                m_BondComparisonFlag->GetFirstParameter()->GetValue()
              );
        }
        else
        {
          bond_data = chemistry::ConfigurationalBondTypeData::e_FuzzyBondOrderAmideOrAromaticWithRingness;
        }

        // set up substructure alignment object
        chemistry::FragmentAlignToScaffold ats
        (
          atom_data,
          bond_data,
          size_t( 3)
        );

        // sample conformations
        chemistry::SampleConformations sample_confs;
        sample_confs.TryRead( m_SampleConfsFlag->GetFirstParameter()->GetValue(), util::GetLogger());

        // molalign
        chemistry::ConformationComparisonPsiField molalign_rigid;
        molalign_rigid.TryRead( m_MolAlignRigidFlag->GetFirstParameter()->GetValue(), util::GetLogger());
        chemistry::ConformationComparisonPsiFlexField molalign_flex;
        molalign_flex.TryRead( m_MolAlignFlexFlag->GetFirstParameter()->GetValue(), util::GetLogger());

        // get the local sampling preferences
        storage::Vector< size_t> local_sampling_preferences;
        if( m_LocalSamplingPrefsFlag->GetFlag())
        {
          for
          (
              size_t i( 0), sz( m_LocalSamplingPrefsFlag->GetParameterList().GetSize());
              i < sz;
              ++i
          )
          {
            local_sampling_preferences.PushBack
            (
              m_LocalSamplingPrefsFlag->GetParameterList()( i)->GetNumericalValue< size_t>()
            );
          }
        }
        else
        {
          // only sample bond angles/lengths as default
          local_sampling_preferences = storage::Vector< size_t>( size_t( 4), size_t( 0));
          local_sampling_preferences( 2) = size_t( 1);
        }

        // set alignment refinement options
        bool refine_alignment
        (
          m_RefineAlignmentFlag->GetFlag() ? true : false
        );

        // set whether to refine input pose only
        bool skip_initial_alignment
        (
          m_SkipInitialAlignmentFlag->GetFlag() ? true : false
        );

        // set comparer for alignments
        util::Implementation< chemistry::ConformationComparisonInterface> comparer
        (
          m_ConformerComparerFlag->GetFlag() ?
          m_ConformerComparerFlag->GetFirstParameter()->GetValue() :
          m_ConformerComparerFlag->GetFirstParameter()->GetDefaultValue()
        );

        // set optimization goal
        opti::Tracker< chemistry::FragmentComplete, double> opti_goal( opti::e_SmallerIsBetter);
        if( m_LargerScoreBetter->GetFlag())
        {
          opti_goal = opti::e_LargerIsBetter;
        }

      //////////////////////////////
      //  Prepare thread manager  //
      //////////////////////////////

        // Start track time
        util::Stopwatch threadmanager_timer( "Molecule Fitting", util::Time( 1, 0), util::Message::e_Standard, true, false);
        threadmanager_timer.Start();

        // Run multi-threading molecule fitting (alignment, docking, etc.)
        ThreadManager thread_manager
        (
          sp_inputs_v,
          m_OutputFilenameFlag->GetFirstParameter()->GetValue(),
          sched::GetNumberCPUs(),
          sp_scaffold_pool_v,
          property_scorer,
          mdl_property,
          receptor_filename,
          m_RoutineFlag->GetFirstParameter()->GetValue(),
          m_PoseSamplingIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_PoseSamplingCyclesFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_RefinementIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_MetropolisTemperatureFlag->GetFirstParameter()->GetNumericalValue< float>(),
          ats,
          molalign_rigid,
          molalign_flex,
          m_AlignmentSolutionsPerScaffoldFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          opti_goal,
          sample_confs,
          local_sampling_preferences,
          comparer,
          refine_alignment,
          skip_initial_alignment
        );

        // End track time
        threadmanager_timer.Stop();
        return 0;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      // static application instance
      static const ApplicationType MoleculeFit_Instance;

    }; // MoleculeFit

    //! @brief standard constructor
    MoleculeFit::MoleculeFit() :

      /////////////////////
      ///               ///
      ///  Generic I/O  ///
      ///               ///
      /////////////////////
      m_InputFilenamesFlag
      (
        new command::FlagStatic
        (
          "input_filenames",
          "molecules to fit",
          command::Parameter
          (
            "input_filenames",
            "filename(s) for input SDF",
            ""
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_filename",
          "file to which the molecules will be output after fitting",
          command::Parameter
          (
            "output_filename",
            "filename for output SDF of molecules",
            ""
          )
        )
      ),
      m_RoutineFlag
      (
        new command::FlagDynamic
        (
          "routine", "the molecule fit routine to perform",
          command::Parameter
          (
            "routine",
            "0 - MolAlign - performs property-based small molecule flexible alignment of the input molecules to the scaffold molecules and returns the pose of the single best alignment \n"
            "1 - MolAlign with structure-based refinement - performs MolAlign routine followed by local refinement of protein-ligand interactions with PLC-DNN and MCM-guided perturbations \n"
            "2 - Maximum common substructure alignment - enumerates MCS alignments to scaffolds with a conformational ensemble; chooses the best pose by MolAlign score \n"
            "3 - Pose-sensitive maximum common substructure alignment - enumerates MCS alignments to scaffolds with a conformational ensemble; chooses the best pose by PLC-DNN score \n"
            "4 - Protein-ligand refinement - performs refinement of protein-ligand complex using MCM-guided perturbations and PLC-DNN for pose-dependent scoring \n"
            "5 - Protein-ligand docking - docks a flexible small molecule to a rigid protein receptor using MCM-guided perturbations and PLC-DNN for pose-dependent scoring",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "0", "1", "2", "3", "4", "5")),
            "2"
          )
        )
      ),

      /////////////////////
      ///               ///
      ///  Alignment    ///
      ///               ///
      /////////////////////

      m_ScaffoldFragmentsFlag
      (
        new command::FlagStatic
        (
          "scaffold_fragments",
          "molecules against which to align the input molecules",
          command::Parameter
          (
            "scaffold_fragments",
            "filename(s) for input SDF of scaffold molecules",
            ""
          )
        )
      ),

      m_AtomComparisonFlag
      (
        new command::FlagStatic
        (
          "atom_comparison_type",
          "atom data compared to determine whether or not atoms match",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable
            (
              chemistry::ConformationGraphConverter::AtomComparisonTypeEnum()
            ),
            "ElementType"
          )
        )
      ),
      m_BondComparisonFlag
      (
        new command::FlagStatic
        (
          "bond_comparison_type",
          "bond data compared to determine whether or not bonds match",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable
            (
              chemistry::ConfigurationalBondTypeData::DataEnum()
            ),
            "FuzzyBondOrderAmideOrAromaticWithRingness"
          )
        )
      ),
      m_MolAlignRigidFlag
      (
        new command::FlagStatic
        (
          "molalign_rigid",
          "settings for small molecule rigid-conformer property-based alignment",
          command::Parameter
          (
            "",
            "",
            "("
            "properties="
            "Combine("
            "Define(ChargeNegative=Multiply(Less(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
            "Define(ChargePositive=Multiply(Greater(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
            "Define(HbondAcceptorsStrict=Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))),"
            "Define(ENegOffset=Subtract(lhs=Atom_ElectroNegativity,rhs=Constant(2.5))),"
            "Define(IsENeg=Greater(lhs=ENegOffset,rhs=Constant(0))),"
            "Define(PolarTernary=Subtract(lhs=Multiply(Add(HbondAcceptorsStrict,Atom_HbondDonors),Constant(2)),rhs=Constant(1))),"
            "Define(Atom_Hydrophobic=Not(DescriptorSum(Abs(Partial(indices(1,2,5,8),2DASign(property=PolarTernary,steps=3)))))),"
            "Atom_SigmaCharge,HbondAcceptorsStrict,Atom_HbondDonors,Atom_Polarizability,Atom_AromaticityAxes,Multiply(IsENeg,ENegOffset),Atom_Hydrophobic,Atom_VDWVolume),"
            "property_weights(5,3.5,7,1.71,2.41,2.41,2.41,3.0,2.0,0.714),"
            "max_atom_distance_criterion=1,"
            "optimizing_weights=0,"
            "anchor_weight=0.01,"
            "linear_mismatch_penalty=0.01,"
            "heavy_mismatch_penalty_fraction=0.6,"
            "heavy_mismatch_penalty=2,"
            "iterations=400,"
            "max_unimproved_steps=160,"
            "number_rigid_trajectories=5,"
            "number_outputs=1,"
            "align_to_scaffold=0,"
            "initial_rand_rotation=0,"
      "exclusion_indices_a="",exclusion_indices_b="",pose_tolerance=0.125,"
            "pose_score_threshold=2,"
            "flip_prob=0.06,"
            "big_rot_prob=0.06,"
            "bond_swap_prob=0.06,"
            "bond_align_prob=0.3,"
            "bond_align3_prob=0.2,"
            "bond_align_inf_prob=0.1,"
            "small_rot_prob=0.12,"
            "small_trans_prob=0.12,"
            "conf_swap_prob=0.06"
            ")"
          )
        )
      ),
      m_MolAlignFlexFlag
      (
        new command::FlagStatic
        (
          "molalign_flex",
          "settings for small molecule flexible property-based alignment",
          command::Parameter
          (
            "",
            "",
            "("
            "properties="
            "Combine("
            "Define(ChargeNegative=Multiply(Less(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
            "Define(ChargePositive=Multiply(Greater(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
            "Define(HbondAcceptorsStrict=Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))),"
            "Define(ENegOffset=Subtract(lhs=Atom_ElectroNegativity,rhs=Constant(2.5))),"
            "Define(IsENeg=Greater(lhs=ENegOffset,rhs=Constant(0))),"
            "Define(PolarTernary=Subtract(lhs=Multiply(Add(HbondAcceptorsStrict,Atom_HbondDonors),Constant(2)),rhs=Constant(1))),"
            "Define(Atom_Hydrophobic=Not(DescriptorSum(Abs(Partial(indices(1,2,5,8),2DASign(property=PolarTernary,steps=3)))))),"
            "Atom_SigmaCharge,HbondAcceptorsStrict,Atom_HbondDonors,Atom_Polarizability,Atom_AromaticityAxes,Multiply(IsENeg,ENegOffset),Atom_Hydrophobic,Atom_VDWVolume),"
            "property_weights(5,3.5,7,1.71,2.41,2.41,2.41,3.0,2.0,0.714),"
            "max_atom_distance_criterion=1,"
            "optimizing_weights=0,"
            "anchor_weight=0.01,"
            "linear_mismatch_penalty=0.01,"
            "heavy_mismatch_penalty_fraction=0.6,"
            "heavy_mismatch_penalty=2,"
            "iterations=400,"
            "max_unimproved_steps=160,"
            "number_rigid_trajectories=5,"
            "number_outputs=1,"
            "align_to_scaffold=0,"
            "initial_rand_rotation=0,"
      "exclusion_indices_a="",exclusion_indices_b="",pose_tolerance=0.125,"
            "pose_score_threshold=2,"
            "flip_prob=0.06,"
            "big_rot_prob=0.06,"
            "bond_swap_prob=0.06,"
            "bond_align_prob=0.3,"
            "bond_align3_prob=0.2,"
            "bond_align_inf_prob=0.1,"
            "small_rot_prob=0.12,"
            "small_trans_prob=0.12,"
            "conf_swap_prob=0.06,"
            "rigid_mol_b=1,"
            "conformer_pairs=100,"
            "top_pairs=3,"
            "filter_iterations=600,"
            "filter_limit=240,"
            "refinement_iterations=200,"
            "refinement_limit=80,"
            "number_flexible_trajectories=5,"
            "fraction_filtered_initially=0.25,"
            "fraction_filtered_iteratively=0.50"
            ")"
          )
        )
      ),

//      m_MolAlignIterationsFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_iterations",
//          "number of iterations of in external loop of property-based alignment",
//          command::Parameter
//          (
//            "iterations",
//            "number of molalign iterations",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "200"
//          )
//        )
//      ),
//      m_MolAlignMaxUnimprovedFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_max_unimproved",
//          "termination criterion - stops alignment early after this number consecutive iterations "
//          "in external loop of property-based alignment fails to improve score",
//          command::Parameter
//          (
//            "max_unimproved",
//            "max number of unimproved molalign iterations",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "200"
//          )
//        )
//      ),
//      m_MolAlignFilterIterationsFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_filter_iterations",
//          "number of iterations of in repeating filter loop of property-based alignment",
//          command::Parameter
//          (
//            "refinement_iterations",
//            "number of molalign filter iterations",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "150"
//          )
//        )
//      ),
//      m_MolAlignFilterMaxUnimprovedFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_filter_max_unimproved",
//          "termination criterion - stops alignment early after this number consecutive iterations "
//          "in repeating filter loop of property-based alignment fails to improve score",
//          command::Parameter
//          (
//            "filter_max_unimproved",
//            "max number of unimproved molalign refinement iterations",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "150"
//          )
//        )
//      ),
//      m_MolAlignRefinementIterationsFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_refinement_iterations",
//          "number of iterations of in final refinement of property-based alignment",
//          command::Parameter
//          (
//            "refinement_iterations",
//            "number of molalign refinement iterations",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "80"
//          )
//        )
//      ),
//      m_MolAlignRefinementMaxUnimprovedFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_refinement_max_unimproved",
//          "termination criterion - stops alignment early after this number consecutive iterations "
//          "in final refinement of property-based alignment fails to improve score",
//          command::Parameter
//          (
//            "refinement_max_unimproved",
//            "max number of unimproved molalign refinement iterations",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "80"
//          )
//        )
//      ),
//      m_MolAlignRigidTrajFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_rigid_traj",
//          "number of independent molalign rigid alignment trajectories with linearly "
//          "scaling variable max atom distances",
//          command::Parameter
//          (
//            "molalign_rigid_traj",
//            "number of molalign rigid alignment trajectories",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "3"
//          )
//        )
//      ),
//      m_MolAlignFlexTrajFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_flex_traj",
//          "number of independent molalign flexible alignment trajectories with linearly "
//          "scaling variable max atom distances",
//          command::Parameter
//          (
//            "molalign_flex_traj",
//            "number of molalign flex alignment trajectories",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "4"
//          )
//        )
//      ),
//      m_MolAlignMaxAtomDistanceUpperFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_max_atom_dist_upper",
//          "upper limit on the maximum atom distance for molalign mutually matching atom pair detection; "
//          "default value should be acceptable",
//          command::Parameter
//          (
//            "molalign_max_atom_dist_upper",
//            "upper limit max atom distance",
//            command::ParameterCheckRanged< size_t>( 0.0, std::numeric_limits< float>::max()),
//            "1.15"
//          )
//        )
//      ),
//      m_MolAlignMaxAtomDistanceLowerFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_max_atom_dist_lower",
//          "lower limit on the maximum atom distance for molalign mutually matching atom pair detection; "
//          "default value should be acceptable",
//          command::Parameter
//          (
//            "molalign_max_atom_dist_lower",
//            "lower limit max atom distance",
//            command::ParameterCheckRanged< size_t>( 0.0, std::numeric_limits< float>::max()),
//            "0.70"
//          )
//        )
//      ),
//      m_MolAlignPoseToleranceFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_pose_tolerance",
//          "conformation comparison (of type SymmetryRealSpaceRMSD) tolerance determining whether a pose is unique",
//          command::Parameter
//          (
//            "tolerance",
//            "pose discrimination tolerance",
//            command::ParameterCheckRanged< size_t>( 0.0, std::numeric_limits< float>::max()),
//            "0.125"
//          )
//        )
//      ),
//      m_MolAlignPoseThresholdFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_pose_threshold",
//          "upper limit on acceptable alignment score determining whether a pose is a solution",
//          command::Parameter
//          (
//            "threshold",
//            "pose acceptance threshold",
//            command::ParameterCheckRanged< size_t>( 0.0, std::numeric_limits< float>::max()),
//            "2.0"
//          )
//        )
//      ),
//      m_MolAlignConfPairsFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_conf_pairs",
//          "number of conformer pairs with which the alignment simulation is initialized",
//          command::Parameter
//          (
//            "conformer pairs",
//            "pairs of conformers of the molecules being aligned",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "250"
//          )
//        )
//      ),
//      m_MolAlignBestPairsFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_best_pairs",
//          "number of aligned poses to incluide in the final refinement for comparison",
//          command::Parameter
//          (
//            "best pairs",
//            "set of alignment pairs compared at the end",
//            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
//            "1"
//          )
//        )
//      ),
//      m_MolAlignFractionFilterInitialFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_filter_fraction_initial",
//          "fraction of initial conformer pair alignments from the first round of flexible alignments "
//          "that will move on to the iterative inner loop round of alignments",
//          command::Parameter
//          (
//            "filter fraction",
//            "fraction of initial conformer pair alignments to keep after external alignment loop",
//            command::ParameterCheckRanged< size_t>( 0.0, 1.0),
//            "0.25"
//          )
//        )
//      ),
//    m_MolAlignFractionFilterIterativeFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_filter_fraction_iterative",
//          "fraction of conformer pair alignments that are retained after each round of filter iterations "
//          "during flexible alignment; continues until 'molalign_best_pairs' candidate alignments remain",
//          command::Parameter
//          (
//            "filter fraction",
//            "fraction of conformer pair alignments to keep after each internal alignment loop",
//            command::ParameterCheckRanged< size_t>( 0.0, 1.0),
//            "0.50"
//          )
//        )
//      ),
//      m_MolAlignRigidScaffoldFlag
//      (
//        new command::FlagStatic
//        (
//          "molalign_rigid_scaffold",
//          "keep scaffold molecule rigid (do not sample conformers)"
//        )
//      ),
      m_RefineAlignmentFlag
      (
        new command::FlagStatic
        (
          "refine_alignment",
          "perform an alignment with local conformers beginning with the pose obtained from an alignment"
        )
      ),
      m_SkipInitialAlignmentFlag
      (
        new command::FlagStatic
        (
          "skip_initial_alignment",
          "perform an alignment with local conformers beginning with the input molecule pose(s)"
        )
      ),
      m_ConformerComparerFlag
      (
        new command::FlagStatic
        (
          "alignment_score",
          "score minimized during small molecule substructure-based alignment; "
          "default is PropertyFieldDistance using the MolAlign weights",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< chemistry::ConformationComparisonInterface>()),
            "PropertyFieldDistance"
          )
        )
      ),

      /////////////////////
      ///               ///
      ///  SampleConfs  ///
      ///               ///
      /////////////////////
      m_SampleConfsFlag
      (
        new command::FlagStatic
        (
          "sample_confs",
          "settings for small molecule conformer generation",
          command::Parameter
          (
            "",
            "",
            "("
            "conformation_comparer=SymmetryRMSD,"
            "tolerance=0.25,"
            "generate_3D=0,"
            "cluster=true,"
            "max_iterations=2000,"
            "max_conformations=200,"
            "change_chirality=0"
            ")"
          )
        )
      ),
      m_LocalSamplingPrefsFlag
      (
        new command::FlagDynamic
        (
          "local_sampling_preferences",
          "restrictions to impose on conformational sampling when generating local conformational ensembles",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "dihedrals",
              "sample dihedrals?",
              "0"
            ),
            command::Parameter
            (
              "rings",
              "sample rings?",
              "0"
            ),
            command::Parameter
            (
              "bonds",
              "sample bond lengths/angles?",
              "1"
            ),
            command::Parameter
            (
              "chirality",
              "sample chirality? Can also be controlled with 'sample_confs' flag; generally recommend 'false' "
              "unless you truly do not know the chirality and/or have a very specific use-case.",
              "0"
            )
          )
        )
      ),
      m_RemoveSampleByPartsFlag
      (
        new command::FlagStatic
        (
          "remove_sample_by_parts",
          "remove any SampleByParts MDL property labels from all molecules prior to conformer generation"
        )
      ),

      /////////////////////
      ///               ///
      /// P-L Interface ///
      ///               ///
      /////////////////////

      m_PoseDependentScoringFunctionFlag
      (
        new command::FlagDynamic
        (
          "pose_dependent_scoring_function",
          "the scoring function to use for protein-ligand interactions",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "function",
              "the scoring function implementation to use for protein-ligand interactions",
              command::ParameterCheckSerializable
              (
                chemistry::ScoreFunctionGeneric()
              )
            ),
            command::Parameter
            (
              "mdl_property",
              "MDL property specifying the PDB filename for the receptor; "
              "needs to match MDL property name used to train the corresponding machine learning model; "
              "if no mdl property is required then pass an empty string"
            ),
            command::Parameter
            (
              "receptor_filename",
              "PDB filename specifying the protein with which to score interactions; "
              "if no receptor is required then pass an empty string"
            )
          )
        )
      ),
      m_PoseSamplingIterationsFlag
      (
        new command::FlagStatic
        (
          "pose_sampling_iterations",
          "number of iterations of pose sampling to perform;"
          " NOTE - for local pose-dependent docking, the number of iterations that will be performed"
          " is set by 'refinement_pose_sampling_iterations' instead",
          command::Parameter
          (
            "iterations",
            "number of sampling iterations",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "20"
          )
        )
      ),
      m_PoseSamplingCyclesFlag
      (
        new command::FlagStatic
        (
          "n_cycles",
          "number of cycles of MCM sampling to do during pose optimization",
          command::Parameter
          (
            "cycles",
            "number of MCM runs",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "5"
          )
        )
      ),
      m_RefinementIterationsFlag
      (
        new command::FlagStatic
        (
          "refinement_pose_sampling_iterations",
          "number of iterations of pose sampling to perform during refinement;"
          " sampling restricted to local conformers and small - medium perturbations; "
          " NOTE - for local pose-dependent docking, this flag sets the number of iterations that"
          " will be done, not 'pose_sampling_iterations'",
          command::Parameter
          (
            "refinement_iterations",
            "number of sampling iterations",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "20"
          )
        )
      ),
      m_AlignmentSolutionsPerScaffoldFlag
      (
        new command::FlagStatic
        (
          "alignment_solutions",
          "max number of alignment solutions per scaffold to move forward toward pose-dependent refinement; "
          "only applicable if pose-dependent scoring refinement is enabled",
          command::Parameter
          (
            "alignment_solutions",
            "max number of alignment solutions per scaffold",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "1"
          )
        )
      ),
      m_MetropolisTemperatureFlag
      (
        new command::FlagStatic
        (
          "temperature", "flag for the temperature used in the Metropolis criterion;"
          " units match the units of the score function",
          command::Parameter
          (
            "temperature", "temperature for MCM evaluation",
            command::ParameterCheckRanged< float>( 0.00001, std::numeric_limits< float>::max()), "0.6"
          )
        )
      ),
      m_LargerScoreBetter
      (
        new command::FlagStatic
        (
          "larger_score_better",
          "sets all pose-dependent scoring objects to increase the score during optimization; "
          "default behavior attempts to decrease the score during optimization"
        )
      )
    {
    }

    const ApplicationType MoleculeFit::MoleculeFit_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeFit(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
