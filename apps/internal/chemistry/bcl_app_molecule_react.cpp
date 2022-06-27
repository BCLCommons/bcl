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
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
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
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_react.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_stochastic_pose_optimizer.h"
#include "chemistry/bcl_chemistry_reaction_search.h"
#include "chemistry/bcl_chemistry_reaction_worker.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_score_function_generic.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
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
#include "io/bcl_io_directory_entry.h"
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
    //! @class MoleculeReact
    //! @brief Application for performing chemical reactions to generate new molecules
    //!
    //! @author brownbp1
    //! @date 11/21/2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeReact : public InterfaceRelease
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
      //! @date 11/21/2021
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class ThreadManager : public util::ObjectInterface
      {
      ///////////
      // Data //
      ///////////

      private:

        // For input
        const size_t                                        m_Threads;                    // Number of threads
        const size_t                                        m_NumberInputFragments;       // Number of input fragments
        const size_t                                        m_NumberRepeatCycles;         // Number of repeat cycles to perform
        const size_t                                        m_NumberTotalJobs;            // Total number of jobs
        size_t                                              m_NumberCompletedJobs;        // Total number of completed jobs
        size_t                                              m_CurrentMolIndex;            // Index of the molecule to fit
        size_t                                              m_CurrentCycleRepeat;         // Index of repeat cycle
        storage::Vector< chemistry::FragmentComplete>       m_Molecules;                  // The molecules which have been built and are ready for output
        chemistry::SampleConformations                      m_SampleConfs;                // Sample conformations object

        // For output
        io::OFStream                                        m_OutputStream;               // Output file to write molecules to
        sched::Mutex                                        m_Mutex;                      // Lock for updating Workers
        bool                                                m_Corina;                     // Corina 3D conformer?

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //!
        //! @class Worker
        //! @brief Runs the threads, though threads are probably only useful for Exhaustive mode when using multiple
        //! starting fragments. Otherwise, the slow step is mutex locked and you may as well just use one thread.
        //!
        //! @author brownbp1
        //! @date 11/21/2021
        //!
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        struct Worker
        {

          // Basic data for general use
          util::SiPtr< ThreadManager>                                   m_ThreadManager;             // Pointer to the parent thread manager; needed to update Worker
          std::string                                                   m_WorkerRoutine;             // The fit routine to perform
          util::ShPtr< storage::Vector< chemistry::FragmentComplete>>   m_InputFragments;            // Input molecules
          size_t                                                        m_CurrentWorkerMolIndex;     // Index of the molecule to fit
          size_t                                                        m_CurrentWorkerRepeatCycle;  // Index of repeat cycle

          // Reaction
          chemistry::FragmentReact                                      m_React;                     // The reaction operation class, used for modifying structures
          chemistry::ReactionSearch                                     m_ReactionSearch;            // The reaction operation class, used for modifying structures
          chemistry::ReactionWorker                                     m_ReactionWorker;            // Heavy lifting for the reaction
          util::ShPtr< chemistry::FragmentEnsemble>                     m_Reagents;                  // The reagents to use in the reactions
          std::string                                                   m_ReactionsDirectory;        // The directory containing the RXN files
          bool                                                          m_Initialized;               // Whether the reaction has been prepared
          storage::Vector< size_t>                                      m_TargetReactantPositions;   // The reactant positions available to input fragments
          bool                                                          m_LigandBased;               // Do not preserve any coordinate information during reaction
          bool                                                          m_FixGeometry;               // Fix bad geometry atoms
          bool                                                          m_FixRingGeometry;           // Fix bad ring geometry atoms
          size_t                                                        m_AdjacentAtoms;             // Adjacent atoms in any direction to add to mobile selection
          chemistry::FragmentMapConformer                               m_Cleaner;                   // For final 3D conformer generation
          chemistry::SampleConformations                                m_WorkerSampleConfs;         // Sample conformations object

        //////////////
        // operator //
        //////////////

          //! @brief Run the thread to perform the chosen fit operations
          // TODO: double check that success rate around 50% for routine Random is actually expected and not a borked feature of app
          void RunThread()
          {
            // run threads
            do
            {
              // Initialization is required
              if( !m_Initialized)
              {
                PrepareReaction
                (
                  ( *m_InputFragments)( m_CurrentWorkerMolIndex),
                  *m_Reagents,
                  m_ReactionsDirectory
                );
              }

              // random
              if( m_WorkerRoutine == "Random")
              {
                // react
                BCL_MessageStd( "Performing a random reaction!");
                auto products( m_React.ReactRandom( ( *m_InputFragments)( m_CurrentWorkerMolIndex)));

                // add products to molecule vector and increment index
                util::ShPtr< chemistry::FragmentComplete> clean_prod;
                for
                (
                    auto prod_itr( products.Second().Begin()), prod_itr_end( products.Second().End());
                    prod_itr != prod_itr_end;
                    ++prod_itr
                )
                {
                  // corina pose-independent conformer
                  if( m_ThreadManager->GetIsCorina())
                  {
                    BCL_MessageStd( "Getting corina conformer!");
                    clean_prod = util::ShPtr< chemistry::FragmentComplete>( GetCorina3DCoordinates( *prod_itr));
                  }
                  // bcl pose-independent conformer
                  else if( m_React.GetReactionWorker().GetIsProductConformerArbitrary())
                  {
                    BCL_MessageStd( "Getting pose-independent BCL conformer!");
                    Clean( *prod_itr);
                    chemistry::FragmentEnsemble prod_confs( GenerateConformers( *prod_itr));
                    if( prod_confs.GetSize())
                    {
                      clean_prod = util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete( prod_confs.GetMolecules().FirstElement()));
                    }
                  }
                  // keep pose from ReactionWorker, which tries to preserve coordinate information for the
                  // primary reactant (unless ligand-based is specified, in which case it operates more or less as
                  // legacy mode from Alex), but that goes to the else-if block above this
                  else
                  {
                    BCL_MessageStd( "Getting pose-dependent BCL conformer!");
                    if( prod_itr->GetSize() && !prod_itr->HasBadGeometry())
                    {
                      clean_prod = util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete( *prod_itr));
                    }
                  }
                  if( clean_prod.IsDefined())
                  {
                    m_ThreadManager->m_Mutex.Lock();
                    m_ThreadManager->AddMolecule( *clean_prod);
                    m_ThreadManager->m_Mutex.Unlock();
                  }
                  else
                  {
                    BCL_MessageStd( "Could not add molecule to final products - failed to generate valid 3D conformer");
                  }
                }
                m_ThreadManager->m_Mutex.Lock();
                m_ThreadManager->UpdateWorker( *this);
                m_ThreadManager->m_Mutex.Unlock();
              }
              // enumerative
              else if( m_WorkerRoutine == "Exhaustive")
              {
                // react
                BCL_MessageStd( "Performing exhaustive reactions!");
                auto products( m_React.ReactExhaustive( ( *m_InputFragments)( m_CurrentWorkerMolIndex)));

                // add products to molecule vector and increment index
                for
                (
                    auto rxn_itr( products.First().Begin()), rxn_itr_end( products.First().End());
                    rxn_itr != rxn_itr_end;
                    ++rxn_itr
                )
                {
                  util::ShPtr< chemistry::FragmentComplete> clean_prod;
                  for
                  (
                      auto prod_itr( rxn_itr->Second().Begin()), prod_itr_end( rxn_itr->Second().End());
                      prod_itr != prod_itr_end;
                      ++prod_itr
                  )
                  {
                    // corina pose-independent conformer
                    if( m_ThreadManager->GetIsCorina())
                    {
                      BCL_MessageStd( "Getting corina conformer!");
                      clean_prod = util::ShPtr< chemistry::FragmentComplete>( GetCorina3DCoordinates( *prod_itr));
                    }
                    // bcl pose-independent conformer
                    else if( m_React.GetReactionWorker().GetIsProductConformerArbitrary())
                    {
                      BCL_MessageStd( "Getting pose-independent BCL conformer!");
                      Clean( *prod_itr);
                      chemistry::FragmentEnsemble prod_confs( GenerateConformers( *prod_itr));
                      if( prod_confs.GetSize())
                      {
                        clean_prod = util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete( prod_confs.GetMolecules().FirstElement()));
                      }
                    }
                    // keep pose from ReactionWorker, which tries to preserve coordinate information for the
                    // primary reactant (unless ligand-based is specified, in which case it operates more or less as
                    // legacy mode from Alex), but that goes to the else-if block above this
                    else
                    {
                      BCL_MessageStd( "Getting pose-dependent BCL conformer!");
                      if( prod_itr->GetSize() && !prod_itr->HasBadGeometry())
                      {
                        clean_prod = util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete( *prod_itr));
                      }
                    }
                    if( clean_prod.IsDefined())
                    {
                      m_ThreadManager->m_Mutex.Lock();
                      m_ThreadManager->AddMolecule( *clean_prod);
                      m_ThreadManager->m_Mutex.Unlock();
                    }
                    else
                    {
                      BCL_MessageStd( "Could not add molecule to final products - failed to generate valid 3D conformer");
                    }
                  }
                }
                m_ThreadManager->m_Mutex.Lock();
                m_ThreadManager->UpdateWorker( *this);
                m_ThreadManager->m_Mutex.Unlock();
              }
              // bail
              else
              {
                BCL_Assert( false, "No valid routine selected. Exiting...");
              }

              // update progress
              if
              (
                  m_ThreadManager->m_NumberCompletedJobs <
                  m_ThreadManager->GetNumberInputFragments() * ( m_ThreadManager->GetNumberRepeatCycles() + 1)
              )
              {
                util::GetLogger().LogStatus
                (
                  "% complete: " +
                  util::Format().FFP( 2)( float( m_ThreadManager->m_NumberCompletedJobs) / float( m_ThreadManager->m_NumberTotalJobs) * 100)
                );
              }
            } while
              (
                  m_ThreadManager->m_NumberCompletedJobs <
                  m_ThreadManager->GetNumberInputFragments() * ( m_ThreadManager->GetNumberRepeatCycles() + 1)
              );
          } // RunThread()

        ////////////////
        // operations //
        ////////////////

          //! @brief Creates a reaction search object for the current molecule
          //! @STARTNG_MOLECULE starting molecule with which to react reagents
          //! @REAGENTS reagents to react
          //! @REACTIONS_DIRECTORY directory containing RXN files
          void PrepareReaction
          (
            const chemistry::FragmentComplete &STARTING_MOLECULE,
            const chemistry::FragmentEnsemble &REAGENTS,
            const std::string &REACTIONS_DIRECTORY
          )
          {
            // add the starting reagent to the pool of reagents
            chemistry::FragmentEnsemble reagents( storage::List< chemistry::FragmentComplete>( 1, STARTING_MOLECULE));
            reagents.Append( REAGENTS);

            // prepare and return reaction search object
            m_ReactionSearch = chemistry::ReactionSearch( reagents, REACTIONS_DIRECTORY);
            m_ReactionSearch.Initialize();

            // prepare the react object and corresponding reaction worker
            m_React = chemistry::FragmentReact( m_ReactionWorker, m_ReactionSearch);
            if( !m_LigandBased && !m_ThreadManager->GetIsCorina())
            {
              BCL_MessageStd( "Setting pose-dependent options");
              m_React.GetReactionWorker().SetProductConformerArbitrary( m_LigandBased);
              m_React.GetReactionWorker().SetCorrectGeometry( m_FixGeometry);
              m_React.GetReactionWorker().SetCorrectNonReferenceRingGeometry( m_FixRingGeometry);
              m_React.GetReactionWorker().SetAdditionalAdjacentAtoms( m_AdjacentAtoms);
              BCL_MessageStd( "Ligand-based: " + util::Format()( m_LigandBased ? "true" : "false"));
              BCL_MessageStd( "Fix bad geometry: " + util::Format()( m_FixGeometry ? "true" : "false"));
              BCL_MessageStd("Fix bad ring geometry: " + util::Format()( m_FixRingGeometry ? "true" : "false"));
              BCL_MessageStd("Extend adjacent atoms: " + util::Format()( m_AdjacentAtoms));
            }
            m_Initialized = true;
          }

          void Clean( chemistry::FragmentComplete &MOLECULE)
          {
            // clean properties and add all atoms
            MOLECULE.RemoveProperty( "SampleByParts");
            MOLECULE.RemoveH();

            // clean at atom_vector level
            auto atom_vec( MOLECULE.GetAtomVector());
            chemistry::AtomsCompleteStandardizer standardizer( atom_vec, "", true);
            standardizer.SetConjugationOfBondTypes( atom_vec);
            chemistry::BondIsometryHandler::AddIsometryInformation( atom_vec, true);
            chemistry::StereocentersHandler::AddChiralityFromConformation( atom_vec);

            // reconstruct our clean mol
            MOLECULE = chemistry::FragmentComplete( atom_vec, MOLECULE.GetName());
            MOLECULE.SaturateWithH();
          }

          //! @brief Generate small molecule conformational ensemble
          //! @param MOLECULE molecule of interest
          chemistry::FragmentEnsemble GenerateConformers
          (
            const chemistry::FragmentComplete &MOLECULE
          )
          {
            // generate conformers
            chemistry::FragmentEnsemble ens;
            ens = m_WorkerSampleConfs( MOLECULE).First();

            if( ens.GetSize())
            {
              return ens;
            }
            else
            {
              // empty
              return chemistry::FragmentEnsemble();
            }
          }

          //! @brief system call to corina to generate 3d conformers
          //! @param MOLECULE the small molecule for which 3d coordinates will be generated
          //! @return the 3d generated small molecule or empty molecule if corina was unavailable
          util::ShPtr< chemistry::FragmentComplete> GetCorina3DCoordinates
          (
            const chemistry::ConformationInterface &MOLECULE
          )
          {

              #if !defined(__MINGW32__) && !defined(__WIN32__) && !defined(_WIN32)

            BCL_MessageVrb( "GetCorina3DCoordinates");

            if( MOLECULE.GetNumberAtoms() == 0)
            {
              BCL_MessageVrb( "GetCorina3DCoordinates: empty molecule as input");
              return util::ShPtr< chemistry::FragmentComplete>();
            }

            time_t cur_time;
            time( &cur_time);
            std::string file_basename( util::Format()( &MOLECULE) + util::Format()( cur_time));

            io::DirectoryEntry filename2d( "/tmp/" + file_basename + "_gen2D.sdf");
            io::DirectoryEntry filename3d( "/tmp/" + file_basename + "_gen3D.sdf");

            io::OFStream out;
            if( !io::File::TryOpenOFStream( out, filename2d.GetFullName()))
            {
              BCL_MessageCrt( "unable to open " + filename2d.GetFullName() + " for writing");
              return util::ShPtr< chemistry::FragmentComplete>();
            }
            MOLECULE.WriteMDL( out);
            io::File::CloseClearFStream( out);
            // generate meaningful coordinates for the ensemble of generated constitutions using corina until somebody implements code to do that into the bcl!!
            const std::string command( "corina -dwh " + filename2d.GetFullName() + " " + filename3d.GetFullName());
            const int error( system( command.c_str()));

            if( error != 0)
            {
              BCL_MessageCrt( "unable to execute command: " + command + " with error: " + util::Format()( error));
            }

            io::IFStream input;
            io::File::MustOpenIFStream( input, filename3d.GetFullName());
            util::ShPtr< chemistry::FragmentComplete> fragment
            (
              new chemistry::FragmentComplete( sdf::FragmentFactory::MakeFragment( input, sdf::e_Saturate))
            );
            io::File::CloseClearFStream( input);
            filename2d.Remove();
            filename3d.Remove();
            if( fragment.IsDefined() && !fragment->HasBadGeometry() && !fragment->HasNonGasteigerAtomTypes())
            {
              return fragment;
            }

            return util::ShPtr< chemistry::FragmentComplete>();

              #else // !defined(__MINGW32__) && !defined(__WIN32__) && !defined(_WIN32)
            return util::ShPtr< chemistry::FragmentComplete>();
              #endif
          }

        }; // struct Worker

      ////////////////
      // operations //
      ////////////////

        //! @brief Updates worker threads
        //! @param WORKER worker thread whose index is incremented
        void UpdateWorker( ThreadManager::Worker &WORKER)
        {
          // repeat the current index until we hit requested number
          if( m_CurrentCycleRepeat >= m_NumberRepeatCycles && m_CurrentMolIndex + 1 < m_NumberInputFragments)
          {
            // increase worker molecule index
            ++m_CurrentMolIndex;
            WORKER.m_CurrentWorkerMolIndex = m_CurrentMolIndex;
            WORKER.m_Initialized = false;

            // reset cycle count
            WORKER.m_CurrentWorkerRepeatCycle = 0;
            m_CurrentCycleRepeat = 0;
          }
          // increase the cycle count
          else if( m_CurrentCycleRepeat < m_NumberRepeatCycles)
          {
            WORKER.m_CurrentWorkerMolIndex = m_CurrentMolIndex;
            ++m_CurrentCycleRepeat;
            WORKER.m_CurrentWorkerRepeatCycle = m_CurrentCycleRepeat;
          }

          // either way increment the total job count
          ++m_NumberCompletedJobs;

        } // UpdateWorker

      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! brief constructor
        ThreadManager
        (
          const util::ShPtr< storage::Vector< chemistry::FragmentComplete> >   &INPUT_FRAGMENTS,       // Input fragments
          const util::ShPtr< chemistry::FragmentEnsemble>                     &REAGENTS,              // Input reagents
          const std::string                                                   &REACTIONS_DIRECTORY,   // Reactions directory
          const std::string                                                   &ROUTINE,               // Molecule react routine
          const size_t                                                         REPEATS,               // Reaction routine repeats
          const storage::Vector< size_t>                                      &TARGET_REACTION_POS,   // Target reactant reaction positions
          const bool                                                           LIGAND_BASED,          // If true, do not preserve coordinate info during reaction
          const chemistry::SampleConformations                                &SAMPLE_CONFS,          // Sample conformers object
          const bool                                                           FIX_GEOMETRY,          // Fix any bad atoms during reactions
          const bool                                                           FIX_RING_GEOMETRY,     // Fix any bad ring atoms during reaction
          const size_t                                                         EXTEND_ADJACENT,       // Number of atoms to extend from mobile into fixed region
          const std::string                                                   &OUTPUT_FILENAME,       // Output SDF filename
          const bool                                                           CORINA,                // Corina 3D conformer
          const size_t                                                         NUMBER_THREADS         // Number of threads (from scheduler)
        ) :
          m_Threads( std::min( NUMBER_THREADS, INPUT_FRAGMENTS->GetSize() * ( REPEATS + 1))),
          m_NumberInputFragments( INPUT_FRAGMENTS->GetSize()),
          m_NumberRepeatCycles( REPEATS),
          m_NumberTotalJobs( ( REPEATS + 1) * INPUT_FRAGMENTS->GetSize()),
          m_NumberCompletedJobs( 0),
          m_CurrentMolIndex
          (
            REPEATS < INPUT_FRAGMENTS->GetSize() ?
            std::min( NUMBER_THREADS, INPUT_FRAGMENTS->GetSize()) - 1 :
            0
          ),
          m_CurrentCycleRepeat( 0),
          m_SampleConfs( SAMPLE_CONFS),
          m_Corina( CORINA)
        {
          // do not allow people to do multiple cycles of enumeration - it is wasteful
          if( ROUTINE == "Exhaustive" && REPEATS)
          {
            BCL_ExitWithoutCallstack(
              "ERROR: Running the 'Exhaustive' routine with multiple repeats is not allowed; "
              "please either change to a stochastic routine or set repeats to 0.",
              -1
            );
          }

          // prepare output filestream
          io::DirectoryEntry entry( OUTPUT_FILENAME);
          if( entry.DoesExist())
          {
            entry.Remove();
          }
          io::File::MustOpenOFStream( m_OutputStream, OUTPUT_FILENAME, std::ios::app);

          // initialize one worker per thread
          std::vector< Worker> workers( m_Threads);
          size_t current_mol_index( 0), current_cycle_index( 0);
          for
          (
              std::vector< Worker>::iterator itr( workers.begin()), end( workers.end());
              itr != end;
              ++itr
          )
          {
            // setup worker
            Worker &worker_ref( *itr);
            worker_ref.m_ThreadManager                    = this;                          // Copy of pointer to thread manager
            worker_ref.m_InputFragments                   = INPUT_FRAGMENTS;               // Input molecules
            worker_ref.m_CurrentWorkerMolIndex            = current_mol_index;             // Current molecule index
            worker_ref.m_CurrentWorkerRepeatCycle         = current_mol_index;             // Current molecule index
            worker_ref.m_WorkerRoutine                    = ROUTINE;                       // Set the routine for the worker
            worker_ref.m_Reagents                         = REAGENTS;                      // Reagents for the reaction
            worker_ref.m_ReactionsDirectory               = REACTIONS_DIRECTORY;           // Reactions
            worker_ref.m_Initialized                      = false;
            worker_ref.m_TargetReactantPositions          = TARGET_REACTION_POS;
            worker_ref.m_LigandBased                      = LIGAND_BASED;
            worker_ref.m_WorkerSampleConfs                = SAMPLE_CONFS;                  // Conformer generator
            worker_ref.m_FixGeometry                      = FIX_GEOMETRY;
            worker_ref.m_FixRingGeometry                  = FIX_RING_GEOMETRY;
            worker_ref.m_AdjacentAtoms                    = EXTEND_ADJACENT;

            // manage molecule vs. cycle threads
            REPEATS < INPUT_FRAGMENTS->GetSize() ?
            ++current_mol_index :
            ++current_cycle_index;

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

        // Return the number of threads being used
        size_t GetNumberThreads()
        {
          return m_Threads;
        }

        // Return the number of threads being used
        size_t GetNumberInputFragments()
        {
          return m_NumberInputFragments;
        }

        // Return the number of threads being used
        size_t GetNumberRepeatCycles()
        {
          return m_NumberRepeatCycles;
        }

        // Return molecules
        storage::Vector< chemistry::FragmentComplete> &GetMolecules()
        {
          return m_Molecules;
        }

        // Return yay or nay on Corina conformer
        bool GetIsCorina()
        {
          return m_Corina;
        }

      //////////////////////
      // helper functions //
      //////////////////////

        // Add fit molecules to output vector
        void AddMolecule
        (
          const chemistry::FragmentComplete &MOLECULE
        )
        {
          MOLECULE.WriteMDL( m_OutputStream);
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

      //! flag for input molecules to be base reagents
      util::ShPtr< command::FlagInterface> m_StartingFragmentsFlag;

      //! the molecules to be reacted with the input molecules
      util::ShPtr< command::FlagInterface> m_ReagentsFlag;

      //! the directory containing the reaction files
      util::ShPtr< command::FlagInterface> m_ReactionsFlag;

      //! flag to specify what type of reaction routine to perform
      util::ShPtr< command::FlagInterface> m_RoutineFlag;

      //! flag to specify the number of times to repeat a reaction routine for a starting fragment
      util::ShPtr< command::FlagInterface> m_RepeatsFlag;

      //! flag to specify which reactant positions are allowed by the starting fragments
      util::ShPtr< command::FlagInterface> m_TargetReactantPositionsFlag;

      //! flag to indicate that we do not need to preserve coordinate information at all during reaction
      util::ShPtr< command::FlagInterface> m_LigandBasedFlag;

      //! flag to provide sample conformations object
      util::ShPtr< command::FlagInterface> m_SampleConfsFlag;

      //! flag to fix bad geometry atoms anywhere in product molecules
      util::ShPtr< command::FlagInterface> m_FixGeometryFlag;

      //! flag to fix bad geometry atoms in rings of product molecules
      util::ShPtr< command::FlagInterface> m_FixRingGeometryFlag;

      //! flag to enable conformational sampling of atoms in mobile group adjacent to non-mobile atoms
      util::ShPtr< command::FlagInterface> m_ExtendAdjacentAtomsFlag;

      //! flag for defining output filename,
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! flag to enable corina conformer generation (if installed)
      util::ShPtr< command::FlagInterface> m_CorinaFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeReact();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MoleculeReact *Clone() const
      {
        return new MoleculeReact( *this);
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
        return "Performs chemical reactions to produce new molecules.";
      }

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const
      {
        static std::string s_read_me =
            "React performs chemical reactions to produce new molecules. Currently, intramolecular and partial intramolecular reactions not supported. Multiple threads "
            "can be used to perform reactions off of multiple input starting fragments simultaneously."
            "Random - Perform a random reaction on each input fragment \n"
            "Exhaustive - Perform all available reactions on each input fragment \n";
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

        sp_cmd->AddFlag( m_StartingFragmentsFlag);
        sp_cmd->AddFlag( m_ReagentsFlag);
        sp_cmd->AddFlag( m_ReactionsFlag);
        sp_cmd->AddFlag( m_RoutineFlag);
        sp_cmd->AddFlag( m_RepeatsFlag);
        sp_cmd->AddFlag( m_TargetReactantPositionsFlag);
        sp_cmd->AddFlag( m_LigandBasedFlag);
        sp_cmd->AddFlag( m_SampleConfsFlag);
        sp_cmd->AddFlag( m_FixGeometryFlag);
        sp_cmd->AddFlag( m_FixRingGeometryFlag);
        sp_cmd->AddFlag( m_ExtendAdjacentAtomsFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_CorinaFlag);

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

        // before we go any further make sure they passed the required flags
        BCL_Assert( m_StartingFragmentsFlag->GetFlag(), "Must provide starting fragments!");
        BCL_Assert( m_ReagentsFlag->GetFlag(), "Must provide reagents!");

        // open input stream for molecules
        io::IFStream input;

        // read in staring fragments
        storage::Vector< std::string> filenames( m_StartingFragmentsFlag->GetStringList());
        util::ShPtr< chemistry::FragmentEnsemble> sp_inputs( new chemistry::FragmentEnsemble);
        for
        (
            storage::Vector< std::string>::const_iterator
            itr( filenames.Begin()), itr_end( filenames.End());
            itr != itr_end;
            ++itr
        )
        {
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
          inputs_v.PushBack( *itr);
        }
        util::ShPtr< storage::Vector< chemistry::FragmentComplete >> sp_inputs_v( new storage::Vector< chemistry::FragmentComplete>( inputs_v));

        // read in reagents
        storage::Vector< std::string> reagents_filenames( m_ReagentsFlag->GetStringList());
        util::ShPtr< chemistry::FragmentEnsemble> sp_inputs_reagents( new chemistry::FragmentEnsemble);
        for
        (
            storage::Vector< std::string>::const_iterator
            itr( reagents_filenames.Begin()), itr_end( reagents_filenames.End());
            itr != itr_end;
            ++itr
        )
        {
          io::File::MustOpenIFStream( input, *itr);
          sp_inputs_reagents->ReadMoreFromMdl( input, sdf::e_Maintain);
          io::File::CloseClearFStream( input);
        }

      /////////////////////////
      // parse the arguments //
      /////////////////////////

        // sample conformations
        chemistry::SampleConformations sample_confs;
        sample_confs.TryRead( m_SampleConfsFlag->GetFirstParameter()->GetValue(), util::GetLogger());

      //////////////////////////////
      //  Prepare thread manager  //
      //////////////////////////////

        // Start track time
        util::Stopwatch threadmanager_timer( "Reactions", util::Time( 1, 0), util::Message::e_Standard, true, false);
        threadmanager_timer.Start();

        // Run multi-threading molecule fitting (alignment, docking, etc.)
        ThreadManager thread_manager
        (
          sp_inputs_v,
          sp_inputs_reagents,
          m_ReactionsFlag->GetFirstParameter()->GetValue(),
          m_RoutineFlag->GetFirstParameter()->GetValue(),
          m_RepeatsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_TargetReactantPositionsFlag->GetNumericalList< size_t>(),
          m_LigandBasedFlag->GetFlag(),
          sample_confs,
          m_FixGeometryFlag->GetFlag(),
          m_FixRingGeometryFlag->GetFlag(),
          m_ExtendAdjacentAtomsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_OutputFilenameFlag->GetFirstParameter()->GetValue(),
          m_CorinaFlag->GetFlag(),
          sched::GetNumberCPUs()
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
      static const ApplicationType MoleculeReact_Instance;

    }; // MoleculeReact

    //! @brief standard constructor
    MoleculeReact::MoleculeReact() :

          m_StartingFragmentsFlag
          (
            new command::FlagStatic
            (
              "starting_fragments",
              "fragments to react with reagents",
              command::Parameter
              (
                "starting_fragments",
                "filename(s) for input SDF",
                ""
              )
            )
          ),
          m_ReagentsFlag
          (
            new command::FlagStatic
            (
              "reagents",
              "reagents to react with starting fragments",
              command::Parameter
              (
                "reagents",
                "filename(s) for input reagent SDF",
                ""
              )
            )
          ),
          m_ReactionsFlag
          (
            new command::FlagStatic
            (
              "reactions",
              "directory containing reaction files",
              command::Parameter
              (
                "reactions",
                "directory name for input RXN files",
                chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "functional_reactions"
              )
            )
          ),
          m_RoutineFlag
          (
            new command::FlagDynamic
            (
              "routine", "the reaction routine to perform",
              command::Parameter
              (
                "routine",
                "Random - Perform a random reaction on each starting fragment \n"
                "Exhaustive - Perform all reactions on each starting fragment \n",
                command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "Random", "Exhaustive")),
                "Random"
              )
            )
          ),
          m_RepeatsFlag
          (
            new command::FlagStatic
            (
              "repeats",
              "the number of times to repeat a reaction routine for a given 'starting_fragment'; "
              "this flag is incompatible with the 'Exhaustive' routine; intended to be "
              "used with stochastic routines",
              command::Parameter
              (
                "repeats",
                "",
                command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
                "0"
              )
            )
          ),
          m_TargetReactantPositionsFlag
          (
            new command::FlagDynamic
            (
              "allowed_reactant_positions",
              "reactant position indices (0-indexed) indicating which parts of the "
              "reaction(s) the target molecule is allowed to match; i.e., if you only "
              "want the target fragment to be considered a potential candidate for the "
              "second position of a reaction, pass '1'; default is to test all positions",
              command::Parameter
              (
                "allowed_reactant_positions",
                "defaults to all reactant positions",
                command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
                ""
              )
            )
          ),
          m_LigandBasedFlag
          (
            new command::FlagStatic
            (
              "ligand_based",
              "overrides 3D conformer settings to just produce an arbitrary conformer without preserving spatial information; "
              "by default all of the mutates derived from FragmentMutateInterface attempt to preserve as much coordinate "
              "information as possible while still yielding high quality 3D conformers. In ligand-based drug discovery we "
              "frequently do not care about real space coordinates. Because preserving the real space information is "
              "considerably more expensive than generating a de novo conformer with complex reaction-based mutates and "
              "the resulting pose is less likely to be biologically relevant (because usually with reactions we have "
              "multiple smaller fragments that we start with and frequently lose atoms along the way), it may be desirable to "
              "simply generate a random 3D conformer."
            )
          ),
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
                "max_iterations=1000,"
                "max_conformations=1,"
                "change_chirality=0"
                ")"
              )
            )
          ),
          m_FixGeometryFlag
          (
            new command::FlagStatic
            (
              "fix_geometry",
              "pose-dependent; "
              "if 3D conformer matters, fix atoms with bad geometry even if they are in reference structure"
            )
          ),
          m_FixRingGeometryFlag
          (
            new command::FlagStatic
            (
              "fix_ring_geometry",
              "pose-dependent; "
              "if 3D conformer matters, add all ring atoms from non-reference scaffolds to mobile selection"
            )
          ),
          m_ExtendAdjacentAtomsFlag
          (
            new command::FlagStatic
            (
              "extend_adjacent_atoms",
              "pose-dependent; "
              "include adjacent atoms out this many bonds from any perturbed atom when generating a new 3D conformer",
              command::Parameter
              (
                "adjacent_atoms",
                "",
                command::ParameterCheckRanged< size_t>( 0, 4),
                "0"
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
          m_CorinaFlag
          (
            new command::FlagStatic
            (
              "corina",
              "make a system call to the external program Corina to make the final 3D conformer of each product"
            )
          )
    {
    }

    const ApplicationType MoleculeReact::MoleculeReact_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeReact(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
