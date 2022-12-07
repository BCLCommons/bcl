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
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_grow.h"
#include "chemistry/bcl_chemistry_fragment_mutate_add_med_chem.h"
#include "chemistry/bcl_chemistry_fragment_mutate_alchemy.h"
#include "chemistry/bcl_chemistry_fragment_mutate_cyclize.h"
#include "chemistry/bcl_chemistry_fragment_mutate_extend_with_linker.h"
#include "chemistry/bcl_chemistry_fragment_mutate_fluorinate.h"
#include "chemistry/bcl_chemistry_fragment_mutate_halogenate.h"
#include "chemistry/bcl_chemistry_fragment_mutate_mcm.h"
#include "chemistry/bcl_chemistry_fragment_mutate_remove_atom.h"
#include "chemistry/bcl_chemistry_fragment_mutate_remove_bond.h"
#include "chemistry/bcl_chemistry_fragment_mutate_ring_swap.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_pick_atom_random.h"
#include "chemistry/bcl_chemistry_pick_fragment_random.h"
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
#include "io/bcl_io_file.h"
#include "math/bcl_math_const_function.h"
#include "math/bcl_math_mutate_combine.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "math/bcl_math_mutate_repeat.h"
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
#include "random/bcl_random_uniform_distribution.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf_mdl_handler.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FocusedLibraryDesign
    //! @brief Application for generating libraries for synthesis using QSAR models and a MCM structure generator
    //!
    //! @author brownbp1, mendenjl, loweew, geanesar
    //! @date 05/09/2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FocusedLibraryDesign :
      public InterfaceRelease
    {

    private:

      // ThreadManager needs access to private nested classes
      friend class ThreadManager;
      friend class Worker;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ThreadManager
      //! @brief manages threads for multithreaded structure generation
      //!
      //! @author mendenjl, geanesar, brownbp1
      //! @date Nov 7, 2013
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class ThreadManager :
        public util::ObjectInterface
      {

      private:

      ///////////
      // Data //
      ///////////

          const size_t                                            m_NumberOfMoleculesRequested; // Number of molecules to build
          size_t                                                  m_NumberOfMoleculesBuilt; // Number of molecules already built
          const size_t                                            m_NumberMCIterations; // Number of iterations in the MC approximator
          const size_t                                            m_NumberMCUnimproved; // Number of allowed consecutive unimproved MC iterations
          const size_t                                            m_NumberMCSkipped; // Number of allowed skipped MC iterations
          const float                                             m_MetropolisTemperature; // Tempterature during Metropolis criterion evaluation
          const size_t                                            m_Threads; // Number of threads
          chemistry::FragmentEnsemble                             m_Molecules; // The molecules which have been built and are ready for output
          chemistry::ConstitutionSet                              m_UniqueConsts; // The unique molecules which have been built
          chemistry::ConfigurationSet                             m_UniqueConfigs; // The unique molecules which have been built
          io::OFStream                                            m_OutputStream; // Output file to write molecules to
          const std::string                                       m_DrugLikenessType; // type of druglikeness filter to use for skipping MCM steps
          const float                                             m_VDWScoreCutoff; // internal VDW score cutoff for 3D conformer (used to check for mols with reasonable substitutions)
          const util::Implementation< chemistry::FragmentSplitInterface> m_SplitImplementation; // splitter to use when making fragments for internal MCM optimization
          const std::string                                       m_PoseDependentMDLProperty; // enable pose-dependent scoring with the receptor indicated by this property
          const std::string                                       m_PoseDependentResolveClashes; // resolve clashes between ligand and receptor
          const size_t                                            m_MaxSequentialMutates;
          const float                                             m_RingSwapProb;
          const float                                             m_CyclizeProb;
          const float                                             m_AlchemyProb;
          const float                                             m_RemoveAtomProb;
          const float                                             m_RemoveBondProb;
          const float                                             m_AddMedChemProb;
          const float                                             m_FluorinateProb;
          const float                                             m_HalogenateProb;
          const float                                             m_ExtendWithLinkerProb;
          sched::Mutex                                            m_Mutex; // Lock for updating Workers

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //!
          //! @class Worker
          //! @brief runs the threads for Worker - builds molecules using metropolis monte-carlo routines
          //!
          //! @author brownbp1, mendenjl, geanesar
          //! @date May 09, 2020
          //!
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          struct Worker
          {
            // Rotamer library to use - read in at Main()
            util::ShPtr< chemistry::FragmentComplete>                                                m_StartFragment; // Base fragment to use
//            util::ShPtr< chemistry::FragmentComplete>                                                m_MutableFragment; // mutable fragment in base fragment
            chemistry::FragmentEnsemble                                                              m_MutableFragment; // mutable fragment in base fragment
            storage::Vector< size_t>                                                                 m_MutableAtomIndices; // mutable atoms in base fragment
            descriptor::CheminfoProperty                                                             m_PropertyScorer; // Set objective function with property instead of model
            size_t                                                                                   m_NumberMCIterations; // Number of MC iterations
            size_t                                                                                   m_NumberMCUnimproved; // Number of allowed consecutive unimproved MC iterations
            size_t                                                                                   m_NumberMCSkipped; // Number of allowed consecutive unimproved MC iterations
            float                                                                                    m_MetropolisTemperature; // Temperature during Metropolis criterion evaluation
            opti::Tracker< chemistry::FragmentComplete, double>                                      m_OptiGoal;
            bool                                                                                     m_SaveAllAcceptedImproved;
            std::string                                                                              m_ConformationComparer; // Conformation comparer
            util::ShPtr< math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double> > m_Score; // Objective function
            util::ShPtr< math::MutateInterface< chemistry::FragmentComplete> >                       m_Mutate; // Grow molecules from scaffold
            bool                                                                                     m_Corina; // enables corina conformers during cleaning
            util::SiPtr< ThreadManager>                                                              m_ThreadManager; // Pointer to the thread manager, needed so Worker can be updated

           // Builds and score the molecule
            void RunThread()
            {
              util::ShPtr< storage::Pair< chemistry::FragmentComplete, double> > last_accepted;
              do
              {
                // create the temperature control
                util::ShPtr< mc::TemperatureInterface> sp_temperature( new mc::TemperatureDefault( float( m_MetropolisTemperature)));

                // create the metropolis
                mc::Metropolis< double> metropolis( sp_temperature, true);

                // create the termination criterion
                opti::CriterionCombine< chemistry::FragmentComplete, double> criterion_combine;

                // insert termination criteria that depends on the total number of MC iterations
                opti::CriterionNumberIterations< chemistry::FragmentComplete, double> maximum_number_iterations( m_NumberMCIterations);
                criterion_combine.InsertCriteria( maximum_number_iterations);

                // insert termination criteria that depends on the total number of unimproved MC iterations
                opti::CriterionUnimproved< chemistry::FragmentComplete, double> maximum_number_unimproved_iterations( m_NumberMCUnimproved);
                criterion_combine.InsertCriteria( maximum_number_unimproved_iterations);

                // insert termination criteria that depends on the total number of skipped MC iterations
                opti::CriterionSkippedSteps< chemistry::FragmentComplete, double> maximum_number_skipped_iterations( m_NumberMCSkipped);
                criterion_combine.InsertCriteria( maximum_number_skipped_iterations);

                // Set up MC method
                mc::Approximator< chemistry::FragmentComplete, double> approximator
                (
                  *m_Score,
                  *m_Mutate,
                  metropolis,
                  criterion_combine,
                  *m_StartFragment,
                  m_OptiGoal
                );

                // assume we start with druglike molecule
                static descriptor::CheminfoProperty bonde( "MoleculeTotalDruglikeBondEnergy");
                double druglike_mol_activity( ( *m_Score)( approximator.GetTracker().GetCurrent()->First()));

                // tell me about the scaffold
                BCL_MessageStd("Scaffold properties");
                BCL_MessageStd( "MolWeight: " + util::Format()( descriptor::GetCheminfoProperties().calc_MolWeight->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "# of HBondAcceptors + HBondDonors: " +
                  util::Format()( descriptor::GetCheminfoProperties().calc_HbondAcceptor->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)
                      + descriptor::GetCheminfoProperties().calc_HbondDonor->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "# of NRotBonds: " + util::Format()( descriptor::GetCheminfoProperties().calc_NRotBond->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "LogP: " + util::Format()( descriptor::GetCheminfoProperties().calc_XLogP->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "Bond energy and atom propensity score: " + util::Format()( bonde->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 3)));
                BCL_MessageStd( "# of F: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsF->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "# of Cl: " + util::Format()(descriptor::GetCheminfoProperties().calc_IsCl->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "# of Br: " + util::Format()(descriptor::GetCheminfoProperties().calc_IsBr->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "# of I: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsI->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "# of Halogens: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsHalogen->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "Complexity : " + util::Format()( descriptor::GetCheminfoProperties().calc_MolComplexity->SumOverObject( approximator.GetTracker().GetCurrent()->First())( 0)));
                BCL_MessageStd( "FLD_Score: " + util::Format()( druglike_mol_activity));

                // run the approximator
                BCL_MessageStd( "MCM BEGIN");
                while( approximator.CanContinue() && approximator.ShouldContinue() && m_ThreadManager->GetNumberMoleculesBuilt() + 1 <= m_ThreadManager->GetNumberMoleculesToBuild())
                {
                  // do next step in approximation and get new molecule
                  approximator.Next();
                  const util::ShPtr< storage::Pair< chemistry::FragmentComplete, double> > &current_mol( approximator.GetTracker().GetCurrent());

                  // Check for undruglike properties of the current molecule
                  if( approximator.GetTracker().GetStatusOfLastStep() == opti::e_Accepted || approximator.GetTracker().GetStatusOfLastStep() == opti::e_Improved)
                  {
                    // save the new molecule
                    last_accepted = current_mol;

                    // tell me about the new mol
                    BCL_MessageStd( "Molecule tracker updated at iteration: " + util::Format()( approximator.GetTracker().GetIteration()));
                    BCL_MessageStd( "MolWeight: " + util::Format()( descriptor::GetCheminfoProperties().calc_MolWeight->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd(
                      "# of HBondAcceptors + HBondDonors: " +
                      util::Format()(
                        descriptor::GetCheminfoProperties().calc_HbondAcceptor->SumOverObject( last_accepted->First())( 0)
                        + descriptor::GetCheminfoProperties().calc_HbondDonor->SumOverObject( last_accepted->First())( 0)
                      )
                    );
                    BCL_MessageStd( "# of NRotBonds: " + util::Format()( descriptor::GetCheminfoProperties().calc_NRotBond->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd( "LogP: " + util::Format()( descriptor::GetCheminfoProperties().calc_XLogP->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd( "Bond energy and atom propensity score: " + util::Format()( bonde->SumOverObject( last_accepted->First())( 3)));
                    BCL_MessageStd( "# of F: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsF->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd( "# of Cl: " + util::Format()(descriptor::GetCheminfoProperties().calc_IsCl->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd( "# of Br: " + util::Format()(descriptor::GetCheminfoProperties().calc_IsBr->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd( "# of I: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsI->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd( "# of Halogens: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsHalogen->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd( "Complexity : " + util::Format()( descriptor::GetCheminfoProperties().calc_MolComplexity->SumOverObject( last_accepted->First())( 0)));
                    BCL_MessageStd( "FLD_Score: " + util::Format()( last_accepted->Second()));
                  }
                  else
                  {
                    BCL_MessageStd( "FLD_Score: " + util::Format()( current_mol->Second()));
                    BCL_MessageStd( "MCM Rejected");
                  }
                  // save every accepted/improved step of MCM
                  // hack - add this to approximator at some point
                  if( last_accepted.IsDefined() && m_SaveAllAcceptedImproved)
                  {
                    m_ThreadManager->m_Mutex.Lock();
                    if( m_ThreadManager->GetNumberMoleculesBuilt() + 1 <= m_ThreadManager->GetNumberMoleculesToBuild())
                    {
                      // get best molecule and best score
                      chemistry::FragmentComplete best_mol( last_accepted->First());
                      linal::Vector< double> best_score( 1, last_accepted->Second());
                      best_mol.StoreProperty( "FLD_Score", best_score);

                      // save the final MCM molecule
                      if( m_ThreadManager->CheckUniqueConfiguration( best_mol))
                      {
                        m_ThreadManager->AddMolecule( best_mol);
                        m_ThreadManager->IncreaseMoleculeBuiltCount();
                      }
                    }
                    m_ThreadManager->m_Mutex.Unlock();
                  }
                }
                BCL_MessageStd( "MCM END");

                // save molecules to output, lock it down with a mutex
                if( last_accepted.IsDefined())
                {
                  m_ThreadManager->m_Mutex.Lock();
                  if( m_ThreadManager->GetNumberMoleculesBuilt() + 1 <= m_ThreadManager->GetNumberMoleculesToBuild())
                  {
                    // print approximator endhook
                    // since the PrintEndHook() method is protected and we're sort of
                    // hacking the approximator here, I just copied this.
                    // TODO: add the conditional output accepted/improved molecules to approximator function so that
                    // I can get rid of this hack
                    BCL_MessageStd( "MC Minimization ended");
                    const storage::Vector< size_t> &counts( approximator.GetTracker().GetCounts());
                    const size_t &tot_nr_steps( approximator.GetTracker().GetIteration());
                    const size_t &nr_improved( counts( opti::e_Improved));
                    const size_t &nr_accepted( counts( opti::e_Accepted));
                    const size_t &nr_rejected( counts( opti::e_Rejected));
                    const size_t &nr_skipped( counts( opti::e_Skipped));
                    util::Format format_a, format_b;
                    format_a.W( 5);
                    format_b.W( 5).FFP( 2);
                    BCL_MessageStd( "#MC steps: " + util::Format()( tot_nr_steps));
                    BCL_MessageStd
                    (
                      "#MC steps improved:\t" + format_a( nr_improved) + "\t%" + format_b( 100.0 * nr_improved / tot_nr_steps)
                    );
                    BCL_MessageStd
                    (
                      "#MC steps accepted:\t" + format_a( nr_accepted) + "\t%" + format_b( 100.0 * nr_accepted / tot_nr_steps)
                    );
                    BCL_MessageStd
                    (
                      "#MC steps rejected:\t" + format_a( nr_rejected) + "\t%" + format_b( 100.0 * nr_rejected / tot_nr_steps)
                    );
                    BCL_MessageStd
                    (
                      "#MC steps skipped:\t" + format_a( nr_skipped) + "\t%" + format_b( 100.0 * nr_skipped / tot_nr_steps)
                    );

                    // get best molecule and best score
                    chemistry::FragmentComplete best_mol( last_accepted->First());
                    linal::Vector< double> best_score( 1, last_accepted->Second());
                    best_mol.StoreProperty( "FLD_Score", best_score);

                    // save the final MCM molecule
                    if( m_ThreadManager->CheckUniqueConfiguration( best_mol))
                    {
                      m_ThreadManager->AddMolecule( best_mol);
                      m_ThreadManager->IncreaseMoleculeBuiltCount();
                    }
                  }
                  m_ThreadManager->m_Mutex.Unlock();
                }
              } while( m_ThreadManager->UpdateWorker( *this));
            } // RunThread()

          }; // struct Worker

          // Tests to see if the worker should keep running
          bool UpdateWorker( ThreadManager::Worker &WORKER)
          {
            // Lock structure during the modification
            m_Mutex.Lock();

            if( m_NumberOfMoleculesBuilt >= m_NumberOfMoleculesRequested)
            {
              m_Mutex.Unlock();
              return false;
            }

            m_Mutex.Unlock();
            return true;
          } // UpdateWorker()

      public:

          //! param
          ThreadManager(
            util::ShPtr< chemistry::FragmentComplete>              START_FRAGMENT, // Base fragment to use
//            util::ShPtr< chemistry::FragmentComplete>              MUTABLE_FRAGMENT, // mutable fragment in base fragment
            chemistry::FragmentEnsemble                            MUTABLE_FRAGMENT, // mutable fragment in base fragment
            storage::Vector< size_t>                               MUTABLE_ATOM_INDICES, // mutable atom indices in base fragment
            util::ShPtr< chemistry::FragmentEnsemble>              FRAGMENT_POOL, // Fragments to add to base fragment
            descriptor::CheminfoProperty                           PROPERTY_SCORER, // alternative scorer
            bool                                                   INTERNAL_MCM_OPTI,
            opti::Tracker< chemistry::FragmentComplete, double>    MCM_OPTI_GOAL,
            bool                                                   SAVE_ALL_ACCEPTED_IMPROVED,
            const size_t                                           &NUMBER_OF_MOLECULES, // Number to build
            const size_t                                           &NUMBER_OF_ITERATIONS, // Number of MC iterations
            const size_t                                           &NUMBER_UNIMPROVED_ITERATIONS, // Number of allowed consecutive unimproved MC iterations
            const size_t                                           &NUMBER_SKIPPED_ITERATIONS, // Number of allowed consecutive unimproved MC iterations
            const float                                            &METROPOLIS_TEMPERATURE, // Temperature during Metropolis criterion evaluation
            const size_t                                           &NUMBER_THREADS, // Number of threads (from scheduler)
            const std::string                                      &OUTPUT_FILENAME,
            const std::string                                      &DRUG_LIKENESS_TYPE,
            const float                                            &VDW_SCORE_CUTOFF,
            const util::Implementation< chemistry::FragmentSplitInterface> &SPLIT_IMPLEMENTATION,
            const size_t                                           &MAX_SEQUENTIAL_MUTATES,
            const float                                            &RING_SWAP_PROB,
            const float                                            &CYCLIZE_PROB,
            const float                                            &ALCHEMY_PROB,
            const float                                            &REMOVE_ATOM_PROB,
            const float                                            &REMOVE_BOND_PROB,
            const float                                            &ADD_MEDCHEM_PROB,
            const float                                            &FLUORINATE_PROB,
            const float                                            &HALOGENATE_PROB,
            const float                                            &EXTEND_WITH_LINKER_PROB,
            const std::string                                      &POSE_DEPENDENT_MDL_PROPERTY, // pose-dependent scoring
            const std::string                                      &POSE_DEPENDENT_RESOLVE_CLASHES, // resolve clashes
            const bool                                             &CORINA_CONFS // enables cornina conformers during cleaning
          ) :
            m_NumberOfMoleculesRequested( NUMBER_OF_MOLECULES),
            m_NumberOfMoleculesBuilt( 0),
            m_NumberMCIterations( NUMBER_OF_ITERATIONS),
            m_NumberMCUnimproved( NUMBER_UNIMPROVED_ITERATIONS),
            m_NumberMCSkipped( NUMBER_SKIPPED_ITERATIONS),
            m_MetropolisTemperature( METROPOLIS_TEMPERATURE),
            m_Threads( std::min( NUMBER_THREADS, NUMBER_OF_MOLECULES)),
            m_DrugLikenessType( DRUG_LIKENESS_TYPE),
            m_VDWScoreCutoff( VDW_SCORE_CUTOFF),
            m_SplitImplementation( SPLIT_IMPLEMENTATION),
            m_MaxSequentialMutates( MAX_SEQUENTIAL_MUTATES),
            m_RingSwapProb( RING_SWAP_PROB),
            m_CyclizeProb( CYCLIZE_PROB),
            m_AlchemyProb( ALCHEMY_PROB),
            m_RemoveAtomProb( REMOVE_ATOM_PROB),
            m_RemoveBondProb( REMOVE_BOND_PROB),
            m_AddMedChemProb( ADD_MEDCHEM_PROB),
            m_FluorinateProb( FLUORINATE_PROB),
            m_HalogenateProb( HALOGENATE_PROB),
            m_ExtendWithLinkerProb( EXTEND_WITH_LINKER_PROB),
            m_PoseDependentMDLProperty( POSE_DEPENDENT_MDL_PROPERTY),
            m_PoseDependentResolveClashes( POSE_DEPENDENT_RESOLVE_CLASHES)
          {
            // prepare output filestream
            io::File::MustOpenOFStream( m_OutputStream, OUTPUT_FILENAME);

            // tree search for RingSwap
            util::ShPtr< chemistry::SearchFragmentLibraryFromTree> tree_search
            (
              new chemistry::SearchFragmentLibraryFromTree
              (
                *util::Implementation< chemistry::RotamerLibraryInterface>( chemistry::RotamerLibraryInterface::GetDefault())
              )
            );

            // set up our primary mutater object
            util::ShPtr< math::MutateDecisionNode< chemistry::FragmentComplete> > mutater
            (
              new math::MutateDecisionNode< chemistry::FragmentComplete>()
            );

            // get the starting molecule minus the mutable region for local mutations
             util::ShPtr< chemistry::FragmentComplete> scaffold_fragment( new chemistry::FragmentComplete());
             if( MUTABLE_FRAGMENT.GetMolecules().FirstElement().GetSize() || MUTABLE_ATOM_INDICES.GetSize())
             {
               static chemistry::FragmentTrackMutableAtoms atom_tracker;
               scaffold_fragment =
                   util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete
                     (
                       atom_tracker.GetBaseFragment
                       (
                         *START_FRAGMENT,
                         MUTABLE_FRAGMENT.GetMolecules().FirstElement(),
                         MUTABLE_ATOM_INDICES
                       )
                     ));
               BCL_Assert( scaffold_fragment->GetSize(), "Exiting because of incompatible mutable options");
             }

            // if the internal MCM local optimization option is selected
            if( INTERNAL_MCM_OPTI)
            {
              // POSE-DEPENDENT CONSTRUCTION OF MUTATES //
              if( !POSE_DEPENDENT_MDL_PROPERTY.empty())
              {
                BCL_MessageStd( "Pose-dependent scoring enabled");
                // set clash resolver
                bool clash_resolver;
                POSE_DEPENDENT_RESOLVE_CLASHES == "true" ?
                    clash_resolver = true:
                    clash_resolver = false;
                mutater->AddMutate
                (
                  chemistry::FragmentMutateMCM
                  (
                    MCM_OPTI_GOAL,
                    SPLIT_IMPLEMENTATION,
                    tree_search,
                    FRAGMENT_POOL,
                    m_DrugLikenessType,
                    *START_FRAGMENT,
                    chemistry::FragmentEnsemble( storage::List< chemistry::FragmentComplete>( 1, *START_FRAGMENT)),
                    MUTABLE_ATOM_INDICES,
                    POSE_DEPENDENT_MDL_PROPERTY,
                    PROPERTY_SCORER,
                    clash_resolver,
                    storage::Vector< float>(),
                    CORINA_CONFS,
                    m_MaxSequentialMutates,
                    m_RingSwapProb,
                    m_CyclizeProb,
                    m_AlchemyProb,
                    m_RemoveAtomProb,
                    m_RemoveBondProb,
                    m_AddMedChemProb,
                    m_FluorinateProb,
                    m_HalogenateProb,
                    m_ExtendWithLinkerProb
                  ),
                  1.0
                );
              }
              // POSE-INDEPENDENT CONSTRUCTION OF MUTATES //
              else
              {
                mutater->AddMutate
                (
                  chemistry::FragmentMutateMCM
                  (
                    MCM_OPTI_GOAL,
                    SPLIT_IMPLEMENTATION,
                    tree_search,
                    FRAGMENT_POOL,
                    m_DrugLikenessType,
                    *START_FRAGMENT,
                    chemistry::FragmentEnsemble( storage::List< chemistry::FragmentComplete>( 1, *START_FRAGMENT)),
                    MUTABLE_ATOM_INDICES,
                    PROPERTY_SCORER,
                    CORINA_CONFS,
                    m_MaxSequentialMutates,
                    m_RingSwapProb,
                    m_CyclizeProb,
                    m_AlchemyProb,
                    m_RemoveAtomProb,
                    m_RemoveBondProb,
                    m_AddMedChemProb,
                    m_FluorinateProb,
                    m_HalogenateProb,
                    m_ExtendWithLinkerProb
                  ),
                  1.0
                );
              }
            }
            // otherwise, just add the mutates and let them fly
            else
            {
              // POSE-DEPENDENT CONSTRUCTION OF MUTATES //
              chemistry::FragmentEnsemble scaffold_ens( storage::List< chemistry::FragmentComplete>( 1, *scaffold_fragment));
              if( !POSE_DEPENDENT_MDL_PROPERTY.empty())
              {
                BCL_MessageStd( "Pose-dependent scoring enabled");
                // set clash resolver
                bool clash_resolver;
                POSE_DEPENDENT_RESOLVE_CLASHES == "true" ?
                    clash_resolver = true:
                    clash_resolver = false;
                mutater->AddMutate( chemistry::FragmentMutateRingSwap( tree_search, m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS, true, false, 0.1, true, true), m_RingSwapProb);
                mutater->AddMutate( chemistry::FragmentMutateCyclize( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS), m_CyclizeProb);
                mutater->AddMutate( chemistry::FragmentMutateAlchemy( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS), m_AlchemyProb);
                mutater->AddMutate( chemistry::FragmentMutateRemoveAtom( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS), m_RemoveAtomProb);
                mutater->AddMutate( chemistry::FragmentMutateRemoveBond( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS), m_RemoveBondProb);
                mutater->AddMutate( chemistry::FragmentMutateExtendWithLinker( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS), m_ExtendWithLinkerProb);
                mutater->AddMutate( chemistry::FragmentMutateAddMedChem( FRAGMENT_POOL, m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS), m_AddMedChemProb);
                mutater->AddMutate( chemistry::FragmentMutateFluorinate( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS), m_FluorinateProb);
                mutater->AddMutate( chemistry::FragmentMutateHalogenate( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, POSE_DEPENDENT_MDL_PROPERTY, PROPERTY_SCORER, clash_resolver, storage::Vector< float>(), CORINA_CONFS), m_HalogenateProb);
              }
              // POSE-INDEPENDENT CONSTRUCTION OF MUTATES //
              else
              {
                BCL_MessageStd( "Pose-independent scoring enabled");
                mutater->AddMutate( chemistry::FragmentMutateRingSwap( tree_search, m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS, true, false, 0.1, true, true), m_RingSwapProb);
                mutater->AddMutate( chemistry::FragmentMutateCyclize( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS), m_CyclizeProb);
                mutater->AddMutate( chemistry::FragmentMutateAlchemy( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS), m_AlchemyProb);
                mutater->AddMutate( chemistry::FragmentMutateRemoveAtom( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS), m_RemoveAtomProb);
                mutater->AddMutate( chemistry::FragmentMutateRemoveBond( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS), m_RemoveBondProb);
                mutater->AddMutate( chemistry::FragmentMutateExtendWithLinker( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS), m_ExtendWithLinkerProb);
                mutater->AddMutate( chemistry::FragmentMutateAddMedChem( FRAGMENT_POOL, m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS), m_AddMedChemProb);
                mutater->AddMutate( chemistry::FragmentMutateFluorinate( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS), m_FluorinateProb);
                mutater->AddMutate( chemistry::FragmentMutateHalogenate( m_DrugLikenessType, *START_FRAGMENT, scaffold_ens, MUTABLE_ATOM_INDICES, CORINA_CONFS), m_HalogenateProb);
              }
            }

            // set up sequential mutate to perform 1 to N mutates in a row prior to scoring (does not bypass druglikeness filtering)
            util::ShPtr< math::MutateInterface< chemistry::FragmentComplete> > mutate_repeater
            (
              new math::MutateRepeat< chemistry::FragmentComplete>
              (
                mutater,
                1,
                m_MaxSequentialMutates
              )
            );

            // Set up workers
            std::vector< Worker> workers( m_Threads);
            for(
                std::vector< Worker>::iterator itr( workers.begin()), end( workers.end());
                itr != end;
                ++itr
            )
            {
              Worker &worker_ref( *itr);
              worker_ref.m_StartFragment                = START_FRAGMENT.HardCopy();
              worker_ref.m_StartFragment->GetCacheMap() = util::ShPtr< descriptor::CacheMap>( new descriptor::CacheMap);
              worker_ref.m_MutableFragment              = MUTABLE_FRAGMENT;
              worker_ref.m_MutableAtomIndices           = MUTABLE_ATOM_INDICES;
              worker_ref.m_NumberMCIterations           = m_NumberMCIterations;
              worker_ref.m_NumberMCUnimproved           = m_NumberMCUnimproved;
              worker_ref.m_NumberMCSkipped              = m_NumberMCSkipped;
              worker_ref.m_MetropolisTemperature        = m_MetropolisTemperature;
              worker_ref.m_OptiGoal                     = MCM_OPTI_GOAL;
              worker_ref.m_SaveAllAcceptedImproved      = SAVE_ALL_ACCEPTED_IMPROVED;
              worker_ref.m_Corina                       = CORINA_CONFS;
              worker_ref.m_ThreadManager                = this;
              worker_ref.m_PropertyScorer               = PROPERTY_SCORER.HardCopy();
              worker_ref.m_Score                        = util::ShPtr< math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double> >
              (
                new chemistry::ScoreFunctionGeneric( worker_ref.m_PropertyScorer)
              );
              if( INTERNAL_MCM_OPTI)
              {
                worker_ref.m_Mutate = mutater;
              }
              else
              {
                worker_ref.m_Mutate = mutate_repeater;
              }
            }

            // Allocate space for jobs
            util::ShPtrVector< sched::JobInterface> jobs;
            jobs.AllocateMemory( m_Threads);

            const size_t group_id( 1);
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

          // Increase the number of molecules that have been built
          void IncreaseMoleculeBuiltCount()
          {
            BCL_MessageVrb( "Number of molecules built: " + util::Format()( m_NumberOfMoleculesBuilt + 1));
            m_NumberOfMoleculesBuilt++;
          }

          size_t GetNumberMoleculesBuilt()
          {
            return m_NumberOfMoleculesBuilt;
          }

          size_t GetNumberMoleculesToBuild()
          {
            return m_NumberOfMoleculesRequested;
          }

          size_t GetNumberMCIterations()
          {
            return m_NumberMCIterations;
          }

          size_t GetNumberMCUnimproved()
          {
            return m_NumberMCUnimproved;
          }

          size_t GetNumberMCSkipped()
          {
            return m_NumberMCSkipped;
          }

          bool CheckUniqueConstitution( const chemistry::FragmentComplete &MOLECULE)
          {
            bool unique( m_UniqueConsts.Insert( chemistry::FragmentConstitutionShared( MOLECULE)).second);
            BCL_MessageStd( "Number in ConstitutionSet: " + util::Format()( m_UniqueConsts.GetSize()));
            BCL_MessageVrb( "Unique? : " + util::Format()( unique));
            return unique;
          }

          bool CheckUniqueConfiguration( const chemistry::FragmentComplete &MOLECULE)
          {
            bool unique( m_UniqueConfigs.Insert( chemistry::FragmentConfigurationShared( MOLECULE)).second);
            BCL_MessageStd( "Number in ConfigurationSet: " + util::Format()( m_UniqueConfigs.GetSize()));
            BCL_MessageVrb( "Unique? : " + util::Format()( unique));
            return unique;
          }

          void AddMolecule( const chemistry::FragmentComplete &MOLECULE)
          {
            m_Molecules.PushBack( MOLECULE);
            MOLECULE.WriteMDL( m_OutputStream);
          }

          // Return FragmentEnsemble of the generated molecules
          chemistry::FragmentEnsemble &GetMolecules()
          {
            return m_Molecules;
          }

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

      //! flag to control number of molecules to be generated
      util::ShPtr< command::FlagInterface> m_NumberMoleculesFlag;

      //! flag to control the number of MC iterations in molecule optimization
      util::ShPtr< command::FlagInterface> m_NumberIterationsFlag;

      //! flag to control the number of maximum allowed consecutive unimproved MC iterations
      util::ShPtr< command::FlagInterface> m_NumberUnimprovedFlag;

      //! flag to control the number of maximum allowed skipped MC iterations
      util::ShPtr< command::FlagInterface> m_NumberSkippedFlag;

      //! flag to control the temperature for the Metropolis criterion
      util::ShPtr< command::FlagInterface> m_MetropolisTemperatureFlag;

      //! flag to control input base fragment
      util::ShPtr< command::FlagInterface> m_StartFragmentFlag;

      //! flag to control input mutable fragment within base fragment
      util::ShPtr< command::FlagInterface> m_MutableFragmentFlag;

      //! flag to control input mutable atoms within base fragment
      util::ShPtr< command::FlagInterface> m_MutableAtomsFlag;

      //! flag for defining output filename,
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! flag for defining input fragments
      util::ShPtr< command::FlagInterface> m_GrowFragmentsFlag;

      //! flag for an alternative score function to just the trained model
      util::ShPtr< command::FlagInterface> m_PropertyScoringFunctionFlag;

      //! flag for the druglikeness filter to use
      util::ShPtr< command::FlagInterface> m_DrugLikenessTypeFlag;

      //! flag to split molecules
      util::ShPtr< command::FlagInterface> m_SplitImplementationFlag;

      //! flag to do an internal MCM optimization
      util::ShPtr< command::FlagInterface> m_SimulatedAnnealingFlag;

      //! flag to maximize score istead of minimize
      util::ShPtr< command::FlagInterface> m_LargerIsBetterFlag;

      //! flag to save all molecules accepted or improved by main MCM
      util::ShPtr< command::FlagInterface> m_SaveAllAcceptedImprovedFlag;

      //! flag to use corina to generate starting conformer
      util::ShPtr< command::FlagInterface> m_Corina;

      //! flag to set 3D VDW score cutoff
      util::ShPtr< command::FlagInterface> m_VDWClashCutoffFlag;

      //! flag to enable pose-dependent scoring (default is ligand-based scoring)
      util::ShPtr< command::FlagInterface> m_PoseDependentFlag;

      //! flag controlling the maximum possible number of sequential mutates that can occur between MCM evaluation
      util::ShPtr< command::FlagInterface> m_MaxSequentialMutatesFlag;

      //! flags controling relative probabilities of different mutate objects
      util::ShPtr< command::FlagInterface> m_RingSwapProbFlag;
      util::ShPtr< command::FlagInterface> m_CyclizeProbFlag;
      util::ShPtr< command::FlagInterface> m_AlchemyProbFlag;
      util::ShPtr< command::FlagInterface> m_RemoveAtomProbFlag;
      util::ShPtr< command::FlagInterface> m_RemoveBondProbFlag;
      util::ShPtr< command::FlagInterface> m_AddMedChemProbFlag;
      util::ShPtr< command::FlagInterface> m_FluorinateProbFlag;
      util::ShPtr< command::FlagInterface> m_HalogenateProbFlag;
      util::ShPtr< command::FlagInterface> m_ExtendWithLinkerProbFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      FocusedLibraryDesign();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      FocusedLibraryDesign *Clone() const
      {
        return new FocusedLibraryDesign( *this);
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

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const
      {
        static std::string s_read_me =
          "FocusedLibraryDesign generates a distribution of new molecules by applying alchemical and medicinal chemistry-like transformations"
          " to a starting scaffold or molecule. A quantitative structure-activity relationship (QSAR) model scores generated molecules "
          " and a Monte Carlo - Metropolis (MCM) algorithm samples the distribution of generated molecules based on QSAR score. Optionally,"
          " an internal MCM-simulated annealing (SA) optimization trial can be performed for each individual transformation, localized to a randomly"
          " selected or manually specified substructure(s). After every transformation, the new molecule is evaluated for drug-likeness according to"
          " a user-specified composite metric. If the molecule is deemed non-druglike, then the molecule is excluded from further analysis and the MCM"
          " move is repeated. When the MCM termination criteria are met, the new molecule is output into an SDF file. Optionally, all generated molecules"
          " that are accepted and/or improved by the MCM engine can be output to the SDF file instead of just the final molecule (note that this will"
          " reduce the QSAR score distribution of the final library). The program terminates when the requested number of unique molecules have been"
          " generated.";
        return s_read_me;
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Generates distributions of molecules utilizing alchemical mutations and a "
            "property-based score metric, such as a QSAR model.";
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

        // insert all the flags and params
        sp_cmd->AddFlag( m_NumberMoleculesFlag);
        sp_cmd->AddFlag( m_NumberIterationsFlag);
        sp_cmd->AddFlag( m_NumberUnimprovedFlag);
        sp_cmd->AddFlag( m_NumberSkippedFlag);
        sp_cmd->AddFlag( m_MetropolisTemperatureFlag);
        sp_cmd->AddFlag( m_StartFragmentFlag);
        sp_cmd->AddFlag( m_MutableFragmentFlag);
        sp_cmd->AddFlag( m_MutableAtomsFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_GrowFragmentsFlag);
        sp_cmd->AddFlag( m_PropertyScoringFunctionFlag);
        sp_cmd->AddFlag( m_DrugLikenessTypeFlag);
        sp_cmd->AddFlag( m_SplitImplementationFlag);
        sp_cmd->AddFlag( m_SimulatedAnnealingFlag);
        sp_cmd->AddFlag( m_LargerIsBetterFlag);
        sp_cmd->AddFlag( m_SaveAllAcceptedImprovedFlag);
        sp_cmd->AddFlag( m_Corina);
        sp_cmd->AddFlag( m_VDWClashCutoffFlag);
        sp_cmd->AddFlag( m_PoseDependentFlag);
        sp_cmd->AddFlag( m_MaxSequentialMutatesFlag);
        sp_cmd->AddFlag( m_RingSwapProbFlag);
        sp_cmd->AddFlag( m_CyclizeProbFlag);
        sp_cmd->AddFlag( m_AlchemyProbFlag);
        sp_cmd->AddFlag( m_RemoveAtomProbFlag);
        sp_cmd->AddFlag( m_RemoveBondProbFlag);
        sp_cmd->AddFlag( m_AddMedChemProbFlag);
        sp_cmd->AddFlag( m_FluorinateProbFlag);
        sp_cmd->AddFlag( m_HalogenateProbFlag);
        sp_cmd->AddFlag( m_ExtendWithLinkerProbFlag);

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

        // setup the base fragment
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_StartFragmentFlag->GetFirstParameter()->GetValue());

        // Needs to be wrapped in a ShPtr so it can be passed to ThreadManager
        util::ShPtr< chemistry::FragmentComplete> sp_startfragment
        (
          new chemistry::FragmentComplete( sdf::FragmentFactory::MakeFragment( input, sdf::e_Maintain))
        );
        io::File::CloseClearFStream( input);

        // setup the mutable fragment
        util::ShPtr< chemistry::FragmentComplete> sp_mutablefragment( new chemistry::FragmentComplete());
        chemistry::FragmentEnsemble mutable_fragments;
        if( m_MutableFragmentFlag->GetFlag())
        {
          io::File::MustOpenIFStream( input, m_MutableFragmentFlag->GetFirstParameter()->GetValue());

          // Needs to be wrapped in a ShPtr so it can be passed to ThreadManager
          chemistry::FragmentComplete frag( sdf::FragmentFactory::MakeFragment( input, sdf::e_Maintain));
          sp_mutablefragment = util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete( frag));
          mutable_fragments = chemistry::FragmentEnsemble( storage::List< chemistry::FragmentComplete>( 1, *sp_mutablefragment));
          // message indicating using mutable fragment
          BCL_MessageStd(
            "Mutating substructure atoms specified in the file '" +
            util::Format()( m_MutableFragmentFlag->GetFirstParameter()->GetValue()) + "'"
          );
        }
        io::File::CloseClearFStream( input);

        // setup the mutable atom indices
        storage::Vector< size_t> mutable_atom_indices;
        if( m_MutableAtomsFlag->GetFlag())
        {
          // convert the functionalization points to numeric values
          storage::Vector< std::string> fxnl_pts_str( m_MutableAtomsFlag->GetStringList());

          // output mutable indices to terminal
          std::string mutable_indices_message;
          for
          (
              auto itr( fxnl_pts_str.Begin()), itr_end( fxnl_pts_str.End());
              itr != itr_end;
              ++itr
          )
          {
            mutable_indices_message.append(util::Format()( *itr));
            mutable_indices_message.append(",");
          }
          BCL_MessageStd( "Mutable atom indices: " + util::Format()( mutable_indices_message));

          storage::Set< size_t> fxnl_pts_set;
          for( size_t i( 0), l( fxnl_pts_str.GetSize()); i < l; ++i)
          {
            size_t point;
            if( !util::TryConvertFromString( point, fxnl_pts_str( i), util::GetLogger()))
            {
              BCL_MessageStd( "Could not parse \"" + fxnl_pts_str( i) + "\" as a number");
              continue;
            }
            if( point < sp_startfragment->GetSize())
            {
              fxnl_pts_set.Insert( point);
            }
            else
            {
              BCL_MessageStd
              (
                "Warning: specified point \"" + util::Format()( point) + "\""
                " has an index greater than the number of atoms in the molecule, not using this point"
              );
            }
          }
          mutable_atom_indices = storage::Vector< size_t>( fxnl_pts_set.Begin(), fxnl_pts_set.End());
        }

        // try to read cheminfo property scorer
        descriptor::CheminfoProperty property_scorer;
        if( m_PropertyScoringFunctionFlag->GetFlag())
        {
          property_scorer = m_PropertyScoringFunctionFlag->GetFirstParameter()->GetValue();
        }
        else
        {
          // flat energy surface
          property_scorer = "Constant(0.0)";
        }

        // read internal mcm opti flag
        bool internal_mcm_opti( false);
        if( m_SimulatedAnnealingFlag->GetFlag())
        {
          internal_mcm_opti = true;
        }

        // set MCM optimization goal
        opti::Tracker< chemistry::FragmentComplete, double> mcm_opti_goal( opti::e_SmallerIsBetter);
        if( m_LargerIsBetterFlag->GetFlag())
        {
          mcm_opti_goal = opti::e_LargerIsBetter;
        }

        // set save options for intermediate molecules
        bool save_all_accepted_improved( false);
        if( m_SaveAllAcceptedImprovedFlag->GetFlag())
        {
          save_all_accepted_improved = true;
        }

        // set corina conformers
        bool corina_confs( false);
        if( m_Corina->GetFlag())
        {
          corina_confs = true;
        }

        // get splitter
        util::Implementation< chemistry::FragmentSplitInterface> splitter;
        if( m_SplitImplementationFlag->GetFlag())
        {
          splitter = m_SplitImplementationFlag->GetFirstParameter()->GetValue();
        }

      /////////////////////////
      // parse the arguments //
      /////////////////////////

        // get all filename for grow fragments
        const storage::Vector< std::string> filenames( m_GrowFragmentsFlag->GetStringList());

        // creating ShPtr of growfragments
        util::ShPtr< chemistry::FragmentEnsemble> sp_fragment_pool( new chemistry::FragmentEnsemble);

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
          sp_fragment_pool->ReadMoreFromMdl( input, sdf::e_Maintain);
          io::File::CloseClearFStream( input);
        }

      /////////////////////////////
      // Prepare rotamer library //
      /////////////////////////////

        // Start track time
        util::Stopwatch threadmanager_timer( "Molecule Building", util::Time( 1, 0), util::Message::e_Standard, true, false);
        threadmanager_timer.Start();

        // Build the molecules using metropolis monte-carlo
        if( m_PoseDependentFlag->GetFlag())
        {
          ThreadManager thread_manager
          (
            sp_startfragment,
//            sp_mutablefragment,
            mutable_fragments,
            mutable_atom_indices,
            sp_fragment_pool,
            property_scorer,
            internal_mcm_opti,
            mcm_opti_goal,
            save_all_accepted_improved,
            m_NumberMoleculesFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_NumberIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_NumberUnimprovedFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_NumberSkippedFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_MetropolisTemperatureFlag->GetFirstParameter()->GetNumericalValue< float>(),
            sched::GetNumberCPUs(),
            m_OutputFilenameFlag->GetFirstParameter()->GetValue(),
            m_DrugLikenessTypeFlag->GetFirstParameter()->GetValue(),
            m_VDWClashCutoffFlag->GetFirstParameter()->GetNumericalValue< float>(),
            splitter,
            m_MaxSequentialMutatesFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_RingSwapProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_CyclizeProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_AlchemyProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_RemoveAtomProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_RemoveBondProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_AddMedChemProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_FluorinateProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_HalogenateProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_ExtendWithLinkerProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_PoseDependentFlag->GetParameterList()( 0)->GetValue(),
            m_PoseDependentFlag->GetParameterList()( 1)->GetValue(),
            corina_confs
          );
        }
        else
        {
          ThreadManager thread_manager
          (
            sp_startfragment,
//            sp_mutablefragment,
            mutable_fragments,
            mutable_atom_indices,
            sp_fragment_pool,
            property_scorer,
            internal_mcm_opti,
            mcm_opti_goal,
            save_all_accepted_improved,
            m_NumberMoleculesFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_NumberIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_NumberUnimprovedFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_NumberSkippedFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_MetropolisTemperatureFlag->GetFirstParameter()->GetNumericalValue< float>(),
            sched::GetNumberCPUs(),
            m_OutputFilenameFlag->GetFirstParameter()->GetValue(),
            m_DrugLikenessTypeFlag->GetFirstParameter()->GetValue(),
            m_VDWClashCutoffFlag->GetFirstParameter()->GetNumericalValue< float>(),
            splitter,
            m_MaxSequentialMutatesFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
            m_RingSwapProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_CyclizeProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_AlchemyProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_RemoveAtomProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_RemoveBondProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_AddMedChemProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_FluorinateProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_HalogenateProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            m_ExtendWithLinkerProbFlag->GetFirstParameter()->GetNumericalValue< float>(),
            std::string(),
            std::string(),
            corina_confs
          );
        }

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

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType FocusedLibraryDesign_Instance;

    }; // FocusedLibraryDesign

    //! @brief standard constructor
    FocusedLibraryDesign::FocusedLibraryDesign() :
      m_NumberMoleculesFlag
      (
        new command::FlagStatic
        (
          "number_molecules", "flag for number of molecules to generate",
          command::Parameter
          (
            "number_molecules", "total number of molecules",
            command::ParameterCheckRanged< int>( 1, std::numeric_limits< int>::max()), "10"
          )
        )
      ),
      m_NumberIterationsFlag
      (
        new command::FlagStatic
        (
          "number_iterations", "flag for number of MC iterations",
          command::Parameter
          (
            "number_iterations", "maximum number of MC iterations",
            command::ParameterCheckRanged< int>( 1, std::numeric_limits< int>::max()), "100"
          )
        )
      ),
      m_NumberUnimprovedFlag
      (
        new command::FlagStatic
        (
          "number_unimproved", "flag for number of maximum allowed consecutive unimproved MC iterations",
          command::Parameter
          (
            "number_unimproved", "maximum number of allowed consecutive unimproved MC iterations",
            command::ParameterCheckRanged< int>( 1, std::numeric_limits< int>::max()), "100"
          )
        )
      ),
      m_NumberSkippedFlag
      (
        new command::FlagStatic
        (
          "number_skipped", "flag for number of maximum allowed skipped MC iterations",
          command::Parameter
          (
            "number_skipped", "maximum number of allowed skipped MC iterations",
            command::ParameterCheckRanged< int>( 1, std::numeric_limits< int>::max()), "100"
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
            command::ParameterCheckRanged< float>( 1.0, std::numeric_limits< float>::max()), "1.0"
          )
        )
      ),
      m_PoseDependentFlag
      (
        new command::FlagDynamic
        (
          "pose_dependent_scoring", "enables pose-dependent scoring",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "MDL_property",
              "MDL property specifying the PDB filename for the receptor; "
              "needs to match MDL property name used to train the corresponding machine learning model"
            ),
            command::Parameter
            (
              "resolve_clashes", "resolve clashes between the protein and ligand",
              command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "true", "false"))
            )
          )
        )
      ),
      m_MaxSequentialMutatesFlag
      (
        new command::FlagStatic
        (
          "max_sequential_mutates", "flag for the maximum number of mutates that can occur between MCM evaluations; "
              "a number between 1 and <max_sequential_mutates> will be randomly selected with uniform probability",
          command::Parameter
          (
            "max_sequential_mutates", "perform multiple mutates in a row prior to druglikeness filtering and MCM evaluation",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()), "1"
          )
        )
      ),
      m_RingSwapProbFlag
      (
        new command::FlagStatic
        (
          "mutate_ringswap_prob", "flag for the relative probability of performing a ringswap mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_ringswap", "replace a single atom or whole ring with a new ring structure from an internal ring library; "
                "substituents on the altered ring are randomly assigned a position on the new ring",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_CyclizeProbFlag
      (
        new command::FlagStatic
        (
          "mutate_cyclize_prob", "flag for the relative probability of performing a cyclize mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_cyclize", "connect non-ring moieties with other non-ring or ring moieties",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_AlchemyProbFlag
      (
        new command::FlagStatic
        (
          "mutate_alchemy_prob", "flag for the relative probability of performing an alchemy mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_alchemy", "transform one element type into another preserving charge and optimizing the bond type to the new element; "
                "restricted to H, C, N, O, S, F, Cl, Br; "
                "the probability of a transformation is based on the prevalence of each element type in a sample of druglike molecules",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_RemoveAtomProbFlag
      (
        new command::FlagStatic
        (
          "mutate_remove_atom_prob", "flag for the relative probability of performing an atom removal mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_remove_atom", "remove an atom from the molecule; if the "
                "molecule is split into multiple fragments, keep only the largest fragment",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_RemoveBondProbFlag
      (
        new command::FlagStatic
        (
          "mutate_remove_bond_prob", "flag for the relative probability of performing a bond removal mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_remove_bond", "remove a bond from the molecule; preserving only the "
                "largest fragment after the molecule is split",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_AddMedChemProbFlag
      (
        new command::FlagStatic
        (
          "mutate_add_medchem_prob", "flag for the relative probability of performing an add medchem group mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_add_medchem", "append medchem-like functional groups to the molecule from an internal library",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_FluorinateProbFlag
      (
        new command::FlagStatic
        (
          "mutate_fluorinate_prob", "flag for the relative probability of performing a fluorinate mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_fluorinate", "strips a carbon atom of any bonded hydrogen atoms and replaces them with fluorine atoms; "
                "equal probability to strip all fluorine atoms from a carbon and replace with hydrogen atoms",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_HalogenateProbFlag
      (
        new command::FlagStatic
        (
          "mutate_halogenate_prob", "flag for the relative probability of performing a halogenate mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_halogenate", "add an F, Cl, Br, or I atom to an aromatic ring system; "
                "the probability of a transformation is based on the relative prevalence of each "
                "halogen type in a sample of druglike molecules; "
                "F - 25%, Cl - 60%, Br - 10%, I - 5%",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_ExtendWithLinkerProbFlag
      (
        new command::FlagStatic
        (
          "mutate_extendwithlinker_prob", "flag for the relative probability of performing an extendwithlinker mutation during molecule optimization; "
              "automatically rescaled between 0 and 1 with any other mutates",
          command::Parameter
          (
            "mutate_extendwithlinker", "split current molecule into two fragments and re-connect them with a linker consisting of "
                "either a ring, alkyl/methoxy/ethoxy chain, or single element; alternatively, beginning with a ring, create a linker"
                " to a new ring system",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()), "0.1"
          )
        )
      ),
      m_StartFragmentFlag
      (
        new command::FlagStatic
        (
          "start_fragment", "filename for input starting fragment",
          command::Parameter
          (
            "fragment_filename", "filename for input sdf of molecules", ""
          )
        )
      ),
      m_MutableFragmentFlag
      (
        new command::FlagStatic
        (
          "mutable_fragment", "filename for fragment in base fragment that can be mutated",
          command::Parameter
          (
            "mutable_fragment_filename", "if no filename is supplied, defaults to all fragments in base fragment are mutable;"
            " note that this option cannot be used simultaneously with 'mutable_atoms'",
            ""
          )
        )
      ),
      m_MutableAtomsFlag
      (
        new command::FlagDynamic
        (
          "mutable_atoms", "filename for 0-indexed atom indices that can be mutated",
          command::Parameter
          (
            "mutable_atom_index", "if no filename is supplied, defaults to all atoms are mutable;"
            " note that this option cannot be used simultaneously with 'mutable_fragment'",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
            ""
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_filename", "flag selecting the output file name",
          command::Parameter
          (
            "output_filename_param", "filename for output sdf of molecules"
          )
        )
      ),
      m_GrowFragmentsFlag
      (
        new command::FlagStatic
        (
          "grow_fragments",
          "files containing fragments to append to the molecule",
          command::Parameter
          (
            "grow fragments filename",
            "name of file containing grow fragments",
            command::ParameterCheckFileExistence(),
            chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + ( "bcl_buildfrag_0.sdf.gz")
          )
        )
      ),
      m_PropertyScoringFunctionFlag
      (
        new command::FlagDynamic
        (
          "scoring_function",
          "the scoring function to use",
          command::Parameter
          (
            "function",
            "the scoring function implementation to use",
            command::ParameterCheckSerializable
            (
              chemistry::ScoreFunctionGeneric()
            )
          )
        )
      ),
      m_DrugLikenessTypeFlag
      (
        new command::FlagStatic
        (
          "druglikeness_type",
          "the type of druglikeness filter to use to determine when a molecule is skipped by the Monte Carlo algorithm",
          command::Parameter
          (
            "type",
            "the type of druglikeness to use",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create
              (
                "IsConstitutionDruglike",
                "IsConstitutionDruglikeAndHitlike",
                "None"
              )
            ),
            "IsConstitutionDruglike"
          )
        )
      ),
      m_SplitImplementationFlag
      (
        new command::FlagStatic
        (
          "split_implementation",
          "method to split molecules if performing internal MCM optimization",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< chemistry::FragmentSplitInterface>()),
            "Rings"
          )
        )
      ),
      m_SimulatedAnnealingFlag
      (
        new command::FlagDynamic
        (
          "internal_mcm_simulated_annealing",
          "Perform MCM simulated annealing (SA) optimization trials for randomly chosen or manually specified"
          "substructures of molecule. The output from each MCM-SA trial is passed to the main MCM "
          "engine for scoring and decision-making. Repeats until main MCM reaches termination criteria. "
          "This is an alternative to performing general MCM optimization of the whole molecule or a single specific fragment "
          "or set of atoms using the standard mutates with one MCM engine at a fixed temperature. In principle, it allows "
          "the main MCM to evaluate optimized substructure changes; however, it is computationally more demanding. \n"
          "The SA engine dynamically adjusts the temperature according to the ratio of accepted steps. It tries the match the"
          "current ratio of accepted steps to the one calculated between given start ratio and end ratio and using the"
          "number of steps as a linear predictor. It only tries to adjust the temperature every Nth step where N is"
          "specified by the user. If the ratio is lower than the expected one, it increases the temperature, otherwise it"
          "decreases the temperature. For cases where the actual ratio is close to expected ratio ( within 0.1) it"
          "adjusts the temperature by multiplying/diving by a small coefficient, while for cases where the difference"
          "is larger, the temperature is adjusted by multiplying/diving by a larger coefficient.",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "temp_accept_start",
              "fraction of MCM moves to be accepted initially",
              command::ParameterCheckRanged< float>( 0.0, 1.0),
              "0.90"
            ),
            command::Parameter
            (
              "temp_accept_end",
              "fraction of MCM moves to be accepted by the end",
              command::ParameterCheckRanged< float>( 0.0, 1.0),
              "0.10"
            ),
            command::Parameter
            (
              "initial_temp",
              "starting temperature",
              command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()),
              "1.0"
            ),
            command::Parameter
            (
              "steps_per_update",
              "number of steps between ",
              command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
              "10"
            )
          )
        )
      ),
      m_LargerIsBetterFlag
      (
        new command::FlagStatic
        (
          "larger_score_better",
          "sets all MCM objects to increase the score during optimization; "
          "default behavior attempts to decrease the score during optimization"
        )
      ),
      m_SaveAllAcceptedImprovedFlag
      (
        new command::FlagStatic
        (
          "save_all_accepted_improved",
          "save all molecules that are accepted or improved according to the main (outer) MCM; "
          "overall distribution of molecule scores will be skewed worse, but intermediate "
          "structures will be available for analysis"
        )
      ),
      m_Corina
      (
        new command::FlagStatic
        (
          "corina",
          "make a system call to Corina to make the starting conformer during molecule cleaning;"
          "this means that if only 1 conformer is desired (i.e. pose-independent scoring) it will be the corina default conformer,"
          "while if multiple conformers are desired (i.e. refinement phase of pose-dependent scoring) there will be no effect; "
          "this option is meant primarily to allow backward compatibility for QSAR models generated with Corina conformers"
        )
      ),
      m_VDWClashCutoffFlag
      (
        new command::FlagStatic
        (
          "conf_vdw_cutoff",
          "maximum Van der Waals score of a valid 3D conformer",
          command::Parameter
          (
            "vdw_score",
            "Internal Van der Waals score of the molecule conformer normalized by the number of atoms in the molecule; "
            "computed with MoleculeVDWScore; "
            "note that ligand-based methods do not perform much conformer optimization (for efficiency purposes), and "
            "may have high VDW scores even if the molecule is consistitutionally acceptable (especially if generated "
            "with Corina). Thus, for ligand-based design tasks we recommend increasing the VDW score cutoff to at least "
            "2.0 - 5.0, though higher can also be appropriate. For structure-based pose-dependent optimization lower "
            "cutoff scores can be used, and this may depend on the number of iterations used for pose optimization "
            "after the generation of each molecule.",
            command::ParameterCheckRanged< float>( 0.0, 1000.0), "5.0"
          )
        )
      )
    {
    }

    const ApplicationType FocusedLibraryDesign::FocusedLibraryDesign_Instance
    (
      GetAppGroups().AddAppToGroup( new FocusedLibraryDesign(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
