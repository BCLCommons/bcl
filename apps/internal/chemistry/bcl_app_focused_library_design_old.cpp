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

// include header for this application
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "internal/chemistry/bcl_app_focused_library_design_old.h"
#include "io/bcl_io_file.h"
#include "opti/bcl_opti_criterion_function.h"
#include "opti/bcl_opti_criterion_skipped_steps.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace app
  {

    namespace
    {
      class FragmentMutateReact :
        public math::MutateInterface< chemistry::FragmentComplete>
      {
      private:

        std::string m_ReactionFilename;
        std::string m_ReactantFilename;

        chemistry::FragmentReact m_Reactor;

      public:
        FragmentMutateReact()
        {
        }

        FragmentMutateReact *Clone() const
        {
          return new FragmentMutateReact( *this);
        }

        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

        const std::string &GetAlias() const
        {
          static const std::string s_name( "React");
          return s_name;
        }

        math::MutateResult< chemistry::FragmentComplete> operator()( const chemistry::FragmentComplete &ARGUMENT) const
        {

          math::MutateResult< chemistry::FragmentComplete> ret_mutate( util::CloneToShPtr( ARGUMENT), *this);

          storage::Pair< util::SiPtr< const chemistry::ReactionComplete>, chemistry::FragmentEnsemble> res_pair
          (
            m_Reactor.ReactRandom( ARGUMENT)
          );
          chemistry::FragmentEnsemble &mols( res_pair.Second());
          size_t n_mols( mols.GetSize());

          if( !res_pair.First().IsDefined() || !n_mols || !m_Reactor.IsGood())
          {
            return ret_mutate;
          }

          chemistry::FragmentEnsemble::const_iterator itr( mols.Begin());
          std::advance
          (
            itr,
            n_mols > 1 ? random::GetGlobalRandom().Random< size_t>( 0, n_mols - 1) : 0
          );

          util::ShPtr< chemistry::FragmentComplete> cleaned_mol( chemistry::FragmentEvolveBase::FinalizeMolecule( *itr));
          std::string rxn_desc( res_pair.First()->GetDescription());
          for( size_t i( 0), l( rxn_desc.length()); i < l; ++i)
          {
            if( rxn_desc[ i] == '\n')
            {
              rxn_desc[ i] = ' ';
            }
          }

          if( !cleaned_mol.IsDefined())
          {
            return ret_mutate;
          }

          cleaned_mol->StoreProperty( "FocusedLibraryDesignOldRxn", rxn_desc);

          ret_mutate = math::MutateResult< chemistry::FragmentComplete>( cleaned_mol, *this);
          return ret_mutate;
        }

        bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
        {
          m_Reactor = chemistry::FragmentReact( chemistry::ReactionSearch( m_ReactantFilename, m_ReactionFilename));
          return true;
        }

        io::Serializer GetSerializer() const
        {
          io::Serializer member_data;
          member_data.SetClassDescription( "randomly mutates a molecule using a set of reactions");
          member_data.AddInitializer
          (
            "reaction filename",
            "set the file to pull reactions from",
            io::Serialization::GetAgent( &m_ReactionFilename)
          );
          member_data.AddInitializer
          (
            "reactant filename",
            "set the file to pull reactants from",
            io::Serialization::GetAgent( &m_ReactantFilename)
          );
          return member_data;
        }

      };

      class CriterionMoleculeTooBig :
        public opti::CriterionInterface< chemistry::FragmentComplete, float>
      {

        //! @brief copies this class
        //! @return a pointer to a copy of this class
        CriterionMoleculeTooBig *Clone() const
        {
          return new CriterionMoleculeTooBig( *this);
        }

        //! @brief Get the name of this class
        //! @return the name of this class
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( this);
        }

        //! @brief get the class name when used in a dynamic context
        //! @return the class name when used in a dynamic context
        const std::string &GetAlias() const
        {
          static const std::string s_alias( "MoleculeTooBig");
          return s_alias;
        }

        bool CriteriaMet( const opti::Tracker< chemistry::FragmentComplete, float> &TRACKER) const
        {
          descriptor::CheminfoProperty weight( "Weight");
          float mol_wt( weight->SumOverObject( TRACKER.GetCurrent()->First())( 0));
          if( mol_wt > 700)
          {
            return true;
          }
          return false;
        }

        io::Serializer GetSerializer() const
        {
          return io::Serializer();
        }

      };
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief Checks if a worker should be updated/continue
    //! @param WORKER reference to the worker that should be updated
    //! @return true if the worker should be updated and continue, false otherwise
    bool FocusedLibraryDesignOld::ThreadManager::WorkerShouldContinue( ThreadManager::Worker &WORKER)
    {
      // Lock structure during the modification
      m_WorkerMutex.Lock();

      bool continue_working( true);
      if( m_WorkersShouldStop)
      {
        BCL_MessageDbg( "Thread " + util::Format()( WORKER.m_ThreadNumber) + " should exit");
        continue_working = false;
      }
      else
      {
        BCL_MessageDbg( "Thread " + util::Format()( WORKER.m_ThreadNumber) + " should continue");
      }

      m_WorkerMutex.Unlock();
      return continue_working;
    } // WorkerShouldContinue()

    //! @brief Constructor
    FocusedLibraryDesignOld::ThreadManager::ThreadManager
    (
      const chemistry::FragmentEnsemble &INITIAL_MOLS,
      const mc::Approximator< chemistry::FragmentComplete, float> &APPROXIMATOR,
      const float &SCORE_CUTOFF,
      const size_t &REQUESTED_MOLECULES,
      const size_t &NUMBER_THREADS,
      const size_t &MOL_WRITE_SIZE,
      const std::string &OUTPUT_FILENAME
    ) :
      m_NumMolsRequested( REQUESTED_MOLECULES),
      m_NumMolsBuilt( 0),
      m_NumThreads( NUMBER_THREADS),
      m_Molecules(),
      m_ConfigurationSet(),
      m_WorkersShouldStop( false),
      m_NumMolsRunningTotal( 0),
      m_ScoreCutoff( SCORE_CUTOFF),
      m_LastMolWritten( m_Molecules.Begin()),
      m_NumMolsWritten( 0),
      m_OutputFilename( OUTPUT_FILENAME)
    {

      // initialize the output file by clearing it
      io::OFStream out( m_OutputFilename.c_str(), std::ios::out);
      io::File::CloseClearFStream( out);

      // for each molecule in the initial mols ensemble run the FLD algorithm
      size_t mol_no( 0);
      for
      (
        chemistry::FragmentEnsemble::const_iterator itr_init( INITIAL_MOLS.Begin()), itr_init_end( INITIAL_MOLS.End());
        itr_init != itr_init_end;
        ++itr_init, ++mol_no
      )
      {
        BCL_MessageStd( "Generating derivatives of molecule " + util::Format()( mol_no));

        // Reset everything for a new run
        m_Molecules = chemistry::FragmentEnsemble();
        m_LastMolWritten = m_Molecules.Begin();
        m_NumMolsBuilt = 0;
        m_ConfigurationSet = chemistry::ConfigurationSet();
        m_WorkersShouldStop = false;
        m_NumMolsWritten = 0;

        // Set up workers
        size_t thread_id( 0);
        std::vector< Worker> workers( m_NumThreads);
        for
        (
          std::vector< Worker>::iterator itr( workers.begin()), itr_end( workers.end());
          itr != itr_end;
          ++itr, ++thread_id
        )
        {
          Worker &worker_ref( *itr);
          worker_ref.m_ThreadManager        = this;
          worker_ref.m_StartMolecule        = *itr_init;
          worker_ref.m_StartMoleculeID      = mol_no;
          worker_ref.m_Approximator         = util::ShPtr< opti::ApproximatorModularBase< chemistry::FragmentComplete, float> >
                                              (
                                                new mc::Approximator< chemistry::FragmentComplete, float>( APPROXIMATOR)
                                              );
          worker_ref.m_ThreadNumber         = thread_id;
        }

        if( m_NumThreads > 1)
        {
          // Allocate space for jobs
          util::ShPtrVector< sched::JobInterface> jobs;
          jobs.AllocateMemory( m_NumThreads);

          const size_t group_id( 1);

          // Run jobs
          for( size_t proc_number( 0); proc_number < m_NumThreads; ++proc_number)
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
            BCL_MessageStd( "Submitted thread " + util::Format()( proc_number));
          }

          // Join all the jobs
          for( size_t proc_number( 0); proc_number < m_NumThreads; ++proc_number)
          {
            sched::GetScheduler().Join( jobs( proc_number));
          }
        }
        else if( m_NumThreads == 1)
        {
          BCL_MessageStd( "Running job");
          workers[ 0].RunThread();
        }

        // force a write to the output file
        WriteMols( true);
      }

    } // ThreadManager()

    //! @brief Add a molecule to the ensemble, and write it to the output file.  This is also the point to check
    //! @param MOLECULE Molecule that should be added
    //! @return True if addition/writing was done
    bool FocusedLibraryDesignOld::ThreadManager::ThreadAddMolecule
    (
      const chemistry::FragmentComplete &MOLECULE,
      const size_t &START_MOL_ID,
      const float &SCORE,
      const size_t &TOTAL_MOLS_BUILT
    )
    {
      // If the molecule has a good geometry, an acceptable score, and is unique, write it to file
      m_WorkerMutex.Lock();
      bool saved_mol( false);

      m_NumMolsRunningTotal += TOTAL_MOLS_BUILT;

      chemistry::FragmentComplete mutable_mol( MOLECULE);
      mutable_mol.SaturateWithH();
      mutable_mol.Translate( -mutable_mol.GetCenter());
      mutable_mol.StoreProperty( "FocusedLibraryDesignOldScore", util::Format()( SCORE));
      mutable_mol.StoreProperty( "FocusedLibraryDesignOldInitialIndex", util::Format()( START_MOL_ID));

      if
      (
        !mutable_mol.HasBadGeometry()
//        && SCORE >= m_ScoreCutoff
        && m_ConfigurationSet.Insert( chemistry::FragmentConfigurationShared( mutable_mol)).second
      )
      {
        BCL_MessageStd( "Generated a unique molecule");
        m_Molecules.PushBack( mutable_mol);

        ++m_NumMolsBuilt;
        std::cerr << "Generated a new molecule.  Total: " << m_NumMolsBuilt << " generated" << std::endl;

        saved_mol = true;
      }

      if( m_NumMolsBuilt >= m_NumMolsRequested)
      {
        m_WorkersShouldStop = true;
      }

      // write mols (if necessary)
      WriteMols();

      m_WorkerMutex.Unlock();
      return saved_mol;
    } // AddMolecule

    void FocusedLibraryDesignOld::ThreadManager::WriteMols( const bool &FLUSH) const
    {
      // if write is forced, or enough mols have been generated, write to file
      if( FLUSH || ( m_NumMolsBuilt - m_NumMolsWritten) >= m_NumMolsBeforeWrite)
      {
        io::OFStream out( m_OutputFilename.c_str(), std::ios::app); // append to the output file

        chemistry::FragmentEnsemble::const_iterator itr_mol( m_NumMolsWritten == 0 ? m_Molecules.Begin() : m_LastMolWritten),
          itr_mol_end( m_Molecules.End());

        // move to the first unwritten molecule
        if( m_NumMolsWritten && itr_mol != itr_mol_end)
        {
          ++itr_mol;
        }

        for( ; itr_mol != itr_mol_end; ++itr_mol)
        {
          itr_mol->WriteMDL( out);
          m_LastMolWritten = itr_mol;
          ++m_NumMolsWritten;
        }
      }
    }

    //! @brief Fetches the molecule ensemble
    //! @return Reference to FragmentEnsemble of the generated molecules
    const chemistry::FragmentEnsemble &FocusedLibraryDesignOld::ThreadManager::GetMolecules() const
    {
      return m_Molecules;
    }
    //! @brief initializes the command object for this application
    //! @return a ShPtr to a Command containing all of this applications parameters
    util::ShPtr< command::Command> FocusedLibraryDesignOld::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // insert all the flags and params

      // flags for input/output
      sp_cmd->AddFlag( m_InitialMoleculesFlag);
      sp_cmd->AddFlag( m_MutateFunctionFlag);
      sp_cmd->AddFlag( m_OutputFilenameFlag);
      sp_cmd->AddFlag( m_NumberMoleculesFlag);
      sp_cmd->AddFlag( m_ScoreCutoffFlag);
      sp_cmd->AddFlag( m_ScoringFunctionFlag);
      sp_cmd->AddFlag( m_FunctionalizationPointsFlag);
      sp_cmd->AddFlag( m_FunctionalGroupsFileFlag);
      sp_cmd->AddFlag( m_OutputInitialFlag);

      // default flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

    void FocusedLibraryDesignOld::ThreadManager::Worker::RunThread()
    {
      do
      {
        BCL_MessageVrb( "Thread " + util::Format()( m_ThreadNumber) + ": optimizing molecule");

        // Run the MC optimizer (generate molecules)
        m_Approximator->GetTracker().Reset();
        m_Approximator->GetTracker().SetPhase( opti::e_Start);
        m_Approximator->GetTracker().Track
        (
          util::ShPtr< storage::Pair< chemistry::FragmentComplete, float> >
          (
            new storage::Pair< chemistry::FragmentComplete, float>
            (
              m_StartMolecule,
              float( 0.0)
            )
          )
        );
        m_Approximator->Approximate();
        storage::Pair< chemistry::FragmentComplete, float> mutation( *m_Approximator->GetTracker().GetBest());
        //storage::Pair< chemistry::FragmentComplete, float> mutation_cur( *m_Approximator->GetTracker().GetCurrent());

        BCL_MessageVrb( "Thread " + util::Format()( m_ThreadNumber) + ": done optimizing molecule");

        m_ThreadManager->ThreadAddMolecule( mutation.First(), m_StartMoleculeID, mutation.Second(), m_Approximator->GetTracker().GetIteration());
        //m_ThreadManager->ThreadAddMolecule( mutation_cur.First(), m_StartMoleculeID, mutation_cur.Second(), m_Approximator->GetTracker().GetIteration());
      }
      while( m_ThreadManager->WorkerShouldContinue( *this)); // Update until WorkerShouldContinue returns false

      BCL_MessageVrb( "Thread " + util::Format()( m_ThreadNumber) + " exiting");
    } // RunThread()

    void FocusedLibraryDesignOld::EnumerateSites
    (
      const chemistry::AtomVector< chemistry::AtomComplete> &BASE,
      const storage::Vector< size_t> &SITES,
      const util::SiPtrVector< const chemistry::AtomVector< chemistry::AtomComplete> > &FRAGS,
      storage::Vector< chemistry::AtomVector< chemistry::AtomComplete> > &ACCUM
    ) const
    {
      if( BASE.IsEmpty() || SITES.IsEmpty() || FRAGS.IsEmpty())
      {
        return;
      }

      size_t site( SITES.LastElement());

      storage::Vector< size_t> new_sites( SITES);
      new_sites.PopBack(); // remove the last element from the vector

      // add each fragment to the chosen site
      storage::Vector< chemistry::AtomVector< chemistry::AtomComplete> > tempv( FRAGS.GetSize(), BASE);

      for( size_t i( 0), last( tempv.GetSize()); i < last; ++i)
      {
        if( !FRAGS( i).IsDefined())
        {
          continue;
        }
        AddSubToMolecule( tempv( i), *FRAGS( i), site);

        // also add fragments to any remaining sites
        if( !new_sites.IsEmpty())
        {
          EnumerateSites( tempv( i), new_sites, FRAGS, tempv);
        }
      }

      ACCUM.Append( tempv);

      /*
      for( size_t frag_no( 0), last_frag( FRAGS.GetSize()); frag_no < last_frag; ++frag_no)
      {
        util::SiPtrVector< const chemistry::AtomVector< chemistry::AtomComplete> > new_frags( FRAGS);
        new_frags.RemoveElements( frag_no, 1);
        chemistry::AtomVector< chemistry::AtomComplete> new_base( BASE);
        AddSubToMolecule( new_base, *FRAGS( frag_no), SITES( 0));
        ACCUM.PushBack( new_base);
        EnumerateSites( new_base, new_sites, new_frags, ACCUM);
      }
      */
    }

    void FocusedLibraryDesignOld::AddSubToMolecule
    (
      chemistry::AtomVector< chemistry::AtomComplete> &PARENT_ATOMS,
      const chemistry::AtomVector< chemistry::AtomComplete> &SUB_ATOMS,
      const int &SITE
    ) const
    {
      // find open valences.  hopefully nobody will provide more than one of these or we'll have to write more code...
      const chemistry::AtomVector< chemistry::AtomComplete> &grp_atms( SUB_ATOMS);
      size_t attach_atom( grp_atms.GetSize());
      // determine which atom in the functional gets attached to the input structure
      for( size_t g_atom_no( 0); g_atom_no < grp_atms.GetSize(); ++g_atom_no)
      {
        size_t n_val( grp_atms( g_atom_no).GetNumberValenceBonds());
        if( n_val && util::IsDefined( n_val))
        {
          attach_atom = g_atom_no; // attach at the first open valence that is found
          break;
        }
      }

      if( attach_atom >= grp_atms.GetSize()) // no attachment could be found
      {
        return;
      }

      // new bond to the functional group
      storage::Vector< sdf::BondInfo> bonds;
      bonds.PushBack
      (
        sdf::BondInfo
        (
          SITE,
          attach_atom + PARENT_ATOMS.GetSize(),
          chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond
        )
      );

      PARENT_ATOMS.AddAtomsWithConnectivity( grp_atms, bonds); // connect the functional group atoms to the original mol
    }

    //! @brief the Main function
    //! @return 0 for success
    int FocusedLibraryDesignOld::Main() const
    {

      chemistry::FragmentEnsemble input_mols;
      io::IFStream in_fs;
      io::File::MustOpenIFStream( in_fs, m_InitialMoleculesFlag->GetFirstParameter()->GetValue());
      input_mols.ReadMoreFromMdl( in_fs);
      io::File::CloseClearFStream( in_fs);
      //input_mols.SaturateWithH();

      if( input_mols.GetSize() > 1)
      {
        BCL_MessageCrt( "WARNING: Only the first molecule from the input ensemble will be used for molecule design");
        BCL_MessageCrt( "         If you want all of them used, re-run this program with each molecule in a different file");
      }

      BCL_Assert( !input_mols.IsEmpty(), "No molecules could be read from input.");

      // copy the first molecule
      chemistry::FragmentEnsemble initial_mols;
      initial_mols.PushBack( *input_mols.Begin());
      size_t n_atoms( initial_mols.Begin()->GetNumberAtoms());
      BCL_Assert( n_atoms, "Input molecule is empty, exitting...");

      storage::Vector< size_t> fxnl_pts; // functionalized points
      if( m_FunctionalizationPointsFlag->GetFlag())
      {
        // ensure that functionalization sites and a set of functional groups are both provided
        BCL_Assert
        (
          m_FunctionalGroupsFileFlag->GetFlag(),
          "Using " + m_FunctionalizationPointsFlag->GetName() + " requires "
            + m_FunctionalGroupsFileFlag->GetName() + " to be specified too"
        );

        // convert the functionalization points to numeric values
        storage::Vector< std::string> fxnl_pts_str( m_FunctionalizationPointsFlag->GetStringList());
        storage::Set< size_t> fxnl_pts_set;

        for( size_t i( 0), l( fxnl_pts_str.GetSize()); i < l; ++i)
        {
          size_t point;
          if( !util::TryConvertFromString( point, fxnl_pts_str( i), util::GetLogger()))
          {
            BCL_MessageStd( "Could not parse \"" + fxnl_pts_str( i) + "\" as a number");
            continue;
          }
          if( point < n_atoms)
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

        fxnl_pts = storage::Vector< size_t>( fxnl_pts_set.Begin(), fxnl_pts_set.End());

      }

      if( m_FunctionalGroupsFileFlag->GetFlag())
      {
        BCL_Assert
        (
          input_mols.GetSize() == 1,
          "If functionalizations are specified, only one molecule may be given in the input file."
        );

        // Try to determine functionalization points based on lack of hydrogens
        if( fxnl_pts.IsEmpty())
        {
          chemistry::CollectorValence valence_collector;
          util::SiPtrList< const chemistry::AtomConformationalInterface> open_atoms( valence_collector.Collect( *initial_mols.Begin()));
          for
          (
            util::SiPtrList< const chemistry::AtomConformationalInterface>::const_iterator itr_atom( open_atoms.Begin()), itr_atom_end( open_atoms.End());
            itr_atom != itr_atom_end;
            ++itr_atom
          )
          {
            fxnl_pts.PushBack( initial_mols.Begin()->GetAtomIndex( **itr_atom));
          }
        }

        BCL_Assert
        (
          !fxnl_pts.IsEmpty(),
          "Could not determine where to functionalize the input molecule.  Please either specify it with the "
          + m_FunctionalizationPointsFlag->GetName() + " flag or remove hydrogens from the input molecule"
        );

        // sort functionalization points so we can remove hydrogens sensibly
        fxnl_pts.Sort( std::less< size_t>());

        // read reactive functional groups
        io::IFStream fxn_file;
        io::File::MustOpenIFStream( fxn_file, m_FunctionalGroupsFileFlag->GetFirstParameter()->GetValue());
        chemistry::FragmentEnsemble fxnl_grp( fxn_file);
        io::File::CloseClearFStream( fxn_file);

        BCL_Assert( !fxnl_grp.IsEmpty(), "No functional groups could be read from given file");

        // grab the first molecule and then reset the initial_mols ensemble
        // we will add functionalized version of the molecule to initial_mols when we're done
        chemistry::FragmentComplete in_mol( *initial_mols.Begin());
        initial_mols = chemistry::FragmentEnsemble();
        chemistry::AtomVector< chemistry::AtomComplete> mol_atoms( in_mol.GetAtomVector());
        size_t n_atoms( in_mol.GetNumberAtoms());

        // now remove hydrogens and determine how much each functional point value needs to be offset to compensate
        size_t tot_h( 0);
        for( size_t ano( 0), cur_pt( 0); ano < n_atoms && cur_pt < fxnl_pts.GetSize(); ++ano)
        {
          // here we have counted up the number of hydrogens that occurred before our functional point
          // remove this quanity from our functional point value, then move to the next functional point
          if( ano == fxnl_pts( cur_pt))
          {
            fxnl_pts( cur_pt) -= tot_h;
            ++cur_pt;
            continue;
          }

          if( mol_atoms( ano).GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
          {
            ++tot_h;
          }
        }

        // remove hydrogens
        in_mol.RemoveH();
        mol_atoms = in_mol.GetAtomVector();

        if( fxnl_pts.GetSize() > 2)
        {
          BCL_MessageStd( "Only 2 sites are handled now, discarding all but the first two");
          fxnl_pts.Resize( 2);
        }

        storage::Vector< chemistry::AtomVector< chemistry::AtomComplete> > grp_atoms;
        util::SiPtrVector< const chemistry::AtomVector< chemistry::AtomComplete> > grp_atoms_ptrs;
        for
        (
          chemistry::FragmentEnsemble::const_iterator itr_grp( fxnl_grp.Begin()), itr_grp_end( fxnl_grp.End());
          itr_grp != itr_grp_end;
          ++itr_grp
        )
        {
          grp_atoms.PushBack( itr_grp->GetAtomVector());
        }

        for( size_t i( 0), last( grp_atoms.GetSize()); i < last; ++i)
        {
          grp_atoms_ptrs.PushBack( &( grp_atoms( i)));
        }

        // loop over sites to make single substitutions, then add fragments to each other grow point
        storage::Vector< chemistry::AtomVector< chemistry::AtomComplete> > accum_v;
        storage::Vector< size_t>::const_iterator itr_beg( fxnl_pts.Begin()), itr_end( fxnl_pts.End());
        for
        (
          storage::Vector< size_t>::const_iterator itr_pt( itr_beg), itr_pt_end( itr_end);
          itr_pt != itr_pt_end;
          ++itr_pt
        )
        {
          storage::Vector< size_t> single_pt( size_t( 1), *itr_pt);
          storage::Vector< size_t> other_pts( itr_beg, itr_pt);
          storage::Vector< size_t>::const_iterator itr_next( itr_pt);
          ++itr_next;
          for( ; itr_next != itr_end; ++itr_next)
          {
            other_pts.PushBack( *itr_next);
          }
          //EnumerateSites( mol_atoms, fxnl_pts, grp_atoms_ptrs, accum_v);
          EnumerateSites( mol_atoms, single_pt, grp_atoms_ptrs, accum_v);
          EnumerateSites( accum_v.LastElement(), other_pts, grp_atoms_ptrs, accum_v);
        }

        chemistry::ConfigurationSet const_set;
        std::string mol_name;
        for( size_t v( 0); v < accum_v.GetSize(); ++v)
        {
          chemistry::AtomsCompleteStandardizer atoms_standardizer( accum_v( v), mol_name, true);
          chemistry::BondIsometryHandler::AddIsometryInformation( accum_v( v), true);
          chemistry::StereocentersHandler::AddChiralityFromConformation( accum_v( v));

          chemistry::FragmentComplete new_mol( accum_v( v), "");
          util::ShPtr< chemistry::FragmentComplete> finalized_mol( chemistry::FragmentEvolveBase::FinalizeMolecule( new_mol));
          if( finalized_mol.IsDefined() && const_set.Insert( chemistry::FragmentConfigurationShared( *finalized_mol)).second)
          {
            initial_mols.PushBack( *finalized_mol);
          }
        }
      }

      initial_mols.SaturateWithH();
      if( m_OutputInitialFlag->GetFlag())
      {
        io::OFStream out_init;
        io::File::MustOpenOFStream( out_init, m_OutputInitialFlag->GetFirstParameter()->GetValue());
        initial_mols.WriteMDL( out_init);
        io::File::CloseClearFStream( out_init);
      }

    ///////////////////////
    // Molecule building //
    ///////////////////////

      // Time the molecule building
      util::Stopwatch threadmanager_timer( "Molecule Building", util::Time( 1, 0), util::Message::e_Standard, true, false);
      threadmanager_timer.Start();

      FitnessFunctionRaw obj_function;
      obj_function.AssertRead( m_ScoringFunctionFlag->GetFirstParameter()->GetValue());

      FragmentMutateReact mutate_function;
      mutate_function.AssertRead( m_MutateFunctionFlag->GetFirstParameter()->GetValue());

      opti::CriterionCombine< chemistry::FragmentComplete, float> criteria;
      criteria.InsertCriteria( opti::CriterionNumberIterations< chemistry::FragmentComplete, float>( 10));
      criteria.InsertCriteria( CriterionMoleculeTooBig());

      util::ShPtr< mc::TemperatureInterface> temperature( new mc::TemperatureDefault( 100));
      mc::Metropolis< float> metropolis( temperature);

      mc::Approximator< chemistry::FragmentComplete, float> approximator
      (
        obj_function,
        mutate_function,
        metropolis,
        criteria
      );

      size_t num_mols_requested( m_NumberMoleculesFlag->GetFirstParameter()->GetNumericalValue< float>());

      // Build the molecules using metropolis Monte-Carlo
      ThreadManager thread_manager
      (
        initial_mols,
        approximator,
        m_ScoreCutoffFlag->GetFirstParameter()->GetNumericalValue< float>(),
        num_mols_requested,    // Number of molecules to output
        sched::GetNumberCPUs(),                                                     // Number of CPUs to use
        std::max( size_t( 1), size_t( num_mols_requested / 10)),
        m_OutputFilenameFlag->GetFirstParameter()->GetValue() // Output filename
      );

      threadmanager_timer.Stop();

      //BCL_MessageStd( "Generated a total of " + util::Format()( thread_manager.GetNumberGeneratedMolecules()) + " molecules while running");

      // Fetch the list of molecules that were generated
      chemistry::FragmentEnsemble generated_molecules
      (
        thread_manager.GetMolecules()
      );

      return 0;
    } // Main

    //! @brief standard constructor
    FocusedLibraryDesignOld::FocusedLibraryDesignOld() :
      m_InitialMoleculesFlag
      (
        new command::FlagStatic
        (
          "initial_molecules", "file containing the starting pieces",
          command::Parameter
          (
            "filename",
            "filename containing a molecule or set of molecules to begin optimization from.  "
            "a full optimization run is performed on each of the molecules, "
            "i.e. the specified number of molecules will be generated for every initial molecule provided",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_MutateFunctionFlag
      (
        new command::FlagDynamic
        (
          "grow_method",
          "the mutation method for the algorithm",
          command::Parameter
          (
            "method",
            "the method used to modify molecules",
            command::ParameterCheckSerializable( FragmentMutateReact())
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output",
          "prefix of where to output molecules",
          command::Parameter
          (
            "output", "filename for where to write generated molecules"
          )
        )
      ),
      m_NumberMoleculesFlag
      (
        new command::FlagStatic
        (
          "number_molecules",
          "flag for number of molecules to generate",
          command::Parameter
          (
            "number_molecules", "total number of molecules",
            command::ParameterCheckRanged< int>( 0, std::numeric_limits< int>::max()), "50"
          )
        )
      ),
      m_ScoreCutoffFlag
      (
        new command::FlagStatic
        (
          "score_cutoff",
          "if a molecule's score is below this value, it is not kept",
          command::Parameter
          (
            "cutoff",
            "the cutoff value used for filtering undesirable molecules"
          )
        )
      ),
      m_ScoringFunctionFlag
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
              FitnessFunctionRaw()
            )
          )
        )
      ),
      m_FunctionalizationPointsFlag
      (
        new command::FlagDynamic
        (
          "functionalize",
          "which atoms (zero-indexed) of the input molecule should be functionalized with reactive groups",
          command::Parameter
          (
            "atoms",
            "the atom indices to modify during preprocessing",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max())
          )
        )
      ),
      m_FunctionalGroupsFileFlag
      (
        new command::FlagDynamic
        (
          "functional_groups",
          "reactive functional groups to add to the input molecule before starting the optimization",
          command::Parameter
          (
            "groups",
            "file containing functional groups to add to the molecule during preprocessing",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_OutputInitialFlag
      (
        new command::FlagDynamic
        (
          "output_init",
          "filename to output initial molecules to (good if automatic functionalization is done)",
          command::Parameter
          (
            "filename",
            "the filename to output the initial molecules to (sdf)"
          )
        )
      )
    {
    }

    // Construct the static instance of this application, and add it to the ChemInfo group
    const ApplicationType FocusedLibraryDesignOld::FocusedLibraryDesignOld_Instance
    (
      GetAppGroups().AddAppToGroup( new FocusedLibraryDesignOld(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app

  template<>
  const util::SiPtr< const util::ObjectInterface>
    opti::CriterionCombine< chemistry::FragmentComplete, float>::s_Instance
  (
    util::Enumerated< opti::CriterionCombine< chemistry::FragmentComplete, float> >::AddInstance
    (
      new opti::CriterionCombine< chemistry::FragmentComplete, float>()
    )
  );
} // namespace bcl

