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

#ifndef BCL_APP_FOCUSED_LIBRARY_DESIGN_OLD_H_
#define BCL_APP_FOCUSED_LIBRARY_DESIGN_OLD_H_
// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_collector_valence.h"
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_react.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_template_instantiations.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_default.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_stopwatch.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FocusedLibraryDesignOld
    //! @brief Application for generating libraries for synthesis using QSAR models and random structure generator
    //!
    //! @author loweew, geanesar
    //! @date 08/31/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class FocusedLibraryDesignOld :
      public Interface
    {
    public:

    private:

      friend class ThreadManager;
      friend class Worker;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class FitnessFunctionRaw
      //! @brief calculates a molecule's fitness based on a descriptor
      //!
      //! @author geanesar
      //! @date 05/12/2014
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class FitnessFunctionRaw :
        public math::FunctionInterfaceSerializable< chemistry::FragmentComplete, float>
      {

      private:

        // the descriptor to use
        descriptor::CheminfoProperty m_Descriptor;

        // a static instance of this class
        static const util::SiPtr< util::ObjectInterface> s_Instance;

      public:

        //! @brief default constructor
        FitnessFunctionRaw() :
          m_Descriptor()
        {
        }

        //! @brief constructor with parameters
        //! @param DESCRIPTOR the descriptor to use
        FitnessFunctionRaw
        (
          const descriptor::CheminfoProperty &DESCRIPTOR
        ) :
          m_Descriptor( DESCRIPTOR)
        {
        }

        //! @brief copies this class
        //! @return a pointer to a copy of this class
        FitnessFunctionRaw *Clone() const
        {
          return new FitnessFunctionRaw( *this);
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
          static const std::string s_alias( "Descriptor");
          return s_alias;
        }

        //! @brief calculates the descriptor for the molecule
        //! @param FRAGMENT the fragment to calculate the descriptor for
        //! @return the mean value of the descriptor
        float operator()( const chemistry::FragmentComplete &FRAGMENT) const
        {
          /*io::OFStream out( "/tmp/last_mol.sdf", std::ios::out);
          FRAGMENT.WriteMDL( out);
          io::File::CloseClearFStream( out); */
          chemistry::FragmentComplete saturated_frag( FRAGMENT);
          saturated_frag.SaturateWithH();
          linal::Vector< float> res( m_Descriptor->SumOverObject( saturated_frag));
          return res.Sum() / res.GetSize();
        }

        //! @brief gets a serializer for constructing this class in a dynamic context
        //! @return the serializer containing member data
        io::Serializer GetSerializer() const
        {
          io::Serializer member_data;

          member_data.SetClassDescription( "scores molecules using the raw mean output from a descriptor");

          member_data.AddInitializer
          (
            "descriptor",
            "the descriptor to calculate.  if multi-valued, this will return the mean value.",
            io::Serialization::GetAgent( &m_Descriptor)
          );

          return member_data;
        }

      }; // FitnessFunctionRaw

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ThreadManager
      //! @brief Manages threads for multithreaded structure generation
      //!
      //! @author geanesar
      //! @date Nov 7, 2013
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class ThreadManager :
        public util::ObjectInterface
      {
      private:

      ///////////
      // Data //
      ///////////

        //! Number of molecules that should be made
        const size_t m_NumMolsRequested;

        //! Number of molecules that have been built
        size_t m_NumMolsBuilt;

        //! Number of threads that should be used
        const size_t m_NumThreads;

        //! Ensemble of (unique) molecules that have been made
        chemistry::FragmentEnsemble m_Molecules;

        //! ConfigurationSet used to determine if two structures are configurationally identical
        chemistry::ConfigurationSet m_ConfigurationSet;

        //! A variable used to check if workers (i.e. threads) should stop working
        bool m_WorkersShouldStop;

        //! Mutex used for adding molecules to the molecule list
        mutable sched::Mutex m_WorkerMutex;

        //! The number of molecules that were generated in total (sum of all monte-carlo steps)
        size_t m_NumMolsRunningTotal;

        //! The score cutoff to accept molecules
        float m_ScoreCutoff;

        //! how many mols should be generated before writing to output
        size_t m_NumMolsBeforeWrite;

        //! last molecule that was written to a file
        mutable chemistry::FragmentEnsemble::const_iterator m_LastMolWritten;

        //! index of the last molecule written
        mutable size_t m_NumMolsWritten;

        //! the filename to output molecules to
        std::string m_OutputFilename;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //!
        //! @class Worker
        //! @brief Helper struct used to run threads for the ThreadManager, builds molecules using Monte-Carlo routines
        //!
        //! @author geanesar
        //! @date Nov 7, 2013
        //!
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        struct Worker
        {

          util::SiPtr< ThreadManager> m_ThreadManager; //!< pointer to the thread manager that controls this
          chemistry::FragmentComplete m_StartMolecule;
          size_t m_StartMoleculeID;
          util::ShPtr< opti::ApproximatorModularBase< chemistry::FragmentComplete, float> > m_Approximator; //!< the approximator that does the work
          size_t m_ThreadNumber; // Numerical ID of the thread

          //! @brief run the thread
          void RunThread();
        }; // stuct Worker

        //! @brief Checks if a worker should be updated/continue
        //! @param WORKER reference to the worker that should be updated
        //! @return true if the worker should be updated and continue, false otherwise
        bool WorkerShouldContinue( ThreadManager::Worker &WORKER);

      public:

        //! @brief Constructor
        ThreadManager
        (
          const chemistry::FragmentEnsemble &INITIAL_MOLS,
          const mc::Approximator< chemistry::FragmentComplete, float> &APPROXIMATOR,
          const float &SCORE_CUTOFF,
          const size_t &REQUESTED_MOLECULES,
          const size_t &NUMBER_THREADS,
          const size_t &MOL_WRITE_SIZE,
          const std::string &OUTPUT_FILENAME
        );

        //! @brief clone function, should NEVER be used
        ThreadManager *Clone() const
        {
          BCL_Exit( "ThreadManager cannot be cloned.", -1);
          return NULL;
        }

        //! @brief Get class identifier string
        //! @return class string identifier (i.e. class name)
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

        //! @brief Add a molecule to the ensemble, and write it to the output file.  This is also the point to check
        //! @param MOLECULE Molecule that should be added
        //! @param SCORE the molecule's score
        //! @return True if addition/writing was done
        bool ThreadAddMolecule
        (
          const chemistry::FragmentComplete &MOLECULE,
          const size_t &START_MOL_ID,
          const float &SCORE,
          const size_t &TOTAL_MOLS_BUILT
        );

        //! @brief Fetches the molecule ensemble
        //! @return Reference to FragmentEnsemble of the generated molecules
        const chemistry::FragmentEnsemble &GetMolecules() const;

        //! @brief write molecules if it is time to do so
        //! @param FLUSH force a write of everything in m_Molecules to the output file
        void WriteMols( const bool &FLUSH = false) const;

      protected:

        //! @brief Reads from an input stream - reads nothing
        //! @return The same input stream
        std::istream &Read( std::istream &ISTREAM)
        {
          return ISTREAM;
        }

        //! @brief Writes to an output stream - writes nothing
        //! @return The same output stream
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
        {
          return OSTREAM;
        }

      }; // class ThreadManager

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_InitialMoleculesFlag;
      util::ShPtr< command::FlagInterface> m_MutateFunctionFlag;
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;
      util::ShPtr< command::FlagInterface> m_NumberMoleculesFlag;
      util::ShPtr< command::FlagInterface> m_ScoreCutoffFlag;
      util::ShPtr< command::FlagInterface> m_ScoringFunctionFlag;
      util::ShPtr< command::FlagInterface> m_FunctionalizationPointsFlag;
      util::ShPtr< command::FlagInterface> m_FunctionalGroupsFileFlag;
      util::ShPtr< command::FlagInterface> m_OutputInitialFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief Default constructor
      FocusedLibraryDesignOld();

    //////////////////////
    // helper functions //
    //////////////////////

      void EnumerateSites
      (
        const chemistry::AtomVector< chemistry::AtomComplete> &BASE,
        const storage::Vector< size_t> &SITES,
        const util::SiPtrVector< const chemistry::AtomVector< chemistry::AtomComplete> > &FRAGS,
        storage::Vector< chemistry::AtomVector< chemistry::AtomComplete> > &ACCUM
      ) const;

      //! @brief adds a substituent to a molecule at a particular site
      void AddSubToMolecule
      (
        chemistry::AtomVector< chemistry::AtomComplete> &PARENT_ATOMS,
        const chemistry::AtomVector< chemistry::AtomComplete> &SUB_ATOMS,
        const int &SITE
      ) const;

    public:

      //! @brief Clone function
      //! @return pointer to new FocusedLibraryDesignOld
      FocusedLibraryDesignOld *Clone() const
      {
        return new FocusedLibraryDesignOld( *this);
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

      //! @brief initializes the command object for this application
      //! @return a ShPtr to a Command containing all of this applications parameters
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return 0 for success
      int Main() const;

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

      // Static instance of FocusedLibraryDesignOld
      static const ApplicationType FocusedLibraryDesignOld_Instance;

    }; // class FocusedLibraryDesignOld

  } // namespace app
} // namespace bcl

#endif // BCL_APP_FOCUSED_LIBRARY_DESIGN_OLD_H_
