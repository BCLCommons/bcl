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
#include "chemistry/bcl_chemistry_collector_valence.h"
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_grow.h"
#include "chemistry/bcl_chemistry_pick_atom_random.h"
#include "chemistry/bcl_chemistry_pick_fragment_random.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_combine.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_template_instantiations.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_default.h"
#include "model/bcl_model_retrieve_interface.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_function.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf_mdl_handler.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FocusedLibraryDesignRecombinationOld
    //! @brief Application for generating libraries for synthesis using QSAR models and random structure generator
    //!
    //! @author loweew, geanesar, brownbp1
    //! @date 08/31/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FocusedLibraryDesignRecombinationOld :
      public Interface
    {

    private:

      // ThreadManager needs access to private nested classes
      friend class ThreadManager;
      friend class Worker;

      //! @brief generates the 3d conformation from the 2d conformation using internal bcl routines
      //! @param MOLECULE the small molecule for which 3d coordinates will be generated
      //! @param CONFORMATION_SAMPLER SampleConformations instance with which to generate conformations of MOLECULE
      //! @return a conformation of MOLECULE, or empty structure if a good conformation couldn't be found
      static chemistry::FragmentComplete GetBCL3DCoordinates
      (
        const chemistry::ConformationInterface &MOLECULE,
        const chemistry::SampleConformations &CONFORMATION_SAMPLER
      )
      {
        BCL_MessageDbg( "MOLECULE has " + util::Format()( MOLECULE.GetNumberAtoms()) + " atoms");
        BCL_MessageVrb( "Generating conformation...");
        chemistry::FragmentEnsemble conf_ensemble( CONFORMATION_SAMPLER( MOLECULE).First());
        const storage::List< chemistry::FragmentComplete> &ConformationList( conf_ensemble.GetMolecules());
        BCL_MessageDbg( "After fetching conformation list");

        // Return the highest scoring conformer if one exists
        if( ConformationList.GetSize() > 0)
        {
          BCL_MessageVrb( "GetBCL3DCoordinates: Generated 3D structure successfully");
          return ConformationList.FirstElement();
        }

        // Otherwise return the starting structure
        else
        {
          // Bailout
          BCL_MessageStd( "GetBCL3DCoordinates: Could not generate a 3D structure!");
          return MOLECULE;
        }
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class HasValenceElectrons
      //! @brief wrapper class checking for valence electrons
      //! @details This class is a wrapper class around chemistry::CollectorValence that checks if any valence
      //! electrons are found
      //!
      //! @author loweew
      //! @date 08/31/2011
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class HasValenceElectrons :
        public util::FunctionInterfaceSerializable< chemistry::FragmentComplete, bool>
      {
      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief copy constructor
        HasValenceElectrons *Clone() const
        {
          return new HasValenceElectrons( *this);
        };

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! @brief the class name as const ref std::string
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

        //! @brief returns the class name of the object behind a pointer or the current object
        //! @return the class name
        const std::string &GetAlias() const
        {
          static const std::string s_alias( "HasValenceElectrons");
          return s_alias;
        }

      ///////////////
      // operators //
      ///////////////

        //! @brief returns whether the provided molecule has any valence electrons
        //! @return whether the provided molecule has any valence electrons
        bool operator()( const chemistry::FragmentComplete &MOLECULE) const
        {
          BCL_MessageDbg
          (
            "# valence electrons in the molecule: "
            + util::Format()( MOLECULE.GetNumberValences())
          );
          BCL_MessageDbg
          (
            "SumFormula of the checked molecule:  " + MOLECULE.GetSumFormula()
          );

          // collect atoms with valence electrons and check if there are any
          return MOLECULE.GetNumberValences() == 0;
        }

        //! @brief return parameters for member data that are set up from the labels
        //! @return parameters for member data that are set up from the labels
        io::Serializer GetSerializer() const
        {
          io::Serializer serializer;
          serializer.SetClassDescription( "Returns true if the molecule has any unsaturated valences");
          return serializer;
        }

      }; // class HasValenceElectrons

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class MoleculeTooBig
      //! @brief helper determining whether a molecule is too large
      //!
      //! This class is a helper for determining whether a molecule has grown too large
      //!
      //! @author loweew
      //! @date 08/31/2011
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class MoleculeTooBig :
        public util::FunctionInterfaceSerializable< chemistry::FragmentComplete, bool>
      {
      //////////
      // data //
      //////////

        double m_MaxWeight;

      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief default constructor
        MoleculeTooBig() :
          m_MaxWeight( 750)
        {
        }

        //! @brief copy constructor
        MoleculeTooBig *Clone() const
        {
          return new MoleculeTooBig( *this);
        }

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! @brief the class name as const ref std::string
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

        //! @brief returns the class name of the object behind a pointer or the current object
        //! @return the class name
        const std::string &GetAlias() const
        {
          static const std::string s_alias( "MaxWeight");
          return s_alias;
        }

      ///////////////
      // operators //
      ///////////////

        //! @brief returns whether the provided molecule is too large
        //! @return whether the provided molecule is too large
        bool operator()( const chemistry::FragmentComplete &MOLECULE) const
        {
          // collect atoms with valence electrons and check if there are any
          return descriptor::GetCheminfoProperties().calc_MolWeight->SumOverObject( MOLECULE)( 0) >= m_MaxWeight;
        }

      //////////////////////
      // input and output //
      //////////////////////

        //! @brief return parameters for member data that are set up from the labels
        //! @return parameters for member data that are set up from the labels
        io::Serializer GetSerializer() const
        {
          io::Serializer serializer;
          serializer.SetClassDescription
          (
            "Maximum weight allowed for the molecule"
          );
          serializer.AddInitializer
          (
            "",
            "",
            io::Serialization::GetAgent( &m_MaxWeight)
          );
          return serializer;
        }
      }; // class MoleculeTooBig

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
        public math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double>
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
        double operator()( const chemistry::FragmentComplete &FRAGMENT) const
        {
          chemistry::FragmentComplete saturated_frag( FRAGMENT);
          saturated_frag.SaturateWithH();
          linal::Vector< double> res( m_Descriptor->SumOverObject( saturated_frag));
          double activity( res.Sum() / res.GetSize());
          return activity;
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

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ObjectiveFunctionPredict
      //! @brief helper class for the scoring object which uses a model::Interface to predict the generated mol activity
      //!
      //! @author loweew, geanesar
      //! @date Aug 31, 2011
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class ObjectiveFunctionPredict :
        public math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double>
      {

      private:

      //////////
      // data //
      //////////

      public:

          //! model
          model::RetrieveInterface::t_Container m_Models;

          //! code
          util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> > m_Codes;

          //! Number of conformations to generate
          size_t m_ConformationNumber;

          //! Bin size to use when generating conformations
          size_t m_BinSize;

          //! Conformation comparer to use
          std::string m_ConformationComparer;

          // Conformation sampler to use
          util::ShPtr< chemistry::SampleConformations> m_SampleConformations;

          // Force the use of a single Corina 3D conformer
          bool m_Force3DCorina;

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

          //! @brief default constructor
          ObjectiveFunctionPredict();

          //! @brief constructor
          ObjectiveFunctionPredict
          (
            model::RetrieveInterface::t_Container MODELS,
            util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> > CODES,
            const size_t CONF_NUM,
            const size_t BIN_SIZE,
            const std::string &CONF_COMPARER
          )
          : m_Models( MODELS),
            m_Codes( CODES),
            m_ConformationNumber( CONF_NUM),
            m_BinSize( BIN_SIZE),
            m_ConformationComparer( CONF_COMPARER),
            m_SampleConformations( new chemistry::SampleConformations()),
            m_Force3DCorina( false)
          {
          }

          //! @brief constructor
          ObjectiveFunctionPredict
          (
            model::RetrieveInterface::t_Container MODELS,
            util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> > CODES,
            const size_t CONF_NUM,
            const size_t BIN_SIZE,
            const std::string &CONF_COMPARER,
            const bool CORINA_BOOL
          )
          : m_Models( MODELS),
            m_Codes( CODES),
            m_ConformationNumber( CONF_NUM),
            m_BinSize( BIN_SIZE),
            m_ConformationComparer( CONF_COMPARER),
            m_SampleConformations( new chemistry::SampleConformations()),
            m_Force3DCorina( CORINA_BOOL)
          {
          }

          //! @brief constructor
          ObjectiveFunctionPredict
          (
            model::RetrieveInterface::t_Container MODELS,
            util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> > CODES,
            util::ShPtr< chemistry::SampleConformations> SAMPLE_CONFORMATIONS,
            const bool CORINA_BOOL
          )
          : m_Models( MODELS),
            m_Codes( CODES),
            m_SampleConformations( SAMPLE_CONFORMATIONS),
            m_Force3DCorina( CORINA_BOOL)
          {
          }

          //! @brief constructor
          ObjectiveFunctionPredict
          (
            model::RetrieveInterface::t_Container MODELS,
            util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> > CODES,
            util::ShPtr< chemistry::SampleConformations> SAMPLE_CONFORMATIONS
          )
          : m_Models( MODELS),
            m_Codes( CODES),
            m_SampleConformations( SAMPLE_CONFORMATIONS),
            m_Force3DCorina( false)
          {
          }

          //! @brief Clone function
          //! @return pointer to new ObjectiveFunctionPredict
          ObjectiveFunctionPredict *Clone() const
          {
            return new ObjectiveFunctionPredict( *this);
          }

      /////////////////
      // data access //
      /////////////////

          //! @brief returns class name
          //! @return the class name as const ref std::string
          const std::string &GetClassIdentifier() const
          {
            return GetStaticClassName( *this);
          }

      ////////////////
      // operations //
      ////////////////

      ///////////////
      // operators //
      ///////////////

          double operator()( const chemistry::FragmentComplete &MOLECULE) const
          {
            BCL_MessageDbg( "scoring...");

            // making sure any free valence electrons are quenched
            //chemistry::FragmentComplete mol( GetCorina3DCoordinates( MOLECULE));

            chemistry::FragmentComplete mol_wo_h( MOLECULE);
            mol_wo_h.RemoveH();

            // Quench free valence electrons with hydrogens
            chemistry::FragmentComplete mol;
//            if( m_Force3DCorina)
//            {
            util::ShPtr< chemistry::FragmentComplete> corina_mol = chemistry::FragmentEvolveBase::GetCorina3DCoordinates( mol_wo_h);
            mol = *corina_mol;
//            }
//            else
//            {
//
//              chemistry::FragmentComplete bcl_mol
//              (
//                GetBCL3DCoordinates
//                (
//                  mol_wo_h,
//                  *m_SampleConformations
//                )
//              );
//              mol = bcl_mol;
//            }

            // Re-saturate
            mol.SaturateWithH();

            if( mol.GetNumberAtoms() < 1)
            {
              BCL_MessageStd( "ObjectiveFunctionPredict: No 3D coordinates; rejecting.");
              return math::GetHighestUnboundedValue< float>();
            }

            double avg_activity( 0);
            size_t model_count( 0);

            util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> >::const_iterator itr_code( m_Codes.Begin());
            // For each model, calculate and predict activity:
            for
            (
                model::RetrieveInterface::t_Container::const_iterator itr_model( m_Models.Begin()), itr_model_end( m_Models.End());
                itr_model != itr_model_end; ++itr_model, ++itr_code
            )
            {
              model::RetrieveInterface::t_ModelPtr model( *itr_model);

              BCL_Assert( itr_code->IsDefined(), "Code is undefined");
              const descriptor::Combine< chemistry::AtomConformationalInterface, float> &code_ref( **itr_code);

              // calculate activity for that molecule
              const linal::Vector< float> result( code_ref.SumOverObject( mol));

              if( result.GetSize() != code_ref.GetSizeOfFeatures())
              {
                BCL_MessageStd( "could not calculate descriptors for generated molecule - skipping model!");
                continue;
              }

              linal::Matrix< float> feature_mat( 1, code_ref.GetSizeOfFeatures());

              std::copy( result.Begin(), result.End(), feature_mat.Begin());

              // predict activity for that molecule
              double activity( model->operator()( model::FeatureDataSet< float>( feature_mat)).GetMatrix()( 0, 0));
              avg_activity += activity;
              ++model_count;
              BCL_MessageVrb( "predicted activity: " + util::Format()( activity));
            }
            BCL_Assert( model_count > 0, "Number of models is zero!");

            avg_activity /= model_count;

            BCL_MessageStd( "Average activity: " + util::Format()( avg_activity));
            return avg_activity;
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
            // end
            return ISTREAM;
          }

          //! @brief write to std::ostream
          //! @param OSTREAM output stream to write to
          //! @param INDENT number of indentations
          //! @return output stream which was written to
          std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
          {
            // end
            return OSTREAM;
          }

      //////////////////////
      // helper functions //
      //////////////////////

      private:

      }; // class ObjectiveFunctionPredict

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ThreadManager
      //! @brief manages threads for multithreaded structure generation
      //!
      //! @author geanesar
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
//          const size_t                                            m_NumberOfIterations; // Number of iterations in the MC approximator
          const size_t                                            m_Threads; // Number of threads
          chemistry::FragmentEnsemble                             m_Molecules; // The molecules which have been built
          util::ShPtr< chemistry::SampleConformations>            m_SampleConformations; // Sample conformation search class
          io::OFStream                                            m_OutputStream; // Output file to write molecules to

          // Lock for updating Workers
          sched::Mutex m_WorkerStateMutex; // Mutex for determining update status
          sched::Mutex m_WorkerMoleculeMutex; // Mutex for adding molecules

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //!
          //! @class Worker
          //! @brief runs the threads for Worker - builds molecules using metropolis monte-carlo routines
          //!
          //! @author geanesar
          //! @date Nov 7, 2013
          //!
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          struct Worker
          {
            // Rotamer library to use - read in at Main()
            util::ShPtr< chemistry::FragmentComplete>                                                m_BaseFragment; // Base fragment to use
            util::ShPtr< chemistry::FragmentEnsemble>                                                m_FragmentPool; // Fragments to expand with
            model::RetrieveInterface::t_Container                                                    m_Models; // Scoring models
            util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> >   m_Codes; // Code (also for scoring)
            descriptor::CheminfoProperty                                                             m_PropertyScorer; // Set objective function with property instead of model
            util::ShPtr< chemistry::SampleConformations>                                             m_SampleConformations; // Conformation sampler
            bool                                                                                     m_Force3DCorina; // Force use of a single Corina 3D conformer
            size_t                                                                                   m_ConformationNumber; // Number of conformations to generate
            size_t                                                                                   m_BinSize; // Bin size
            std::string                                                                              m_ConformationComparer; // Conformation comparer
            chemistry::FragmentEnsemble                                                              m_ThreadMolecules; // Number of molecules generated
            util::ShPtr< math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double> > m_Score; // Objective function
            util::ShPtr< math::MutateInterface< chemistry::FragmentComplete> >                       m_Mutate; // Grow molecules from scaffold
            util::SiPtr< ThreadManager>                                                              m_ThreadManager; // Pointer to the thread manager, needed so Worker can be updated

            // Builds and score the molecule
            void RunThread()
            {
              do
              {
                // Set up MC method

            //////////////////
            // approximator //
            //////////////////

                // create the temperature control
                util::ShPtr< mc::TemperatureInterface> sp_temperature( new mc::TemperatureDefault( 75.0));

                // create the metropolis
                mc::Metropolis< double> metropolis( sp_temperature, true);

                // create the termination criterion
                opti::CriterionCombine< chemistry::FragmentComplete, double> criterion_combine;

                // insert termination criteria that depends on the number of valence electrons
                criterion_combine.InsertCriteria
                (
                  opti::CriterionFunction< chemistry::FragmentComplete, double>( FocusedLibraryDesignRecombinationOld::HasValenceElectrons())
                );

                BCL_MessageDbg( " after create the termination criteria dependent on number of valence electrons");

                // insert termination criteria that depends on the size of the molecule
                criterion_combine.InsertCriteria
                (
                  opti::CriterionFunction< chemistry::FragmentComplete, double>( FocusedLibraryDesignRecombinationOld::MoleculeTooBig())
                );

                BCL_MessageDbg( " after create the termination criteria dependent on size of the molecule");

                // insert termination criteria that depends on the total number of MC iterations
                opti::CriterionNumberIterations< chemistry::FragmentComplete, double> maximum_number_iterations( 100);
                criterion_combine.InsertCriteria( maximum_number_iterations);

                BCL_MessageDbg( " after create the termination criteria dependent on total number of MC iterations");

                mc::Approximator< chemistry::FragmentComplete, double> approximator
                (
                  *m_Score,
                  *m_Mutate,
                  metropolis,
                  criterion_combine,
                  *m_BaseFragment
                );

                BCL_MessageStd( " before approx");

            //////////////////
            // approximator //
            //////////////////

                // run the approximator
                approximator.Approximate();
                BCL_MessageStd( " after run the approximator");

                // creating shared pointer to small molecule of generated molecule
                const chemistry::FragmentComplete gen_mol_orig( approximator.GetTracker().GetBest()->First());
                double activity( approximator.GetTracker().GetBest()->Second());

                // Need to saturate with H for checking - copy it
                chemistry::FragmentComplete gen_mol_2d( gen_mol_orig);

                // Generate final 3D coordinates
//                chemistry::FragmentComplete gen_mol_3d;
//                if( m_Force3DCorina)
//                {
//                util::ShPtr< chemistry::FragmentComplete> corina_mol = chemistry::FragmentEvolveBase::GetCorina3DCoordinates(gen_mol_orig);
//                gen_mol_3d = *corina_mol;
//                }
//                else
//                {
//
//                  chemistry::FragmentComplete bcl_mol
//                  (
//                    GetBCL3DCoordinates
//                    (
//                      gen_mol_orig,
//                      *m_SampleConformations
//                    )
//                  );
//                  gen_mol_3d = bcl_mol;
//                }

//                gen_mol_3d.SaturateWithH();
                gen_mol_2d.SaturateWithH();

                // Make sure the generated molecule is constitutionally the same as the original
                chemistry::ConstitutionGraphConverter graph_maker;
                const graph::ConstGraph< size_t, size_t> graph_2D( graph_maker( chemistry::FragmentConstitutionShared( gen_mol_2d)));
//                const graph::ConstGraph< size_t, size_t> graph_3D( graph_maker( chemistry::FragmentConstitutionShared( gen_mol_3d)));

                graph::SubgraphIsomorphism< size_t, size_t> isomorphism;

                isomorphism.SetSubgraphExternalOwnership( graph_2D);
//                isomorphism.SetGraphExternalOwnership( graph_3D);

                // store score as MDL property on new molecule
                gen_mol_2d.StoreProperty( "FLDR_Score", util::Format()( activity));

//                if( gen_mol_2d.GetNumberAtoms() > 0 && isomorphism.FindIsomorphism())
                if( gen_mol_2d.GetNumberAtoms() > 0)
                {
                  m_ThreadManager->m_WorkerMoleculeMutex.Lock();
                  if( m_ThreadManager->GetNumberMoleculesBuilt() + 1 <= m_ThreadManager->GetNumberMoleculesToBuild())
                  {
                    m_ThreadManager->AddMolecule( gen_mol_2d);
                    m_ThreadManager->IncreaseMoleculeBuiltCount();
                  }

                  m_ThreadManager->m_WorkerMoleculeMutex.Unlock();
                  BCL_MessageStd( "Added new molecule");
                }
                else
                {
                  if( !isomorphism.FindIsomorphism())
                  {
                    BCL_MessageStd( "Constitution of generated 3D conformation differs from the original");
                  }
                  else
                  {
                    BCL_MessageStd( "No 3D coordinates for this molecule");
                  }
                }
              } while( m_ThreadManager->UpdateWorker( *this));
            } // RunThread()

          }; // stuct Worker

          // Tests to see if the worker should keep running
          bool UpdateWorker( ThreadManager::Worker &WORKER)
          {
            // Lock structure during the modification
            m_WorkerStateMutex.Lock();

            if( m_NumberOfMoleculesBuilt >= m_NumberOfMoleculesRequested)
            {
              m_WorkerStateMutex.Unlock();
              return false;
            }

            m_WorkerStateMutex.Unlock();
            return true;
          } // UpdateWorker()

      public:

          //! param
          ThreadManager(
            util::ShPtr< chemistry::FragmentComplete>                   BASE_FRAGMENT, // Base fragment to use
            util::ShPtr< chemistry::FragmentEnsemble>                   FRAGMENT_POOL, // Fragments to add to base fragment
            model::RetrieveInterface::t_Container                       MODELS, // Models
            util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> > CODES, // Input code
            descriptor::CheminfoProperty                           PROPERTY_SCORER, // alternative scorer
            const size_t                                           &CONFORMATION_NUMBER, // Number of conformations
            const size_t                                           &BIN_SIZE, // Bin size
            const std::string                                      &CONFORMATION_COMPARER, // How to compare conformations
            const size_t                                           &NUMBER_OF_MOLECULES, // Number to build
//            const size_t                                           &NUMBER_OF_ITERATIONS, // Number of MC iterations
            const size_t                                           &NUMBER_THREADS, // Number of threads (from scheduler)
            const std::string                                      &OUTPUT_FILENAME
          ) :
            m_NumberOfMoleculesRequested( NUMBER_OF_MOLECULES),
            m_NumberOfMoleculesBuilt( 0),
//            m_NumberOfIterations( NUMBER_OF_ITERATIONS),
            m_Threads( std::min( NUMBER_THREADS, NUMBER_OF_MOLECULES)),
            m_SampleConformations
            (
              new chemistry::SampleConformations()
            )
          {

            io::File::MustOpenOFStream( m_OutputStream, OUTPUT_FILENAME);

            // Set up workers
            std::vector< Worker> workers( m_Threads);
            for(
                std::vector< Worker>::iterator itr( workers.begin()), end( workers.end());
                itr != end;
                ++itr
            )
            {
              Worker &worker_ref( *itr);
              worker_ref.m_BaseFragment         = BASE_FRAGMENT;
              worker_ref.m_FragmentPool         = FRAGMENT_POOL;
              worker_ref.m_Models               = MODELS;
              worker_ref.m_Codes                = CODES;
              worker_ref.m_ConformationNumber   = CONFORMATION_NUMBER;
              worker_ref.m_BinSize              = BIN_SIZE;
              worker_ref.m_ConformationComparer = CONFORMATION_COMPARER;
              worker_ref.m_SampleConformations  = m_SampleConformations;
              worker_ref.m_PropertyScorer       = PROPERTY_SCORER;
              worker_ref.m_ThreadManager        = this;

              BCL_Assert( worker_ref.m_PropertyScorer.IsDefined(), "PropertyScorer is undefined");
              if( !worker_ref.m_PropertyScorer.IsDefined())
              {
                worker_ref.m_Score = util::ShPtr< math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double> >
                (
                  new FocusedLibraryDesignRecombinationOld::ObjectiveFunctionPredict
                  (
                    MODELS,
                    CODES,
                    m_SampleConformations
                  )
                );
              }
              else
              {
                worker_ref.m_Score = util::ShPtr< math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double> >
                (
                  new FocusedLibraryDesignRecombinationOld::FitnessFunctionRaw
                  (
                    PROPERTY_SCORER
                  )
                );
              }
              worker_ref.m_Mutate = util::ShPtr< math::MutateInterface< chemistry::FragmentComplete> >
              (
                new chemistry::FragmentGrow
                (
                  FRAGMENT_POOL,
                  chemistry::CollectorValence(),
                  chemistry::PickAtomRandom(),
                  chemistry::PickFragmentRandom()
                )
              );
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

          int GetNumberMoleculesBuilt()
          {
            return m_NumberOfMoleculesBuilt;
          }

          int GetNumberMoleculesToBuild()
          {
            return m_NumberOfMoleculesRequested;
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

      //input sdf files and output matrix file

      //! flag to control number of molecules to be generated
      util::ShPtr< command::FlagInterface> m_NumberMoleculesFlag;

      //! flag to control the number of MC iterations in molecule optimization
      util::ShPtr< command::FlagInterface> m_NumberIterationsFlag;

      //! flag to control
      util::ShPtr< command::FlagInterface> m_BaseFragmentFlag;

      //! flag for defining output filename,
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! flag for defining input fragments
      util::ShPtr< command::FlagInterface> m_GrowFragmentsFlag;

      //! flag for model id
      util::ShPtr< command::FlagInterface> m_ModelIDFlag;

      //! flag for storing the final trained model
      util::ShPtr< command::FlagInterface> m_FlagModelInterfaceStorage;

      //! flag for an alternative score function to just the trained model
      util::ShPtr< command::FlagInterface> m_PropertyScoringFunctionFlag;

      //! flag for setting the conformation comparer
      util::ShPtr< command::FlagInterface> m_ConformationComparerFlag;

      //! flag for the number of conformations to generate
      util::ShPtr< command::FlagInterface> m_ConformationNumberFlag;

      //! flag for the number of conformations to generate
      util::ShPtr< command::FlagInterface> m_Force3DCorina;

      //! flag for specifying bin size
      util::ShPtr< command::FlagInterface> m_BinSizeFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      FocusedLibraryDesignRecombinationOld();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      FocusedLibraryDesignRecombinationOld *Clone() const
      {
        return new FocusedLibraryDesignRecombinationOld( *this);
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

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_NumberMoleculesFlag);
        sp_cmd->AddFlag( m_NumberIterationsFlag);
        sp_cmd->AddFlag( m_BaseFragmentFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_GrowFragmentsFlag);
        sp_cmd->AddFlag( m_ModelIDFlag);

        // flag for model storage
        sp_cmd->AddFlag( m_FlagModelInterfaceStorage);

        // flag for alternative scoring function
        sp_cmd->AddFlag( m_PropertyScoringFunctionFlag);

        // flags for controlling the 3D conformation builder
        sp_cmd->AddFlag( m_ConformationNumberFlag);
        sp_cmd->AddFlag( m_ConformationComparerFlag);
        sp_cmd->AddFlag( m_BinSizeFlag);
        sp_cmd->AddFlag( m_Force3DCorina);
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
        io::File::MustOpenIFStream( input, m_BaseFragmentFlag->GetFirstParameter()->GetValue());

        // Needs to be wrapped in a ShPtr so it can be passed to ThreadManager
        util::ShPtr< chemistry::FragmentComplete> sp_basefragment
        (
          new chemistry::FragmentComplete( sdf::FragmentFactory::MakeFragment( input, sdf::e_Maintain))
        );
        io::File::CloseClearFStream( input);

        // try to read cheminfo property scorer
        descriptor::CheminfoProperty property_scorer;
        if( m_PropertyScoringFunctionFlag->GetFlag())
        {
          property_scorer = m_PropertyScoringFunctionFlag->GetFirstParameter()->GetValue();
        }

        // storage for the final model::Interface and descriptor
        linal::Matrix< float> feature_mat;

        // store models from storage
        storage::Vector< std::string> model_ids;
        model::RetrieveInterface::t_Container models;
        storage::List< util::ObjectDataLabel> descriptor_model_list;
        util::ShPtrList< descriptor::Combine< chemistry::AtomConformationalInterface, float> > code_inputs;

        if( m_FlagModelInterfaceStorage->GetFlag())
        {

          // storage for the final model::Interface
          util::Implementation< model::RetrieveInterface> model_storage
          (
            m_FlagModelInterfaceStorage->GetFirstParameter()->GetValue()
          );

          //        BCL_Assert
          //        (
          //          model_storage.IsDefined(),
          //          "Model storage was undefined. see help for information about how to initialize model storage."
          //        );
          //        BCL_MessageDbg( "Model storage initialized ... ");

          // Read specified models, if set, or all models if not
          if( m_ModelIDFlag->GetFlag())
          {
            for
            (
                util::ShPtrVector< command::ParameterInterface>::const_iterator
                itr_model( m_ModelIDFlag->GetParameterList().Begin()),
                itr_model_end( m_ModelIDFlag->GetParameterList().End());
                itr_model != itr_model_end;
                ++itr_model
            )
            {
              const std::string id_model( ( *itr_model)->GetValue());

              // store model id
              model_ids.PushBack( id_model);

              // store models from storage
              models.PushBack( model_storage->Retrieve( id_model));

              // store descriptor set from storage
              descriptor_model_list.PushBack( model_storage->RetrieveDescriptorSet( id_model));

              code_inputs.PushBack
              (
                util::ShPtr< descriptor::Combine< chemistry::AtomConformationalInterface, float> >( new descriptor::Combine< chemistry::AtomConformationalInterface, float>)
              );

              BCL_MessageStd( "Model id: " + id_model + " retrieved.");
            }
          }
          else
          {
            // store models from storage
            model_ids = model_storage->GetAllKeys();
            models = model_storage->RetrieveEnsemble();
            descriptor_model_list = model_storage->RetrieveEnsembleDescriptors();

            for
            (
                storage::List< util::ObjectDataLabel>::const_iterator itr_descriptors( descriptor_model_list.Begin()),
                itr_descriptors_end( descriptor_model_list.End());
                itr_descriptors != itr_descriptors_end;
                ++itr_descriptors
            )
            {
              code_inputs.PushBack
              (
                util::ShPtr< descriptor::Combine< chemistry::AtomConformationalInterface, float> >( new descriptor::Combine< chemistry::AtomConformationalInterface, float>)
              );
            }

            BCL_MessageStd( "Read models from database");
            BCL_MessageStd( "Read " + util::Format()( models.GetSize()) + " models");
          }

          // read in the descriptors, set codes accordingly
          storage::List< util::ShPtr< descriptor::Combine< chemistry::AtomConformationalInterface, float> > >::iterator itr_code = code_inputs.Begin();
          for
          (
              storage::List< util::ObjectDataLabel>::const_iterator itr_descriptors( descriptor_model_list.Begin()),
              itr_descriptors_end( descriptor_model_list.End());
              itr_descriptors != itr_descriptors_end;
              ++itr_descriptors, ++itr_code
          )
          {
            code_inputs.PushBack
            (
              util::ShPtr< descriptor::Combine< chemistry::AtomConformationalInterface, float> >( new descriptor::Combine< chemistry::AtomConformationalInterface, float>)
            );
            ( *itr_code)->AssertRead( *itr_descriptors);
          }
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
        ThreadManager thread_manager
        (
          sp_basefragment,
          sp_fragment_pool,
          models,
          code_inputs,
          property_scorer,
          m_ConformationNumberFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_BinSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_ConformationComparerFlag->GetFirstParameter()->GetValue(),
          m_NumberMoleculesFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
//          m_NumberIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          sched::GetNumberCPUs(),
          m_OutputFilenameFlag->GetFirstParameter()->GetValue()
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

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType FocusedLibraryDesignRecombinationOld_Instance;

    }; // FocusedLibraryDesignRecombinationOld

    //! @brief standard constructor
    FocusedLibraryDesignRecombinationOld::FocusedLibraryDesignRecombinationOld() :
      m_NumberMoleculesFlag
      (
        new command::FlagStatic
        (
          "number_molecules", "flag for number of molecules to generate",
          command::Parameter
          (
            "number_molecules", "total number of molecules",
            command::ParameterCheckRanged< int>( 0, std::numeric_limits< int>::max()), "10"
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
            command::ParameterCheckRanged< int>( 0, std::numeric_limits< int>::max()), "10"
          )
        )
      ),
      m_BaseFragmentFlag
      (
        new command::FlagStatic
        (
          "base_fragment", "filename for input starting fragment",
          command::Parameter
          (
            "fragment_filename", "filename for input sdf of molecules", ""
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
        new command::FlagDynamic
        (
          "grow_fragments",
          "files containing fragments to append to the molecule",
          command::Parameter
          (
            "grow fragments filename",
            "name of file containing grow fragments",
            command::ParameterCheckFileExistence()
          ),
          1
        )
      ),
      m_ModelIDFlag
      (
        new command::FlagDynamic
        (
          "model_ids",
          "flag for model ids",
          command::Parameter
          (
            "model_ids",
            "model ids"
          ),
          0,
          std::numeric_limits< size_t>::max()
        )
      ),
      m_FlagModelInterfaceStorage
      (
        new command::FlagDynamic
        (
          "storage_model",
          "choice of model storage",
          command::Parameter
          (
            "model_storage",
            "choice of model storage",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveInterface>())
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
              FitnessFunctionRaw()
            )
          )
        )
      ),
      m_ConformationComparerFlag
      (
        new command::FlagStatic
        (
          "conformation_comparer",
          "method to compare conformers with",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< chemistry::ConformationComparisonInterface>())
          )
        )
      ),
      m_ConformationNumberFlag
      (
        new command::FlagStatic
        (
          "conformation_number",
          "the number of conformations to build",
          command::Parameter
          (
            "conf_num",
            "the number of conformations",
            command::ParameterCheckRanged< unsigned int>( 0, std::numeric_limits< unsigned int>::max()),
            "1"
          )
        )
      ),
      m_Force3DCorina
      (
        new command::FlagStatic
        (
          "corina_conf",
          "generate a single conformer with Corina instead of the BCL",
          command::Parameter
          (
            "corina_conf",
            "generate a single conformer with Corina instead of the BCL",
            command::ParameterCheckRanged< bool>( false, true),
            "false"
          )
        )
      ),
      m_BinSizeFlag
      (
        new command::FlagStatic
        (
          "bin_size",
          "size of the bin to use",
          command::Parameter(
            "bin_size",
            "size of the bin",
            command::ParameterCheckRanged< unsigned int>( 0, std::numeric_limits< unsigned int>::max()),
            "30"
          )
        )
      )
    {
    }

    const ApplicationType FocusedLibraryDesignRecombinationOld::FocusedLibraryDesignRecombinationOld_Instance
    (
      GetAppGroups().AddAppToGroup( new FocusedLibraryDesignRecombinationOld(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
