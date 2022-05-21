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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_quality_batch.h"
#include "assemble/bcl_assemble_sse_compare_type.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "coord/bcl_coord_move_rotate_defined.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "find/bcl_find_template_instantiations.h"
#include "fold/bcl_fold_default_mutates.h"
#include "fold/bcl_fold_mutate_protein_model_sse.h"
#include "fold/bcl_fold_mutate_protein_model_sse_remove.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap.h"
#include "fold/bcl_fold_mutate_sse_bend_random.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_repeat.h"
#include "math/bcl_math_template_instantiations.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "storage/bcl_storage_template_instantiations.h"
namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinPerturb
    //! @brief This class is the application for folding amino acid sequences into tertiary structures.
    //!
    //! @author woetzen
    //! @date Mar 24, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinPerturb :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! input pdb file flag
      util::ShPtr< command::FlagInterface> m_PdbFilenameFlag;

      //! idealize flag
      util::ShPtr< command::FlagInterface> m_IdealizeFlag;

      //! starting "mutation" flag
      util::ShPtr< command::FlagInterface> m_StartMutationFilenameParam;

      //! filename where mutate object
      util::ShPtr< command::FlagInterface> m_PerturbationFilenameParam;

      //! number of repeats
      util::ShPtr< command::FlagInterface> m_NumberRepeatsFlag;

      //! prefix to be used for writing files
      util::ShPtr< command::FlagInterface> m_FlagFilePrefix;

      //! aa sequences to be folded
      mutable util::ShPtr< assemble::ProteinModel> m_StartModel;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinPerturb();

    public:

      //! @brief Clone function
      //! @return pointer to new ProteinPerturb
      ProteinPerturb *Clone() const
      {
        return new ProteinPerturb( *this);
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

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>( size_t( 1), "PerturbProtein");
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief create default perturbation function for protein models
      util::ShPtr< math::MutatePerturbation< assemble::ProteinModel> > GetDefaultPerturbation() const
      {
        fold::MutateTree mutate_tree;
        fold::DefaultMutates::GetInstance().InitializeMutates();
        fold::DefaultMutates::GetInstance().InitializeMutateTree( mutate_tree);
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> > mutates
        (
          mutate_tree.ConstructMutates()
        );

        // create the combined mutate
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> > combined_mutate
        (
          new math::MutateRepeat< assemble::ProteinModel>
          (
            mutates,
            1,
            GetStartModel().GetNumberSSEs()
          )
        );

        // create the perturbation
        const size_t nr_global_repeats( 3);
        util::ShPtr< math::MutatePerturbation< assemble::ProteinModel> > default_perturbation
        (
          new math::MutatePerturbation< assemble::ProteinModel>()
        );
        default_perturbation->InsertPerturbation( combined_mutate, nr_global_repeats);

        // end
        return default_perturbation;
      }

      //! @brief create perturbation function for protein models
      util::ShPtr< math::MutatePerturbation< assemble::ProteinModel> > GetPerturbation() const
      {
        util::ShPtr< math::MutatePerturbation< assemble::ProteinModel> > mutate;

        if
        (
             m_PerturbationFilenameParam->GetFlag()
          && m_PerturbationFilenameParam->GetFirstParameter()->GetWasSetInCommandLine()
        )
        {
          io::IFStream read;

          if( io::File::TryOpenIFStream( read, m_PerturbationFilenameParam->GetFirstParameter()->GetValue()))
          {
            read >> mutate;
            io::File::CloseClearFStream( read);
          }
          else
          {
            io::OFStream write;
            io::File::MustOpenOFStream( write, m_PerturbationFilenameParam->GetFirstParameter()->GetValue());
            mutate = GetDefaultPerturbation();
            write << mutate;
            io::File::CloseClearFStream( write);
          }
        }
        else
        {
          mutate = GetDefaultPerturbation();
        }
        return mutate;
      }

      //! @brief return the aasequences that are supposed to be folded - initializes them once from given fasta
      const assemble::ProteinModel &GetStartModel() const
      {
        if( !m_StartModel.IsDefined())
        {
          // open pdb file
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_PdbFilenameFlag->GetFirstParameter()->GetValue());

          pdb::Handler handler( read);
          m_StartModel = util::ShPtr< assemble::ProteinModel>( pdb::Factory().ProteinModelFromPDB( handler).Clone());

          // close stream
          io::File::CloseClearFStream( read);

          // idealize
          if( m_IdealizeFlag->GetFlag())
          {
            m_StartModel->SetToIdealConformation();
          }

          // is start mutation was given, apply
          if( m_StartMutationFilenameParam->GetFlag())
          {
            util::ShPtr< math::MutateInterface< assemble::ProteinModel> > start_mutate;

            // open mutate file
            io::IFStream read;
            io::File::MustOpenIFStream( read, m_StartMutationFilenameParam->GetFirstParameter()->GetValue());

            read >> start_mutate;
            m_StartModel = start_mutate->operator()( *m_StartModel).GetArgument();

            // close stream
            io::File::CloseClearFStream( read);
          }
        }

        // return
        return *m_StartModel;
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        //input pdb file parameter
        sp_cmd->AddFlag( m_PdbFilenameFlag);

        // idealize flag
        sp_cmd->AddFlag( m_IdealizeFlag);

        // start mutation
        sp_cmd->AddFlag( m_StartMutationFilenameParam);

        // flag for mutate objects filename
        sp_cmd->AddFlag( m_PerturbationFilenameParam);

        //input m_QualityList
        sp_cmd->AddFlag( quality::Measures::GetFlagQualityMeasures());

        // flag for number of repeats
        sp_cmd->AddFlag( m_NumberRepeatsFlag);

        // add file prefix flag
        sp_cmd->AddFlag( m_FlagFilePrefix);

        // add pdb flags
        sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());
        sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
        sp_cmd->AddFlag( pdb::Factory::GetFlagConvertToNaturalAAType());
        sp_cmd->AddFlag( pdb::Factory::GetFlagWritePDBResID());
        sp_cmd->AddFlag( pdb::Factory::GetFlagWritePDBAtomID());
        sp_cmd->AddFlag( pdb::Factory::GetFlagDSSP());
        sp_cmd->AddFlag( pdb::Factory::GetFlagWriteZeroCoordinatesForUndefinedAminoAcids());

        // add default bcl parameters
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const
      {
        // explicit template instantiation for those classes do not work in vs - this guarantees
        // that the s_Instances are initialized so that ShPtr to those classes can be read
        math::MutateCombine< assemble::ProteinModel>::s_Instance.IsDefined();
        math::MutateDecisionNode< assemble::ProteinModel>::s_Instance.IsDefined();
        math::MutateRepeat< assemble::ProteinModel>::s_Instance.IsDefined();

        util::ShPtr< math::MutatePerturbation< assemble::ProteinModel> > perturbation( GetPerturbation());

        // create factory
        pdb::Factory factory;

        util::Format number_format;
        number_format.W( 5).R().Fill( '0');
        size_t iteration_nr( 0);

        perturbation->SetStartArgument( GetStartModel());

        // can calculate the quality mesarues for a given protein model
        const assemble::QualityBatch qualities
        (
          quality::Measures::GetCommandLineQualityMeasures(),
          biol::GetAtomTypes().GetBackBoneAtomTypes()
        );

        // write start model
        {
          // write the current model to file
          const std::string out_filename
          (
            m_FlagFilePrefix->GetFirstParameter()->GetValue()
            + "/perturbed_" + number_format( iteration_nr) + ".pdb"
          );
          io::OFStream write;
          io::File::MustOpenOFStream( write, out_filename);
          factory.WriteModelToPDB( GetStartModel(), write);

          // write measures
          qualities.ConstructTable( GetStartModel()).WriteFormatted( write);

          io::File::CloseClearFStream( write);
        }

        // get the number repeats
        const size_t total_number_repeats( m_NumberRepeatsFlag->GetFirstParameter()->GetNumericalValue< size_t>());

        for( size_t number_repeat( 0); number_repeat != total_number_repeats; ++number_repeat)
        {
          while( !perturbation->IsFinished())
          {
            ++iteration_nr;

            // perturb
            perturbation->Perturb();

            // get current
            const assemble::ProteinModel &current( perturbation->GetCurrent());

            // write the current model to file
            const std::string out_filename
            (
                m_FlagFilePrefix->GetFirstParameter()->GetValue()
              + "/perturbed_" + number_format( iteration_nr) + ".pdb"
            );
            io::OFStream write;
            io::File::MustOpenOFStream( write, out_filename);
            factory.WriteModelToPDB( current, write);

            // write measures
            qualities.ConstructTable( current).WriteFormatted( write);

            // clear stream
            io::File::CloseClearFStream( write);
          }
          perturbation->Reset();
        }
        // end
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

    private:

      static const ApplicationType ProteinPerturb_Instance;

    }; // ProteinPerturb

    // default constructor
    ProteinPerturb::ProteinPerturb() :
      m_PdbFilenameFlag
      (
        new command::FlagStatic
        (
          pdb::GetDefaultFileExtension(),
          "\t\tpdb file of starting structure to be perturbed",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename for template pdb to which rmsd will be calculated",
            command::ParameterCheckExtension( std::string( ".pdb")),
            ""
          )
        )
      ),
      m_IdealizeFlag
      (
        new command::FlagStatic
        (
          "idealize",
          "idealize sses"
        )
      ),
      m_StartMutationFilenameParam
      (
        new command::FlagStatic
        (
          "start_mutate",
          "\t\tfile with starting mutation",
          command::Parameter
          (
            "starting_mutation_filename",
            "\tfilename for starting mutation object",
            ""
          )
        )
      ),
      m_PerturbationFilenameParam
      (
        new command::FlagStatic
        (
          "perturb",
          "\t\tfile with perturbation object - if it does not exist, default will be written to it",
          command::Parameter
          (
            "perturbation_filename",
            "\tfilename for input perturbation object",
            ""
          )
        )
      ),
      m_NumberRepeatsFlag
      (
        new command::FlagStatic
        (
          "number_repeats",
          "\t\tnumber of repeats of the perturbation ( only makes sense with random-based mutates",
          command::Parameter
          (
            "number_of_the_repeats",
            "\tnumber of repeats",
            "1"
          )
        )
      ),
      m_FlagFilePrefix
      (
        new command::FlagStatic
        (
          "prefix",
          "file prefix to be used for writing files",
          command::Parameter
          (
            "file_prefix",
            "\twill be prepended to each file - can contain a path",
            "."
          )
        )
      ),
      m_StartModel()
    {
    }

    const ApplicationType ProteinPerturb::ProteinPerturb_Instance
    (
      GetAppGroups().AddAppToGroup( new ProteinPerturb(), GetAppGroups().e_Protein)
    );

  } // namespace app
} // namespace bcl
