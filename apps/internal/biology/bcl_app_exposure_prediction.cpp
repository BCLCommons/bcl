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

// include the interface for all apps
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_with_cache_storage_file.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "biol/bcl_biol_exposure_prediction.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdlib.h>

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ExposurePrediction
    //! @brief Reads in a sequence and writes out the predicted exposure values.
    //!
    //! @author weinerbe, lib14
    //! @date Oct 3, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ExposurePrediction :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

        //! fasta file
        util::ShPtr< command::FlagInterface> m_FastaFileFlag;

        //! exposure types
        util::ShPtr< command::FlagInterface> m_ExposureTypesFlag;

        //! output prefix
        util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ExposurePrediction();

      //! @brief Clone function
      //! @return pointer to new ExposurePrediction
      ExposurePrediction *Clone() const
      {
        return new ExposurePrediction( *this);
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

      // instantiate enumerator for ExposurePrediction class
      static const ApplicationType ExposurePrediction_Instance;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! Main
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const
      {
        static const std::string s_readme;
        return s_readme;
      }

    ///////////////
    // operators //
    ///////////////

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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class ExposurePrediction

      //! @brief default constructor
      ExposurePrediction::ExposurePrediction() :
          m_FastaFileFlag
          (
            new command::FlagStatic
            (
              "fasta",
              "fasta file of the protein for which to predict residue exposure",
              command::Parameter
              (
                "fasta file",
                "amino acid sequence file in fasta format"
              )
            )
          ),
          m_ExposureTypesFlag
          (
            new command::FlagDynamic
            (
              "exposure_types",
              "specify what types of exposure to predict: protomeric_cn, oligomeric_cn, "
              "protomeric_rsa, oligomeric_rsa",
              command::Parameter
              (
                "exposure types",
                "protomeric_cn, oligomeric_cn, protomeric_rsa, oligomeric_rsa"
              ),
              1, // mininum number of parameters
              4 // maximum number of parameters
            )
          ),
          m_OutputPrefixFlag
          (
            new command::FlagStatic
            (
              "output_prefix",
              "prefix of the output file(s)",
              command::Parameter
              (
                "output prefix",
                "prefix for the output file(s)",
                m_FastaFileFlag->GetFirstParameter()->GetValue().substr
                (
                  0,
                  m_FastaFileFlag->GetFirstParameter()->GetDefaultValue().rfind( '.')
                )
              )
            )
          )
          {
            // nothing else to do
          }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> ExposurePrediction::InitializeCommand() const
      {
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // add member flags
        sp_cmd->AddFlag( m_FastaFileFlag);
        sp_cmd->AddFlag( m_ExposureTypesFlag);
        sp_cmd->AddFlag( m_OutputPrefixFlag);

        // add default bcl flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

    //! Main
    int ExposurePrediction::Main() const
    {
      // read fasta
      const std::string fasta_filename( m_FastaFileFlag->GetFirstParameter()->GetValue());
      assemble::ProteinModelWithCache model
      (
        assemble::ProteinWithCacheStorageFile::GenerateProteinModelFromFile
        (
          fasta_filename,
          false,
          biol::GetAAClasses().e_AA
        ),
        false
      );

      // determine the prefix
      const size_t pos_period( fasta_filename.rfind( "."));
      const std::string prefix( fasta_filename.substr( 0, pos_period));

      // read ascii, use parameter if given, if not replace ".fasta" with ".ascii" in fasta_filename and try that
      if( !biol::BlastProfileHandler::TryReadProfileForProteinModel( model, prefix, ".uniref50.ascii5"))
      {
        // if the script for running psiblast does not even exist, terminate the application
        if( !io::DirectoryEntry( "/dors/meilerlab/apps/PISCES/opm/scripts").DoesExist())
        {
          BCL_MessageTop
          (
            "Required blast profile file could not be located! Run psiblast on the uniref50 database, "
            "5 iterations with ethreshold/evalue set to 0.01, and save the file to " + prefix + ".uniref50.ascii5"
          );
          BCL_ExitWithoutCallstack( "Required blast profile not found", -1);
        }
        else
        {
          std::string command
          (
            "/dors/meilerlab/apps/PISCES/opm/scripts/psiblast.db_itr_eval.local.allcores.sh "
            + fasta_filename + " uniref50nofrag 5 0.01"
          );
          BCL_MessageStd
          (
            "Required blast profile file could not be located. Running psiblast on the uniref50 database, "
            "5 iterations with ethreshold/evalue set to 0.01, and saving the file to " + prefix + ".uniref50.ascii5 . Command:\n"
            + command
          );
          // run the script for running psiblast
          system( command.c_str());
          BCL_Assert
          (
            biol::BlastProfileHandler::TryReadProfileForProteinModel( model, prefix, ".uniref50.ascii5"),
            "Could not generate pssm for fasta file. See " + prefix + ".uniref50.ascii5.log for command and errors"
          );
        }
      }

      // initialize output
      io::OFStream write;

      // parameters for the exposure type flag
      util::ShPtrVector< command::ParameterInterface> par_list( m_ExposureTypesFlag->GetParameterList());

      biol::ExposurePrediction::ExposureType exposure_type;
//      (
//        biol::ExposurePrediction::ExposureType::e_ProtomericContactNumber
//      );

      // iterate over exposure types
      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          par_itr( par_list.Begin()), par_itr_end( par_list.End());
        par_itr != par_itr_end;
        ++par_itr
      )
      {
        // exposure type specified on the command line
        const std::string exposure_type_string( ( *par_itr)->GetValue());

        // get parameter value
        if( exposure_type_string.compare( "protomeric_cn") == 0)
        {
          exposure_type = biol::ExposurePrediction::ExposureType::e_ProtomericContactNumber;
        }
        else if( exposure_type_string.compare( "oligomeric_cn") == 0)
        {
          exposure_type = biol::ExposurePrediction::ExposureType::e_OligomericContactNumber;
        }
        else if( exposure_type_string.compare( "protomeric_rsa") == 0)
        {
          exposure_type = biol::ExposurePrediction::ExposureType::e_ProtomericRSA;
        }
        else if( exposure_type_string.compare( "oligomeric_rsa") == 0)
        {
          exposure_type = biol::ExposurePrediction::ExposureType::e_OligomericRSA;
        }
        else
        {
          BCL_MessageStd
          (
            "Unknown exposure type " + exposure_type_string + " was given. Use one or several of: "
            "protomeric_cn, oligomeric_cn, protomeric_rsa, oligomeric_rsa.\n"
            + "Predicting protomeric_cn exposure as the default."
          )
        }

        // output file name
        const std::string output_name
        (
          m_OutputPrefixFlag->GetFirstParameter()->GetValue() +
            prefix + biol::ExposurePrediction::GetFileExtension( exposure_type)
        );

        // calculate predictions
        biol::ExposurePrediction::Calculate( model, exposure_type);

        // output the predictions
        io::File::MustOpenOFStream( write, output_name);
        biol::ExposurePrediction::WritePredictions( write, *model.GetChains().FirstElement()->GetSequence());
        io::File::CloseClearFStream( write);
      }

      return 0;
    }

    const ApplicationType ExposurePrediction::ExposurePrediction_Instance
    (
      GetAppGroups().AddAppToGroup( new ExposurePrediction(), GetAppGroups().e_BioInfo)
    );

  } // namespace app
} // namespace bcl
