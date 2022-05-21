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
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_sas_transformation.h"
#include "score/bcl_score_restraint_saxs.h"
#include "score/bcl_score_sas_type.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SaxsStatistics
    //! @brief generates statistics for protein models score by SAXS data
    //! @details compares pairwise scores for protein models and SAXS data in order to calculate statistics and other
    //!          required info for ROC curves etc
    //!
    //! @author weinerbe, putnamdk
    //! @date Jun 14, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SaxsStatistics :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! list of pdbs for doing statistics over
      util::ShPtr< command::FlagInterface> m_PDBListFlag;

      //! list of saxs data files
      util::ShPtr< command::FlagInterface> m_SaxsListFlag;

      //! output file name
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

      //! SAXS Debye implementation
      util::ShPtr< command::FlagInterface> m_DebyeImplementationFlag;

      //! Use Error Flag
      util::ShPtr< command::FlagInterface> m_UseErrors;

      //! Use Derivative
      util::ShPtr< command::FlagInterface> m_UseDerivative;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SaxsStatistics();

      //! @brief Clone function
      //! @return pointer to new SaxsStatistics
      SaxsStatistics *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

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
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType SaxsStatistics_Instance;

    }; // class SaxsStatistics

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SaxsStatistics::SaxsStatistics() :
      m_PDBListFlag
      (
        new command::FlagStatic
        (
          "pdb_list",
          "\tlist of pdbs to do statistics over",
          command::Parameter
          (
            "pdb_list",
            "\tname of file which is the list of pdbs to do statistics over",
            ""
          )
        )
      ),
      m_SaxsListFlag
      (
        new command::FlagStatic
        (
          "saxs_list",
          "\tlist of SAXS data files to do statistics over",
          command::Parameter
          (
            "saxs_list",
            "\tname of file which is the list of SAXS data files to do statistics over",
            ""
          )
        )
      ),
      m_OutputFileFlag
      (
        new command::FlagStatic
        (
          "output",
          "\toutput file name",
          command::Parameter
          (
            "output",
            "\toutput file name",
            "saxs.data"
          )
        )
      ),
      m_DebyeImplementationFlag
      (
        new command::FlagStatic
        (
          "implementation",
          "\tSAXS Debye implementation",
          command::Parameter
          (
            "implementation",
            "the object label for the saxs debye calculation enabling C++ or OpenCL implementations",
            command::ParameterCheckSerializable( util::Implementation< restraint::SasDebyeInterface>()),
            "SasDebye(consider loops=0,analytic=0, excluded volume=1, hydration shell=0, approximate_sidechains=0)"
          )
        )
      ),
      m_UseErrors
      (
        new command::FlagStatic
        (
          "use_errors",
          "pass this flag if you want to use errors in analysis"
        )
      ),
      m_UseDerivative
      (
        new command::FlagStatic
        (
          "use_derivative",
          "pass this flag if you want to use the derivative of the sas profiles"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new SaxsStatistics
    SaxsStatistics *SaxsStatistics::Clone() const
    {
      return new SaxsStatistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SaxsStatistics::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> SaxsStatistics::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add member flags
      sp_cmd->AddFlag( m_PDBListFlag);
      sp_cmd->AddFlag( m_SaxsListFlag);
      sp_cmd->AddFlag( m_OutputFileFlag);
      sp_cmd->AddFlag( m_DebyeImplementationFlag);
      sp_cmd->AddFlag( m_UseErrors);
      sp_cmd->AddFlag( m_UseDerivative);

      // add PDB factory flags
      sp_cmd->AddFlag( pdb::Factory::GetFlagConvertToNaturalAAType());
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // add default bcl flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int SaxsStatistics::Main() const
    {
      // get the list of pdb files
      const storage::Vector< std::string> pdb_list
      (
        command::StringVectorFromFilenameParameter( *( m_PDBListFlag->GetFirstParameter()))
      );

      // get the list of SAXS data files
      const storage::Vector< std::string> saxs_list
      (
        command::StringVectorFromFilenameParameter( *( m_SaxsListFlag->GetFirstParameter()))
      );

      BCL_Assert
      (
        pdb_list.GetSize() == saxs_list.GetSize(),
        "Number of PDB files does not match number of SAXS data files"
      );

      // initialize vector to hold models and data
      storage::Vector< assemble::ProteinModel> protein_models;
      protein_models.AllocateMemory( pdb_list.GetSize());
      util::ShPtrVector< restraint::SasScatteringData> saxs_data;
      saxs_data.AllocateMemory( saxs_list.GetSize());
      const pdb::Factory factory;
      size_t empty_model( 0);

      // iterate over the the lists
      for
      (
        storage::Vector< std::string>::const_iterator pdb_itr( pdb_list.Begin()), pdb_itr_end( pdb_list.End()),
          saxs_itr( saxs_list.Begin()), saxs_itr_end( saxs_list.End());
        pdb_itr != pdb_itr_end && saxs_itr != saxs_itr_end; ++pdb_itr, ++saxs_itr
      )
      {
        // read the model
        const assemble::ProteinModel protein_model( factory.ProteinModelFromPDBFilename( *pdb_itr));

        // only add the model if SSE's exist for the model
        if( protein_model.GetSSEs().IsEmpty())
        {
          ++empty_model;
        }
        else
        {
          protein_models.PushBack( protein_model);

          // read in the saxs data
          io::IFStream read;
          io::File::MustOpenIFStream( read, *saxs_itr);
          util::ShPtr< restraint::SasScatteringData> this_saxs_data( new restraint::SasScatteringData());
          this_saxs_data->ReadFromDataFile( read);
          io::File::CloseClearFStream( read);
          saxs_data.PushBack( this_saxs_data);
        }
      }

      // initialize saxs scoring functions
      util::Implementation< restraint::SasDebyeInterface>
        saxs_debye( m_DebyeImplementationFlag->GetFirstParameter()->GetValue());

      // initialize scoring matrix with only valid models
      size_t valid_model( pdb_list.GetSize() - empty_model);

      linal::Matrix< float> scores( math::Sqr( valid_model), 2);
      size_t counter( 0);

      // iterate over the protein models
      for
      (
        storage::Vector< assemble::ProteinModel>::iterator protein_itr( protein_models.Begin()),
          protein_itr_end( protein_models.End());
        protein_itr != protein_itr_end; ++protein_itr
      )
      {
        // iterate over the saxs data
        for
        (
          util::ShPtrVector< restraint::SasScatteringData>::const_iterator saxs_itr( saxs_data.Begin()),
            saxs_itr_end( saxs_data.End());
          saxs_itr != saxs_itr_end; ++saxs_itr, ++counter
        )
        {
          saxs_debye->SetExperimentalData( *saxs_itr);

          //Create container to hold experimental and computed SAXS profiles and store results in container
          util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_experimental_and_calculated_data
          (
            new restraint::SasExperimentalAndCalculatedData( saxs_debye->operator()( *protein_itr))
          );

          //Transform the data
          storage::Vector< restraint::SasTransformation::TransformationTypeEnum> data_transform;
          data_transform.PushBack( restraint::SasTransformation::e_NormalizeData);
          data_transform.PushBack( restraint::SasTransformation::e_Log10Profiles);
          data_transform.PushBack( restraint::SasTransformation::e_ScaleCalculatedProfile);

          if( m_UseDerivative->GetFlag())
          {
            data_transform.PushBack( restraint::SasTransformation::e_DerivativeProfiles);
          }

          restraint::SasExperimentalAndCalculatedData transformed_data
          (
            restraint::SasTransformation
            (
              data_transform,           // Transformations to use
              false,                    // Print Transformations
              m_UseErrors->GetFlag(),   // Use Errors
              1.0                       // Set Y Scale
            )( *sp_experimental_and_calculated_data)
          );

          // Score the transformed data
          double sas_score( score::SasType( m_UseErrors->GetFlag(), score::SasType::e_chi)( transformed_data));

          // Convert score to float
          float sas_score_float( sas_score);

          //const score::RestraintSaxs saxs_score( *saxs_debye, score::SaxsType());

          // add the score to the matrix
          scores( counter, 0) = protein_itr - protein_models.Begin() == saxs_itr - saxs_data.Begin() ? 0 : 100;
          //scores( counter, 1) = saxs_score( *protein_itr);
          scores( counter, 1) = sas_score_float;

          if( counter % 100 == 0)
          {
            BCL_MessageStd
            (
              "line: " +
              util::Format()( counter)
            );
          }

        }
      }

      // write out the file
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_OutputFileFlag->GetFirstParameter()->GetValue());
      write << scores;
      io::File::CloseClearFStream( write);

      size_t total_model( pdb_list.GetSize());
      {
        BCL_MessageStd
        (
          "Total pdbs in list: " +
          util::Format()( total_model)
        );
      }

      {
        BCL_MessageStd
        (
          "number of invalid models: " +
          util::Format()( empty_model)
        );
      }

      {
        BCL_MessageStd
        (
          "number of scored models: " +
          util::Format()( valid_model)
        );
      }

      // end
      return 0;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SaxsStatistics::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SaxsStatistics::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    const ApplicationType SaxsStatistics::SaxsStatistics_Instance
    (
      GetAppGroups().AddAppToGroup( new SaxsStatistics(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl
