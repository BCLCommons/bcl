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
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "fold/bcl_fold_default_flags.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_cone_model.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CalculateConeModelParameters
    //! @brief calculates the parameters associated with the cone model
    //! @details Calculates the SL-CB->SL max angle, SLeffective->CB->CA angle, and SLeffective->CB distance for
    //!          ensemble(s) of proteins. See bcl::restraint::ConeModel for references with more details on the cone
    //!          model.
    //!
    //! @author alexanns
    //! @date Feb 27, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CalculateConeModelParameters :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! list of protein pdbs for which the values of the cone model parameters can be calculated
      util::ShPtr< command::FlagStatic> m_PDBList;
      //! filename of file containing a list of pdbs
      util::ShPtr< command::ParameterInterface> m_PDBListFilename;

      //! file containing a list of pdb lists so the cone model parameters will be calculated for multiple mutants
      util::ShPtr< command::FlagStatic> m_ListOfPDBLists;
      //! filename of file containing a list of pdb lists
      util::ShPtr< command::ParameterInterface> m_ListOfPDBListsFilename;

      //! flag for specifying the histogram that will hold the SL->CB->SL maximum angle statistics
      util::ShPtr< command::FlagStatic> m_SLCBSLMaxAngleHistogram;
      util::ShPtr< command::ParameterInterface> m_SLCBSLMaxAngleHistogramMinimumValue;
      util::ShPtr< command::ParameterInterface> m_SLCBSLMaxAngleHistogramBinSize;
      util::ShPtr< command::ParameterInterface> m_SLCBSLMaxAngleHistogramNumberOfBins;
      util::ShPtr< command::ParameterInterface> m_SLCBSLMaxAngleHistogramOutputFilename;

      //! flag for specifying the histogram that will hold the SLeffective->CB->CA angle statistics
      util::ShPtr< command::FlagStatic> m_SLeffectiveCBCAAngleHistogram;
      util::ShPtr< command::ParameterInterface> m_SLeffectiveCBCAAngleHistogramMinimumValue;
      util::ShPtr< command::ParameterInterface> m_SLeffectiveCBCAAngleHistogramBinSize;
      util::ShPtr< command::ParameterInterface> m_SLeffectiveCBCAAngleHistogramNumberOfBins;
      util::ShPtr< command::ParameterInterface> m_SLeffectiveCBCAAngleHistogramOutputFilename;

      //! flag for specifying the histogram that will hold the SLeffective->CB distance statistics
      util::ShPtr< command::FlagStatic> m_SLeffectiveCBDistanceHistogram;
      util::ShPtr< command::ParameterInterface> m_SLeffectiveCBDistanceHistogramMinimumValue;
      util::ShPtr< command::ParameterInterface> m_SLeffectiveCBDistanceHistogramBinSize;
      util::ShPtr< command::ParameterInterface> m_SLeffectiveCBDistanceHistogramNumberOfBins;
      util::ShPtr< command::ParameterInterface> m_SLeffectiveCBDistanceHistogramOutputFilename;

      //! flag for specifying the output filename of the file that will contain the parameters as calculated for
      //! every list of pdbs
      util::ShPtr< command::FlagStatic> m_TableOutput;
      util::ShPtr< command::ParameterInterface> m_TableOutputFilename;

      //! flag indicating that histogram specifications were given in degrees and statistics should be done in degrees
      util::ShPtr< command::FlagStatic> m_UseDegrees;

      //! flag indicating that angles should be added to histogram with weight cos( ANGLE)
      util::ShPtr< command::FlagStatic> m_NormalizeAngles;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CalculateConeModelParameters();

    public:

      //! @brief Clone function
      //! @return pointer to new FitEPRDistribution
      CalculateConeModelParameters *Clone() const
      {
        return new CalculateConeModelParameters( *this);
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

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

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

      //! @brief creates a protein ensemble from a vector of pdb filenames
      //! @param PDB_FILENAMES the vector of pdb filenames from which the ensemble will be created
      //! @return protein ensemble which has the protein models that have been created out of PDB_FILENAMES
      assemble::ProteinEnsemble GetProteinEnsemble( const storage::Vector< std::string> &PDB_FILENAMES) const;

      //! @brief provides the angle in either the provided radians or converts it to degrees depending on m_UseDegrees
      //! @param RADIAN angle which might be converted to degrees if desired by user according to m_UseDegrees
      //! @return double which is either the provided RADIAN angle or the same angle converted to degrees
      double GetDegreeOrRadian( const double RADIAN) const;

      //! @brief determines the weight that should be given to the provided angle (in radians)
      //!        this is the weight that is given to the angle when it is added to the histogram
      //! @param RADIAN the angle in radians whose weight will be determined
      //! @return double which is the weight of the angle provided in terms of how much it adds to the histogram counts
      double GetAngleWeight( const double RADIAN) const;

    private:

      //! single instance of that class
      static const ApplicationType s_Instance;

    }; // class CalculateConeModelParameters

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> CalculateConeModelParameters::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // list of pdbs
      sp_cmd->AddFlag( m_PDBList);

      // list of pdb lists
      sp_cmd->AddFlag( m_ListOfPDBLists);

      // flag for specifying amino acid class type
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // prefix flag
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagPrefix());

      // flag for specifying the histogram that will hold the SL->CB->SL maximum angle statistics
      sp_cmd->AddFlag( m_SLCBSLMaxAngleHistogram);

      // flag for specifying the histogram that will hold the SLeffective->CB->CA angle statistics
      sp_cmd->AddFlag( m_SLeffectiveCBCAAngleHistogram);

      // flag for specifying the histogram that will hold the SLeffective->CB->CA angle statistics
      sp_cmd->AddFlag( m_SLeffectiveCBDistanceHistogram);

      // output filename of the file that will contain the parameters as calculated for every list of pdbs
      sp_cmd->AddFlag( m_TableOutput);

      // flag indicating that histogram specifications were given in degrees and statistics should be done in degrees
      sp_cmd->AddFlag( m_UseDegrees);

      // flag indicating that angles should be added to histogram with weight cos( ANGLE)
      sp_cmd->AddFlag( m_NormalizeAngles);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return command
      return sp_cmd;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &CalculateConeModelParameters::GetReadMe() const
    {
      static const std::string s_readme;
      return s_readme;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int CalculateConeModelParameters::Main() const
    {
      // will hold the pdb list filenames that are provided - sometimes just one, sometimes more
      storage::Vector< std::string> pdb_list_filenames;

      // just a single pdb list filename was provided
      if( m_PDBList->GetFlag())
      {
        // add the filename containing the list of pdbs to pdb_lists
        pdb_list_filenames.PushBack( m_PDBListFilename->GetValue());
      }
      else if( m_ListOfPDBLists->GetFlag()) //< list of pdb list files was given
      {
        // read all the pdb list filenames from the provided file and set pdb_lists to the list of pdb list filenames
        io::IFStream read;
        const std::string &list_of_pdb_lists_filename( m_ListOfPDBListsFilename->GetValue());
        io::File::MustOpenIFStream( read, list_of_pdb_lists_filename);
        pdb_list_filenames = util::StringLineListFromIStream( read);
        io::File::CloseClearFStream( read);
      }

      // table to hold the results
      const storage::TableHeader header
      (
        storage::Vector< std::string>::Create
        (
          "sl_cb_max_angle", "cos( sl_cb_max_angle)", "SLeffective_cb_ca_angle", "cos( SLeffective_cb_ca_angle)",
          "SLeffective_cb_distance"
        )
      );
      storage::Table< double> results_table( header);

      // data sets to hold the statistics of the cone model parameters as calculated over the pdb lists
      math::RunningAverageSD< double> sl_cb_sl_max_angle_statistics;
      math::RunningAverageSD< double> sleffective_cb_ca_angle_statistics;
      math::RunningAverageSD< double> sleffective_cb_distance_statistics;
      math::RunningAverageSD< double> cos_sl_cb_sl_max_angle_statistics;
      math::RunningAverageSD< double> cos_sleffective_cb_ca_angle_statistics;

      // histograms to hold the frequencies with which a parameter value is observed
      math::Histogram sl_cb_sl_max_angle_histogram
      (
        m_SLCBSLMaxAngleHistogramMinimumValue->GetNumericalValue< double>(),
        m_SLCBSLMaxAngleHistogramBinSize->GetNumericalValue< double>(),
        m_SLCBSLMaxAngleHistogramNumberOfBins->GetNumericalValue< double>()
      );
      math::Histogram sleffective_cb_ca_angle_histogram
      (
        m_SLeffectiveCBCAAngleHistogramMinimumValue->GetNumericalValue< double>(),
        m_SLeffectiveCBCAAngleHistogramBinSize->GetNumericalValue< double>(),
        m_SLeffectiveCBCAAngleHistogramNumberOfBins->GetNumericalValue< double>()
      );
      math::Histogram sl_cb_sl_max_angle_histogram_cos
      (
        -1.0,
         0.2,
          10
      );
      math::Histogram sleffective_cb_ca_angle_histogram_cos
      (
        -1.0,
         0.2,
          10
      );
      math::Histogram sleffective_cb_distance_histogram
      (
        m_SLeffectiveCBDistanceHistogramMinimumValue->GetNumericalValue< double>(),
        m_SLeffectiveCBDistanceHistogramBinSize->GetNumericalValue< double>(),
        m_SLeffectiveCBDistanceHistogramNumberOfBins->GetNumericalValue< double>()
      );

      // iterate through the pdb list filenames
      for
      (
        storage::Vector< std::string>::const_iterator
          pdb_list_filename_itr( pdb_list_filenames.Begin()), pdb_list_filename_itr_end( pdb_list_filenames.End());
        pdb_list_filename_itr != pdb_list_filename_itr_end;
        ++pdb_list_filename_itr
      )
      {
        // get the list of pdbs contained in the list file
        io::IFStream read;
        io::File::MustOpenIFStream( read, *pdb_list_filename_itr);
        storage::Vector< std::string> pdb_filenames( util::StringLineListFromIStream( read));
        io::File::CloseClearFStream( read);

        // get the ensemble from the current list of pdbs
        assemble::ProteinEnsemble protein_ensemble( GetProteinEnsemble( pdb_filenames));

        // get the cone model parameters as calculated from the current ensemble
        const double sl_cb_sl_max_angle( restraint::ConeModel::SLCBSLMaxAngle( protein_ensemble));
        BCL_MessageDbg( "sl_cb_sl_max_angle" + util::Format()( sl_cb_sl_max_angle));
        const double sleffective_cb_ca_angle( restraint::ConeModel::SLeffectiveCBCAAngle( protein_ensemble));
        const double sleffective_cb_distance( restraint::ConeModel::SLeffectiveCBDistance( protein_ensemble));

        // add the current parameter values to the statistics objects
        sl_cb_sl_max_angle_statistics += GetDegreeOrRadian( sl_cb_sl_max_angle);
        sleffective_cb_ca_angle_statistics +=  GetDegreeOrRadian( sleffective_cb_ca_angle);
        sleffective_cb_distance_statistics +=  sleffective_cb_distance;
        cos_sl_cb_sl_max_angle_statistics +=  cos( sl_cb_sl_max_angle);
        cos_sleffective_cb_ca_angle_statistics += cos( sleffective_cb_ca_angle);

        // add the current parameter values to the histograms
        // normalize these two by cosine by giving the value a weight of the cosine of the value
        // if the flag to use degrees is set then convert the radians into degrees
        sl_cb_sl_max_angle_histogram.PushBack
        (
          GetDegreeOrRadian( sl_cb_sl_max_angle), GetAngleWeight( sl_cb_sl_max_angle)
        );
        sleffective_cb_ca_angle_histogram.PushBack
        (
          GetDegreeOrRadian( sleffective_cb_ca_angle), GetAngleWeight( sleffective_cb_ca_angle)
        );
        sl_cb_sl_max_angle_histogram_cos.PushBack( cos( sl_cb_sl_max_angle));
        sleffective_cb_ca_angle_histogram_cos.PushBack( cos( sleffective_cb_ca_angle));
        // always has weight of one
        sleffective_cb_distance_histogram.PushBack( sleffective_cb_distance);

        // add the current parameter values to the table of output information
        results_table.InsertRow
        (
          *pdb_list_filename_itr,
          storage::Vector< double>::Create
          (
            GetDegreeOrRadian( sl_cb_sl_max_angle),
            cos( sl_cb_sl_max_angle),
            GetDegreeOrRadian( sleffective_cb_ca_angle),
            cos( sleffective_cb_ca_angle),
            sleffective_cb_distance
          )
        );
      } //< end iteration through list of pdb list filenames

      // add the mean values to the table
      results_table.InsertRow
      (
        "mean",
        storage::Vector< double>::Create
        (
              sl_cb_sl_max_angle_statistics.GetAverage(),
          cos_sl_cb_sl_max_angle_statistics.GetAverage(),
              sleffective_cb_ca_angle_statistics.GetAverage(),
          cos_sleffective_cb_ca_angle_statistics.GetAverage(),
              sleffective_cb_distance_statistics.GetAverage()
        )
      );

      // add the standard deviation values to the table
      results_table.InsertRow
      (
        "standard_deviation",
        storage::Vector< double>::Create
        (
              sl_cb_sl_max_angle_statistics.GetStandardDeviation(),
          cos_sl_cb_sl_max_angle_statistics.GetStandardDeviation(),
              sleffective_cb_ca_angle_statistics.GetStandardDeviation(),
          cos_sleffective_cb_ca_angle_statistics.GetStandardDeviation(),
              sleffective_cb_distance_statistics.GetStandardDeviation()
        )
      );

      // write the table information to the table output file specified by user
      {
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_TableOutputFilename->GetValue());
        results_table.WriteFormatted( write);
        io::File::CloseClearFStream( write);
      }

      // write the histograms to their files
      { // SLCBSLMaxAngleHistogram
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_SLCBSLMaxAngleHistogramOutputFilename->GetValue());
        sl_cb_sl_max_angle_histogram.WriteHorizontally( write);
        io::File::CloseClearFStream( write);
      }
      { // SLCBSLMaxAngleHistogram cosine
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_SLCBSLMaxAngleHistogramOutputFilename->GetValue() + ".cos");
        sl_cb_sl_max_angle_histogram_cos.WriteHorizontally( write);
        io::File::CloseClearFStream( write);
      }
      { // SLeffectiveCBCAAngleHistogram
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_SLeffectiveCBCAAngleHistogramOutputFilename->GetValue());
        sleffective_cb_ca_angle_histogram.WriteHorizontally( write);
        io::File::CloseClearFStream( write);
      }
      { // SLeffectiveCBCAAngleHistogram cosine
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_SLeffectiveCBCAAngleHistogramOutputFilename->GetValue() + ".cos");
        sleffective_cb_ca_angle_histogram_cos.WriteHorizontally( write);
        io::File::CloseClearFStream( write);
      }
      { // SLeffectiveCBDistanceHistogram
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_SLeffectiveCBDistanceHistogramOutputFilename->GetValue());
        sleffective_cb_distance_histogram.WriteHorizontally( write);
        io::File::CloseClearFStream( write);
      }

      return 0;
    }

    //! @brief creates a protein ensemble from a vector of pdb filenames
    //! @param PDB_FILENAMES the vector of pdb filenames from which the ensemble will be created
    //! @return protein ensemble which has the protein models that have been created out of PDB_FILENAMES
    assemble::ProteinEnsemble
    CalculateConeModelParameters::GetProteinEnsemble( const storage::Vector< std::string> &PDB_FILENAMES) const
    {
      // ensemble to hold models created from pdb list
      assemble::ProteinEnsemble ensemble;

      // factory to create protein models
      pdb::Factory factory;

      // iterate through the pdb filenames and create protein models from them and add them to "ensemble"
      for
      (
        storage::Vector< std::string>::const_iterator
          pdb_itr( PDB_FILENAMES.Begin()), pdb_itr_end( PDB_FILENAMES.End());
        pdb_itr != pdb_itr_end; ++pdb_itr
      )
      {
        ensemble.InsertElement
        (
          util::ShPtr< assemble::ProteinModel>( factory.ProteinModelFromPDBFilename( *pdb_itr).Clone())
        );
      }

      BCL_Assert( !ensemble.IsEmpty(), "ensemble is empty");

      // return the ensemble containing all of the models
      return ensemble;
    }

    //! @brief provides the angle in either the provided radians or converts it to degrees depending on m_UseDegrees
    //! @param RADIAN angle which might be converted to degrees if desired by user according to m_UseDegrees
    //! @return double which is either the provided RADIAN angle or the same angle converted to degrees
    double CalculateConeModelParameters::GetDegreeOrRadian( const double RADIAN) const
    {
      return m_UseDegrees->GetFlag() ? math::Angle::Degree( RADIAN) : RADIAN;
    }

    //! @brief determines the weight that should be given to the provided angle (in radians)
    //!        this is the weight that is given to the angle when it is added to the histogram
    //! @param RADIAN the angle in radians whose weight will be determined
    //! @return double which is the weight of the angle provided in terms of how much it adds to the histogram counts
    double CalculateConeModelParameters::GetAngleWeight( const double RADIAN) const
    {
      return m_NormalizeAngles->GetFlag() ? cos( RADIAN) : 1.0;
    }

    //! @brief default constructor
    CalculateConeModelParameters::CalculateConeModelParameters() :
      m_PDBList
      (
        new command::FlagStatic
        (
          "pdb_list",
          "\tlist of protein pdbs for which the values of the cone model parameters can be calculated. Typically this should be an ensemble of models for a one spin label mutant (both in the sense that the mutant should have a single spin label and the ensemble should have the spin label at the same sequence position.)"
        )
      ),
      m_PDBListFilename
      (
        new command::Parameter
        (
          "filename",
          "\tfilename of file containing a list of pdbs",
          "single_mutant_pdbs.ls"
        )
      ),
      m_ListOfPDBLists
      (
        new command::FlagStatic
        (
          "list_of_pdb_lists",
          "\tfile containing a list of pdb lists so the cone model parameters will be calculated for single mutants and each list is an ensemble with the spin label at a different position. This is useful for doing statistics on the cone model parameters over multiple sequence positions and/or proteins"
        )
      ),
      m_ListOfPDBListsFilename
      (
        new command::Parameter
        (
          "filename",
          "\tfilename of file containing a list of pdb lists",
          "single_mutant_pdb_lists.ls"
        )
      ),
      m_SLCBSLMaxAngleHistogram
      (
        new command::FlagStatic
        (
          "sl_cb_sl_max_angle",
          "\tflag for specifying the histogram that will hold the SL->CB->SL maximum angle statistics"
        )
      ),
      m_SLCBSLMaxAngleHistogramMinimumValue
      (
        new command::Parameter
        (
          "min_value",
          "\tminimum value the histogram will hold",
          util::Format()( 0)
        )
      ),
      m_SLCBSLMaxAngleHistogramBinSize
      (
        new command::Parameter
        (
          "bin_size",
          "\tthe size of the bins that make up the histogram",
          util::Format()( math::Angle::Radian( 20.0))
        )
      ),
      m_SLCBSLMaxAngleHistogramNumberOfBins
      (
        new command::Parameter
        (
          "number_of_bins",
          "\tnumber of bins in the histogram",
          util::Format()( 9)
        )
      ),
      m_SLCBSLMaxAngleHistogramOutputFilename
      (
        new command::Parameter
        (
          "output_filename",
          "\tfilename that the histogram will be written to",
          "SLCBSLMaxAngle.bcl_histogram"
        )
      ),
      m_SLeffectiveCBCAAngleHistogram
      (
        new command::FlagStatic
        (
          "sleffective_cb_ca_angle",
          "\tflag for specifying the histogram that will hold the SLeffective->CB->CA angle statistics"
        )
      ),
      m_SLeffectiveCBCAAngleHistogramMinimumValue
      (
        new command::Parameter
        (
          "min_value",
          "\tminimum value the histogram will hold",
          util::Format()( 0)
        )
      ),
      m_SLeffectiveCBCAAngleHistogramBinSize
      (
        new command::Parameter
        (
          "bin_size",
          "\tthe size of the bins that make up the histogram",
          util::Format()( math::Angle::Radian( 20.0))
        )
      ),
      m_SLeffectiveCBCAAngleHistogramNumberOfBins
      (
        new command::Parameter
        (
          "number_of_bins",
          "\tnumber of bins in the histogram",
          util::Format()( 9)
        )
      ),
      m_SLeffectiveCBCAAngleHistogramOutputFilename
      (
        new command::Parameter
        (
          "output_filename",
          "\tfilename that the histogram will be written to",
          "SLeffectiveCBCAAngle.bcl_histogram"
        )
      ),
      m_SLeffectiveCBDistanceHistogram
      (
        new command::FlagStatic
        (
          "sleffective_cb_distance",
          "\tflag for specifying the histogram that will hold the SLeffective->CB distance statistics"
        )
      ),
      m_SLeffectiveCBDistanceHistogramMinimumValue
      (
        new command::Parameter
        (
          "min_value",
          "\tminimum value the histogram will hold",
          util::Format()( 0)
        )
      ),
      m_SLeffectiveCBDistanceHistogramBinSize
      (
        new command::Parameter
        (
          "bin_size",
          "\tthe size of the bins that make up the histogram",
          util::Format()( 1.0)
        )
      ),
      m_SLeffectiveCBDistanceHistogramNumberOfBins
      (
        new command::Parameter
        (
          "number_of_bins",
          "\tnumber of bins in the histogram",
          util::Format()( 9)
        )
      ),
      m_SLeffectiveCBDistanceHistogramOutputFilename
      (
        new command::Parameter
        (
          "output_filename",
          "\tfilename that the histogram will be written to",
          "SLeffectiveCBDistance.bcl_histogram"
        )
      ),
      m_TableOutput
      (
        new command::FlagStatic
        (
          "parameter_table_filename",
          "\toutput filename of the file that will contain the parameters as calculated for every list of pdbs"
        )
      ),
      m_TableOutputFilename
      (
        new command::Parameter
        (
          "filename",
          "\toutput filename for the table",
          "cone_model_parameters.table"
        )
      ),
      m_UseDegrees
      (
        new command::FlagStatic
        (
          "use_degrees",
          "\tindicates that histogram specifications are given in degrees and statistics should be calculated in degrees"
        )
      ),
      m_NormalizeAngles
      (
        new command::FlagStatic
        (
          "normalize_angles_by_cosine",
          "\tangles will be added to histogram with weight of cos( ANGLE)"
        )
      )
    {
      // attach parameters to flags
      m_PDBList->PushBack( m_PDBListFilename);
      m_ListOfPDBLists->PushBack( m_ListOfPDBListsFilename);
      m_SLCBSLMaxAngleHistogram->PushBack( m_SLCBSLMaxAngleHistogramMinimumValue);
      m_SLCBSLMaxAngleHistogram->PushBack( m_SLCBSLMaxAngleHistogramBinSize);
      m_SLCBSLMaxAngleHistogram->PushBack( m_SLCBSLMaxAngleHistogramNumberOfBins);
      m_SLCBSLMaxAngleHistogram->PushBack( m_SLCBSLMaxAngleHistogramOutputFilename);
      m_SLeffectiveCBCAAngleHistogram->PushBack( m_SLeffectiveCBCAAngleHistogramMinimumValue);
      m_SLeffectiveCBCAAngleHistogram->PushBack( m_SLeffectiveCBCAAngleHistogramBinSize);
      m_SLeffectiveCBCAAngleHistogram->PushBack( m_SLeffectiveCBCAAngleHistogramNumberOfBins);
      m_SLeffectiveCBCAAngleHistogram->PushBack( m_SLeffectiveCBCAAngleHistogramOutputFilename);
      m_SLeffectiveCBDistanceHistogram->PushBack( m_SLeffectiveCBDistanceHistogramMinimumValue);
      m_SLeffectiveCBDistanceHistogram->PushBack( m_SLeffectiveCBDistanceHistogramBinSize);
      m_SLeffectiveCBDistanceHistogram->PushBack( m_SLeffectiveCBDistanceHistogramNumberOfBins);
      m_SLeffectiveCBDistanceHistogram->PushBack( m_SLeffectiveCBDistanceHistogramOutputFilename);
      m_TableOutput->PushBack( m_TableOutputFilename);
    }

    const ApplicationType CalculateConeModelParameters::s_Instance
    (
      GetAppGroups().AddAppToGroup( new CalculateConeModelParameters(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl
