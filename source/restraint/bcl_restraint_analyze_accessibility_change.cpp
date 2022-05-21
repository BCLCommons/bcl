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
#include "restraint/bcl_restraint_analyze_accessibility_change.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "biol/bcl_biol_aa_classes.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialization_via_static_functions.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"
#include "restraint/bcl_restraint_accessibility_aa_assignment.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAccessibilityChange::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAccessibilityChange())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAccessibilityChange::AnalyzeAccessibilityChange() :
      m_OutFilePostFix( ".AnalyzeAccessibilityChange"),
      m_StartEnsemble(),
      m_ExposureMethod( util::CloneToShPtr( assemble::AANeighborCount())),
      m_ExperimentalExposures(),
      m_MeanMinCutoff(),
      m_MeanMaxCutoff(),
      m_ZScoreMinCutoff(),
      m_ZScoreMaxCutoff(),
      m_PymolOuputFilename(),
      m_EnsembleRepresentativeIndex(),
      m_EnsembleRepresentativeFromStartEnsemble(),
      m_GradientMin(),
      m_GradientMax(),
      m_DirectRelation(),
      m_SetNCRange(),
      m_NCRangeMin(),
      m_NCRangeMax()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAccessibilityChange
    AnalyzeAccessibilityChange *AnalyzeAccessibilityChange::Clone() const
    {
      return new AnalyzeAccessibilityChange( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAccessibilityChange::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAccessibilityChange::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAccessibilityChange::GetAlias() const
    {
      static const std::string s_Name( "AccessibilityChange");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAccessibilityChange::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // to hold, over all profiles, the average of the average exposure difference calculated across all the residues
      // in the profile
      math::RunningAverageSD< double> average_mean_diff_calculated;
      // ranges for normalizing the calculated exposures in to the range of the experimental data
      const storage::VectorND< 2, math::Range< double> > mean_range_exp_range
      (
        GetCalculatedAndExperimentalRanges( ENSEMBLE, average_mean_diff_calculated)
      );
      const math::Range< double> &mean_range( mean_range_exp_range.First());
      const math::Range< double> &exp_range( mean_range_exp_range.Second());

      const double all_mean_change_mean( average_mean_diff_calculated.GetAverage()),
                   all_mean_change_stddev( average_mean_diff_calculated.GetStandardDeviation());

      BCL_MessageDbg( "all_mean_change_mean "   + util::Format()( all_mean_change_mean));
      BCL_MessageDbg( "all_mean_change_stddev " + util::Format()( all_mean_change_stddev));

      // to hold the analysis text
      std::string analysis;

      // open write to pymol output script
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_PymolOuputFilename);

      // iterator to the structure for visualization in pymol, by default comes from the end ensemble
      assemble::ProteinEnsemble::const_iterator represent_itr( ENSEMBLE.Begin()), represent_itr_end( ENSEMBLE.End());

      // true if the representative structure should come from the start ensemble
      if( m_EnsembleRepresentativeFromStartEnsemble)
      {
        // set iterator to start ensemble
        represent_itr     = m_StartEnsemble.Begin();
        represent_itr_end = m_StartEnsemble.End();
      }

      // move the iterator to the correct structure
      storage::AdvanceIterator( represent_itr, represent_itr_end, m_EnsembleRepresentativeIndex);

      // get the name of the pdb the structure came from
      util::ShPtr< util::Wrapper< std::string> > pdb_filename
      (
        ( *represent_itr)->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );

      // write intitial pymol commands
      {
        std::string pdb( "dummy.pdb");

        if( pdb_filename.IsDefined())
        {
          pdb = std::string( pdb_filename->GetData());
        }

        write <<  "load " << pdb << ",access_model\n";
        write <<  "cmd.show_as(\"cartoon\"   ,\"access_model\")\n";
        write <<  "color grey70, access_model\n";
        write << "alter access_model, q=" + util::Format()( 0.0) << '\n';
      }

      // string for selection of residues with experimental information
      std::string experimental_residues( "select exp_resi, ");

      // vectors for experimental and accessibility values for calculating the spearman correlation
      storage::Vector< double> exp_vals;
      storage::Vector< double> access_vals;

      // iterate through the data
      for
      (
        storage::List< AccessibilityProfile>::const_iterator profile_itr( m_ExperimentalExposures.Begin()),
        profile_itr_end( m_ExperimentalExposures.End());
        profile_itr != profile_itr_end;
        ++profile_itr
      )
      {
        // to hold, for a profile, the average exposure difference calculated across all the residues in the profile
        math::RunningAverageSD< double> average_profile_diff;

        // to hold the experimental data associated with this profile
        double exp_data;

        // iterate through the profile
        for
        (
          storage::List< AccessibilityAA>::const_iterator aa_itr( profile_itr->GetAccessibilities().Begin()),
          aa_itr_end( profile_itr->GetAccessibilities().End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // reference to current data
          const AccessibilityAA &current_data( *aa_itr);

          // iterate over the start ensemble
          for
          (
            assemble::ProteinEnsemble::const_iterator
              start_ensemble_itr( m_StartEnsemble.Begin()), start_ensemble_itr_end( m_StartEnsemble.End());
            start_ensemble_itr != start_ensemble_itr_end;
            ++start_ensemble_itr
          )
          {
            // the accessibility assignment for the starting ensemble
            const AccessibilityAAAssignment assignment_start( current_data.GenerateAssignment( **start_ensemble_itr));

            // iterate over the end ensemble
            for
            (
              assemble::ProteinEnsemble::const_iterator
               end_ensemble_itr( ENSEMBLE.Begin()), end_ensemble_itr_end( ENSEMBLE.End());
              end_ensemble_itr != end_ensemble_itr_end;
              ++end_ensemble_itr
            )
            {
              // the accessibility assignment for the ending ensemble
              AccessibilityAAAssignment assignment_end( current_data.GenerateAssignment( **end_ensemble_itr));

              // right now the data is just stored as oxygen accessibility
              // this is reset with each iteration but all exp data is the same across the entire profile
              exp_data = assignment_end.GetAccessibilityByEnvironment
              (
                AccessibilityAA::e_Oxygen
              ).Second();

              // add the current exposure difference
              average_profile_diff += assignment_end.GetExposureValue() - assignment_start.GetExposureValue();
            }
          } // start ensemble
        } // profile aas

        // get the profile and locators for the beginning and ending residues
        const storage::List< AccessibilityAA> &profile( profile_itr->GetAccessibilities());
        // casts from LocatorInterface to LocatorAA
        util::ShPtr< assemble::LocatorAA> start_locator( profile.FirstElement().GetAA());
        util::ShPtr< assemble::LocatorAA> end_locator( profile.LastElement().GetAA());

        // statistical values related to the accessibility profile
        const double stddev( average_profile_diff.GetStandardDeviation());
        const double mean( average_profile_diff.GetAverage());
        double zscore( stddev == 0 ? util::GetUndefinedDouble() : mean / stddev);

        // true if the profile meets the user specified cutoffs
        if
        (
          math::Absolute( mean) >= m_MeanMinCutoff && math::Absolute( mean) <= m_MeanMaxCutoff &&
          (
            (
              util::IsDefined( zscore) && math::Absolute( zscore) >= m_ZScoreMinCutoff &&
              math::Absolute( zscore) <= m_ZScoreMaxCutoff
            )
            ||
            ( !util::IsDefined( zscore))
          )
        )
        {
          double overall_zscore( ( mean - all_mean_change_mean) / all_mean_change_stddev);
          if( !m_DirectRelation)
          {
            overall_zscore = ( -mean - all_mean_change_mean) / all_mean_change_stddev;
          }

          // analysis file
          analysis += "|" + start_locator->GetIdentification() + "| |";
          analysis += end_locator->GetIdentification() + "|\texp_data:\t";
          analysis += util::Format()( exp_data) + "\tmodel_mean_change:\t";
          analysis += util::Format()( mean) + "\tmodel_stdev:\t";
          analysis += util::Format()( stddev) + "\tzscore:\t";
          analysis += util::Format()( zscore) + "\tre-ranged_mean:\t";
          analysis += util::Format()( exp_range.Rescale( mean, mean_range));
          analysis += util::Format()( "\toverall_zscore_exp_related\t");
          analysis += util::Format()( overall_zscore) + "\n";

          // for spearman correlation
          access_vals.PushBack( mean);
          exp_vals.PushBack( exp_data);
          BCL_MessageDbg( "mean " + util::Format()( mean) + " expval " + util::Format()( exp_data));

          // pymol file
          // iterate over the residue locators for the current profile
          for
          (
            storage::List< AccessibilityAA>::const_iterator aa_itr( profile_itr->GetAccessibilities().Begin()),
            aa_itr_end( profile_itr->GetAccessibilities().End());
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            // need to have a "+" if not the first residue for the selection
            if( aa_itr != profile_itr->GetAccessibilities().Begin())
            {
              experimental_residues += "+";
            }

            // reference to data
            const AccessibilityAA &current_data( *aa_itr);

            // try to cast the locator to a LocatorAA
            util::ShPtr< assemble::LocatorAA> current_locator( current_data.GetAA());

            // message and continue if the locator is not defined
            if( !current_locator.IsDefined())
            {
              BCL_MessageCrt
              (
                "locator " + util::Format()( current_data.GetAA()) +
                " could not be cast to a LocatorAA"
              );
              continue;
            }

            // set b factor to calculated exposure
//            write << "alter access_model and chain " << current_locator->GetLocatorChain().GetChainID()
//                  << " and resi " + util::Format()( current_locator->GetAAID()) + ", b=" +
//                     util::Format()( exp_range.Rescale( mean, mean_range))
//                  << '\n';

            // show sticks
//            write << "cmd.show( \"sticks\", \"access_model and chain "
//                  << current_locator->GetLocatorChain().GetChainID()
//                  << " and resi " + util::Format()( current_locator->GetAAID()) << "\")\n";

            // set q value to experimental data
            write << "alter access_model and chain " << current_locator->GetLocatorChain().GetChainID()
                  << " and resi " + util::Format()( current_locator->GetAAID()) + " and name CA, q=" +
                     util::Format()( exp_data - overall_zscore)
                  << '\n';

            // add residue to experimental residues selection
            experimental_residues +=
            (
                std::string( "(access_model and chain ")
              + util::Format()( current_locator->GetLocatorChain().GetChainID())
              + std::string( " and resi ") + util::Format()( current_locator->GetAAID()) + " )"
            );
          } // pymol file : iterate over the residue locators for the current profile
        } // if residue meets criteria
      } // profiles

      // calculate spearman correlation
      const double correlation
      (
        math::Statistics::CorrelationSpearman( access_vals.Begin(), access_vals.End(), exp_vals.Begin(), exp_vals.End())
      );

      BCL_MessageStd( "CorrelationSpearman " + util::Format()( correlation));

      // write experimental residues selection
      write << experimental_residues << "\n";

      // write the ramp information
//      write <<  "ramp_new grad, access_model, [" << m_GradientMin << ", " << ( m_GradientMax + m_GradientMin) / 2.0
//            << ", "
//            << m_GradientMax << "],[red,green,blue]\n";

      // write the spectrum information
      write <<  "spectrum q, rainbow_rev, exp_resi, minimum=" << m_GradientMin << "," << "maximum=" << m_GradientMax
            << "\n";

      // write spectrum for experimental values
//      std::string experimental_gradient( m_DirectRelation ? "blue_white_red" : "red_white_blue");
//      write <<  "spectrum q, " << experimental_gradient << ", exp_resi and name CA, minimum=" << m_GradientMin << ","
//            << "maximum=" << m_GradientMax << "\n";

      // hide hydrogens and close pymol script file
      write << "cmd.hide(\"(all and hydro)\")\n";
      io::File::CloseClearFStream( write);

      // return the analysis text
      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAccessibilityChange::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_OutFilePostFix,                          ISTREAM);
      io::Serialize::Read( m_StartEnsemble,                           ISTREAM);
      io::Serialize::Read( m_ExposureMethod,                          ISTREAM);
      io::Serialize::Read( m_ExperimentalExposures,                   ISTREAM);
      io::Serialize::Read( m_MeanMinCutoff,                           ISTREAM);
      io::Serialize::Read( m_MeanMaxCutoff,                           ISTREAM);
      io::Serialize::Read( m_ZScoreMinCutoff,                         ISTREAM);
      io::Serialize::Read( m_ZScoreMaxCutoff,                         ISTREAM);
      io::Serialize::Read( m_PymolOuputFilename,                      ISTREAM);
      io::Serialize::Read( m_EnsembleRepresentativeIndex,             ISTREAM);
      io::Serialize::Read( m_EnsembleRepresentativeFromStartEnsemble, ISTREAM);
      io::Serialize::Read( m_GradientMin,                             ISTREAM);
      io::Serialize::Read( m_GradientMax,                             ISTREAM);
      io::Serialize::Read( m_DirectRelation,                          ISTREAM);
      io::Serialize::Read( m_SetNCRange,                              ISTREAM);
      io::Serialize::Read( m_NCRangeMin,                              ISTREAM);
      io::Serialize::Read( m_NCRangeMax,                              ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAccessibilityChange::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_OutFilePostFix,                          OSTREAM, INDENT);
      io::Serialize::Write( m_StartEnsemble,                           OSTREAM, INDENT);
      io::Serialize::Write( m_ExposureMethod,                          OSTREAM, INDENT);
      io::Serialize::Write( m_ExperimentalExposures,                   OSTREAM, INDENT);
      io::Serialize::Write( m_MeanMinCutoff,                           OSTREAM, INDENT);
      io::Serialize::Write( m_MeanMaxCutoff,                           OSTREAM, INDENT);
      io::Serialize::Write( m_ZScoreMinCutoff,                         OSTREAM, INDENT);
      io::Serialize::Write( m_ZScoreMaxCutoff,                         OSTREAM, INDENT);
      io::Serialize::Write( m_PymolOuputFilename,                      OSTREAM, INDENT);
      io::Serialize::Write( m_EnsembleRepresentativeIndex,             OSTREAM, INDENT);
      io::Serialize::Write( m_EnsembleRepresentativeFromStartEnsemble, OSTREAM, INDENT);
      io::Serialize::Write( m_GradientMin,                             OSTREAM, INDENT);
      io::Serialize::Write( m_GradientMax,                             OSTREAM, INDENT);
      io::Serialize::Write( m_DirectRelation,                          OSTREAM, INDENT);
      io::Serialize::Write( m_SetNCRange,                              OSTREAM, INDENT);
      io::Serialize::Write( m_NCRangeMin,                              OSTREAM, INDENT);
      io::Serialize::Write( m_NCRangeMax,                              OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAccessibilityChange::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates accessibility changes as (EndEnsembleAccessibility) - (BeginEnsembleAccessibility). Takes two "
        "ensembles of models. Calculates the mean, stddev, zscore of exposure change between them"
        " for residues that have experimental data. The experimental data is provided by an input file. Cutoffs "
        "(min and max) can be given for a residue to be outputted as of interest. Outputs a text file with "
        "these numbers and the information about the residues involved. Also, outputs a pymol script file "
        "which colors the backbone of residues according to the experimental accessibility change and colors the "
        "side chains of residues according to the calculated accessibility change."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAccessibilityChange"
      );
      parameters.AddInitializer
      (
        "start_ensemble_filename",
        "the name of the file containing the list of pdbs for the starting ensemble",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< assemble::ProteinEnsemble>
          (
            command::ParameterCheckFileExistence(),
            ( &EnsembleAsFilename),
            ( &EnsembleFromFilename),
            &m_StartEnsemble
          )
        ),
        "start_ensemble_pdbs.ls"
      );

      parameters.AddInitializer
      (
        "experimental_data_filename",
        "the name of the file containing the list of experimental data\nShould have the format\n"
        "'A' 33  'A' 38  2\n'A' 41 'A' 52  2\n'A' 53  'A' 81  1\n"
        "Where -2 = large accessibility decrease; -1 = decrease; 0 = no change; 1 = increase; 2 = large "
        "accessibility increase",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< storage::List< AccessibilityProfile> >
          (
            command::ParameterCheckFileExistence(),
            ( &ExposureDataAsFilename),
            ( &ExposureDataFromFilename),
            &m_ExperimentalExposures
          )
        ),
        "exposure_data.cst"
      );

      parameters.AddInitializer
      (
        "mean_min_cutoff",
        "The absolute value (inclusive) of the minimum mean calculated accessibility change for a residue to be "
        "considered at all",
        io::Serialization::GetAgent( &m_MeanMinCutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "mean_max_cutoff",
        "The absolute value (inclusive) of the maximum mean calculated accessibility change for a residue to be "
        "considered at all",
        io::Serialization::GetAgent( &m_MeanMaxCutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "zscore_min_cutoff",
        "The absolute value (inclusive) of the minimum zscore calculated for an accessibility change for a "
        "residue to be considered at all",
        io::Serialization::GetAgent( &m_ZScoreMinCutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "zscore_max_cutoff",
        "The absolute value (inclusive) of the maximum zscore calculated for an accessibility change for a "
        "residue to be considered at all",
        io::Serialization::GetAgent( &m_ZScoreMaxCutoff),
        "0"
      );

      parameters.AddInitializer
      (
        "pymol_output_filename",
        "The filename of the pymol script that will be outputted showing the accessibility changes",
        io::Serialization::GetAgent( &m_PymolOuputFilename),
        "accessibilities.pml"
      );

      parameters.AddInitializer
      (
        "ensemble_representative_index",
        "The index of the model in the desired ensemble list that should be used as the representative in the "
        "pymol script. First model = index 0",
        io::Serialization::GetAgent( &m_EnsembleRepresentativeIndex),
        "0"
      );

      parameters.AddInitializer
      (
        "ensemble_representative_from_start_ensemble",
        "The model that should be used as the representative in the pymol script should come from the start "
        "ensemble. 1 = true; 0=false. If false representative will be taken from the end ensemble",
        io::Serialization::GetAgent( &m_EnsembleRepresentativeFromStartEnsemble),
        "0"
      );

      parameters.AddInitializer
      (
        "gradient_min",
        "The minimum value for the color gradient in pymol",
        io::Serialization::GetAgent( &m_GradientMin),
        "0"
      );

      parameters.AddInitializer
      (
        "gradient_max",
        "The maximum value for the color gradient in pymol",
        io::Serialization::GetAgent( &m_GradientMax),
        "0"
      );

      parameters.AddInitializer
      (
        "direct_relation",
        "Indicates the experimental data and exposures calculated are directly related. If true, a larger "
        "experimental value indicates the neighbor count should be larger (experiment measures buriedness). If"
        "false, a larger experimental value indicates the neigbor count should be smaller (experiment measures"
        "exposure).1=true;0=false",
        io::Serialization::GetAgent( &m_DirectRelation),
        "0"
      );

      parameters.AddInitializer
      (
        "set_nc_range",
        "Indicates that the range calculated from the NC's should not be determined automatically. Automatically"
        "the range is calculated as between the negative and postive value of the most extremem absolute value."
        "If this flag is set to true, the next parameter will be used as the range of NCs used to transform"
        "their values into the range of the experimental values. 1=true;0=false",
        io::Serialization::GetAgent( &m_SetNCRange),
        "0"
      );

      parameters.AddInitializer
      (
        "nc_range_min",
        "The minimimum value for the range of NCs (see parameter set_nc_range)",
        io::Serialization::GetAgent( &m_NCRangeMin),
        "0"
      );

      parameters.AddInitializer
      (
        "nc_range_max",
        "The maximum value for the range of NCs (see parameter set_nc_range)",
        io::Serialization::GetAgent( &m_NCRangeMax),
        "0"
      );

      return parameters;
    }

    //! @brief returns dummy name for ensemble
    //! @param ENSEMBLE for which a name will be created
    //! @return string which is the dummy name of ensemble
    std::string AnalyzeAccessibilityChange::EnsembleAsFilename
    (
      const assemble::ProteinEnsemble &ENSEMBLE
    )
    {
      return "start_ensemble_pdbs.ls";
    }

    //! @brief create ensemble from filename
    //! @param ENSEMBLE ensemble to setup
    //! @param NAME string name of file which will be used to create the ensemble
    //! @param ERR_STREAM stream to write out errors to
    //! @return ensemble created from the filename
    bool AnalyzeAccessibilityChange::EnsembleFromFilename( assemble::ProteinEnsemble &ENSEMBLE, const std::string &NAME, std::ostream &ERR_STREAM)
    {
      ENSEMBLE = assemble::ProteinEnsemble( NAME, 0, biol::GetAAClasses().e_AAComplete);
      return true;
    }

    //! @brief returns dummy name for exposure data file
    //! @return string which could be the name of the file the data comes from
    std::string AnalyzeAccessibilityChange::ExposureDataAsFilename
    (
      const storage::List< AccessibilityProfile> &DATA
    )
    {
      return "exposure_data.cst";
    }

    //! @brief reads in exposure data from a file given the filename
    //! @param PROFILES a list of accessibility profiles to set
    //! @param NAME the name of the file the data will be read from
    //! @param ERR_STREAM the stream any error will be written to
    //! @return true on success
    bool AnalyzeAccessibilityChange::ExposureDataFromFilename
    (
      storage::List< AccessibilityProfile> &PROFILES,
      const std::string &NAME,
      std::ostream &ERR_STREAM
    )
    {
      PROFILES.Reset();
      // open the data file
      io::IFStream read;
      io::File::MustOpenIFStream( read, NAME);

      // read in data file
      while( !read.eof() && read.peek() != std::istream::traits_type::eof())
      {
        // read in the chain and residue numbers
        char chain_a, chain_b;
        int aa_id_a, aa_id_b;
        io::Serialize::Read( chain_a, read);
        io::Serialize::Read( aa_id_a, read);
        io::Serialize::Read( chain_b, read);
        io::Serialize::Read( aa_id_b, read);

        // read in the exposure value
        double exposure;
        io::Serialize::Read( exposure, read);

        // create profile
        storage::List< AccessibilityAA> current_profile;
        // iterate through the residue range to make the profile
        for( int current_aa_id( aa_id_a); current_aa_id <= aa_id_b; ++current_aa_id)
        {
          // accessibility restraint
          AccessibilityAA restraint
          (
            storage::Map< AccessibilityAA::EnvironmentEnum, double>
            (
              storage::Map< AccessibilityAA::EnvironmentEnum, double>::Create
              (
                std::make_pair( AccessibilityAA::EnvironmentEnum( AccessibilityAA::e_Oxygen), exposure)
              )
            ),
            util::CloneToShPtr( assemble::LocatorAA( chain_a, current_aa_id, true)),
            util::ShPtr< assemble::AAExposureInterface>( new assemble::AANeighborCount())
          );

          // add current restraint
          current_profile.PushBack( restraint);
        }

        // create profile and add it to data
        const AccessibilityProfile profile( current_profile);
        PROFILES.PushBack( profile);

        // read end of line character
        char tmp;
        read.get( tmp);
      }

      BCL_MessageStd( "read in " + util::Format()( PROFILES.GetSize()) + " profiles");

      // return the data
      return true;
    }

    //! @brief determines the range of the calculated exposures and experimental exposures
    //! @param MEAN_DIFF_STATS to hold, over all profiles, the average of the average exposure difference calculated
    //!                        across all the residues in the profiles
    //! @return vector nd with two ranges, one for calculate the other for experimental exposures
    storage::VectorND< 2, math::Range< double> > AnalyzeAccessibilityChange::GetCalculatedAndExperimentalRanges
    (
      const assemble::ProteinEnsemble &ENSEMBLE, math::RunningAverageSD< double> &MEAN_DIFF_STATS
    ) const
    {
      // vectors to hold the values of calculated and experimental exposures
      storage::Vector< double> profile_mean;
      storage::Vector< double> profile_expval;

      // to hold, over all profiles, the average of the average exposure difference calculated across all the residues
      // in the profile
      MEAN_DIFF_STATS.Reset();

      // iterate through the data
      for
      (
        storage::List< AccessibilityProfile>::const_iterator profile_itr( m_ExperimentalExposures.Begin()),
        profile_itr_end( m_ExperimentalExposures.End());
        profile_itr != profile_itr_end;
        ++profile_itr
      )
      {
        // to hold, for a profile, the average exposure difference calculated across all the residues in the profile
        math::RunningAverageSD< double> average_profile_diff;

        // to hold the experimental data associated with this profile
        double exp_data;

        // iterate through the profile
        for
        (
          storage::List< AccessibilityAA>::const_iterator aa_itr( profile_itr->GetAccessibilities().Begin()),
          aa_itr_end( profile_itr->GetAccessibilities().End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // reference to current data
          const AccessibilityAA &current_data( *aa_itr);

          // iterate over the start ensemble
          for
          (
            assemble::ProteinEnsemble::const_iterator
              start_ensemble_itr( m_StartEnsemble.Begin()), start_ensemble_itr_end( m_StartEnsemble.End());
            start_ensemble_itr != start_ensemble_itr_end;
            ++start_ensemble_itr
          )
          {
            // the accessibility assignment for the starting ensemble
            const AccessibilityAAAssignment assignment_start( current_data.GenerateAssignment( **start_ensemble_itr));

            // iterate over the end ensemble
            for
            (
              assemble::ProteinEnsemble::const_iterator
               end_ensemble_itr( ENSEMBLE.Begin()), end_ensemble_itr_end( ENSEMBLE.End());
              end_ensemble_itr != end_ensemble_itr_end;
              ++end_ensemble_itr
            )
            {
              // the accessibility assignment for the ending ensemble
              const AccessibilityAAAssignment assignment_end( current_data.GenerateAssignment( **end_ensemble_itr));

              // right now the data is just stored as oxygen accessibility
              // this is reset with each iteration but all exp data is the same across the entire profile
              exp_data = assignment_end.GetAccessibilityByEnvironment
              (
                AccessibilityAA::e_Oxygen
              ).Second();

              // add in the current exposure difference to the statistics calculator
              average_profile_diff += assignment_end.GetExposureValue() - assignment_start.GetExposureValue();
            }
          } // start ensemble
        } // profile aas

        // get the mean and add it and the experimental data to their vectors
        const double mean( average_profile_diff.GetAverage());
        profile_mean.PushBack( mean);
        profile_expval.PushBack( exp_data);

        // add the mean change for the current profile into the overall statistics mean
        if( !m_DirectRelation)
        {
          MEAN_DIFF_STATS += -mean;
        }
        else
        {
          MEAN_DIFF_STATS += mean;
        }
      } // profiles

      // get range of means
      const storage::Vector< double> &mean_values( profile_mean);
      const double mean_min( math::Absolute( math::Statistics::MinimumValue( mean_values.Begin(), mean_values.End())));
      const double mean_max( math::Absolute( math::Statistics::MaximumValue( mean_values.Begin(), mean_values.End())));
      const double mean_range_value( mean_min > mean_max ? mean_min : mean_max);
      math::Range< double> mean_range( -mean_range_value, mean_range_value);
      if( m_SetNCRange)
      {
        mean_range = math::Range< double>( m_NCRangeMin, m_NCRangeMax);
      }

      // get range of exp_data
      const storage::Vector< double> &exp_values( profile_expval);
      const double exp_min( math::Absolute( math::Statistics::MinimumValue( exp_values.Begin(), exp_values.End())));
      const double exp_max( math::Absolute( math::Statistics::MaximumValue( exp_values.Begin(), exp_values.End())));
      const double exp_range_value( exp_min > exp_max ? exp_min : exp_max);
      const math::Range< double> exp_range( -exp_range_value, exp_range_value);
      BCL_MessageStd
      (
        "calculated accessibilities will be transformed from range\n"
        + util::Format()( mean_range) + "\ninto range\n" + util::Format()( exp_range)
      );

      // holds the calculated and experimental ranges
      storage::VectorND< 2, math::Range< double> > ranges( mean_range, exp_range);

      // return ranges
      return ranges;
    }

  } // namespace restraint
} // namespace bcl
