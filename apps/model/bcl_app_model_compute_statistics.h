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

#ifndef BCL_APP_MODEL_COMPUTE_STATISTICS_H_
#define BCL_APP_MODEL_COMPUTE_STATISTICS_H_

// include the namespace header

// include forward headers
#include "io/bcl_io.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "model/bcl_model.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "command/bcl_command_flag_interface.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ModelComputeStatistics
    //! @brief Application for computing statistics rmsd, correlation, r^2 for predicted outputs and
    //!        ROC curve related measures. If different model evaluations for the same dataset are provided
    //!        a jury evaluation can be conducted.
    //!
    //! @see @link example_app_model_compute_statistics.cpp @endlink
    //! @author loweew, butkiem1
    //! @date 10/28/2011, 12/01/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ModelComputeStatistics :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! input file names
      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

      //! allows limiting the number of file combinations returned by considering only combinations with N elements
      util::ShPtr< command::FlagInterface> m_MaxFilesPerConsensusFlag;

      //! table name
      util::ShPtr< command::FlagInterface> m_TableNameFlag;

      //! sort by
      util::ShPtr< command::FlagInterface> m_SortByFlag;

      //! control which measures will be plotted
      //! Mutable to allow setting based off the objective_function flag
      mutable util::ShPtr< command::FlagInterface> m_PlotXFlag;
      mutable util::ShPtr< command::FlagInterface> m_PlotLogXFlag;
      mutable util::ShPtr< command::FlagInterface> m_PlotYFlag;

      //! Flag to prevent automatically deciding what to plot based on other options
      util::ShPtr< command::FlagInterface> m_NoPlotFlag;

      //! take log
      util::ShPtr< command::FlagInterface> m_TakeLog10Flag;

      //! correlation plots desired
      //! Mutable to allow setting based off the objective_function flag
      mutable util::ShPtr< command::FlagInterface> m_CorrelationFlag;

      //! potency cutoff when plotting roc curves
      //! Mutable to allow setting based off the objective_function flag
      mutable util::ShPtr< command::FlagInterface> m_PotencyCutOffFlag;

      //! indicates actives are those predicted below the potency cutoff
      //! Mutable to allow setting based off the objective_function flag
      mutable util::ShPtr< command::FlagInterface> m_ActivesBelowCutoffFlag;

      //! output file names
      mutable util::ShPtr< command::FlagInterface> m_OutputDirectory;

      //! input objective function labels which will be evaluated
      util::ShPtr< command::FlagInterface> m_ObjectiveFunction;

      //! filename to store evaluated objective function values
      util::ShPtr< command::FlagInterface> m_EvaluateObjectiveFunctionFilename;
      
      //! the image format to output
      util::ShPtr< command::FlagInterface> m_ImageFormatFlag;

      //! Objective functions
      mutable storage::Vector< util::Implementation< model::ObjectiveFunctionInterface> > m_ObjectiveFunctions;

      //! Test whether gnuplot is available
      mutable bool m_GnuplotAvailable;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      ModelComputeStatistics();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      ModelComputeStatistics *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const;

    private:

      //! @brief create the raw data output to create a file that can be parsed by gnuplot
      //! @param TABLE table of interest that will be filled up with results data
      //! @param PRED predicted data
      //! @param EXP experimental data
      //! @param NAMES base name to create appropriate output filenames
      void CreateDesiredPlot
      (
        storage::Table< float> &TABLE,
        const linal::Vector< float> &PRED,
        const linal::Vector< float> &EXP,
        const std::string &NAMES
      ) const;

      //! @brief create the raw data output to create a file that can be parsed by gnuplot
      //! @param TABLE table of interest that will be filled up with results data
      //! @param PRED predicted data matrix
      //! @param EXP experimental data matrix
      //! @param NAMES base name to create appropriate output filenames
      void CreateDesiredPlots
      (
        storage::Table< float> &TABLE,
        const linal::Matrix< float> &PRED,
        const linal::Matrix< float> &EXP,
        const std::string &NAMES
      ) const;

      //! @brief generate the table header for the final results table
      //! @param CORRELATION flag whether a correlation plot will be considered
      //! @return table header
      storage::TableHeader GetTableHeader( bool CORRELATION) const;

      //! @brief gets the string for the terminal type to use for gnuplot graphs
      //! @return a string with the gnuplot terminal type
      std::string GetGnuplotTerminalType() const;

      //! @brief get the extension for the given image type
      //! @return the three-letter extension for the image output type
      std::string GetImageFormatExtension() const;
        
      //! @brief write the correlation plot header for a gnuplot file to output stream
      //! @param OSTREAM output stream of interest
      //! @param OUTPUT_FILENAME output file name OSTREAM will write to
      void WriteCorrelationHeader( io::OFStream &OSTREAM, std::string &OUTPUT_FILENAME) const;

      //! @brief write the correlation plot legend for a gnuplot file to output stream
      //! @param OSTREAM output stream of interest
      //! @param RMSD rmsd value that is written on the plot
      //! @param RSQAURED r^2 value that is written on the plot
      void WriteCorrelationLegend( io::OFStream &OSTREAM, double RMSD, double RSQAURED) const;

      //! @brief write the correlation plot data to a filestream that is readable by gnuplot
      //! @param OSTREAM output stream of interest
      //! @param EXP experimental data
      //! @param PRED predicted data
      void WriteFormattedCorrelationInput
      (
        io::OFStream &OSTREAM,
        const linal::Vector< float> &EXP,
        const linal::Vector< float> &PRED
      ) const;

      //! @brief write ROC curve plot header for a gnuplot file to output stream
      //! @param OSTREAM output stream of interest
      //! @param OUTPUT_FILENAME output file name OSTREAM will write to
      //! @param ENRICHMENT enrichment value that is written on the plot
      //! @param AUC area under the curve value that is written on the plot
      void WriteROCHeader
      (
        io::OFStream &OSTREAM,
        const std::string &OUTPUT_FILENAME,
        const float ENRICHMENT, const float AUC
      ) const;

      //! @brief write ROC curve plot header for a gnuplot file with all contingency matrix measures to output stream
      //! @param OSTREAM output stream of interest
      //! @param OUTPUT_FILENAME output file name OSTREAM will write to
      //! @param MEASURES a list of all measures that were outuput, in the order they were output
      void WriteROCCompleteHeader
      (
        io::OFStream &OSTREAM,
        const std::string &OUTPUT_FILENAME,
        const storage::Vector< std::string> &MEASURES
      ) const;

      //! @brief evaluate objective function values for given obj function labels per commandline flag
      //! @param PRED predicted data
      //! @param EXP experimental data
      //! @param NAME base name for writing out
      //! @param IDS any ids associated with the predictions
      void EvaluateObjectiveFunctions
      (
        const linal::Matrix< float> &PRED,
        const model::FeatureDataSet< float> &EXP,
        const std::string &NAME
      ) const;

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType ModelComputeStatistics_Instance;

    }; // ModelComputeStatistics

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MODEL_COMPUTE_STATISTICS_H_
