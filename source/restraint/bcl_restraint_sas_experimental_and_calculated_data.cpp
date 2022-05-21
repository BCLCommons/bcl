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
#include "restraint/bcl_restraint_sas_experimental_and_calculated_data.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math.h"
#include "restraint/bcl_restraint_sas_analysis.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasExperimentalAndCalculatedData::s_Instance
    (
      GetObjectInstances().AddInstance( new SasExperimentalAndCalculatedData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasExperimentalAndCalculatedData::SasExperimentalAndCalculatedData() :
      m_ExperimentalData(),
      m_CalculatedData()
    {
    }

    //! @brief constructor from given input data
    SasExperimentalAndCalculatedData::SasExperimentalAndCalculatedData
    (
      const SasScatteringData &EXPERIMENTAL,
      const SasScatteringData &CALCULATED
    ) :
      m_ExperimentalData( EXPERIMENTAL),
      m_CalculatedData( CALCULATED)
    {
      BCL_Assert
      (
        EXPERIMENTAL.GetScatteringData().GetSize() == CALCULATED.GetScatteringData().GetSize(),
        "Experimental and calculated SAXS data must have the same size"
      );
    }

    //! @brief Clone function
    //! @return pointer to new SasExperimentalAndCalculatedData
    SasExperimentalAndCalculatedData *SasExperimentalAndCalculatedData::Clone() const
    {
      return new SasExperimentalAndCalculatedData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasExperimentalAndCalculatedData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief scale the calculated and experimental data
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    void SasExperimentalAndCalculatedData::ScaleData( const double &SCALING_FACTOR)
    {
      // Scale Experimental and Calculated Data
      m_ExperimentalData = SasAnalysis::ScaleData( m_ExperimentalData, SCALING_FACTOR);
      m_CalculatedData   = SasAnalysis::ScaleData( m_CalculatedData, SCALING_FACTOR);
    }

    //! @brief scale the calculated intensity to align the experimental and calculated curves
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    void SasExperimentalAndCalculatedData::ScaleCalculatedData( const double &SCALING_FACTOR)
    {
      m_CalculatedData = SasAnalysis::ScaleData( m_CalculatedData, SCALING_FACTOR);
    }

    //! @brief Set the Experimental Error Values of current object to the values of passed object
    //! @param DATA_OBJECT - the object with the error values to copy
    void SasExperimentalAndCalculatedData::SetExperimentalData( const SasScatteringData &DATA_OBJECT)
    {
      m_ExperimentalData = DATA_OBJECT;
    }

    //! @brief move the data to all positive numbers while maintaining identical morphology
    void SasExperimentalAndCalculatedData::SlideData( const double &Y_MIN)
    {
      m_ExperimentalData = SasAnalysis::SlideData( m_ExperimentalData, Y_MIN);
      m_CalculatedData = SasAnalysis::SlideData( m_CalculatedData, Y_MIN);
    }

    //! @brief Set the scale for the experimental and calculated intensities to the specified boundary
    void SasExperimentalAndCalculatedData::SetYScale( const double &Y_MAX)
    {

      double max_intensity( SasAnalysis::MaxIntensity( m_ExperimentalData));

      double min_intensity
      (
        std::min
        (
          SasAnalysis::MinIntensity( m_ExperimentalData), SasAnalysis::MinIntensity( m_CalculatedData)
        )
      );

      if( min_intensity < 0)
      {
        SlideData( min_intensity);
        max_intensity = SasAnalysis::MaxIntensity( m_ExperimentalData);
      }
      const double scale_factor( Y_MAX / max_intensity);
      ScaleData( scale_factor);
    }

    //! @brief Normalize the experimental data set by its largest intensity and then Normalize the calculated data
    //! @brief by its largest intensity value
    void SasExperimentalAndCalculatedData::NormalizeData()
    {
      // Normalize the 2 data sets
      m_ExperimentalData = SasAnalysis::NormalizeData( m_ExperimentalData);
      m_CalculatedData =   SasAnalysis::NormalizeData( m_CalculatedData);
    }

    //! @brief take the data to log base 10
    void SasExperimentalAndCalculatedData::Log10()
    {
      m_ExperimentalData = SasAnalysis::Log10( m_ExperimentalData);
      m_CalculatedData   = SasAnalysis::Log10( m_CalculatedData);
    }

    //! @brief transform the log10 data to the absolute scale
    void SasExperimentalAndCalculatedData::LogtoAbsolute()
    {
      m_ExperimentalData = SasAnalysis::LogtoAbsolute( m_ExperimentalData);
      m_CalculatedData = SasAnalysis::LogtoAbsolute( m_CalculatedData);
    }

    //! @brief take the derivative of the data
    void SasExperimentalAndCalculatedData::Derivative()
    {
      m_ExperimentalData = SasAnalysis::Derivative( m_ExperimentalData);
      m_CalculatedData =   SasAnalysis::Derivative( m_CalculatedData);
    }

    //! @brief compute the chi score between the experimental and calculated curves
    //! @return the chi score between the experimental and calculated curves
    double SasExperimentalAndCalculatedData::ComputeScoringFunction( bool USE_ERRORS) const
    {
      double sum( 0.0);
      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( m_ExperimentalData.Begin()),
          cal_data_itr( m_CalculatedData.Begin()),
          exp_data_itr_end( m_ExperimentalData.End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        if( USE_ERRORS)
        {
          sum += math::Sqr( ( exp_data_itr->GetIntensity() - cal_data_itr->GetIntensity()) / exp_data_itr->GetError());
        }
        else
        {
          sum += math::Sqr( exp_data_itr->GetIntensity() - cal_data_itr->GetIntensity());
        }
      }

      return math::Sqrt( sum / double( GetExperimentalData().GetScatteringData().GetSize()));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &SasExperimentalAndCalculatedData::WriteToGnuplot( std::ostream &OSTREAM) const
    {
      OSTREAM << "Q_Value Experimental_Intensity Experimental_Error Computed_Intensity" << '\n';

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( GetExperimentalData().Begin()),
          cal_data_itr( GetCalculatedData().Begin()),
          exp_data_itr_end( GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // write the data to the ofstream
        OSTREAM << exp_data_itr->GetQvalue()    << ' '
                << exp_data_itr->GetIntensity() << ' '
                << exp_data_itr->GetError()     << ' '
                << cal_data_itr->GetIntensity() << '\n';
      }
      return OSTREAM;
    }

    //! @brief write to filename in three columns : Q-value, calculated intensity, experimental error
    //! @param filename to write to
    std::ostream &SasExperimentalAndCalculatedData::WriteToGnomeFormat( std::ostream &OSTREAM) const
    {
      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( GetExperimentalData().Begin()),
          cal_data_itr( GetCalculatedData().Begin()),
          exp_data_itr_end( GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // write the data to the ofstream
        OSTREAM << cal_data_itr->GetQvalue()    << ' '
                << cal_data_itr->GetIntensity() << ' '
                << exp_data_itr->GetError()     << '\n';
      }
      return OSTREAM;
    }

    //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param filename to write to
    void SasExperimentalAndCalculatedData::WriteToGnomeFileName( const std::string &FILENAME) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToGnomeFormat( write);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param filename to write to
    void SasExperimentalAndCalculatedData::WriteToGnuplotFileName( const std::string &FILENAME) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToGnuplot( write);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief read in the member data from a pre-formatted file
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedData::ReadFromDataFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( line.GetSize());

      // ignore empty lines at the first of the file
      while( line_size < 1)
      {
        std::getline( ISTREAM, read_line);
        line = util::SplitString( read_line, " ");
        line_size = line.GetSize();
      }

      bool is_correct_file( false);

      // verify bcl file type
      if(
          line( 0) == "Q_Value" &&
          line( 1) == "Experimental_Intensity" &&
          line( 2) == "Experimental_Error" &&
          line( 3) == "Computed_Intensity"
        )
      {
        is_correct_file = true;
      }

      BCL_Assert( is_correct_file == true, "Incorrect Input File Type" + util::Format()( line));

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double exp_intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::ConvertStringToNumericalValue< double>( split_line( 2)));
        double cal_intensity( util::ConvertStringToNumericalValue< double>( split_line( 3)));

        m_ExperimentalData.PushBackScattering( SasScatteringPoint( q, exp_intensity, error));
        m_CalculatedData.PushBackScattering( SasScatteringPoint( q, cal_intensity, 0.0));
      }

      size_t experimental_size( m_ExperimentalData.GetScatteringSize());
      size_t calculated_size( m_CalculatedData.GetScatteringSize());

      BCL_Assert( experimental_size != 0 && experimental_size == calculated_size, " Data was incorrectly read");

      return ISTREAM;
    }

    //! @brief read fit file format from CRYSOL
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedData::ReadFromCrysolFitFile
    ( std::istream &ISTREAM, const double &FIRST_EXPERIMENTAL_POINT)
    {

      // build a string to hold the line information
      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( line.GetSize());

      // ignore empty lines at the first of the file
      while( line_size < 1)
      {
        std::getline( ISTREAM, read_line);
        line = util::SplitString( read_line, " ");
        line_size = line.GetSize();
      }

      // Switch to control when to start reading data
      bool use_data( false);

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double exp_intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::GetUndefined< double>());
        double cal_intensity( util::ConvertStringToNumericalValue< double>( split_line( 2)));

        // Do not use extrapolated data from Crysol. Only use the data provided by experiment
        // Use within tolerance function to account for numerical drift between platforms
        if( math::EqualWithinAbsoluteTolerance( q, FIRST_EXPERIMENTAL_POINT, 0.001))
        {
          use_data = true;
        }

        if( use_data)
        {
          m_ExperimentalData.PushBackScattering( SasScatteringPoint( q, exp_intensity, error));
          m_CalculatedData.PushBackScattering( SasScatteringPoint( q, cal_intensity, error));
        }
      }

      size_t experimental_size( m_ExperimentalData.GetScatteringSize());
      size_t calculated_size( m_CalculatedData.GetScatteringSize());

      BCL_Assert( experimental_size != 0 && experimental_size == calculated_size, " Data was incorrectly read");
      return ISTREAM;
    }

    //! @brief read fit file format from CRYSOL
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedData::ReadFromFoxsFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( line.GetSize());

      // ignore empty lines at the first of the file
      while( line_size < 1)
      {
        std::getline( ISTREAM, read_line);
        line = util::SplitString( read_line, " ");
        line_size = line.GetSize();
      }

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        std::string first_value( split_line( 0));

        // Do not use header information from foxs. Only use the data provided by experiment
        if( first_value.compare( "#") == 0)
        {
          continue;
        }

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double exp_intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::GetUndefined< double>());
        double cal_intensity( util::ConvertStringToNumericalValue< double>( split_line( 2)));

        m_ExperimentalData.PushBackScattering( SasScatteringPoint( q, exp_intensity, error));
        m_CalculatedData.PushBackScattering( SasScatteringPoint( q, cal_intensity, error));
      }

      size_t experimental_size( m_ExperimentalData.GetScatteringSize());
      size_t calculated_size( m_CalculatedData.GetScatteringSize());

      BCL_Assert( experimental_size != 0 && experimental_size == calculated_size, " Data was incorrectly read");

      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExperimentalData, ISTREAM);
      io::Serialize::Read( m_CalculatedData, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasExperimentalAndCalculatedData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExperimentalData, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CalculatedData, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
