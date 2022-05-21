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
#include "restraint/bcl_restraint_sas_scattering_data.h"

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
    const util::SiPtr< const util::ObjectInterface> SasScatteringData::s_Instance
    (
      util::Enumerated< HandlerBase< SasScatteringData> >::AddInstance
      (
        new SasScatteringData()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasScatteringData::SasScatteringData( const std::string &EXTENSION) :
      HandlerBase< SasScatteringData>( EXTENSION),
      m_Data()
    {
    }

    //! @brief constructor from given input data
    SasScatteringData::SasScatteringData
    (
      const storage::Vector< SasScatteringPoint> &INIT_DATA
    ) :
      HandlerBase< SasScatteringData>( ".saxs"),
      m_Data( INIT_DATA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasScatteringData
    SasScatteringData *SasScatteringData::Clone() const
    {
      return new SasScatteringData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasScatteringData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &SasScatteringData::GetAlias() const
    {
      // format is automatically determined
      static const std::string s_name( "Detect");
      return s_name;
    }

    //! @brief pushback function to add object to Dataset vector
    //! @param DATAPOINT_OBJECT DataPoint values Q, I, and Error
    void SasScatteringData::PushBackScattering( const SasScatteringPoint &DATAPOINT_OBJECT)
    {
      m_Data.PushBack( DATAPOINT_OBJECT);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief preallocate Scattering Data memory
    //! @param SIZE size to preallocate
    void SasScatteringData::AllocateScatteringMemory( const size_t &SIZE)
    {
      m_Data.AllocateMemory( SIZE);
    }

    //! @brief bool test for error values
    //! @return true if error is defined for all values of the dataset
    const bool SasScatteringData::IsErrorDefined() const
    {
      // iterate over m_Data
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator itr( m_Data.Begin()), itr_end( m_Data.End());
         itr != itr_end; ++itr
      )
      {
        // if error is defined do not terminate
        if( !itr->IsErrorDefined())
        {
          // if error is not defined return false immediately
          return false;
        }
      }

      // error is defined for all cases, return true
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read experimental data from BCL
    //! @param ISTREAM input data stream
    void SasScatteringData::ReadBCLProcessedData( std::istream &ISTREAM)
    {
      // String to hold line data
      std::string read_line;

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::ConvertStringToNumericalValue< double>( split_line( 2)));

        // push back the data
        m_Data.PushBack( SasScatteringPoint( q, intensity, error));
      } // close while loop
    } // close function

    //! @brief read experimental data from Crysol. See Crysol documentation for column definitions
    //! @param ISTREAM input data stream
    void SasScatteringData::ReadCrysolData( std::istream &ISTREAM)
    {
      // String to hold line data
      std::string read_line;

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double q( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));

        m_Data.PushBack( SasScatteringPoint( q, intensity, SasAnalysis::ComputeError( q, intensity)));
      } // close while loop

    } // close function

    //! @brief read experimental data
    //! @param ISTREAM input data stream
    void SasScatteringData::ReadExperimentalData( std::istream &ISTREAM)
    {
      //String to hold line data
      std::string read_line;
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double intensity( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        if( intensity > 0)
        {
          // push back the data
          m_Data.PushBack
          (
            SasScatteringPoint
            (
              util::ConvertStringToNumericalValue< double>( split_line( 0)),    // Q-value
              util::ConvertStringToNumericalValue< double>( split_line( 1)),    // I-value for intensity in solution
              util::ConvertStringToNumericalValue< double>( split_line( 2))     // error
            )
          );
        }
        else
        {
           double value( util::ConvertStringToNumericalValue< double>( split_line( 1)));
           BCL_MessageStd( "Improper intensity value: " + util::Format()( value));
           BCL_Exit( "Intensity values cannot be zero or negative please examine input experimental data: ", -1);
        }
      } // close while loop
    } // close function

    //! @brief read experimental data from Gnom
    //! @param ISTREAM input data stream
    void SasScatteringData::ReadGnomData( std::istream &ISTREAM, size_t &INTENSITY_COLUMN, bool USE_EXTRAPOLATION)
    {

      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> Line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( Line.GetSize());

      // paramaters to read through gnom filetype
      bool readFlag( false);

      // while the end of the file is not reached
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));
        if
        (
          split_line.GetSize() >= 7 &&
          split_line( 0) == "S"     &&
          split_line( 1) == "J"     &&
          split_line( 2) == "EXP"   &&
          split_line( 3) == "ERROR" &&
          split_line( 4) == "J"     &&
          split_line( 5) == "REG"   &&
          split_line( 6) == "I"
        )
        {
          readFlag = true;
        }
        else if( readFlag)
        {
          // get size of the vector
          line_size = split_line.GetSize();

          if( line_size > 1 && !util::IsNumerical( split_line( 0)))
          {
            readFlag = false;
          }
          else if( line_size == 2 && USE_EXTRAPOLATION)
          {
            // push back the extrapolated data
             m_Data.PushBack
             (
               SasScatteringPoint
               (
                 util::ConvertStringToNumericalValue< double>( split_line( 0)),  // Q-value
                 util::ConvertStringToNumericalValue< double>( split_line( 1)),  // I-value
                 util::GetUndefined< double>()                                   // error
               )
             );
          }
          else if( line_size == 5)
          {
            // push back the data
            m_Data.PushBack
            (
              SasScatteringPoint
              (
                util::ConvertStringToNumericalValue< double>( split_line( 0)),                 // Q-value
                util::ConvertStringToNumericalValue< double>( split_line( INTENSITY_COLUMN)),  // I-value
                util::ConvertStringToNumericalValue< double>( split_line( 2))                  // error
              )
            );
          }
        }
      } // close while loop
    } // close the function

    //! @brief reads in the member data from a formatted file containing 3 columns:
    //! @brief scattering angle q ( 4*Pi*sin(theta)/lambda) where lambda is measured in Angstroms, I is the intensity
    //! @brief at a given a value, E is the experimental error.
    //! @brief Crysol generated files must used the paramater /dro 0.0.  The algorithm does not support adding the
    //! @brief hydration layer around the molecule.
    //! @param ISTREAM input stream
    //! @param FORMAT the file format to use for reading
    //! @return istream which was read from
    std::istream &SasScatteringData::ReadFromDataFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      enum Types { bcl_processed_data, crysol, experimental, gnom, unknown};
      Types filetype( unknown);

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> Line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( Line.GetSize());

      // ignore empty lines at the first of the file
      while( line_size <= 1)
      {
        std::getline( ISTREAM, read_line);
        Line = util::SplitString( read_line, " ");
        line_size = Line.GetSize();
      }

      // verify bclfile type
      if(
          line_size == 7 &&
          Line( 0) == "BCL" &&
          Line( 1) == "SCATTERING" &&
          Line( 2) == "INPUT" &&
          Line( 3) == "PARAMETERS:" &&
          Line( 4) == "QValue" &&
          Line( 5) == "Intensity" &&
          Line( 6) == "Error"
        )
      {
        BCL_MessageStd( "bcl_processed_data");
        filetype = bcl_processed_data;
      }

      else if(
          line_size == 4 &&
          Line( 0) == "Q_Value" &&
          Line( 1) == "Experimental_Intensity" &&
          Line( 2) == "Experimental_Error" &&
          Line( 3) == "Computed_Intensity"
        )
      {
        BCL_MessageStd( "bcl_processed_data");
        filetype = bcl_processed_data;
      }

      // verify crysol file type
      else if( Line( 0) == "Dif/Atom/Shape/Bord")
      {
        filetype = crysol;
      }

      // Verify the file is a GNOM file
      else if( line_size >= 5 && Line( 1) == "G" && Line( 2) == "N" && Line( 3) == "O" && Line( 4) == "M")
      {
        filetype = gnom;
      }

      // experimental data file type
      else if( line_size == 3)
      {
        filetype = experimental;
        double intensity( util::ConvertStringToNumericalValue< double>( Line( 1)));
        if( intensity > 0)
        {
          // push back the data
          m_Data.PushBack
          (
            SasScatteringPoint
            (
               util::ConvertStringToNumericalValue< double>( Line( 0)),    // Q-value
               util::ConvertStringToNumericalValue< double>( Line( 1)),    // I-value for intensity in solution
               util::ConvertStringToNumericalValue< double>( Line( 2))     // error
            )
          );
        }
        else
        {
          BCL_MessageStd( "Skipped negative or zero intensity value");
        }

      }
      else
      {
        BCL_Exit( "Could not read file; last read line was: " + read_line, -1);
      }

      size_t experimental_data_column( 1);
      switch( filetype)
      {
        case bcl_processed_data:
          ReadBCLProcessedData( ISTREAM);
          break;
        case crysol:
          ReadCrysolData( ISTREAM);
          break;
        case experimental:
          ReadExperimentalData( ISTREAM);
          break;
        case gnom:
          ReadGnomData( ISTREAM, experimental_data_column);
          break;
        case unknown:
        default:
          break;
      }

      size_t size( 0);

      // make sure there are q values for analysis
      size = m_Data.GetSize();

      BCL_Assert( size != 0, "The number of Q values are: " + util::Format()( size));

      // return the stream
      return ISTREAM;
    }

    //! @brief reads saxs restraints from an input stream
    //! @brief input stream to read the restraints from
    //! @return the read in restraints
    SasScatteringData SasScatteringData::ReadRestraints( std::istream &ISTREAM) const
    {
      SasScatteringData data;
      data.ReadFromDataFile( ISTREAM);
      return data;
    }

    //! @brief reads the computed SAXS profile from the fit density curve from gnom
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasScatteringData::ReadFitFromGnom( std::istream &ISTREAM)
    {
      size_t fit_data_column( 4);
      bool use_extrapolation( true);
      ReadGnomData( ISTREAM, fit_data_column, use_extrapolation);

      size_t size( 0);
      size = m_Data.GetSize();
      BCL_Assert( size != 0, "The number of Q values are: " + util::Format()( size));

      // return the stream
      return ISTREAM;
    }

    //! @brief writes out the member data from a formatted file containing 3 columns
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &SasScatteringData::WriteToDataFile( std::ostream &OSTREAM, bool HEADER) const
    {
      // initialize format
      const util::Format format;

      if( HEADER)
      {
        OSTREAM << format( "BCL SCATTERING INPUT PARAMETERS: QValue Intensity Error") << '\n';
      }

      // iterate over m_Data
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator itr( m_Data.Begin()), itr_end( m_Data.End());
         itr != itr_end; ++itr
      )
      {
        // write the data
        OSTREAM << format( itr->GetQvalue()) << ' '
                << format( itr->GetIntensity()) << ' '
                << format( itr->GetError()) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param filename to write to
    void SasScatteringData::WriteToDataFileName( const std::string &FILENAME, const bool &HEADER) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToDataFile( write, HEADER);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasScatteringData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasScatteringData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
