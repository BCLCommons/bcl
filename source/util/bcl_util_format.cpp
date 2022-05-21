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
#include "util/bcl_util_format.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Format::Format() :
      m_Precision( GetUndefined< size_t>()),
      m_Width( GetUndefined< size_t>()),
      m_SingleSpace( false),
      m_ForceWidth( false),
      m_Fill( ' '),
      m_Format(),
      m_Initialized( false)
    {
    }

    //! copy function
    Format *Format::Clone() const
    {
      return new Format( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &Format::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief turn SingleSpace on and off
    //! @return Format reference to this
    Format &Format::S()
    {
      m_SingleSpace = !m_SingleSpace;
      return *this;
    }

    //! @brief turn fixed width on and off
    //! @return Format reference to this
    Format &Format::ForceW()
    {
      m_ForceWidth = !m_ForceWidth;
      return *this;
    }

    //! @brief set left aligned output
    //! @return Format reference to this
    Format &Format::L()
    {
      m_Format |= std::ios::left;
      m_Initialized = true;
      return *this;
    }

    //! @brief set right aligned output
    //! @return Format reference to this
    Format &Format::R()
    {
      m_Format |= std::ios::right;
      m_Initialized = true;
      return *this;
    }

    //! @brief set width
    //! @param WIDTH width of formatted numerical value
    //! @return Format reference to this
    Format &Format::W( const size_t WIDTH)
    {
      m_Width = WIDTH;
      return *this;
    }

    //! @brief set scientific format, gets width and precision as input
    //! @param PREC the precision - number of digits after '.'
    //! @return Format reference to this
    Format &Format::SFP( const size_t PREC)
    {
      m_Precision = PREC;
      m_Format |= std::ios::scientific;
      m_Initialized = true;
      return *this;
    }

    //! @brief set floating point precision, gets precision as input
    //! @param PREC the precision - number of digits after '.'
    //! @return Format reference to this
    Format &Format::FFP( const size_t PREC)
    {
      m_Precision = PREC;
      m_Format |= std::ios::fixed;
      m_Initialized = true;
      return *this;
    }

    //! @brief set filling character - will be used to fill unused width
    //! @param FILL the character to fill with - default ' '
    //! @return Format reference to this
    Format &Format::Fill( const char FILL)
    {
      m_Fill = FILL;
      if( FILL == '0')
      {
        m_Format |= std::ios_base::internal;
      }
      m_Initialized = true;
      return *this;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief apply flags on stream
    //! @param IOSTREAM the stream to set the flags on
    //! @return std::ios reference to the IOSTREAM
    std::ios &Format::SetFlags( std::ios &IOSTREAM) const
    {
      if( IsDefined( m_Precision))
      {
        IOSTREAM.precision( m_Precision);
      }

      if( IsDefined( m_Width))
      {
        IOSTREAM.width( m_Width);
      }

      if( m_Fill != ' ')
      {
        IOSTREAM.fill( m_Fill);
      }

      if( m_Initialized)
      {
        IOSTREAM.setf( m_Format);
      }

      return IOSTREAM;
    }

    //! @brief remove flags on stream
    //! @param IOSTREAM the stream to disable the flags on
    //! @return std::ios reference to the IOSTREAM
    std::ios &Format::UnsetFlags( std::ios &IOSTREAM) const
    {
      // reset all formats
      IOSTREAM.unsetf( std::ios::floatfield);

      // end
      return IOSTREAM;
    }

    //! output operator() for float type (checks for nans and writes them accordingly)
    std::string Format::operator ()( const float &VALUE) const
    {
      std::stringstream stream; // make a stream

      SetupStream( stream); // set it up with any flags that are set

      if( !IsDefined( VALUE))
      {
        static std::string s_nan( "nan");
        // write nan, since the actual undefined value is machine-dependent
        stream << s_nan; // write out the value
      }
      else
      {
        stream << VALUE;
      }

      // post-process the string in the stream and return it
      return Postprocess( stream.str());
    }

    //! output operator() for double type (checks for nans and writes them accordingly)
    std::string Format::operator ()( const double &VALUE) const
    {
      std::stringstream stream; // make a stream

      SetupStream( stream); // set it up with any flags that are set

      if( !IsDefined( VALUE))
      {
        static std::string s_nan( "nan");
        // write nan, since the actual undefined value is machine-dependent
        stream << s_nan; // write out the value
      }
      else
      {
        stream << VALUE;
      }

      // post-process the string in the stream and return it
      return Postprocess( stream.str());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Format::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Precision, ISTREAM);
      io::Serialize::Read( m_Width, ISTREAM);
      io::Serialize::Read( m_SingleSpace, ISTREAM);
      io::Serialize::Read( m_ForceWidth, ISTREAM);
      io::Serialize::Read( m_Fill, ISTREAM);
      io::Serialize::Read( m_Initialized, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &Format::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Precision,   OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_Width,       OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_SingleSpace, OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_ForceWidth,  OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_Fill,        OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_Initialized, OSTREAM, INDENT) << ' ';

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief setup the stream with the given flags and add the single space, if necessary
    //! @param STREAM the stream to set properties on
    void Format::SetupStream( std::stringstream &STREAM) const
    {
      // output stream
      if( m_SingleSpace)
      {
        STREAM << ' ';
      }

      SetFlags( STREAM);
    }

    //! @brief truncates the string if force width is set
    //! @details will truncate the string if width was forced
    //! @param STRING the string to process
    //! @return the string in the stream after it was post-processed
    std::string Format::Postprocess( std::string STRING) const
    {
      // fix length
      if( m_ForceWidth)
      {
        // resize the stream to width
        STRING.resize( m_Width + size_t( m_SingleSpace), m_Fill);
      }

      return STRING;
    }

  } // namespace util
} // namespace bcl
