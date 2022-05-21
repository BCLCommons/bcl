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

#ifndef BCL_UTIL_FORMAT_H_
#define BCL_UTIL_FORMAT_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <sstream>

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Format
    //! @brief This class allows simple access to the rather complex formatted output functions of C++.
    //! @details The input value is converted into a formated output string. the Format object is
    //! defined in a way that allows the usage in a cout command, e.g. \n
    //!
    //! bcl::Format().F(7,3)(12); \n
    //!
    //! but also its value-seperated definition so that it can be reused and even passed as function
    //! argument, e.g: \n
    //!
    //! bcl::Format fmt; \n
    //! fmt.F(7,3); \n
    //! cout << fmt( 12); \n
    //! Write( fmt); \n
    //!
    //! @see @link example_util_format.cpp @endlink
    //! @author meilerj
    //! @remarks reviewed on Dec. 04, 2010 woetzen, alexanns
    //! @date Jul 31, 2004
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Format
    {

    private:

    //////////
    // data //
    //////////

      size_t                  m_Precision;     //!< precision of floating point and scientific formats
      size_t                  m_Width;         //!< width of output for all formats
      bool                    m_SingleSpace;   //!< flags whether output is preceded by a single space
      bool                    m_ForceWidth;    //!< flags whether width should be enforced
      char                    m_Fill;          //!< character used for filling
      std::ios_base::fmtflags m_Format;        //!< format flags to be applied
      bool                    m_Initialized;   //!< indicates if format object was initialized

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Format();

      //! copy function
      Format *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! returns class name
      //! the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the width
      //! @return the width
      size_t GetWidth() const
      {
        return m_Width;
      }

    ///////////////
    // operators //
    ///////////////

      //! output operator(), gets value T1
      template< typename t_DataType>
      std::string operator ()( const t_DataType &VALUE) const
      {
        std::stringstream stream; // make a stream

        SetupStream( stream); // set it up with any flags that are set

        stream << VALUE; // write out the value

        // post-process the string in the stream and return it
        return Postprocess( stream.str());
      }

      //! output operator() for float type (checks for nans and writes them accordingly)
      std::string operator ()( const float &VALUE) const;

      //! output operator() for double type (checks for nans and writes them accordingly)
      std::string operator ()( const double &VALUE) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief turn SingleSpace on and off
      //! @return Format reference to this
      Format &S();

      //! @brief turn fixed width on and off
      //! @return Format reference to this
      Format &ForceW();

      //! @brief set left aligned output
      //! @return Format reference to this
      Format &L();

      //! @brief set right aligned output
      //! @return Format reference to this
      Format &R();

      //! @brief set width
      //! @param WIDTH width of formatted numerical value
      //! @return Format reference to this
      Format &W( const size_t WIDTH);

      //! @brief set scientific format, gets width and precision as input
      //! @param PREC the precision - number of digits after '.'
      //! @return Format reference to this
      Format &SFP( const size_t PREC);

      //! @brief set floating point precision, gets precision as input
      //! @param PREC the precision - number of digits after '.'
      //! @return Format reference to this
      Format &FFP( const size_t PREC);

      //! @brief set filling character - will be used to fill unused width
      //! @param FILL the character to fill with - default ' '
      //! @return Format reference to this
      Format &Fill( const char FILL = ' ');

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief apply flags on stream
      //! @param IOSTREAM the stream to set the flags on
      //! @return std::ios reference to the IOSTREAM
      std::ios &SetFlags( std::ios &IOSTREAM) const;

      //! @brief remove flags on stream
      //! @param IOSTREAM the stream to disable the flags on
      //! @return std::ios reference to the IOSTREAM
      std::ios &UnsetFlags( std::ios &IOSTREAM) const;

    //////////////////////
    // input and output //
    //////////////////////

      // Read/Write must not be protected, otherwise there is no way to write util::Format

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief setup the stream with the given flags and add the single space, if necessary
      //! @param STREAM the stream to set properties on
      void SetupStream( std::stringstream &STREAM) const;

      //! @brief truncates the string if force width is set
      //! @details will truncate the string if width was forced and replaces new line char with the desired header
      //! @param STRING the string to process
      //! @return the string in the stream after it was post-processed
      std::string Postprocess( std::string STRING) const;

    }; // class Format

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_FORMAT_H_
