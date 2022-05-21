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

#ifndef BCL_IO_VALIDATION_RESULT_H_
#define BCL_IO_VALIDATION_RESULT_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

namespace bcl
{
  namespace io
  {
    //! Result of validation
    enum ValidationResultType
    {
      e_Allowed,  //!< Parameter value is valid
      e_Help,     //!< Parameter value is a request for help
      e_Invalid,  //!< Parameter value is invalid
      e_Complete  //!< Parameter value is valid and has already been read in completely.
                  //! e_Complete can be used for objects that handle their own setting internally
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ValidationResult
    //! @brief enum wrapper indicating the result of validating, setting, or reading something that may be a help request
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Mar 07, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ValidationResult
    {
    public:

    //////////
    // data //
    //////////

    private:

      ValidationResultType m_Type; //!< Actual type

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! constructor from bool; appropriate when the parameter is never help
      explicit ValidationResult( const bool &VALID);

      //! constructor from Type
      ValidationResult( const ValidationResultType &TYPE);

    /////////////////
    // data access //
    /////////////////

      //! allow implicit conversion to bool
      operator bool() const
      {
        return m_Type == e_Allowed || m_Type == e_Complete;
      }

      //! @brief Test whether the type is allowed
      //! @return true if type == e_Allowed
      bool IsAllowed() const
      {
        return m_Type == e_Allowed;
      }

      //! @brief Test whether the type is invalid
      //! @return true if type == e_Invalid
      bool IsInvalid() const
      {
        return m_Type == e_Invalid;
      }

      //! @brief Test whether type is help
      //! @return true if type == e_Help
      bool IsHelp() const
      {
        return m_Type == e_Help;
      }

      //! @brief Test whether type is complete
      //! @return true if type == e_Complete
      bool IsComplete() const
      {
        return m_Type == e_Complete;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief equality operator
      bool operator ==( const ValidationResult &RESULT) const
      {
        return RESULT.m_Type == m_Type;
      }

      //! @brief inequality operator
      bool operator !=( const ValidationResult &RESULT) const
      {
        return RESULT.m_Type != m_Type;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get the string the user can use to indicate help is desired for a particular parameter
      //! @return help
      static const std::string &GetHelpString();

    };
  } // namespace io
} // namespace bcl

#endif // BCL_IO_VALIDATION_RESULT_H_
