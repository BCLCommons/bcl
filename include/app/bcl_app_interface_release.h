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

#ifndef BCL_APP_INTERFACE_RELEASE_H_
#define BCL_APP_INTERFACE_RELEASE_H_

// include the namespace header
#include "bcl_app.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_app_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class InterfaceRelease
    //! @brief interface for all release applications
    //! @remarks every application that will be released should derive from this interface
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Mar 5, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API InterfaceRelease :
      public Interface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Interface
      virtual InterfaceRelease *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      //! @brief returns readme information
      //! @return string containing information about application
      virtual const std::string &GetReadMe() const = 0;

      //! @brief return flag whether application is release application. By default the return value is false.
      //!        Release applications have to overwrite the method to return true and indicate the release status.
      //! @return true if application is release application, false otherwise
      bool IsReleaseApplication() const
      {
        // relase applications have to overwrite this method to change return value!
        return true;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      virtual std::string GetDescription() const = 0;

    }; // class InterfaceRelease

  } // namespace app
} // namespace bcl

#endif // BCL_APP_INTERFACE_RELEASE_H_
