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

#ifndef BCL_UTIL_RUNTIME_ENVIRONMENT_INTERFACE_H_
#define BCL_UTIL_RUNTIME_ENVIRONMENT_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RuntimeEnvironmentInterface
    //! @brief is the bcl environment interface
    //! @details An environment interface defines an environment for applications to do I/O and other operations
    //! in a runtime environment independent way.
    //!
    //! @see RuntimeEnvironments
    //! @see @link example_util_runtime_environment_interface.cpp @endlink
    //!
    //! @author woetzen, karakam, whitebc
    //! @date 04.02.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RuntimeEnvironmentInterface :
      public ObjectInterface
    {
    /////////////
    // friends //
    /////////////

      friend class RuntimeEnvironments; //!< allow RuntimeEnvironments to call protected functions

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RuntimeEnvironmentInterface()
      {
      }

      //! @brief Clone function
      //! @return new Pointer to a copy of the actual object behind the pointer
      virtual RuntimeEnvironmentInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief is a stand alone application
      //! @return true is the application is stand alone, false if it runs in some kind of environment
      virtual bool IsStandAlone() const = 0;

      //! @brief access to the process number or slot number
      //! @return the process or slot, this app is running in
      virtual int GetProcessNumber() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief Resolves the provided file name to an absolute path name
      //! @param FILE_NAME file name to resolve
      //! @return string with the resolved filename, empty if it could not be resolved
      virtual std::string ResolveFileName( const std::string &FILE_NAME) const = 0;

      //! @brief Get the OS file path separator
      //! @return the path separator as a std::string
      const std::string &GetPathSeperator() const;

      //! @brief report progress
      //! @param FRACTION_PROGRESS progress between [0,1]
      //! @return true if successful reported
      virtual bool ReportProgress( const double FRACTION_PROGRESS) const;

      //! @brief add name to the process number
      //! @param NAME name to be added to process number
      //! @return NAME + processnumber as string
      std::string AddProcessNumberToName( const std::string &NAME) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Initializes the environment
      //! @return A boolean representing whether the initialization was a success or failure.
      virtual bool Initialize() const = 0;

      //! @brief Finalizes the environment
      //! @param STATUS indicates if an error occured
      //! @return A boolean representing whether finalize was a success or failure.
      virtual bool Finalize( const int STATUS) const = 0;

    }; // class RuntimeEnvironmentInterface

    // Forward declaration of function GetRuntimeEnvironment
    BCL_API const RuntimeEnvironmentInterface &GetRuntimeEnvironment();

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_RUNTIME_ENVIRONMENT_INTERFACE_H_
