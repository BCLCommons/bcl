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

#ifndef BCL_UTIL_RUNTIME_ENVIRONMENT_DEFAULT_H_
#define BCL_UTIL_RUNTIME_ENVIRONMENT_DEFAULT_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_runtime_environment_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RuntimeEnvironmentDefault
    //! @brief represents a default runtime environment.
    //! @details This means a standard desktop or server running Linux, Windows or Mac OS.
    //! This is opposed to say BOINC or MPI environments.
    //!
    //! @see RuntimeEnvironments
    //! @see @link example_util_runtime_environment_default.cpp @endlink
    //! @author whitebc
    //! @date 07.17.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RuntimeEnvironmentDefault :
      public RuntimeEnvironmentInterface
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RuntimeEnvironmentDefault();

      //! @brief Clone function
      //! @return new Pointer to a copy of the actual object behind the pointer
      RuntimeEnvironmentDefault *Clone() const
      {
        return new RuntimeEnvironmentDefault( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @return the class name if this function is overwritten
      const std::string &GetClassIdentifier() const;

      //! @brief is a stand alone application
      //! @return true is the application is stand alone, false if it runs in some kind of environment
      bool IsStandAlone() const
      {
        return true;
      }

      //! @brief access to the process number or slot number
      //! @return the process or slot, this app is running in
      int GetProcessNumber() const;

    ////////////////
    // operations //
    ////////////////

      //! @copydoc RuntimeEnvironmentInterface::ResolveFileName()
      std::string ResolveFileName( const std::string &FILE_NAME) const;

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
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief Initializes the environment
      //! @return A boolean representing whether the initialization was a success or failure.
      bool Initialize() const;

      //! @brief Finalizes the environment
      //! @param STATUS indicates if an error occured
      //! @return A boolean representing whether finalize was a success or failure.
      bool Finalize( const int STATUS) const;

    }; // class RuntimeEnvironmentDefault

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_RUNTIME_ENVIRONMENT_DEFAULT_H_
