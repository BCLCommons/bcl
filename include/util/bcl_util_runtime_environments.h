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

#ifndef BCL_UTIL_RUNTIME_ENVIRONMENTS_H_
#define BCL_UTIL_RUNTIME_ENVIRONMENTS_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_enumerate.h"
#include "bcl_util_runtime_environment_interface.h"
#include "bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    //! @brief for RuntimeEnvironment
    typedef Enum< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments> RuntimeEnvironment;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RuntimeEnvironments
    //! @brief manages the system-wide list of enumerated runtime environments
    //!
    //! @see @link example_util_runtime_environments.cpp @endlink
    //! @author whitebc, woetzen
    //! @date 07.18.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RuntimeEnvironments :
      public Enumerate< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments>
    {
    public:

      friend class Enumerate< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments>;

    //////////
    // data //
    //////////

      RuntimeEnvironment e_Default; //!< default environment

      //! the current runtime environment
      ShPtr< RuntimeEnvironmentInterface> m_CurrentRuntimeEnvironment;

      //! @brief command line flag to be used to set RuntimeEnvironment over the command line
      //! @return ShPtr to a FlagInterface which is used to set RuntimeEnvironment
      static ShPtr< command::FlagInterface> &GetFlagRuntimeEnvironment();

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RuntimeEnvironments();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Initialize the runtime environment from the command line flag
      static void UpdateCurrentRuntimeEnvironmentFromCommandLineFlag();

      //! @brief access to the current environment
      //! @return one and only reference to one of the environments
      const RuntimeEnvironmentInterface &GetCurrentEnvironment() const;

    }; // class RuntimeEnvironments

    //! @brief returns RuntimeEnvironments
    //! @return RuntimeEnvironments
    BCL_API RuntimeEnvironments &GetRuntimeEnvironments();

    //! @brief returns active RuntimeEnvironment
    //! @return active RuntimeEnvironment
    BCL_API const RuntimeEnvironmentInterface &GetRuntimeEnvironment();

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments>;

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_RUNTIME_ENVIRONMENTS_H_
