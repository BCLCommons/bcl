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

#ifndef BCL_COMMAND_FWD_HH_
#define BCL_COMMAND_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the command namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace command
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AppDefaultFlags;
    class Command;
    class CommandLineWriter;
    class CommandState;
    class FlagDynamic;
    class FlagInterface;
    class FlagSeparator;
    class FlagStatic;
    class FlagStaticAndDynamic;
    class Guesser;
    class Parameter;
    class ParameterCheckAllowed;
    class ParameterCheckAllowedNonConst;
    class ParameterCheckDefault;
    class ParameterCheckExtension;
    class ParameterCheckExtensionsFileExistence;
    class ParameterCheckFileExistence;
    class ParameterCheckFileInSearchPath;
    class ParameterCheckInterface;
    class ParameterCheckOr;
    class ParameterCheckSerializable;
    class ParameterInterface;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_Enumerate>
    class ParameterCheckEnumerate;

    template< typename t_DataType>
    class ParameterCheckRanged;

  //////////////
  // typedefs //
  //////////////

  } // namespace command
} // namespace bcl

#endif // BCL_COMMAND_FWD_HH_
