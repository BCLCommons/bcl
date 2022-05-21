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

#ifndef BCL_APP_FWD_HH_
#define BCL_APP_FWD_HH_

// include the dependency file for this header
#include "bcl_app.depends.fwd.hh"

// This file contains forward declarations for the app namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace app
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class Apps;
    class GroupHandler;
    class Groups;
    class Interface;
    class InterfaceRelease;

  //////////////////////
  // template classes //
  //////////////////////

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< GroupHandler, Groups>          ApplicationGroup;
    typedef util::Enum< util::ShPtr< Interface>, Apps> ApplicationType;

  } // namespace app
} // namespace bcl

#endif // BCL_APP_FWD_HH_
