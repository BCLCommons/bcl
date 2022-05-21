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

#ifndef BCL_SIGNAL_FWD_HH_
#define BCL_SIGNAL_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the signal namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace signal
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class ConnectionInterface0;
    class Signal0;
    class SignalBase;
    class SignalBase0;
    class Slots;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_DestinationType>
    class Connection0;

    template< typename t_DestinationType, typename t_ArgType1>
    class Connection1;

    template< typename t_DestinationType, typename t_ArgType1, typename t_ArgType2>
    class Connection2;

    template< typename t_DestinationType, typename t_ArgType1, typename t_ArgType2, typename t_ArgType3>
    class Connection3;

    template< typename t_DestinationType, typename t_ArgType1, typename t_ArgType2, typename t_ArgType3, typename t_ArgType4>
    class Connection4;

    template< typename t_ArgType1>
    class ConnectionInterface1;

    template< typename t_ArgType1, typename t_ArgType2>
    class ConnectionInterface2;

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3>
    class ConnectionInterface3;

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3, typename t_ArgType4>
    class ConnectionInterface4;

    template< typename t_ArgType1>
    class Signal1;

    template< typename t_ArgType1, typename t_ArgType2>
    class Signal2;

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3>
    class Signal3;

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3, typename t_ArgType4>
    class Signal4;

    template< typename t_ArgType1>
    class SignalBase1;

    template< typename t_ArgType1, typename t_ArgType2>
    class SignalBase2;

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3>
    class SignalBase3;

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3, typename t_ArgType4>
    class SignalBase4;

  //////////////
  // typedefs //
  //////////////

  } // namespace signal
} // namespace bcl

#endif // BCL_SIGNAL_FWD_HH_
