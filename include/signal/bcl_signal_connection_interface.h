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

#ifndef BCL_SIGNAL_CONNECTION_INTERFACE_H_
#define BCL_SIGNAL_CONNECTION_INTERFACE_H_

// include the namespace header
#include "bcl_signal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace signal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConnectionInterface0
    //! @brief represents a ConnectionInterface0 between a Signal0 as Interface
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConnectionInterface0
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtaul destuctor
      virtual ~ConnectionInterface0()
      {
      }

      //! @brief virtual copy constructor
      virtual ConnectionInterface0 *Clone() = 0;

      //! @brief duplicate connection using the new destination slot
      virtual ConnectionInterface0 *Duplicate( Slots *NEW_DESTINATION) = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief access to the destination
      //! @return access to the Slots interface of the destination
      virtual Slots *GetDestination() const = 0;

      //! @brief can this connection be copied
      //! @return true, if connection can be copied, false otherwise
      virtual bool CanBeCopied() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief emit signal through the connection to the destination by calling its connected member function
      virtual void Emit() = 0;

    }; // template class ConnectionInterface0

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConnectionInterface1
    //! @brief represents a ConnectionInterface1 between a Signal1 as Interface
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1>
    class ConnectionInterface1
    {
    public:
      virtual ~ConnectionInterface1< t_ArgType1>()
      {
      }
      virtual ConnectionInterface1< t_ArgType1> *Clone() = 0;
      virtual ConnectionInterface1< t_ArgType1> *Duplicate( Slots *NEW_DESTINATION) = 0;
      virtual Slots *GetDestination() const = 0;
      virtual bool CanBeCopied() const = 0;
      virtual void Emit( t_ArgType1) = 0;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConnectionInterface2
    //! @brief represents a ConnectionInterface2 between a Signal2 as Interface
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2>
    class ConnectionInterface2
    {
    public:
      virtual ~ConnectionInterface2< t_ArgType1, t_ArgType2>()
      {
      }
      virtual ConnectionInterface2< t_ArgType1, t_ArgType2>* Clone() = 0;
      virtual ConnectionInterface2< t_ArgType1, t_ArgType2>* Duplicate( Slots* NEW_DESTINATION) = 0;
      virtual Slots *GetDestination() const = 0;
      virtual bool CanBeCopied() const = 0;
      virtual void Emit( t_ArgType1, t_ArgType2) = 0;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConnectionInterface3
    //! @brief represents a ConnectionInterface3 between a Signal3 as Interface
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3>
    class ConnectionInterface3
    {
    public:
      virtual ~ConnectionInterface3< t_ArgType1, t_ArgType2, t_ArgType3>()
      {
      }
      virtual ConnectionInterface3< t_ArgType1, t_ArgType2, t_ArgType3>* Clone() = 0;
      virtual ConnectionInterface3< t_ArgType1, t_ArgType2, t_ArgType3>* Duplicate( Slots* NEW_DESTINATION) = 0;
      virtual Slots *GetDestination() const = 0;
      virtual bool CanBeCopied() const = 0;
      virtual void Emit( t_ArgType1, t_ArgType2, t_ArgType3) = 0;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConnectionInterface4
    //! @brief represents a ConnectionInterface4 between a Signal4 as Interface
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgType1, typename t_ArgType2, typename t_ArgType3, typename t_ArgType4>
    class ConnectionInterface4
    {
    public:
      virtual ~ConnectionInterface4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>()
      {
      }
      virtual ConnectionInterface4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>* Clone() = 0;
      virtual ConnectionInterface4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>* Duplicate( Slots* NEW_DESTINATION) = 0;
      virtual Slots *GetDestination() const = 0;
      virtual bool CanBeCopied() const = 0;
      virtual void Emit( t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4) = 0;
    };

  } // namespace signal
} // namespace bcl

#endif // BCL_SIGNAL_CONNECTION_INTERFACE_H_
