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

#ifndef BCL_SIGNAL_CONNECTION_H_
#define BCL_SIGNAL_CONNECTION_H_

// include the namespace header
#include "bcl_signal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_signal_connection_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace signal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Connection0
    //! @brief represents a connection between a Signal0 and a destination
    //! @details The connection holds the ptr to the destination and and the pointer to the void member function taking 0 arguments
    //! the Signal can call the emit through the connection to invoke the member function on the destination class
    //!
    //! @tparam t_DestinationType type of Slots derived class and class the member function will be called of
    //!
    //! @see @link example_signal_connection.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DestinationType>
    class Connection0 :
      public ConnectionInterface0
    {
    private:

    //////////
    // data //
    //////////

      //! pointer to destination class
      t_DestinationType *m_PtrObject;

      //! pointer to member function that will be invoked on destination
      void ( t_DestinationType::*m_PtrMemberFunction)();

      //! bool, that indicates, if a connection can be inherited by copies of signals, or should only exist between the
      //! Signal and Destination
      bool m_CopyConnection;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Connection0() :
        m_PtrObject( NULL),
        m_PtrMemberFunction( NULL),
        m_CopyConnection( false)
      {
      }

      //! @brief construct from pointer to destination and ptr to member function
      //! @param DESTINATION ptr to the object this connection is pointing to
      //! @param MEMBER_FUNCTION_PTR the member function, this connection will call on Emit
      //! @param COPY_CONNECTION when true, the connection will be copied to copies of Signals or Slots, when false, they will be not inherited
      Connection0( t_DestinationType *DESTINATION, void ( t_DestinationType::*MEMBER_FUNCTION_PTR)(), const bool COPY_CONNECTION = false) :
        m_PtrObject( DESTINATION),
        m_PtrMemberFunction( MEMBER_FUNCTION_PTR),
        m_CopyConnection( COPY_CONNECTION)
      {
      }

      //! @brief copy
      //! @return ptr to a new Connection0, which is a copy of this
      ConnectionInterface0 *Clone()
      {
        return new Connection0< t_DestinationType>( *this);
      }

      //! @brief Duplicate function, that takes a Pointer to a different destination
      //! @param NEW_DESTINATION pointer to Slots derived class
      //! @return a new connection to the NEW_DESTINATION
      ConnectionInterface0 *Duplicate( Slots *NEW_DESTINATION)
      {
        return new Connection0< t_DestinationType>( ( t_DestinationType *)NEW_DESTINATION, m_PtrMemberFunction);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief access to the destination
      //! @return pointer to the destination
      Slots *GetDestination() const
      {
        return m_PtrObject;
      }

      //! @brief can this connection be copied
      //! @return true, if connection can be copied, false otherwise
      bool CanBeCopied() const
      {
        return m_CopyConnection;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Emit signal to the actual destination
      //! call the member function on the destination object
      void Emit()
      {
        ( m_PtrObject->*m_PtrMemberFunction)();
      }

    }; // template class Connection0

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Connection1
    //! @brief represents a connection between a Signal1 and a destination
    //! @details The connection holds the ptr to the destination and and the pointer to the void member function taking 1 argument
    //! the Signal can call the emit through the connection to invoke the member function on the destination class
    //!
    //! @tparam t_DestinationType type of Slots derived class and class the member function will be called of
    //! @tparam t_ArgType1 type of first argument of member function
    //!
    //! @see Connection0
    //! @see @link example_signal_connection.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DestinationType, typename t_ArgType1>
    class Connection1 :
      public ConnectionInterface1< t_ArgType1>
    {

    private:

      t_DestinationType *m_PtrObject;
      void ( t_DestinationType::*m_PtrMemberFunction)( t_ArgType1);
      bool m_CopyConnection;

    public:

      Connection1() :
        m_PtrObject( NULL),
        m_PtrMemberFunction( NULL),
        m_CopyConnection( false)
      {
      }

      Connection1( t_DestinationType *DESTINATION, void ( t_DestinationType::*MEMBER_FUNCTION_PTR)( t_ArgType1), const bool COPY_CONNECTION = false) :
        m_PtrObject( DESTINATION),
        m_PtrMemberFunction( MEMBER_FUNCTION_PTR),
        m_CopyConnection( COPY_CONNECTION)
      {
      }

      ConnectionInterface1< t_ArgType1> *Clone()
      {
        return new Connection1< t_DestinationType, t_ArgType1>( *this);
      }

      ConnectionInterface1< t_ArgType1> *Duplicate( Slots *NEW_DESTINATION)
      {
        return new Connection1< t_DestinationType, t_ArgType1>( ( t_DestinationType *)NEW_DESTINATION, m_PtrMemberFunction);
      }

      Slots *GetDestination() const
      {
        return m_PtrObject;
      }

      bool CanBeCopied() const
      {
        return m_CopyConnection;
      }

      void Emit( t_ArgType1 ARGUMENT1)
      {
        ( m_PtrObject->*m_PtrMemberFunction)( ARGUMENT1);
      }

    }; // template class Connection1

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Connection2
    //! @brief represents a connection between a Signal2 and a destination
    //! @details The connection holds the ptr to the destination and and the pointer to the void member function taking 2 arguments
    //! the Signal can call the emit through the connection to invoke the member function on the destination class
    //!
    //! @tparam t_DestinationType type of Slots derived class and class the member function will be called of
    //! @tparam t_ArgType1 type of first argument of member function
    //! @tparam t_ArgType2 type of second argument of member function
    //!
    //! @see Connection0
    //! @see @link example_signal_connection.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DestinationType, typename t_ArgType1, typename t_ArgType2>
    class Connection2 :
      public ConnectionInterface2< t_ArgType1, t_ArgType2>
    {

    private:

      t_DestinationType *m_PtrObject;
      void ( t_DestinationType::*m_PtrMemberFunction)( t_ArgType1, t_ArgType2);
      bool m_CopyConnection;

    public:

      Connection2() :
        m_PtrObject( NULL),
        m_PtrMemberFunction( NULL),
        m_CopyConnection( false)
      {
      }

      Connection2( t_DestinationType *DESTINATION, void ( t_DestinationType::*MEMBER_FUNCTION_PTR)( t_ArgType1, t_ArgType2), const bool COPY_CONNECTION = false) :
        m_PtrObject( DESTINATION),
        m_PtrMemberFunction( MEMBER_FUNCTION_PTR),
        m_CopyConnection( COPY_CONNECTION)
      {
      }

      ConnectionInterface2< t_ArgType1, t_ArgType2> *Clone()
      {
        return new Connection2< t_DestinationType, t_ArgType1, t_ArgType2>( *this);
      }

      ConnectionInterface2< t_ArgType1, t_ArgType2> *Duplicate( Slots *NEW_DESTINATION)
      {
        return new Connection2< t_DestinationType, t_ArgType1, t_ArgType2>( ( t_DestinationType *)NEW_DESTINATION, m_PtrMemberFunction);
      }

      Slots *GetDestination() const
      {
        return m_PtrObject;
      }

      bool CanBeCopied() const
      {
        return m_CopyConnection;
      }

      void Emit( t_ArgType1 ARGUMENT1, t_ArgType2 ARGUMENT2)
      {
        ( m_PtrObject->*m_PtrMemberFunction)( ARGUMENT1, ARGUMENT2);
      }

    }; // template class Connection2

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Connection3
    //! @brief represents a connection between a Signal3 and a destination
    //! @details The connection holds the ptr to the destination and and the pointer to the void member function taking 3 arguments
    //! the Signal can call the emit through the connection to invoke the member function on the destination class
    //!
    //! @tparam t_DestinationType type of Slots derived class and class the member function will be called of
    //! @tparam t_ArgType1 type of first argument of member function
    //! @tparam t_ArgType2 type of second argument of member function
    //! @tparam t_ArgType3 type of third argument of member function
    //!
    //! @see Connection0
    //! @see @link example_signal_connection.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DestinationType, typename t_ArgType1, typename t_ArgType2, typename t_ArgType3>
    class Connection3 :
      public ConnectionInterface3< t_ArgType1, t_ArgType2, t_ArgType3>
    {
    private:

      t_DestinationType *m_PtrObject;
      void ( t_DestinationType::*m_PtrMemberFunction)( t_ArgType1, t_ArgType2, t_ArgType3);
      bool m_CopyConnection;

    public:

      Connection3() :
        m_PtrObject( NULL),
        m_PtrMemberFunction( NULL),
        m_CopyConnection( false)
      {
      }

      Connection3( t_DestinationType *DESTINATION, void ( t_DestinationType::*MEMBER_FUNCTION_PTR)( t_ArgType1, t_ArgType2, t_ArgType3), const bool COPY_CONNECTION = false) :
        m_PtrObject( DESTINATION),
        m_PtrMemberFunction( MEMBER_FUNCTION_PTR),
        m_CopyConnection( COPY_CONNECTION)
      {
      }

      ConnectionInterface3< t_ArgType1, t_ArgType2, t_ArgType3> *Clone()
      {
        return new Connection3< t_DestinationType, t_ArgType1, t_ArgType2, t_ArgType3>( *this);
      }

      ConnectionInterface3< t_ArgType1, t_ArgType2, t_ArgType3> *Duplicate( Slots *NEW_DESTINATION)
      {
        return new Connection3< t_DestinationType, t_ArgType1, t_ArgType2, t_ArgType3>( ( t_DestinationType *)NEW_DESTINATION, m_PtrMemberFunction);
      }

      Slots *GetDestination() const
      {
        return m_PtrObject;
      }

      bool CanBeCopied() const
      {
        return m_CopyConnection;
      }

      void Emit( t_ArgType1 ARGUMENT1, t_ArgType2 ARGUMENT2, t_ArgType3 ARGUMENT3)
      {
        ( m_PtrObject->*m_PtrMemberFunction)( ARGUMENT1, ARGUMENT2, ARGUMENT3);
      }

    }; // template class Connection3

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Connection4
    //! @brief represents a connection between a Signal4 and a destination
    //! @details The connection holds the ptr to the destination and and the pointer to the void member function taking 4 arguments
    //! the Signal can call the emit through the connection to invoke the member function on the destination class
    //!
    //! @tparam t_DestinationType type of Slots derived class and class the member function will be called of
    //! @tparam t_ArgType1 type of first argument of member function
    //! @tparam t_ArgType2 type of second argument of member function
    //! @tparam t_ArgType3 type of third argument of member function
    //! @tparam t_ArgType4 type of fourth argument of member function
    //!
    //! @see Connection0
    //! @see @link example_signal_connection.cpp @endlink
    //! @author woetzen
    //! @date 05/10/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DestinationType, typename t_ArgType1, typename t_ArgType2, typename t_ArgType3, typename t_ArgType4>
    class Connection4 :
      public ConnectionInterface4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>
    {
    private:

      t_DestinationType *m_PtrObject;
      void ( t_DestinationType::*m_PtrMemberFunction)( t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4);
      bool m_CopyConnection;

    public:

      Connection4() :
        m_PtrObject( NULL),
        m_PtrMemberFunction( NULL),
        m_CopyConnection( false)
      {
      }

      Connection4( t_DestinationType *DESTINATION, void ( t_DestinationType::*MEMBER_FUNCTION_PTR)( t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4), const bool COPY_CONNECTION = false) :
        m_PtrObject( DESTINATION),
        m_PtrMemberFunction( MEMBER_FUNCTION_PTR),
        m_CopyConnection( COPY_CONNECTION)
      {
      }

      ConnectionInterface4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4> *Clone()
      {
        return new Connection4< t_DestinationType, t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>( *this);
      }

      ConnectionInterface4< t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4> *Duplicate( Slots *NEW_DESTINATION)
      {
        return new Connection4< t_DestinationType, t_ArgType1, t_ArgType2, t_ArgType3, t_ArgType4>( ( t_DestinationType *)NEW_DESTINATION, m_PtrMemberFunction);
      }

      Slots *GetDestination() const
      {
        return m_PtrObject;
      }

      bool CanBeCopied() const
      {
        return m_CopyConnection;
      }

      void Emit( t_ArgType1 ARGUMENT1, t_ArgType2 ARGUMENT2, t_ArgType3 ARGUMENT3, t_ArgType4 ARGUMENT4)
      {
        ( m_PtrObject->*m_PtrMemberFunction)( ARGUMENT1, ARGUMENT2, ARGUMENT3, ARGUMENT4);
      }

    }; // template class Connection4

  } // namespace signal
} // namespace bcl

#endif // BCL_SIGNAL_CONNECTION_H_
