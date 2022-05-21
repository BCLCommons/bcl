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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "coord/bcl_coord_move_combine.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MoveCombine::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveCombine())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveCombine::MoveCombine() :
      m_Moves()
    {
    }

    //! @brief constructor from a list of moves
    //! @param MOVES_LIST List of moves
    MoveCombine::MoveCombine( const util::ShPtrList< MoveInterface> &MOVES_LIST) :
      m_Moves()
    {
      for( auto move_it( MOVES_LIST.Begin()); move_it != MOVES_LIST.End(); ++move_it)
      {
        m_Moves.PushBack( **move_it);
      }
    }

    //! @brief Clone function
    //! @return pointer to new MoveCombine
    MoveCombine *MoveCombine::Clone() const
    {
      return new MoveCombine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoveCombine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveCombine::GetAlias() const
    {
      static const std::string s_name( "MoveCombine");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveCombine::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Combines a list of moves and applies them consecutively.");
      serializer.AddInitializer
      (
        "moves",
        "list of moves to be applied",
        io::Serialization::GetAgent( &m_Moves)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveCombine::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      // iterate over the moves
      for( auto move_itr( m_Moves.Begin()), move_itr_end( m_Moves.End()); move_itr != move_itr_end; ++move_itr)
      {
        // apply the move
        ( *move_itr)->Move( MOVEABLE_OBJECT);
      }

      // end
      return MOVEABLE_OBJECT;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace coord
} // namespace bcl
