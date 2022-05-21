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

#ifndef BCL_COORD_MOVE_COMBINE_H_
#define BCL_COORD_MOVE_COMBINE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_move_interface.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoveCombine
    //! @brief Class that combines a list of moves and applies them consecutively
    //! @details Initialized with a list of moves to apply, its operator applies all of the moves consecutively
    //!
    //! @see @link example_coord_move_combine.cpp @endlink
    //! @author karakam
    //! @date Feb 22, 2011
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoveCombine :
      public MoveInterface
    {

    private:

    //////////
    // data //
    //////////

      //! list of moves
      storage::Vector< util::Implementation< MoveInterface> > m_Moves;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoveCombine();

      //! @brief constructor from a list of moves
      //! @param MOVES_LIST List of moves
      MoveCombine( const util::ShPtrList< MoveInterface> &MOVES_LIST);

      //! @brief Clone function
      //! @return pointer to new MoveCombine
      MoveCombine *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the move list
      //! @return move list
      const util::ShPtrList< MoveInterface> &GetMoves() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief move the MovableInterface derived object
      //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
      //! @return reference to movable object
      MovableInterface &Move( MovableInterface &MOVEABLE_OBJECT) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MoveCombine

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_MOVE_COMBINE_H_ 
