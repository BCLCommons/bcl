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

#ifndef BCL_COORD_MOVE_INTERFACE_H_
#define BCL_COORD_MOVE_INTERFACE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoveInterface
    //! @brief the interface implements a Move function, that takes a reference to a MovabelInterface and moves it
    //! @details MoveInterface and Movable interface provide a framework to move objects in 3D around
    //!
    //! @remarks example unnecessary
    //! @remarks reviewed Nov 6, 2010
    //! @author woetzen
    //! @date Nov 16, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MoveInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to a new MoveInterface
      virtual MoveInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief move the MovableInterface derived object
      //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
      //! @return reference to movable object
      virtual MovableInterface &Move( MovableInterface &MOVEABLE_OBJECT) const = 0;

    }; //class MoveInterface

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_MOVE_INTERFACE_H_
