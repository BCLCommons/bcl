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

#ifndef BCL_COORD_MOVE_TRANSLATE_DEFINED_H_
#define BCL_COORD_MOVE_TRANSLATE_DEFINED_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_move_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoveTranslateDefined
    //! @brief translates a Moveable object by a defined translation vector.
    //!
    //! @see @link example_coord_move_translate_defined.cpp @endlink
    //! @author alexanns, woetzen
    //! @date March 16, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoveTranslateDefined :
      public MoveInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_Translation;    //!< distance for translation in each direction (x, y, z)
      //! flag - if it is set to true, will move objects center to origin before applying translation, and move it back
      bool           m_TranslateInternal;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoveTranslateDefined();

      //! @brief construct from  translation and internal coordinates or not
      //! @param TRANSLATION maximal distance for translation in each direction
      //! @param TRANSLATE_INTERNAL
      MoveTranslateDefined( const linal::Vector3D &TRANSLATION, const bool TRANSLATE_INTERNAL = true);

      //! @brief Clone function
      //! @return pointer to new Move
      MoveTranslateDefined *Clone() const;

      //! @brief destructor
      ~MoveTranslateDefined();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

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
      //! @param MOVEABLE_OBJECT refernce on MovableInterface derived object
      //! @return reference to movable object
      MovableInterface &Move( MovableInterface &MOVEABLE_OBJECT) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MoveTranslateDefined

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_MOVE_TRANSLATE_DEFINED_H_
