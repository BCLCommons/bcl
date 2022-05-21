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

#ifndef BCL_COORD_MOVE_TRANSLATE_EXTERNAL_AXIS_H_
#define BCL_COORD_MOVE_TRANSLATE_EXTERNAL_AXIS_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_axes.h"
#include "bcl_coord_move_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoveTranslateExternalAxis
    //! @brief translates an object towards or away from an external axis
    //! @details translates an object towards or away from an external axis. This is used to move protein subunits
    //!          towards or away from the axis of symmetry in a multimer
    //!
    //! @see @link example_coord_move_translate_external_axis.cpp @endlink
    //! @author weinerbe
    //! @date Oct 7, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MoveTranslateExternalAxis :
      public MoveInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the minimum translation
      double m_MinTranslation;

      //! the maximum translation
      double m_MaxTranslation;

      //! external axis to translate towards or away from
      Axis m_ExternalAxis;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoveTranslateExternalAxis();

      //! @brief construct from translation range and axis
      //! @param MIN_TRANSLATION minimal distance for translation
      //! @param MAX_TRANSLATION maximal distance for translation
      //! @param EXTERNAL_AXIS axis to translate towards or away from
      MoveTranslateExternalAxis
      (
        const double MIN_TRANSLATION,
        const double MAX_TRANSLATION,
        const Axis &EXTERNAL_AXIS
      );

      //! @brief Clone function
      //! @return pointer to new MoveTranslateExternalAxis
      MoveTranslateExternalAxis *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
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
      //! @return reference to movable object
      MovableInterface &Move( MovableInterface &MOVEABLE_OBJECT) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MoveTranslateExternalAxis

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_MOVE_TRANSLATE_EXTERNAL_AXIS_H_ 
