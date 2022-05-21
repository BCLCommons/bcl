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

#ifndef BCL_COORD_MOVE_TRANSLATE_RANDOM_H_
#define BCL_COORD_MOVE_TRANSLATE_RANDOM_H_

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
    //! @class MoveTranslateRandom
    //! @brief TODO: add an general comment to this class
    //!
    //! @see @link example_coord_move_translate_random.cpp @endlink
    //! @author woetzen
    //! @date Feb 1, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoveTranslateRandom :
      public MoveInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_MinTranslation;    //!< minimal distance for translation in each direction (x, y, z)
      linal::Vector3D m_MaxTranslation;    //!< maximal distance for translation in each direction (x, y, z)
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
      MoveTranslateRandom();

      //! @brief construct from undirected maximal translation and rotation
      //! @param MAX_TRANSLATION maximal distance for translation in each direction (x, y, z)
      //! @param TRANSLATE_INTERNAL        maximal distance for translation
      MoveTranslateRandom( const double MAX_TRANSLATION, const bool TRANSLATE_INTERNAL = true);

      //! @brief construct from undirected maximal translation and rotation
      //! @param MAX_TRANSLATIONS maximal distance for translation in each direction (x, y, z)
      //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
      MoveTranslateRandom( const linal::Vector3D &MAX_TRANSLATIONS, const bool TRANSLATE_INTERNAL = true);

      //! @brief construct from undirected maximal translation and rotation
      //! @param MIN_TRANSLATION minimal distance for translation in each direction (x, y, z)
      //! @param MAX_TRANSLATION maximal distance for translation in each direction (x, y, z)
      //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
      MoveTranslateRandom
      (
        const double MIN_TRANSLATION,
        const double MAX_TRANSLATION,
        const bool TRANSLATE_INTERNAL = true
      );

      //! @brief construct from undirected minimal and maximal translations
      //! @param MIN_TRANSLATIONS minimal distances for translation in each direction (x, y, z)
      //! @param MAX_TRANSLATIONS maximal distances for translation in each direction (x, y, z)
      //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
      MoveTranslateRandom
      (
        const linal::Vector3D &MIN_TRANSLATIONS,
        const linal::Vector3D &MAX_TRANSLATIONS,
        const bool TRANSLATE_INTERNAL = true
      );

      //! @brief Clone function
      //! @return pointer to new Move
      MoveTranslateRandom *Clone() const;

      //! @brief destructor
      ~MoveTranslateRandom();

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
      //! @return reference to movable object
      MovableInterface &Move( MovableInterface &MOVEABLE_OBJECT) const;

      //! @brief generate a TransformationMatrix
      //! @return a translation vector
      linal::Vector3D GenerateRandomTranslation() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief generate a TransformationMatrix
      //! @param MIN_TRANSLATIONS
      //! @param MAX_TRANSLATIONS
      //! @return a translation vector
      static linal::Vector3D GenerateRandomTranslation
      (
        const linal::Vector3D &MAX_TRANSLATIONS,
        const linal::Vector3D &MIN_TRANSLATIONS = linal::Vector3D( 0.0)
      );

    }; // class MoveTranslateRandom

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_MOVE_TRANSLATE_RANDOM_H_
