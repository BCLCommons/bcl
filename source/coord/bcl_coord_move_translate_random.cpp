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
#include "coord/bcl_coord_move_translate_random.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveTranslateRandom::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveTranslateRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveTranslateRandom::MoveTranslateRandom()
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_TRANSLATION        maximal distance for translation
    MoveTranslateRandom::MoveTranslateRandom
    (
      const double MAX_TRANSLATION,
      const bool TRANSLATE_INTERNAL
    ) :
      m_MinTranslation(),
      m_MaxTranslation( MAX_TRANSLATION),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_TRANSLATIONS maximal distance for translation in each direction
    //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
    MoveTranslateRandom::MoveTranslateRandom
    (
      const linal::Vector3D &MAX_TRANSLATIONS,
      const bool TRANSLATE_INTERNAL
    ) :
      m_MinTranslation(),
      m_MaxTranslation( MAX_TRANSLATIONS),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_TRANSLATION minimal distance for translation in each direction (x, y, z)
    //! @param MAX_TRANSLATION maximal distance for translation in each direction (x, y, z)
    //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
    MoveTranslateRandom::MoveTranslateRandom
    (
      const double MIN_TRANSLATION,
      const double MAX_TRANSLATION,
      const bool TRANSLATE_INTERNAL
    ) :
      m_MinTranslation( MIN_TRANSLATION),
      m_MaxTranslation( MAX_TRANSLATION),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_TRANSLATIONS minimal distances for translation in each direction (x, y, z)
    //! @param MAX_TRANSLATIONS maximal distances for translation in each direction (x, y, z)
    //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
    MoveTranslateRandom::MoveTranslateRandom
    (
      const linal::Vector3D &MIN_TRANSLATIONS,
      const linal::Vector3D &MAX_TRANSLATIONS,
      const bool TRANSLATE_INTERNAL
    ) :
      m_MinTranslation( MIN_TRANSLATIONS),
      m_MaxTranslation( MAX_TRANSLATIONS),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {

    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveTranslateRandom *MoveTranslateRandom::Clone() const
    {
      return new MoveTranslateRandom( *this);
    }

    //! @brief destructor
    MoveTranslateRandom::~MoveTranslateRandom()
    {}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveTranslateRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveTranslateRandom::GetAlias() const
    {
      static const std::string s_name( "MoveTranslateRandom");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveTranslateRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Applies a random translation.");
      serializer.AddInitializer
      (
        "min translation",
        "minimum translation",
        io::Serialization::GetAgent( &m_MinTranslation)
      );
      serializer.AddInitializer
      (
        "max translation",
        "maximum translation",
        io::Serialization::GetAgent( &m_MaxTranslation)
      );
      serializer.AddInitializer
      (
        "internal",
        "move to origin before applying translation",
        io::Serialization::GetAgent( &m_TranslateInternal),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT refernce on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveTranslateRandom::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      if( m_TranslateInternal)
      {
        // transform with generated translation relative to internal orientation
        MOVEABLE_OBJECT.Translate( GenerateRandomTranslation().Rotate( MOVEABLE_OBJECT.GetOrientation().GetRotation()));
      }
      else
      {
        // transform with generated translation relative to normal coordinatesystem
        MOVEABLE_OBJECT.Translate( GenerateRandomTranslation());
      }

      // end
      return MOVEABLE_OBJECT;
    }

    //! @brief generate a TransformationMatrix
    //! @return a translation vector
    linal::Vector3D MoveTranslateRandom::GenerateRandomTranslation() const
    {
      return GenerateRandomTranslation( m_MaxTranslation, m_MinTranslation);
    }

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
    linal::Vector3D MoveTranslateRandom::GenerateRandomTranslation
    (
      const linal::Vector3D &MAX_TRANSLATIONS,
      const linal::Vector3D &MIN_TRANSLATIONS
    )
    {
      // initialize a translation vector and set to random translation
      linal::Vector3D translation;
      translation.SetRandomTranslation( linal::Vector3D( MAX_TRANSLATIONS - MIN_TRANSLATIONS));

      // now add the min translations, if random translation is negative then subtract the translation otherwise add
      translation.X() += translation.X() < 0.0 ? -MIN_TRANSLATIONS.X() : MIN_TRANSLATIONS.X();
      translation.Y() += translation.Y() < 0.0 ? -MIN_TRANSLATIONS.Y() : MIN_TRANSLATIONS.Y();
      translation.Z() += translation.Z() < 0.0 ? -MIN_TRANSLATIONS.Z() : MIN_TRANSLATIONS.Z();

      // end
      return translation;
    }

  } // namespace coord
} // namespace bcl
