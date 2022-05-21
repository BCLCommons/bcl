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
#include "coord/bcl_coord_move_translate_defined.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "io/bcl_io_serialization.h"
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
    const util::SiPtr< const util::ObjectInterface> MoveTranslateDefined::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveTranslateDefined())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveTranslateDefined::MoveTranslateDefined() :
      m_Translation(),
      m_TranslateInternal( true)
    {
    }

    //! @brief construct from  translation and internal coordinates or not
    //! @param TRANSLATION a translation vector
    //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
    MoveTranslateDefined::MoveTranslateDefined( const linal::Vector3D &TRANSLATION, const bool TRANSLATE_INTERNAL) :
      m_Translation( TRANSLATION),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveTranslateDefined *MoveTranslateDefined::Clone() const
    {
      return new MoveTranslateDefined( *this);
    }

    //! @brief virtual destructor
    MoveTranslateDefined::~MoveTranslateDefined()
    {}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveTranslateDefined::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveTranslateDefined::GetAlias() const
    {
      static const std::string s_name( "MoveTranslateDefined");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveTranslateDefined::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "translate object");
      serializer.AddInitializer
      (
        "translation",
        "translation vector",
        io::Serialization::GetAgent( &m_Translation)
      );
      serializer.AddInitializer
      (
        "internal",
        "move to the origin before applying the translation",
        io::Serialization::GetAgent( &m_TranslateInternal),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveTranslateDefined::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      if( m_TranslateInternal)
      {
        // transform with generated translation relative to internal orientation
        MOVEABLE_OBJECT.Translate( linal::Vector3D( m_Translation).Rotate( MOVEABLE_OBJECT.GetOrientation().GetRotation()));
      }
      else
      {
        // transform with generated translation relative to normal coordinatesystem
        MOVEABLE_OBJECT.Translate( m_Translation);
      }

      // end
      return MOVEABLE_OBJECT;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace coord
} // namespace bcl
