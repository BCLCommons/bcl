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
#include "coord/bcl_coord_move_translate_external_axis.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_range.h"
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MoveTranslateExternalAxis::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveTranslateExternalAxis())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveTranslateExternalAxis::MoveTranslateExternalAxis() :
      m_MinTranslation( util::GetUndefinedDouble()),
      m_MaxTranslation( util::GetUndefinedDouble()),
      m_ExternalAxis( GetAxes().e_Undefined)
    {
    }

    //! @brief construct from translation range and axis
    //! @param MIN_TRANSLATION minimal distance for translation
    //! @param MAX_TRANSLATION maximal distance for translation
    //! @param EXTERNAL_AXIS axis to translate towards or away from
    MoveTranslateExternalAxis::MoveTranslateExternalAxis
    (
      const double MIN_TRANSLATION,
      const double MAX_TRANSLATION,
      const Axis &EXTERNAL_AXIS
    ) :
      m_MinTranslation( MIN_TRANSLATION),
      m_MaxTranslation( MAX_TRANSLATION),
      m_ExternalAxis( EXTERNAL_AXIS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoveTranslateExternalAxis
    MoveTranslateExternalAxis *MoveTranslateExternalAxis::Clone() const
    {
      return new MoveTranslateExternalAxis( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoveTranslateExternalAxis::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveTranslateExternalAxis::GetAlias() const
    {
      static const std::string s_name( "MoveTranslateExternalAxis");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveTranslateExternalAxis::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Translate along an external axis.");
      serializer.AddInitializer
      (
        "min translation",
        "minimum translation distance",
        io::Serialization::GetAgent( &m_MinTranslation)
      );
      serializer.AddInitializer
      (
        "max translation",
        "maximum translation distance",
        io::Serialization::GetAgent( &m_MaxTranslation)
      );
      serializer.AddInitializer
      (
        "external axis",
        "translation axis",
        io::Serialization::GetAgent( &m_ExternalAxis)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveTranslateExternalAxis::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      // get a random translation distance
      double distance
      (
        random::GetGlobalRandom().Double( math::Range< double>( m_MinTranslation, m_MaxTranslation))
      );

      // randomly decide to flip the sign
      if( random::GetGlobalRandom().Boolean())
      {
        distance *= -1.0;
      }

      // initialize translation vector
      linal::Vector3D translate;

      // if the axis is the x-axis
      if( m_ExternalAxis == GetAxes().e_X)
      {
        // get the vector from the x-axis to the object
        translate = linal::Vector3D( 0.0, MOVEABLE_OBJECT.GetCenter().Y(), MOVEABLE_OBJECT.GetCenter().Z());
      }
      // if the axis is the y-axis
      else if( m_ExternalAxis == GetAxes().e_Y)
      {
        // get the vector from the y-axis to the object
        translate = linal::Vector3D( MOVEABLE_OBJECT.GetCenter().X(), 0.0, MOVEABLE_OBJECT.GetCenter().Z());
      }
      // if the axis is the z-axis
      else if( m_ExternalAxis == GetAxes().e_Z)
      {
        // get the vector from the z-axis to the object
        translate = linal::Vector3D( MOVEABLE_OBJECT.GetCenter().X(), MOVEABLE_OBJECT.GetCenter().Y(), 0.0);
      }

      // normalize the vector
      translate.Normalize();

      // translate the object by the chosen distance
      MOVEABLE_OBJECT.Translate( distance * translate);

      // end
      return MOVEABLE_OBJECT;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace coord

} // namespace bcl
