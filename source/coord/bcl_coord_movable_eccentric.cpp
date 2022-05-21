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
#include "coord/bcl_coord_movable_eccentric.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MovableEccentric::s_Instance
    (
      GetObjectInstances().AddInstance( new MovableEccentric())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MovableEccentric::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the geometric center of the hinge
    //! @return the geometric center of the hinge
    linal::Vector3D MovableEccentric::GetCenter() const
    {
      return m_Hinge->GetCenter();
    }

    //! @brief returns the orientation of the hinge
    //! @param AXIS The requested Axis for which orientation will be returned
    //! @return the orientation of the hinge
    linal::Vector3D MovableEccentric::GetAxis( const Axis &AXIS) const
    {
      return m_Hinge->GetAxis( AXIS);
    }

    //! @brief returns the orientation and Position as TransformationMatrix3D
    //! @return the orientation and Position as TransformationMatrix3D
    const math::TransformationMatrix3D MovableEccentric::GetOrientation() const
    {
      return m_Hinge->GetOrientation();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Vector along which the translation will occur
    void MovableEccentric::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_Data->Translate( TRANSLATION);
    }

    //! @brief transform the object by a given TRANSFORMATIONMATRIX3D using m_Hinge as origin
    //! @param TRANSFORMATIONMATRIX3D TransformationMatrix3D that will be applied using m_Hinge as origin
    void MovableEccentric::Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D)
    {
      // apply the transformation to MovableData
      m_Data->Transform( TRANSFORMATIONMATRIX3D);
    }

    //! @brief rotate the object by a given ROTATIONMATRIX3D
    //! @brief ROTATIONMATRIX3D RotationMatrix3D that defines the rotation to be applied
    void MovableEccentric::Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D)
    {
      // apply the rotation to MovableData
      m_Data->Rotate( ROTATIONMATRIX3D);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator =
    //! @param MOVABLE_ECCENTRIC MovableEccentric to be copied
    //! @return MovableEccentric object with members updated to MOVABLE_ECCENTRIC
    MovableEccentric &MovableEccentric::operator =( const MovableEccentric &MOVABLE_ECCENTRIC)
    {
      // update data members
      m_Data = MOVABLE_ECCENTRIC.m_Data;
      m_Hinge = MOVABLE_ECCENTRIC.m_Hinge;

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MovableEccentric::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MovableEccentric::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }
  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace coord
} // namespace bcl
