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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_interface.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "signal/bcl_signal_signal.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometry
    //! @brief SSEGeometryInterface derived class that describes the geometry of an SSE without the sequence information
    //! @details This class provides the geometrical representation of an SSE without storage for sequence information
    //! (amino acids). This way it can be used for describing templates. It has the following information
    //! - SSType and the length of the Z-axis
    //! - a string identifier
    //! - a TransformationMatrix3d that describes the orientation
    //! - ShPtrVector of sub-geometries (so that bent SSEs can be represented)
    //! - RMSD from idealized structure
    //!
    //! @see @link example_assemble_sse_geometry.cpp @endlink
    //! @author weinerbe, karakam
    //! @date Mar 12, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometry :
      public SSEGeometryInterface
    {

    private:

    //////////
    // data //
    //////////

      //! The type of the structure element: SSTypes enumerator
      biol::SSType m_SSType;

      //! Length of the Z axis
      double m_Length;

      //! Identification
      std::string m_Identification;

      //! The orientation of the SSE
      math::TransformationMatrix3D m_Orientation;

      //! Main axis, cached for performance
      mutable coord::LineSegment3D m_MainAxis;

      //! bool: is main axis up-to-date
      mutable bool m_MainAxisUpToDate;

      //! The fragments representing this geometry
      util::ShPtrVector< SSEGeometryInterface> m_Fragments;

      //! signal handler for coordinate changes for SSEGeometryInterface base
      mutable signal::Signal1< const SSEGeometryInterface &> m_GeometryCoordinateChangeSignal;

      //! signal handler for destructor for SSEGeometryInterface base
      mutable signal::Signal1< const SSEGeometryInterface &> m_GeometryDestructorSignal;

      //! amino acid in the center of the sse geometry
      int m_CentralAA;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from SSType
      //! @param SSTYPE SSType of geometry
      //! @param IDENTIFICATION identification of this geometry
      //! @param LENGTH length of the geometry
      SSEGeometry
      (
        const biol::SSType &SSTYPE = biol::GetSSTypes().e_Undefined,
        const std::string &IDENTIFICATION = "",
        const double &LENGTH = util::GetUndefinedDouble()
      );

      //! @brief construct from a geometry interface
      //! @param GEOMETRY Geometry to be stored
      SSEGeometry( const SSEGeometryInterface &GEOMETRY);

      //! @brief construct from a geometry interface and its central amino acid
      //! @param GEOMETRY Geometry to be stored
      //! @param CENTRAL_AA central amino acid of fragment
      SSEGeometry( const SSEGeometryInterface &GEOMETRY, const int &CENTRAL_AA);

      //! @brief copy constructor
      //! @param GEOMETRY SSEGeometry to be copied
      SSEGeometry( const SSEGeometry &GEOMETRY);

      //! @brief Clone function
      //! @return pointer to new SSEGeometry
      SSEGeometry *Clone() const;

      //! @brief destructor
      //! emits destructor signal
      ~SSEGeometry();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the SS type
      //! @return the SS type
      const biol::SSType &GetType() const
      {
        return m_SSType;
      }

      //! @brief returns identification
      //! @return identification
      std::string GetIdentification() const
      {
        return m_Identification;
      }

      //! @brief sets identification
      //! @param IDENTIFICATION identification
      void SetIdentification( const std::string &IDENTIFICATION)
      {
        m_Identification = IDENTIFICATION;
      }

      //! @brief return the length of the Z axis
      //! @return the length of Z axis
      double GetLength() const
      {
        return m_Length;
      }

      //! @brief returns ShPtrVector of SSEGeometries
      //! @return ShPtrVector of SSEGeometries
      const util::ShPtrVector< SSEGeometryInterface> &GetFragments() const
      {
        return m_Fragments;
      }

      //! @brief return central amino acid of sse fragments
      //! @return central amino acid
      const int GetCentralAA() const
      {
        return m_CentralAA;
      }

      //! @brief return SiPtrVector of SSEGeometries
      //! @return SiPtrVector of SSEGeometries
      util::SiPtrVector< const SSEGeometryInterface> GetSSEGeometries() const;

      //! @brief returns SiPtrVector of GeometryInterface that make up this GeometryInterface
      //! @return SiPtrVector of GeometryInterface that make up this GeometryInterface
      util::SiPtrVector< const coord::GeometryInterface> GetGeometries() const;

      //! @brief returns the main axis as a LineSegment3D
      //! @return the main axis as a LineSegment3D
      coord::LineSegment3D GetMainAxis() const;

      //! @brief returns the requested extent
      //! @param AXIS axis of interest
      //! @return the requested extent
      double GetExtent( const coord::Axis &AXIS) const;

      //! @brief returns the radial extent
      //! @return the radial extent
      double GetRadialExtent() const;

      //! @brief returns the orientation of the object
      //! @return the orientation of the object
      const math::TransformationMatrix3D GetOrientation() const
      {
        return m_Orientation;
      }

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      virtual linal::Vector3D GetCenter() const
      {
        return m_Orientation.GetOrigin();
      }

      //! @brief return the orientation of the geometry for a given axis
      //! @param AXIS Axis of interest
      //! @return the orientation of the sse for a given axis
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const
      {
        return m_Orientation.GetAxis( AXIS);
      }

      //! @brief returns the orientation of the geometry
      //! @return the orientation of the geometry
      math::RotationMatrix3D GetRotation() const
      {
        return m_Orientation.GetRotation();
      }

      //! @brief returns whether this SSEGeometryInterface is defined
      //! @return whether this SSEGeometryInterface is defined
      bool IsDefined() const
      {
        return m_Orientation.IsDefined() && m_SSType.IsDefined() && util::IsDefined( m_Length);
      }

      //! @brief access to the coordinate change signal handler for SSEGeometryInterface base
      //! @return const ref to the SignalHandler that emits on coordinate changes
      signal::Signal1< const SSEGeometryInterface &> &GetGeometryCoordinateChangeSignal() const
      {
        return m_GeometryCoordinateChangeSignal;
      }

      //! @brief access to the destructor signal handler for SSEGeometryInterface base
      //! @return const ref to the SignalHandler that emits on destruction
      signal::Signal1< const SSEGeometryInterface &> &GetGeometryDestructorSignal() const
      {
        return m_GeometryDestructorSignal;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION Translation to be applied
      void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief transform the object by a given TransformationMatrix3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D);

      //! @brief rotate the object by a given RotationMatrix3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator = SSEGeometry assign operator
      //! @param SSE_GEOMETRY_RHS SSEGeometry to be assigned to
      //! @return this SSEGeometry after assignment
      SSEGeometry &operator =( const SSEGeometry &SSE_GEOMETRY_RHS);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SSEGeometry

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_H_
