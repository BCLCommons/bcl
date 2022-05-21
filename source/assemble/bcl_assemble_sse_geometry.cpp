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
#include "assemble/bcl_assemble_sse_geometry.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEGeometry::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometry())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from SSType
    //! @param SSTYPE SSType of geometry
    //! @param IDENTIFICATION identification of this geometry
    //! @param LENGTH length of the geometry
    SSEGeometry::SSEGeometry
    (
      const biol::SSType &SSTYPE,
      const std::string &IDENTIFICATION,
      const double &LENGTH
    ) :
      m_SSType( SSTYPE),
      m_Length( LENGTH),
      m_Identification( IDENTIFICATION),
      m_Orientation( math::TransformationMatrix3D()),
      m_MainAxis(),
      m_MainAxisUpToDate( false),
      m_Fragments()
    {
    }

    //! @brief construct from a geometry interface
    //! @param GEOMETRY Geometry to be stored
    SSEGeometry::SSEGeometry( const SSEGeometryInterface &GEOMETRY) :
      m_SSType( GEOMETRY.GetType()),
      m_Length( GEOMETRY.GetLength()),
      m_Identification( GEOMETRY.GetIdentification()),
      m_Orientation( GEOMETRY.GetOrientation()),
      m_MainAxis(),
      m_MainAxisUpToDate( false),
      m_Fragments( GEOMETRY.GetFragments().HardCopy())
    {
    }

    //! @brief construct from a geometry interface and its central amino acid
    //! @param GEOMETRY Geometry to be stored
    //! @param CENTRAL_AA central amino acid
    SSEGeometry::SSEGeometry( const SSEGeometryInterface &GEOMETRY, const int &CENTRAL_AA) :
      m_SSType( GEOMETRY.GetType()),
      m_Length( GEOMETRY.GetLength()),
      m_Identification( GEOMETRY.GetIdentification()),
      m_Orientation( GEOMETRY.GetOrientation()),
      m_MainAxis(),
      m_MainAxisUpToDate( false),
      m_Fragments( GEOMETRY.GetFragments().HardCopy()),
      m_CentralAA( CENTRAL_AA)
    {
    }

    //! @brief copy constructor
    //! @param GEOMETRY SSEGeometry to be copied
    SSEGeometry::SSEGeometry( const SSEGeometry &GEOMETRY) :
      m_SSType( GEOMETRY.m_SSType),
      m_Length( GEOMETRY.m_Length),
      m_Identification( GEOMETRY.m_Identification),
      m_Orientation( GEOMETRY.m_Orientation),
      m_MainAxis( GEOMETRY.m_MainAxis),
      m_MainAxisUpToDate( GEOMETRY.m_MainAxisUpToDate),
      m_Fragments( GEOMETRY.m_Fragments.HardCopy())
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEGeometry
    SSEGeometry *SSEGeometry::Clone() const
    {
      return new SSEGeometry( *this);
    }

    //! @brief destructor
    //! emits destructor signal
    SSEGeometry::~SSEGeometry()
    {
      m_GeometryDestructorSignal.Emit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEGeometry::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return SiPtrVector of SSEGeometries
    //! @return SiPtrVector of SSEGeometries
    util::SiPtrVector< const SSEGeometryInterface> SSEGeometry::GetSSEGeometries() const
    {
      // if the fragments is empty
      if( m_Fragments.IsEmpty())
      {
        return util::SiPtrVector< const SSEGeometryInterface>( 1, this);
      }

      // otherwise return the fragments
      return m_Fragments;
    }

    //! @brief returns SiPtrVector of GeometryInterface that make up this GeometryInterface
    //! @return SiPtrVector of GeometryInterface that make up this GeometryInterface
    util::SiPtrVector< const coord::GeometryInterface> SSEGeometry::GetGeometries() const
    {
      return m_Fragments;
    }

    //! @brief returns the main axis as a LineSegment3D
    //! @return the main axis as a LineSegment3D
    coord::LineSegment3D SSEGeometry::GetMainAxis() const
    {
      if( !m_MainAxisUpToDate)
      {
        m_MainAxisUpToDate = true;
        m_MainAxis = coord::LineSegment3D( BeginOfZ(), EndOfZ());
      }
      return m_MainAxis;
    }

    //! @brief returns the requested extent
    //! @param AXIS axis of interest
    //! @return the requested extent
    double SSEGeometry::GetExtent( const coord::Axis &AXIS) const
    {
      // return the appropriate extent
      if( AXIS == coord::GetAxes().e_X)
      {
        return m_SSType->GetRadialExtent();
      }
      if( AXIS == coord::GetAxes().e_Y)
      {
        return m_SSType->GetRadialExtent();
      }
      if( AXIS == coord::GetAxes().e_Z)
      {
        return m_Length / 2.0;
      }

      // else return undefined
      return util::GetUndefinedDouble();
    }

    //! @brief returns the radial extent
    //! @return the radial extent
    double SSEGeometry::GetRadialExtent() const
    {
      return m_SSType->GetRadialExtent();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Translation to be applied
    void SSEGeometry::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_MainAxisUpToDate = false;
      // translate orientation
      m_Orientation( TRANSLATION);

      // iterate over fragments
      for
      (
        util::ShPtrVector< SSEGeometryInterface>::iterator itr( m_Fragments.Begin()), itr_end( m_Fragments.End());
        itr != itr_end; ++itr
      )
      {
        // translate fragment
        ( *itr)->Translate( TRANSLATION);
      }

      // emit geometry coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);
    }

    //! @brief transform the object by a given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void SSEGeometry::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      m_MainAxisUpToDate = false;
      // transform orientation
      m_Orientation( TRANSFORMATION_MATRIX_3D);

      // iterate over fragments
      for
      (
        util::ShPtrVector< SSEGeometryInterface>::iterator itr( m_Fragments.Begin()), itr_end( m_Fragments.End());
        itr != itr_end; ++itr
      )
      {
        // transform fragment
        ( *itr)->Transform( TRANSFORMATION_MATRIX_3D);
      }

      // emit geometry coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void SSEGeometry::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      m_MainAxisUpToDate = false;
      //transform orientation
      m_Orientation( ROTATION_MATRIX_3D);

      // iterate over fragments
      for
      (
        util::ShPtrVector< SSEGeometryInterface>::iterator itr( m_Fragments.Begin()), itr_end( m_Fragments.End());
        itr != itr_end; ++itr
      )
      {
        // rotate fragment
        ( *itr)->Rotate( ROTATION_MATRIX_3D);
      }

      // emit geometry coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator = SSEGeometry assign operator
    //! @param SSE_GEOMETRY_RHS SSEGeometry to be assigned to
    //! @return this SSEGeometry after assignment
    SSEGeometry &SSEGeometry::operator =( const SSEGeometry &SSE_GEOMETRY_RHS)
    {
      m_SSType = SSE_GEOMETRY_RHS.m_SSType;
      m_Length = SSE_GEOMETRY_RHS.m_Length;
      m_Identification = SSE_GEOMETRY_RHS.m_Identification;
      m_Orientation = SSE_GEOMETRY_RHS.m_Orientation;
      m_Fragments = SSE_GEOMETRY_RHS.m_Fragments.HardCopy();
      m_MainAxisUpToDate = false;

      // emit geometry coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);

      // return this
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometry::Read( std::istream &ISTREAM)
    {
      // read the members
      io::Serialize::Read( m_SSType, ISTREAM);
      io::Serialize::Read( m_Length, ISTREAM);
      io::Serialize::Read( m_Identification, ISTREAM);
      io::Serialize::Read( m_Orientation, ISTREAM);
      io::Serialize::Read( m_Fragments, ISTREAM);

      m_MainAxisUpToDate = false;

      // emit geometry coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometry::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write the members
      io::Serialize::Write( m_SSType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Length, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Identification, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Orientation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Fragments, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
