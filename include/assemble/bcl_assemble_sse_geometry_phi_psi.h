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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PHI_PSI_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PHI_PSI_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry.h"
#include "biol/bcl_biol_aa_sequence_phi_psi.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPhiPsi
    //! @brief Representation of a geometry using an SSEGeometry and an AASequencePhiPsi
    //! @details Contains both SSEGeometry and AASequencePhiPsi information.  Allows for creation of geometries based on
    //!          phi-psi values by fitting a dummy sequence, then creating an SSE.  This class is used to represent
    //!          SSE geometries in fold templates.
    //!
    //! @see @link example_assemble_sse_geometry_phi_psi.cpp @endlink
    //! @author weinerbe
    //! @date Feb 2, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPhiPsi :
      public SSEGeometryInterface
    {

    private:

    //////////
    // data //
    //////////

      //! sse geometry
      util::ShPtr< SSEGeometry> m_Geometry;

      //! sequence phi psi
      util::ShPtr< biol::AASequencePhiPsi> m_PhiPsi;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPhiPsi();

      //! @brief construct from an SSE
      //! @param ORIGINAL_SSE SSE to build from
      SSEGeometryPhiPsi( const SSE &ORIGINAL_SSE);

      //! @brief construct from phi/psi angles
      //! @param PHI_PSI phi/psi angles
      //! @param SS_TYPE SS type to use to build the geometry
      //! @param IDENTIFICATION identification for this geometry
      //! @param SET_GEOMETRY bool whether to set the geometry using the given SS type, if it is coil, or undefined,
      //!        the geometry is not constructed
      SSEGeometryPhiPsi
      (
        const biol::AASequencePhiPsi &PHI_PSI,
        const biol::SSType &SS_TYPE = biol::GetSSTypes().e_Undefined,
        const std::string &IDENTIFICATION = "",
        const bool SET_GEOMETRY = false
      );

      //! @brief copy constructor
      //! @param GEOMETRY SSEGeometryPhiPsi to be copied
      SSEGeometryPhiPsi( const SSEGeometryPhiPsi &GEOMETRY);

      //! @brief Clone function
      //! @return pointer to new SSEGeometryPhiPsi
      SSEGeometryPhiPsi *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the SSEGeometry member
      //! @return the SSEGeometry member
      const util::ShPtr< SSEGeometry> &GetSSEGeometry() const
      {
        return m_Geometry;
      }

      //! @brief gets the AASequencePhiPsi member
      //! @return the AASequencePhiPsi member
      const util::ShPtr< biol::AASequencePhiPsi> GetPhiPsi() const
      {
        return m_PhiPsi;
      }

      //! @brief use the phi/psi information to build the SSEGeometry object
      void SetSSEGeometryUsingPhiPsi();

      //! @brief get SSType
      //! @return SSType
      const biol::SSType &GetType() const
      {
        return m_Geometry->GetType();
      }

      //! @brief returns identification
      //! @return identification
      std::string GetIdentification() const
      {
        return m_Geometry->GetIdentification();
      }

      //! @brief return the length of the Z axis
      //! @return the length of Z axis
      double GetLength() const
      {
        return double( m_PhiPsi->GetAngles().GetSize()) * GetType()->GetRiseInZPerResidue();
      }

      //! @brief returns whether this SSEGeometryInterface is defined
      //! @return whether this SSEGeometryInterface is defined
      bool IsDefined() const
      {
        return m_Geometry->IsDefined() && m_PhiPsi->IsDefined();
      }

      //! @brief return ShPtrVector of sub-SSEGeometries
      //! @return ShPtrVector of sub-SSEGeometries
      const util::ShPtrVector< SSEGeometryInterface> &GetFragments() const
      {
        return m_Geometry->GetFragments();
      }

      //! @brief return central amino acid of sse fragments
      //! @return 0 for phi psi geometries
      const int GetCentralAA() const
      {
        return m_Geometry.IsDefined() ? m_Geometry->GetCentralAA() : 0;
      }

      //! @brief return ShPtrVector of sub-SSEGeometries
      //! @return ShPtrVector of sub-SSEGeometries
      util::SiPtrVector< const SSEGeometryInterface> GetSSEGeometries() const
      {
        return m_Geometry->GetSSEGeometries();
      }

      //! @brief returns SiPtrVector of GeometryInterface that make up this GeometryInterface
      //! @return SiPtrVector of GeometryInterface that make up this GeometryInterface
      util::SiPtrVector< const coord::GeometryInterface> GetGeometries() const
      {
        return m_Geometry->GetGeometries();
      }

      //! @brief returns the main axis as a LineSegment3D
      //! @return the main axis as a LineSegment3D
      coord::LineSegment3D GetMainAxis() const
      {
        return m_Geometry->GetMainAxis();
      }

      //! @brief returns the requested extent
      //! @param AXIS axis of interest
      //! @return the requested extent
      double GetExtent( const coord::Axis &AXIS) const
      {
        return m_Geometry->GetExtent( AXIS);
      }

      //! @brief returns the radial extent
      //! @return the radial extent
      double GetRadialExtent() const
      {
        return m_Geometry->GetRadialExtent();
      }

      //! @brief returns the orientation of the object
      //! @return the orientation of the object
      const math::TransformationMatrix3D GetOrientation() const
      {
        return m_Geometry->GetOrientation();
      }

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      virtual linal::Vector3D GetCenter() const
      {
        return m_Geometry->GetCenter();
      }

      //! @brief return the orientation of the geometry for a given axis
      //! @param AXIS Axis of interest
      //! @return the orientation of the sse for a given axis
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const
      {
        return m_Geometry->GetAxis( AXIS);
      }

      //! @brief returns the orientation of the geometry
      //! @return the orientation of the geometry
      math::RotationMatrix3D GetRotation() const
      {
        return m_Geometry->GetRotation();
      }

      //! @brief access to the coordinate change signal handler
      //! @return const ref to the SignalHandler that emits on coordinate changes
      signal::Signal1< const SSEGeometryInterface &> &GetGeometryCoordinateChangeSignal() const
      {
        return m_Geometry->GetGeometryCoordinateChangeSignal();
      }

      //! @brief access to the destructor signal handler
      //! @return const ref to the SignalHandler that emits on destruction
      signal::Signal1< const SSEGeometryInterface &> &GetGeometryDestructorSignal() const
      {
        return m_Geometry->GetGeometryDestructorSignal();
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

      //! @brief operator = SSEGeometryPhiPsi assign operator
      //! @param SSE_GEOMETRY_PHI_PSI_RHS SSEGeometry to be assigned to
      //! @return this SSEGeometryPhiPsi after assignment
      SSEGeometryPhiPsi &operator =( const SSEGeometryPhiPsi &SSE_GEOMETRY_PHI_PSI_RHS);

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

    }; // class SSEGeometryPhiPsi

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPhiPsiLessThan
    //! @brief This is a function class for comparing two SSEGeometryPhiPsi's in less-than fashion
    //!
    //! @author weinerbe
    //! @date Feb 11, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPhiPsiLessThan
    {
    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if SSE_GEO_A is less than SSE_GEO_B in size
      //! @param SSE_GEO_A first SSEGeometryPhiPsi
      //! @param SSE_GEO_B second SSEGeometryPhiPsi
      //! @return true if SSE_GEO_A is less than SSE_GEO_B in size
      bool operator()( const SSEGeometryPhiPsi &SSE_GEO_A, const SSEGeometryPhiPsi &SSE_GEO_B) const;

      //! @brief return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
      //! @param PTR_SSE_GEO_A first SSEGeometryPhiPsi
      //! @param PTR_SSE_GEO_B second SSEGeometryPhiPsi
      //! @return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
      bool operator()
      (
        const util::PtrInterface< SSEGeometryPhiPsi> &PTR_SSE_GEO_A,
        const util::PtrInterface< SSEGeometryPhiPsi> &PTR_SSE_GEO_B
      ) const;

      //! @brief return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
      //! @param PTR_SSE_GEO_A first SSEGeometryPhiPsi
      //! @param PTR_SSE_GEO_B second SSEGeometryPhiPsi
      //! @return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
      bool operator()
      (
        const util::PtrInterface< const SSEGeometryPhiPsi> &PTR_SSE_GEO_A,
        const util::PtrInterface< const SSEGeometryPhiPsi> &PTR_SSE_GEO_B
      ) const;

    }; // class SSEGeometryLessThan

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PHI_PSI_H_
