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

#ifndef BCL_ASSEMBLE_SSE_H_
#define BCL_ASSEMBLE_SSE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "coord/bcl_coord.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_interface.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_ss_type_data.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "signal/bcl_signal_signal.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSE
    //! @brief This is a class for structure elements of amino acids.
    //! @details It is derived from the AASequence and has its geometric information as members. The geometrical
    //! information includes the orientation, the extent of X,Y and Z axis. It also has a ShPtrVector of sub-geometries
    //! to allow geometrical representation of bent SSEs.
    //! It can be either helix, strand or coil.
    //!
    //! @see @link example_assemble_sse.cpp @endlink
    //! @author woetzen, karakam, staritrd, meilerj
    //! @date 23.04.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSE :
      public biol::AASequence,
      public SSEGeometryInterface
    {

    private:

    //////////
    // data //
    //////////

      //! The type of the structure element: SSTypes enumerator
      biol::SSType m_SSType;

      //! The orientation of the SSE
      math::TransformationMatrix3D m_Orientation;

      //! The extent along the X-axis (or radial extent for helices)
      double m_XExtent;

      //! The extent along the Y-axis
      double m_YExtent;

      //! The list of SSE fragments that make up this SSE
      util::ShPtrVector< SSEGeometryInterface> m_Fragments;

      //! signal handler for coordinate changes
      mutable signal::Signal1< const SSE &> m_CoordinateChangeSignal;

      //! signal handler for coordinate changes for SSEGeometryInterface base
      mutable signal::Signal1< const SSEGeometryInterface &> m_GeometryCoordinateChangeSignal;

      //! signal handler for destructor for SSE
      mutable signal::Signal1< const SSE &> m_DestructorSignal;

      //! signal handler for destructor for SSEGeometryInterface base
      mutable signal::Signal1< const SSEGeometryInterface &> m_GeometryDestructorSignal;

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
      //! @param SS_TYPE specific SSType
      SSE( const biol::SSType &SS_TYPE = biol::GetSSTypes().e_Undefined);

      //! @brief construct by AASequence and SSType
      //! @param SEQUENCE amino acid sequence
      //! @param SS_TYPE specific SSType
      SSE( const biol::AASequence &SEQUENCE, const biol::SSType &SS_TYPE);

      //! @brief copy constructor
      //! @param SSE_TO_COPY
      SSE( const SSE &SSE_TO_COPY);

      //! @brief Clone function
      //! @return pointer to new SSE
      SSE *Clone() const;

      //! @brief virtual hard copy all amino acids and their data
      //! @return SSE with independent hard copied AADatas
      SSE *HardCopy() const;

      //! @brief destructor
      virtual ~SSE();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get identification of this SSE
      //! @return string with identification for this SSE
      std::string GetIdentification() const;

      //! @brief returns ShPtrVector of SSEGeometries
      //! @return ShPtrVector of SSEGeometries
      const util::ShPtrVector< SSEGeometryInterface> &GetFragments() const
      {
        return m_Fragments;
      }

      //! @brief returns central aa for fragments
      //! @return 0 for whole SSEs
      const int GetCentralAA() const
      {
        return GetFirstMember()->GetSeqID() + GetSize() / 2;
      }

      //! @brief returns SiPtrVector of SSEGeometries
      //! @return SiPtrVector of SSEGeometries
      util::SiPtrVector< const SSEGeometryInterface> GetSSEGeometries() const
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
      util::SiPtrVector< const coord::GeometryInterface> GetGeometries() const
      {
        // if the fragments is empty
        if( m_Fragments.IsEmpty())
        {
          return util::SiPtrVector< const coord::GeometryInterface>( 1, this);
        }

        // otherwise return the fragments
        return m_Fragments;
      }

      //! @brief moves origin of SSE to given LOCATION
      //! @param LOCATION new origin
      void SetOrigin( const linal::Vector3D &LOCATION)
      {
        Translate( LOCATION - GetCenter());
      }

      //! @brief returns the main axis as a LineSegment3D
      //! @return the main axis as a LineSegment3D
      coord::LineSegment3D GetMainAxis() const
      {
        return coord::LineSegment3D( BeginOfZ(), EndOfZ());
      }

      //! @brief returns the requested extent
      //! @param AXIS axis of interest
      //! @return the requested extent
      double GetExtent( const coord::Axis &AXIS) const;

      //! @brief sets the extents
      //! @param EXTENT vector containing the extents in x, y, z
      void SetExtents( const linal::Vector3D &EXTENT);

      //! @brief returns the radial extent
      //! @return the radial extent
      double GetRadialExtent() const
      {
        return m_XExtent;
      }

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const
      {
        return m_Orientation.GetOrigin();
      }

      //! @brief returns the orientation of the object
      //! @return the orientation of the object
      const math::TransformationMatrix3D GetOrientation() const
      {
        return m_Orientation;
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

      //! @brief returns the SSType
      //! @return the SSType
      const biol::SSType &GetType() const
      {
        return m_SSType;
      }

      //! @brief Get SSE hash string to aid in identifying similar chains
      std::string GetHashString() const;

      //! @brief sets the SSType
      //! @param SS_TYPE SSType to be set
      void SetType( const biol::SSType &SS_TYPE);

      //! @brief return the length of the Z axis
      //! @return the length of Z axis
      double GetLength() const
      {
        return GetSize() * m_SSType->GetRiseInZPerResidue();
      }

      //! @brief get the longest contiguous predicted TM-segment
      bool IsPredictedTransmembrane() const;

      //! @brief returns whether the sse is defined
      //! @return true SSE is defined
      bool IsDefined() const
      {
        return m_Orientation.IsDefined() && GetData().IsDefined() && m_SSType.IsDefined();
      }

      //! @brief access to the coordinate change signal handler
      //! @return const ref to the SignalHandler that emits on coordinate changes
      signal::Signal1< const SSE &> &GetCoordinateChangeSignal() const
      {
        return m_CoordinateChangeSignal;
      }

      //! @brief access to the coordinate change signal handler for SSEGeometryInterface base
      //! @return const ref to the SignalHandler that emits on coordinate changes
      signal::Signal1< const SSEGeometryInterface &> &GetGeometryCoordinateChangeSignal() const
      {
        return m_GeometryCoordinateChangeSignal;
      }

      //! @brief access to the destructor signal handler
      //! @return const ref to the SignalHandler that emits on destruction
      signal::Signal1< const SSE &> &GetDestructorSignal() const
      {
        return m_DestructorSignal;
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

      //! @brief rotation around ROTATION_AXIS that goes through ROTATION_POINT and given ANGLE
      //! @param ROTATION_POINT point which rotation should pass through
      //! @param ROTATION_AXIS RotationAxis around which SSE will be rotated
      //! @param ANGLE angle for the rotation
      void Rotate( const linal::Vector3D &ROTATION_POINT, const linal::Vector3D &ROTATION_AXIS, const double &ANGLE);

      //! @brief sets the conformation to idealized at origin
      void SetToIdealConformationAtOrigin();

      //! @brief sets the conformation to idealized in place
      void SetToIdealConformationInPlace();

      //! @brief sets all geometries; main geometry and fragment geometries
      void SetGeometry();

      //! @brief sets the fragment geometries
      //! @param FRAGMENT_LENGTH fragment length to use
      void SetFragmentGeometries( size_t FRAGMENT_LENGTH = 0);

      //! @brief returns a storage::Set of SSEs chopped to pieces of minimal size with spacing of one aa
      //! @param SIZE minimal size
      //! @return a storage::Set of SSEs chopped to pieces of minimal size with spacing of one aa
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> Chop( const size_t &SIZE) const;

      //! @brief finds the coordinates for the given amino acid and prepends it to the sequence
      //! @param AMINO_ACID amino acid to be prepended
      //! @param IDEALIZE bool whether to idealize the SSE prior to prepending
      void Prepend( const biol::AABase &AMINO_ACID, const bool IDEALIZE = true);

      //! @brief calculates ideal coordinates for the given amino acid and appends it to the sequence
      //! @param AMINO_ACID amino acid to be appended
      //! @param IDEALIZE bool whether to idealize the SSE prior to appending
      void Append( const biol::AABase &AMINO_ACID, const bool IDEALIZE = true);

      //! @brief calculates ideal coordinates for the amino acids in the given sequence and prepends them to the sequence
      //! @param AA_SEQUENCE sequence to be prepended
      //! @param IDEALIZE bool whether to idealize the SSE prior to prepending
      void PrependSequence( const AASequence &AA_SEQUENCE, const bool IDEALIZE = true);

      //! @brief calculates ideal coordinates for the amino acids in the given sequence and appends them to the sequence
      //! @param AA_SEQUENCE sequence acid to be appended
      //! @param IDEALIZE bool whether to idealize the SSE prior to appending
      void AppendSequence( const AASequence &AA_SEQUENCE, const bool IDEALIZE = true);

      //! @brief fits SSE using this sequence and the passed SSE as a template
      //! @param SSE_TEMPLATE SSE to be used as a template
      //! @return fit SSE
      void FitToSSE( const SSE &SSE_TEMPLATE);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator = SSE : makes HardCopy of the m_Data through the BaseClass operator =
      //! @param SSE_RHS SSE to be assigned to
      //! @return this SSE after assignment
      SSE &operator =( const SSE &SSE_RHS);

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief sets the main geometry to the geometry of an copy SSE that is idealized, does not set fragment geometries
      void SetMainGeometry();

      //! @brief sets the extents for the geometry
      void SetExtents();

      //! @brief transforms only the geometries while leaving the AASequence coords in place
      //! @param TRANSFORMATION_MATRIX_3D transformation to apply to each geometry
      void TransformGeometries( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D);

      //! @brief resets the geometries
      void ResetGeometries();

    }; // class SSE

    //! @brief boolean operator SSE_LHS == SSE_RHS
    //! @param SSE_LHS first SSE
    //! @param SSE_RHS second SSE
    //! @return whether SSE_LHS is equal to SSE_RHS
    BCL_API bool operator ==( const SSE &SSE_LHS, const SSE &SSE_RHS);

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_SSE_H_
