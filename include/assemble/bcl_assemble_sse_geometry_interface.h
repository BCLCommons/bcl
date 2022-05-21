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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_INTERFACE_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_INTERFACE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "coord/bcl_coord_geometry_interface.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "signal/bcl_signal_signal.h"
#include "util/bcl_util_ptr_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryInterface
    //! @brief This is the interface class from which SSE and SSEGeometry is derived
    //! @details This class provides a clean interface to be used in SSE and SSEGeometry. It has functions that allow access
    //! to data members in derived class that are required for SSE Geometry related operations such as packing.
    //!
    //! @remarks example unnecessary
    //! @author karakam, weinerbe
    //! @date Apr 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryInterface :
      public coord::GeometryInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief get SSType
      //! @return SSType
      virtual const biol::SSType &GetType() const = 0;

      //! @brief returns identification
      //! @return identification
      virtual std::string GetIdentification() const = 0;

      //! @brief return the length of the Z axis
      //! @return the length of Z axis
      virtual double GetLength() const = 0;

      //! @brief returns whether this SSEGeometryInterface is defined
      //! @return whether this SSEGeometryInterface is defined
      virtual bool IsDefined() const = 0;

      //! @brief returns the starting point on z axis
      //! @return the starting point on z axis
      linal::Vector3D BeginOfZ() const
      {
        return GetCenter() - GetLength() * GetAxis( coord::GetAxes().e_Z) / double( 2.0);
      }

      //! @brief returns the ending point on z axis
      //! @return the ending point on z axis
      linal::Vector3D EndOfZ() const
      {
        return GetCenter() + GetLength() * GetAxis( coord::GetAxes().e_Z) / double( 2.0);
      }

      //! @brief return ShPtrVector of sub-SSEGeometries
      //! @return ShPtrVector of sub-SSEGeometries
      virtual const util::ShPtrVector< SSEGeometryInterface> &GetFragments() const = 0;

      //! @brief get central amino acid of geometry
      //! @return central amino acid of geometry
      virtual const int GetCentralAA() const = 0;

      //! @brief return ShPtrVector of sub-SSEGeometries
      //! @return ShPtrVector of sub-SSEGeometries
      virtual util::SiPtrVector< const SSEGeometryInterface> GetSSEGeometries() const = 0;

      //! @brief access to the coordinate change signal handler
      //! @return const ref to the SignalHandler that emits on coordinate changes
      virtual signal::Signal1< const SSEGeometryInterface &> &GetGeometryCoordinateChangeSignal() const = 0;

      //! @brief access to the destructor signal handler
      //! @return const ref to the SignalHandler that emits on destruction
      virtual signal::Signal1< const SSEGeometryInterface &> &GetGeometryDestructorSignal() const = 0;

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      virtual linal::Vector3D GetCenter() const
      {
        return GetOrientation().GetOrigin();
      }

      //! @brief return the orientation of the geometry for a given axis
      //! @param AXIS Axis of interest
      //! @return the orientation of the sse for a given axis
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const
      {
        return GetOrientation().GetAxis( AXIS);
      }

      //! @brief returns the orientation of the geometry
      //! @return the orientation of the geometry
      math::RotationMatrix3D GetRotation() const
      {
        return GetOrientation().GetRotation();
      }

    }; // class SSEGeometryInterface

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryInterfaceLessThan
    //! @brief This is a function class for comparing two SSEGeometryInterface derived objects
    //!
    //! @author weinerbe, karakam
    //! @date Jul 22, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryInterfaceLessThan
    {
    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if SSE_GEO_A is less than SSE_GEO_B in size
      //! @param SSE_GEO_A first SSEGeometryInterface
      //! @param SSE_GEO_B second SSEGeometryInterface
      //! @return true if SSE_GEO_A is less than SSE_GEO_B in size
      bool operator()
      (
        const SSEGeometryInterface &SSE_GEO_A,
        const SSEGeometryInterface &SSE_GEO_B
      ) const
      {
        const bool equal_lengths( math::EqualWithinTolerance( SSE_GEO_A.GetLength(), SSE_GEO_B.GetLength()));

        // if A is shorter than B
        if( SSE_GEO_A.GetLength() < SSE_GEO_B.GetLength() && !equal_lengths)
        {
          return true;
        }

        // if A is the same length as B
        if( equal_lengths)
        {
          // if A's type is less than B's type
          if( SSE_GEO_A.GetType() < SSE_GEO_B.GetType())
          {
            return true;
          }

          // if they have the same type
          if( SSE_GEO_A.GetType() == SSE_GEO_B.GetType())
          {
            // check the identification
            if( SSE_GEO_A.GetIdentification() < SSE_GEO_B.GetIdentification())
            {
              return true;
            }
          }
        }

        // for other cases return false
        return false;
      }

      //! @brief return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
      //! @param PTR_SSE_GEO_A first SSEGeometryInterface
      //! @param PTR_SSE_GEO_B second SSEGeometryInterface
      //! @return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
      bool operator()
      (
        const util::PtrInterface< SSEGeometryInterface> &PTR_SSE_GEO_A,
        const util::PtrInterface< SSEGeometryInterface> &PTR_SSE_GEO_B
      ) const
      {
        return operator()( *PTR_SSE_GEO_A, *PTR_SSE_GEO_B);
      }

      //! @brief return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
      //! @param PTR_SSE_GEO_A first SSEGeometryInterface
      //! @param PTR_SSE_GEO_B second SSEGeometryInterface
      //! @return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
      bool operator()
      (
        const util::PtrInterface< const SSEGeometryInterface> &PTR_SSE_GEO_A,
        const util::PtrInterface< const SSEGeometryInterface> &PTR_SSE_GEO_B
      ) const
      {
        return operator()( *PTR_SSE_GEO_A, *PTR_SSE_GEO_B);
      }

    }; // class SSEGeometryInterfaceLessThan

    //! @brief boolean operator SSE_GEOMETRY__LHS == SSE_GEOMETRY_RHS
    //! @param SSE_GEOMETRY_LHS first SSEGeometryInterface derived object
    //! @param SSE_GEOMETRY_RHS second SSEGeometryInterface derived object
    //! @return whether SSE_GEOMETRY_LHS is equal to SSE_GEOMETRY_RHS
    inline bool operator ==
    (
      const SSEGeometryInterface &SSE_GEOMETRY_LHS,
      const SSEGeometryInterface &SSE_GEOMETRY_RHS
    )
    {
      return
      (
        SSE_GEOMETRY_LHS.GetIdentification() == SSE_GEOMETRY_RHS.GetIdentification() &&
        SSE_GEOMETRY_LHS.GetOrientation()    == SSE_GEOMETRY_RHS.GetOrientation() &&
        math::EqualWithinTolerance( SSE_GEOMETRY_LHS.GetLength(), SSE_GEOMETRY_RHS.GetLength()) &&
        SSE_GEOMETRY_LHS.GetType()           == SSE_GEOMETRY_RHS.GetType()
      );
    }

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_INTERFACE_H_ 
