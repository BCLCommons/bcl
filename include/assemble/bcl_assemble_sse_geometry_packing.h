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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_interface.h"
#include "contact/bcl_contact_types.h"
#include "coord/bcl_coord_line_segment_3d.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPacking
    //! @brief This is a class for calculating packing properties of SSE pairs and provide various related information
    //! @details This class describes the packing properties for two given SSEs. It holds the following information
    //! - Pointers to both of the given geometries it was calculated from
    //! - Minimal interface length that was used to calculate it
    //! - Contact type
    //! - Line segment describing the shortest connection and whether this connection was orthogonal or not
    //! - Twist angle between Z-axis of the SSEs and the orientation defined by it ( parallel or anti-parallel)
    //! - relative position and weight
    //! - strand strand pairing weight
    //!
    //! @see @link example_assemble_sse_geometry_packing.cpp @endlink
    //! @author woetzen, karakam
    //! @date 12.02.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPacking :
      public util::ObjectInterface
    {

    public:

      //! enumerator for the orientation of the packing
      enum Orientation { e_Parallel, e_AntiParallel, e_UndefinedOrientation, s_NumberOrientations};

      //! @brief conversion to a string from a Orientation
      //! @param ORIENTATION the type to get a string for
      //! @return a string representing that orientation
      static const std::string &GetOrientationName( const Orientation &ORIENTATION);

      //! @brief wrapper enum for Orientation
      typedef util::WrapperEnum< Orientation, &GetOrientationName, s_NumberOrientations> OrientationEnum;

    private:

    //////////
    // data //
    //////////

      //! pointer to first SSEGeometryInterface
      util::SiPtr< const SSEGeometryInterface> m_FirstSSEGeometry;

      //! pointer to second SSEGeometryInterface
      util::SiPtr< const SSEGeometryInterface> m_SecondSSEGeometry;

      //! contact type
      contact::Type m_ContactType;

      //! line segment representing the shortest connection between the two sses
      coord::LineSegment3D m_ShortestConnection;

      //! distance
      double m_Distance;

      //! boolean whether connection is unique and orthogonal
      bool m_OrthogonalConnection;

      //! weight for packing - if orthogonal 1 - decreasing with the angle of the shortest distance to the sses
      double m_InteractionWeight;

      //! twist angle between sses
      double m_TwistAngle;

      //! Orientation of two SSEs
      OrientationEnum m_Orientation;

      //! relative position to each other
      double m_RelativePosition;

      //! weight resulting from the relative position (angle X-axes( normal) and shortest connection)
      double m_RelativePositionWeight;

      //! weight between two strands in one sheet
      double m_StrandStrandPairingWeight;

      //! minimal length of interface to get an interaction weight of 1.0
      double m_MinimalInterfaceLength;

    public:

      //! @brief static function that returns minimal length of fragment interface to get an interaction weight of 1.0
      static double GetDefaultFragmentMinimalInterfaceLength();

    //////////
    // data //
    //////////

      //! critical distance between two SSEs
      static const double s_CriticalDistance;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPacking();

      //! @brief construct from pair of SSE geometries
      //! @param SSE_GEOMETRY_A first SSE geometry of interest
      //! @param SSE_GEOMETRY_B second SSE geometry of interest
      //! @param MINIMAL_INTERFACE_LENGTH minimal interface length, set to GetDefaultMinimalInterfaceLength by default
      //! @param ALLOW_UNDEFINED_TYPES if true, allow UNDEFINED_HELIX_STRAND and similar as contact types, otherwise,
      //!        choose the closest type to remove ambiguity
      SSEGeometryPacking
      (
        const SSEGeometryInterface &SSE_GEOMETRY_A,
        const SSEGeometryInterface &SSE_GEOMETRY_B,
        const double MINIMAL_INTERFACE_LENGTH = GetDefaultFragmentMinimalInterfaceLength(),
        const bool &ALLOW_UNDEFINED_TYPES = true
      );

      //! @brief virtual copy constructor
      //! @return a new copy of this SSEGeometryPacking
      SSEGeometryPacking *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get util::SiPtr to first SSEGeometryInterface
      //! @return util::SiPtr to first SSEGeometryInterface
      const util::SiPtr< const SSEGeometryInterface> &GetFirstSSEGeometry() const
      {
        return m_FirstSSEGeometry;
      }

      //! @brief get util::SiPtr to first SSEGeometryInterface
      //! @return util::SiPtr to first SSEGeometryInterface
      const util::SiPtr< const SSEGeometryInterface> &GetSecondSSEGeometry() const
      {
        return m_SecondSSEGeometry;
      }

      //! @brief get ContactType
      //! @return ContactType
      const contact::Type &GetContactType() const
      {
        return m_ContactType;
      }

      //! @brief return shortest connection as Line
      //! @return shortest connection as Line
      const coord::LineSegment3D &GetShortestConnection() const
      {
        return m_ShortestConnection;
      }

      //! @brief return orthogonal connection
      //! @return boolean for orthogonal connection
      bool GetOrthogonalConnection() const
      {
        return m_OrthogonalConnection;
      }

      //! @brief get Distance between secondary structure elements
      //! @return Distance between secondary structure elements
      double GetDistance() const
      {
        return m_Distance;
      }

      //! @brief get TwistAngle
      //! @return TwistAngle
      const double &GetTwistAngle() const
      {
        return m_TwistAngle;
      }

      //! @brief return orientation
      //! return Orientation
      Orientation GetOrientation() const
      {
        return m_Orientation;
      }

      //! @brief get Relative position
      //! @return Relative position
      const double &GetRelativePosition() const
      {
        return m_RelativePosition;
      }

      //! @brief get Relative position weight
      //! @return Relative position weight
      const double &GetRelativePositionWeight() const
      {
        return m_RelativePositionWeight;
      }

      //! @brief get strand strand pairing weight
      //! @return strand strand pairing weight
      const double &GetStrandStrandPairingWeight() const
      {
        return m_StrandStrandPairingWeight;
      }

      //! @brief get the packing weight
      //! @return the packing weight
      const double &GetInteractionWeight() const
      {
        return m_InteractionWeight;
      }

      //! @brief get the minimal interface length
      //! @return the minimal interface length
      const double &GetMinimalInterfaceLength() const
      {
        return m_MinimalInterfaceLength;
      }

      //! @brief return if this packing object is defined
      //! @return if this packing object is defined
      bool IsDefined() const
      {
        // the contact type and both SSEGeometry pointers has to be defined for the packing object to be defined
        return m_ContactType.IsDefined() && m_FirstSSEGeometry.IsDefined() && m_SecondSSEGeometry.IsDefined();
      }

      //! @brief function to reverse the order of the geometries stored
      void Reverse();

      //! @brief returns identification string
      //! @return identification string
      std::string GetIdentification() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Calculate all SSEGeometryPacking member information for given SSE geometry pair
      //! @param SSE_GEOMETRY_A first SSE of interest
      //! @param SSE_GEOMETRY_B second SSE of interest
      //! @param ALLOW_UNDEFINED_TYPES ALLOW_UNDEFINED_TYPES if true, allow UNDEFINED_HELIX_STRAND and similar as contact types, otherwise,
      //!        choose the closest type to remove ambiguity
      void Initialize
      (
        const SSEGeometryInterface &SSE_GEOMETRY_A,
        const SSEGeometryInterface &SSE_GEOMETRY_B,
        const bool &ALLOW_UNDEFINED_TYPES = true
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read SSEGeometryPacking from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write SSEGeometryPacking to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief Calculate twist angle between sse pair
      //! @param SSE_GEOMETRY_A first SSE of interest
      //! @param SSE_GEOMETRY_B second SSE of interest
      //! @param SHORTEST_CONNECTION shortest connection to be used
      //! @return twist angle between sse pair
      static
      double
      CalulateTwistAngle
      (
        const SSEGeometryInterface &SSE_GEOMETRY_A,
        const SSEGeometryInterface &SSE_GEOMETRY_B,
        const coord::LineSegment3D &SHORTEST_CONNECTION
      );

      //! @brief calculate the the alpha1 (angle between THIS_SSE and shortest connection)
      //! @param SSE_GEOMETRY SSE geometry of interest
      //! @param SHORTEST_CONNECTION shortest connection to be used in angle calculation
      //! @return the the alpha1 (angle between THIS_SSE and shortest connection)
      static
      double
      AngleSSEConnection
      (
        const SSEGeometryInterface &SSE_GEOMETRY,
        const coord::LineSegment3D &SHORTEST_CONNECTION
      );

      //! @brief Calculate weight of interaction between sse pair
      //! @param SSE_GEOMETRY_A first SSE of interest
      //! @param SSE_GEOMETRY_B second SSE of interest
      //! @param SHORTEST_CONNECTION shortest connection to be used in weight calculation
      //! @return weight of interaction between sse pair
      static
      double
      InteractionWeight
      (
        const SSEGeometryInterface &SSE_GEOMETRY_A,
        const SSEGeometryInterface &SSE_GEOMETRY_B,
        const coord::LineSegment3D &SHORTEST_CONNECTION
      );

      //! @brief static function to calculate the orientation of two SSEs
      //! @param SSE_GEOMETRY_A first SSE of interest
      //! @param SSE_GEOMETRY_B second SSE of interest
      //! @return the orientation of two SSEs
      static OrientationEnum OrientationFromSSEs
      (
        const SSEGeometryInterface &SSE_GEOMETRY_A,
        const SSEGeometryInterface &SSE_GEOMETRY_B
      );

      //! @brief static function to calculate the orientation of two SSEs given the twist angle
      //! @param TWIST_ANGLE twist angle between the z planes of two SSEs
      //! @return the orientation of two SSEs with the given the twist angle
      static Orientation OrientationFromTwistAngle( const double TWIST_ANGLE)
      {
        // return parallel if larger 90 degrees, otherwise return anti-parallel
        return Orientation( TWIST_ANGLE >= ( math::g_Pi / 2));
      }

    }; // class SSEGeometryPacking

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKING_H_
