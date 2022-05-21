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
#include "assemble/bcl_assemble_sse_geometry_packing.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //! @brief conversion to a string from a Orientation
    //! @param ORIENTATION the type to get a string for
    //! @return a string representing that orientation
    const std::string &SSEGeometryPacking::GetOrientationName( const Orientation &ORIENTATION)
    {
      static const std::string s_descriptors[] =
      {
        "PARALLEL",
        "ANTIPARALLEL",
        "UNDEFINED_ORIENTATION",
        GetStaticClassName< Orientation>()
      };

      return s_descriptors[ size_t( ORIENTATION)];
    }

  //////////
  // data //
  //////////

    //! Critical distance between two SSEs
    const double SSEGeometryPacking::s_CriticalDistance( 3.0);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPacking::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPacking())
    );

    //! @brief static function that returns minimal length of fragment interface to get an interaction weight of 1.0
    double SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength()
    {
      // initialize static minimal interface length
      static const double s_default_fragment_minimal_interface_length( 4.0);

      // return
      return s_default_fragment_minimal_interface_length;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPacking::SSEGeometryPacking() :
      m_FirstSSEGeometry(),
      m_SecondSSEGeometry(),
      m_ContactType( contact::GetTypes().e_Undefined),
      m_ShortestConnection(),
      m_Distance( util::GetUndefined< double>()),
      m_OrthogonalConnection( false),
      m_InteractionWeight( util::GetUndefined< double>()),
      m_TwistAngle( util::GetUndefined< double>()),
      m_Orientation( e_UndefinedOrientation),
      m_RelativePosition( util::GetUndefined< double>()),
      m_RelativePositionWeight( util::GetUndefined< double>()),
      m_StrandStrandPairingWeight( util::GetUndefined< double>()),
      m_MinimalInterfaceLength( util::GetUndefined< double>())
    {
    }

    //! @brief construct from pair of SSE geometries
    //! @param SSE_GEOMETRY_A first SSE geometry of interest
    //! @param SSE_GEOMETRY_A second SSE geometry of interest
    //! @param MINIMAL_INTERFACE_LENGTH minimal interface length, set to GetDefaultMinimalInterfaceLength by default
    SSEGeometryPacking::SSEGeometryPacking
    (
      const SSEGeometryInterface &SSE_GEOMETRY_A,
      const SSEGeometryInterface &SSE_GEOMETRY_B,
      const double MINIMAL_INTERFACE_LENGTH,
      const bool &ALLOW_UNDEFINED_TYPES
    ) :
      m_FirstSSEGeometry( SSE_GEOMETRY_A),
      m_SecondSSEGeometry( SSE_GEOMETRY_B),
      m_ShortestConnection(),
      m_OrthogonalConnection( false),
      m_InteractionWeight( util::GetUndefined< double>()),
      m_TwistAngle( util::GetUndefined< double>()),
      m_RelativePosition( util::GetUndefined< double>()),
      m_RelativePositionWeight( util::GetUndefined< double>()),
      m_StrandStrandPairingWeight( util::GetUndefined< double>()),
      m_MinimalInterfaceLength( MINIMAL_INTERFACE_LENGTH)
    {
      Initialize( SSE_GEOMETRY_A, SSE_GEOMETRY_B, ALLOW_UNDEFINED_TYPES);
    }

    //! @brief virtual copy constructor
    //! @return a new copy of this SSEGeometryPacking
    SSEGeometryPacking *SSEGeometryPacking::Clone() const
    {
      return new SSEGeometryPacking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEGeometryPacking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief function to reverse the order of the geometries stored
    void SSEGeometryPacking::Reverse()
    {
      // swap the geometries
      std::swap( m_FirstSSEGeometry, m_SecondSSEGeometry);

      // reverse the contact type
      m_ContactType = contact::GetTypes().Reverse( m_ContactType);

      // reverse the shortest connection
      m_ShortestConnection = m_ShortestConnection.GetReverse();
    }

    //! @brief returns identification string
    //! @return identification string
    std::string SSEGeometryPacking::GetIdentification() const
    {
      // initialize identification string
      const std::string identification
      (
        " GEO_1: " + ( m_FirstSSEGeometry.IsDefined() ? m_FirstSSEGeometry->GetIdentification() : "NULL") +
        "\tGEO_2: " + ( m_SecondSSEGeometry.IsDefined() ? m_SecondSSEGeometry->GetIdentification() : "NULL") +
        "\tco: " + m_ContactType.GetName() +
        "\tdist: " + util::Format()( GetDistance()) +
        "\tortho: " + util::Format()( m_OrthogonalConnection) +
        "\tint_weight: " + util::Format()( m_InteractionWeight) +
        "\tstr_weight: " + util::Format()( m_StrandStrandPairingWeight) +
        "\ttwist:" + util::Format()( m_TwistAngle) +
        "\torient:" + GetOrientationName( m_Orientation) +
        "\trel_pos:" + util::Format()( m_RelativePosition) +
        "\trel_weight:" + util::Format()( m_RelativePositionWeight) +
        "\tint_length:" + util::Format()( m_MinimalInterfaceLength)
      );

      // end
      return identification;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Calculate all SSEGeometryPacking member information for given SSE geometry pair
    //! @param SSE_GEOMETRY_A first SSE of interest
    //! @param SSE_GEOMETRY_B second SSE of interest
    //! @param ALLOW_UNDEFINED_TYPES ALLOW_UNDEFINED_TYPES if true, allow UNDEFINED_HELIX_STRAND and similar as contact types, otherwise,
    //!        choose the closest type to remove ambiguity
    void SSEGeometryPacking::Initialize
    (
      const SSEGeometryInterface &SSE_GEOMETRY_A,
      const SSEGeometryInterface &SSE_GEOMETRY_B,
      const bool &ALLOW_UNDEFINED_TYPES
    )
    {
      // first get the basic contact type
      const contact::Type contact_type( contact::GetTypes().TypeFromSSTypes( SSE_GEOMETRY_A, SSE_GEOMETRY_B));

      // check that the SSEs have defined transformation matrices (bodies)
      if
      (
        SSE_GEOMETRY_A.GetType()->IsStructured() &&
        SSE_GEOMETRY_B.GetType()->IsStructured() &&
        ( !SSE_GEOMETRY_A.IsDefined() || !SSE_GEOMETRY_B.IsDefined())
      )
      {
        // if that's the case set the contact type to undefined
        m_ContactType = contact::GetTypes().e_Undefined;

        // warn user and break
        BCL_MessageCrt
        (
          "The provided sse pair has at least one SSE with invalid body " +
          SSE_GEOMETRY_A.GetIdentification() + " vs " + SSE_GEOMETRY_B.GetIdentification()
        )
        return;
      }

      //determine the shortest connection
      const storage::Pair< coord::LineSegment3D, bool>
        shortest_connection_orthogonal
        (
          coord::ShortestConnectionBetweenGeometries
          (
            SSE_GEOMETRY_A,
            SSE_GEOMETRY_B,
            m_MinimalInterfaceLength
          )
        );
      //set member shortest connection
      m_ShortestConnection = shortest_connection_orthogonal.First();

      // initialize distance
      m_Distance = m_ShortestConnection.GetLength();

      //set member orthogonal connection
      m_OrthogonalConnection = shortest_connection_orthogonal.Second();

      //angle between main axis (z-Axis) - twist angle
      m_TwistAngle = CalulateTwistAngle( SSE_GEOMETRY_A, SSE_GEOMETRY_B, shortest_connection_orthogonal.First());

      // set the orientation
      m_Orientation = OrientationFromSSEs( SSE_GEOMETRY_A, SSE_GEOMETRY_B);

      //weight of interaction
      m_InteractionWeight =
        m_OrthogonalConnection ?
          double( 1.0) :
          InteractionWeight( SSE_GEOMETRY_A, SSE_GEOMETRY_B, shortest_connection_orthogonal.First());

      // switch over contact type
      // if UNKNOWN
      if( contact_type == contact::GetTypes().e_Undefined)
      {
        m_ContactType = contact::GetTypes().e_Undefined;
      }

      // if HELIX_HELIX
      else if( contact_type == contact::GetTypes().HELIX_HELIX)
      {
        m_RelativePosition = 0;
        m_RelativePositionWeight = double( 1.0);
        m_ContactType = contact::GetTypes().HELIX_HELIX;
      }

      // if HELIX_SHEET
      else if( contact_type == contact::GetTypes().HELIX_SHEET)
      {
        m_RelativePosition =
          math::Absolute
          (
            math::g_Pi * 0.5 -
            math::Absolute
            (
              math::g_Pi * 0.5 -
              linal::ProjAngle
              (
                m_ShortestConnection.GetEndPoint(),
                m_ShortestConnection.GetStartPoint(),
                m_ShortestConnection.GetEndPoint() + SSE_GEOMETRY_B.GetAxis( coord::GetAxes().e_X)
              )
            )
          );

        m_RelativePositionWeight = math::WeightBetweenZeroAndPi( 2.0 * m_RelativePosition);

        //determine contact type according to relative position
        //HELIX_SHEET if angle < 30
        if( m_RelativePosition < math::g_Pi / ( ALLOW_UNDEFINED_TYPES ? 6 : 4))
        {
          m_ContactType = contact::GetTypes().HELIX_SHEET;
        }
        //HELIX_STRAND if angle > 60
        else if( !ALLOW_UNDEFINED_TYPES || m_RelativePosition > math::g_Pi / 3)
        {
          m_ContactType = contact::GetTypes().HELIX_STRAND;
        }
        //TWILIGHT_ZONE if 30 < angle < 60
        else
        {
          m_ContactType = contact::GetTypes().UNDEFINED_HELIX_STRAND;
        }
      }

      // if SHEET_HELIX
      else if( contact_type == contact::GetTypes().SHEET_HELIX)
      {
        m_RelativePosition =
          math::Absolute
          (
            math::g_Pi * 0.5 -
            math::Absolute
            (
              math::g_Pi * 0.5 -
              linal::ProjAngle
              (
                m_ShortestConnection.GetStartPoint(),
                m_ShortestConnection.GetEndPoint(),
                m_ShortestConnection.GetStartPoint() + SSE_GEOMETRY_A.GetAxis( coord::GetAxes().e_X)
              )
            )
          );
        m_RelativePositionWeight = math::WeightBetweenZeroAndPi( 2.0 * m_RelativePosition);

        //determine contact type according to relative position
        //SHEET_HELIX if angle < 30
        if( m_RelativePosition < math::g_Pi / ( ALLOW_UNDEFINED_TYPES ? 6 : 4))
        {
          m_ContactType = contact::GetTypes().SHEET_HELIX;
        }
        //STRAND_HELIX if angle > 60
        else if( !ALLOW_UNDEFINED_TYPES || m_RelativePosition > math::g_Pi / 3)
        {
          m_ContactType = contact::GetTypes().STRAND_HELIX;
        }
        //TWILIGHT_ZONE if 30 < angle < 60
        else
        {
          m_ContactType = contact::GetTypes().UNDEFINED_STRAND_HELIX;
        }
      }

      // if STRAND_STRAND
      else if( contact_type == contact::GetTypes().STRAND_STRAND)
      {
        // initialize variables
        const double a
        (
          math::Absolute
          (
            math::g_Pi * 0.5 -
            math::Absolute
            (
              math::g_Pi * 0.5 -
              linal::ProjAngle
              (
                m_ShortestConnection.GetEndPoint(),
                m_ShortestConnection.GetStartPoint(),
                m_ShortestConnection.GetEndPoint() + SSE_GEOMETRY_B.GetAxis( coord::GetAxes().e_X)
              )
            )
          )
        );
        const double b
        (
          math::Absolute
          (
            math::g_Pi * 0.5 -
            math::Absolute
            (
              math::g_Pi * 0.5 -
              linal::ProjAngle
              (
                m_ShortestConnection.GetStartPoint(),
                m_ShortestConnection.GetEndPoint(),
                m_ShortestConnection.GetStartPoint() + SSE_GEOMETRY_A.GetAxis( coord::GetAxes().e_X)
              )
            )
          )
        );

        m_RelativePosition = a + b;

        const double w1( math::WeightBetweenZeroAndPi( 2.0 * a));
        const double w2( math::WeightBetweenZeroAndPi( 2.0 * b));

        m_RelativePositionWeight = w1 * w2;
        m_StrandStrandPairingWeight = ( 1.0 - w1) * ( 1.0 - w2);

        // STRAND_STRAND if angle > 120
        if( m_RelativePosition > 2 * math::g_Pi / ( ALLOW_UNDEFINED_TYPES ? 3 : 4))
        {
          m_ContactType = contact::GetTypes().STRAND_STRAND;
        }
        //SHEET_SHEET if angle < 60
        else if( !ALLOW_UNDEFINED_TYPES || m_RelativePosition < math::g_Pi / 3)
        {
          m_ContactType = contact::GetTypes().SHEET_SHEET;
        }
        //TWILIGHT_ZONE if 60 < angle < 120
        else
        {
          m_ContactType = contact::GetTypes().UNDEFINED_STRAND_STRAND;
        }

      }
      // if the contact type was not matched
      else
      {
        BCL_Exit( "impossible contact type supplied, since TypeFromSSTypes only returns 4 different ones", -1);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read SSEGeometryPacking from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPacking::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ContactType, ISTREAM);
      io::Serialize::Read( m_ShortestConnection, ISTREAM);
      io::Serialize::Read( m_Distance, ISTREAM);
      io::Serialize::Read( m_OrthogonalConnection, ISTREAM);
      io::Serialize::Read( m_InteractionWeight, ISTREAM);
      io::Serialize::Read( m_TwistAngle, ISTREAM);
      io::Serialize::Read( m_Orientation, ISTREAM);
      io::Serialize::Read( m_RelativePosition, ISTREAM);
      io::Serialize::Read( m_RelativePositionWeight, ISTREAM);
      io::Serialize::Read( m_StrandStrandPairingWeight, ISTREAM);
      io::Serialize::Read( m_MinimalInterfaceLength, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write SSEGeometryPacking to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &SSEGeometryPacking::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ContactType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ShortestConnection, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Distance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_OrthogonalConnection, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_InteractionWeight, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TwistAngle, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Orientation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RelativePosition, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RelativePositionWeight, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StrandStrandPairingWeight, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinimalInterfaceLength, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Calculate twist angle between sse pair
    //! @param SSE_GEOMETRY_A first SSE of interest
    //! @param SSE_GEOMETRY_B second SSE of interest
    //! @param SHORTEST_CONNECTION shortest connection to be used
    //! @return twist angle between sse pair
    double
    SSEGeometryPacking::CalulateTwistAngle
    (
      const SSEGeometryInterface &SSE_GEOMETRY_A,
      const SSEGeometryInterface &SSE_GEOMETRY_B,
      const coord::LineSegment3D &SHORTEST_CONNECTION
    )
    {
      return linal::Dihedral
             (
               SHORTEST_CONNECTION.GetStartPoint() + SSE_GEOMETRY_A.GetAxis( coord::GetAxes().e_Z),
               SHORTEST_CONNECTION.GetStartPoint(),
               SHORTEST_CONNECTION.GetEndPoint(),
               SHORTEST_CONNECTION.GetEndPoint()   + SSE_GEOMETRY_B.GetAxis( coord::GetAxes().e_Z)
             );
    }

    //! @brief calculate the the alpha1 (angle between THIS_SSE and shortest connection)
    //! @param SSE_GEOMETRY SSE geometry of interest
    //! @param SHORTEST_CONNECTION shortest connection to be used in angle calculation
    //! @return the the alpha1 (angle between THIS_SSE and shortest connection)
    double SSEGeometryPacking::AngleSSEConnection
    (
      const SSEGeometryInterface &SSE_GEOMETRY,
      const coord::LineSegment3D &SHORTEST_CONNECTION
    )
    {
      return
        math::Absolute
        (
          math::g_Pi * 0.5 - linal::ProjAngle( SSE_GEOMETRY.GetAxis( coord::GetAxes().e_Z), SHORTEST_CONNECTION.GetDirection())
        );
    }

    //! @brief Calculate weight of interaction between sse pair
    //! @param SSE_GEOMETRY_A first SSE of interest
    //! @param SSE_GEOMETRY_B second SSE of interest
    //! @param SHORTEST_CONNECTION shortest connection to be used in weight calculation
    //! @return weight of interaction between sse pair
    double
    SSEGeometryPacking::InteractionWeight
    (
      const SSEGeometryInterface &SSE_GEOMETRY_A,
      const SSEGeometryInterface &SSE_GEOMETRY_B,
      const coord::LineSegment3D &SHORTEST_CONNECTION
    )
    {
      //calculate the deviation of an orthogonal angle between the main axis of each SSE and the shortest connection
      // for a nice pairing where the foot-points are not on the borders of the main axis, they should be 0
      //for a bad pairing (orthogonal) one should be 90 degree
      const double alpha1( AngleSSEConnection( SSE_GEOMETRY_A, SHORTEST_CONNECTION));
      const double alpha2( AngleSSEConnection( SSE_GEOMETRY_B, SHORTEST_CONNECTION));

      //product from two weights, will be maximal if alpha is 0 degree and zero if either is 90 degree
      return math::WeightBetweenZeroAndPi( alpha1 * 2) * math::WeightBetweenZeroAndPi( alpha2 * 2);
    }

    //! @brief static function to calculate the orientation of two SSEs
    //! @param SSE_GEOMETRY_A first SSE of interest
    //! @param SSE_GEOMETRY_B second SSE of interest
    //! @return the orientation of two SSEs
    SSEGeometryPacking::OrientationEnum SSEGeometryPacking::OrientationFromSSEs
    (
      const SSEGeometryInterface &SSE_GEOMETRY_A,
      const SSEGeometryInterface &SSE_GEOMETRY_B
    )
    {
      // check that the SSEs have defined transformation matrices (bodies)
      if
      (
        SSE_GEOMETRY_A.GetType()->IsStructured() &&
        SSE_GEOMETRY_B.GetType()->IsStructured() &&
        ( !SSE_GEOMETRY_A.IsDefined() || !SSE_GEOMETRY_B.IsDefined())
      )
      {
        return OrientationEnum( e_UndefinedOrientation);
      }

      const double angle_start_to_finish
      (
        SSE_GEOMETRY_A.GetCentralAA() < SSE_GEOMETRY_B.GetCentralAA()
        ? linal::ProjAngle( SSE_GEOMETRY_A.BeginOfZ(), SSE_GEOMETRY_A.EndOfZ(), SSE_GEOMETRY_B.BeginOfZ(), SSE_GEOMETRY_B.EndOfZ())
        : linal::ProjAngle( SSE_GEOMETRY_B.BeginOfZ(), SSE_GEOMETRY_B.EndOfZ(), SSE_GEOMETRY_A.BeginOfZ(), SSE_GEOMETRY_A.EndOfZ())
      );

      // set the orientation
      return Orientation( angle_start_to_finish >= math::g_Pi * 0.5);
    }

  } // namespace assemble
} // namespace bcl
