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
#include "score/bcl_score_accessibility_hydrophobic_moment.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_trigonometric_transition.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AccessibilityHydrophobicMoment::s_Instance
    (
      GetObjectInstances().AddInstance( new AccessibilityHydrophobicMoment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AccessibilityHydrophobicMoment::AccessibilityHydrophobicMoment() :
      m_AccessibilityType(),
      m_WindowSizes()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param ENVIRONMENT the type of environment the accessibility was measured in that should be scored
    //! @param WINDOW_SIZES for sstype, number of restraints included in each window the moment will be calculated for
    AccessibilityHydrophobicMoment::AccessibilityHydrophobicMoment
    (
      const restraint::AccessibilityAA::EnvironmentEnum &ENVIRONMENT,
      const storage::Map< biol::SSType, size_t> WINDOW_SIZES
    ) :
      m_AccessibilityType( ENVIRONMENT),
      m_WindowSizes( WINDOW_SIZES)
    {
    }

    //! @brief default constructor
    AccessibilityHydrophobicMoment::Window::Window() :
      m_CalculatedMoment( util::GetUndefinedDouble()),
      m_ExperimentalMoment( util::GetUndefinedDouble()),
      m_AccessibilityAssignments()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param CALCULATED_MOMENT the value of the moment calculated from structure
    //! @param EXPERIMENTAL_MOMENT the experimentally measured moment
    //! @param ACCESSIBILITY_ASSIGNMENTS the list of accessibilities associated with this window
    AccessibilityHydrophobicMoment::Window::Window
    (
      const linal::Vector3D &CALCULATED_MOMENT,
      const linal::Vector3D &EXPERIMENTAL_MOMENT,
      const storage::List< restraint::AccessibilityAAAssignment> &ACCESSIBILITY_ASSIGNMENTS
    ) :
      m_CalculatedMoment( CALCULATED_MOMENT),
      m_ExperimentalMoment( EXPERIMENTAL_MOMENT),
      m_AccessibilityAssignments( ACCESSIBILITY_ASSIGNMENTS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AccessibilityHydrophobicMoment
    AccessibilityHydrophobicMoment *AccessibilityHydrophobicMoment::Clone() const
    {
      return new AccessibilityHydrophobicMoment( *this);
    }

    //! @brief Clone function
    //! @return pointer to new AccessibilityHydrophobicMoment
    AccessibilityHydrophobicMoment::Window *AccessibilityHydrophobicMoment::Window::Clone() const
    {
      return new AccessibilityHydrophobicMoment::Window( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AccessibilityHydrophobicMoment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AccessibilityHydrophobicMoment::Window::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief a string description of this class
    //! @return std::string which gives the data of this class
    std::string AccessibilityHydrophobicMoment::Window::GetIdentification() const
    {
      std::string identification( "StartWindow\n");

      // iterate over the assignments
      for
      (
        storage::List< restraint::AccessibilityAAAssignment>::const_iterator
          assignment_itr( m_AccessibilityAssignments.Begin()), assignment_itr_end( m_AccessibilityAssignments.End());
        assignment_itr != assignment_itr_end;
        ++assignment_itr
      )
      {
        identification +=
        (
          assignment_itr->GetAABase()->GetIdentification() + " exposure " +
          util::Format()( assignment_itr->GetExposureValue()) + " "
        );

        const storage::Map
        <
          restraint::AccessibilityAA::EnvironmentEnum, double
        > &accessibilities( assignment_itr->GetAccessibility());
        for
        (
          storage::Map
          <
            restraint::AccessibilityAA::EnvironmentEnum, double
          >::const_iterator itr( accessibilities.Begin()), itr_end( accessibilities.End());
          itr != itr_end;
          ++itr
        )
        {
          identification += ( util::Format()( itr->second) + " " + itr->first.GetString() + " ");
        }
        identification += '\n';
      }

      identification += ( "calculated_moment " + util::Format()( m_CalculatedMoment) + " length " + util::Format()( m_CalculatedMoment.Norm()) + '\n');
      identification += ( "experiment_moment " + util::Format()( m_ExperimentalMoment) + " length " + util::Format()( m_ExperimentalMoment.Norm()) + '\n');
      identification += "EndWindow\n";

      return identification;
    }

    //! @brief a string description of this class
    //! @return std::string which gives the data of this class
    std::string AccessibilityHydrophobicMoment::Window::GetIdentificationInLine() const
    {
      std::string identification( "|");

      // iterate over the assignments
      for
      (
        storage::List< restraint::AccessibilityAAAssignment>::const_iterator
          assignment_itr( m_AccessibilityAssignments.Begin()), assignment_itr_end( m_AccessibilityAssignments.End());
        assignment_itr != assignment_itr_end;
        ++assignment_itr
      )
      {
        identification += ( '<' + assignment_itr->GetAABase()->GetIdentification() + '>');
      }
      identification += "|";
      return identification;
    }

    //! @brief gives the moment calculated from structure
    //! @return vector 3d which is the moment calculated from structure
    const linal::Vector3D &AccessibilityHydrophobicMoment::Window::GetCalculatedMoment() const
    {
      return m_CalculatedMoment;
    }

    //! @brief gives the moment determined from experiment
    //! @return vector 3d which is the moment determined from experiment
    const linal::Vector3D &AccessibilityHydrophobicMoment::Window::GetExperimentMoment() const
    {
      return m_ExperimentalMoment;
    }

    //! @brief gives the list of accessibilities associated with this window
    //! @return storage::List< restraint::AccessibilityAAAssignment> which are the window's accessibilities
    const storage::List
    <
      restraint::AccessibilityAAAssignment
    > &AccessibilityHydrophobicMoment::Window::GetAccessibilities() const
    {
      return m_AccessibilityAssignments;
    }

  ////////////////
  // operations //
  ////////////////

    storage::List< AccessibilityHydrophobicMoment::Window>
    AccessibilityHydrophobicMoment::CalculateHydrophobicMomentWindows
    (
      const storage::List< restraint::AccessibilityAAAssignment> &ASSIGNMENTS, const size_t WINDOW_SIZE,
      const restraint::AccessibilityAA::EnvironmentEnum &ENVIRONMENT_TYPE
    )
    {
      storage::List< AccessibilityHydrophobicMoment::Window> sse_windows;

      // index into the list
      size_t assignment_index( 0);
      const size_t assignment_list_size( ASSIGNMENTS.GetSize());

      // iterate over all the assignments
      for
      (
        storage::List< restraint::AccessibilityAAAssignment>::const_iterator
          assignment_itr( ASSIGNMENTS.Begin()), assignment_itr_end( ASSIGNMENTS.End());
        assignment_itr != assignment_itr_end && assignment_index + WINDOW_SIZE <= assignment_list_size;
        ++assignment_itr, ++assignment_index
      )
      {
        const AccessibilityHydrophobicMoment::Window window
        (
          CalculateHydrophobicMomentWindow( assignment_itr, assignment_itr_end, WINDOW_SIZE, ENVIRONMENT_TYPE)
        );

        if( !window.GetAccessibilities().IsEmpty())
        {
          sse_windows.PushBack( window);
        }
      }

      return sse_windows;
    }

    AccessibilityHydrophobicMoment::Window AccessibilityHydrophobicMoment::CalculateHydrophobicMomentWindow
    (
      storage::List< restraint::AccessibilityAAAssignment>::const_iterator ITR,
      storage::List< restraint::AccessibilityAAAssignment>::const_iterator ITR_END,
      const size_t WINDOW_SIZE,
      const restraint::AccessibilityAA::EnvironmentEnum &ENVIRONMENT_TYPE
    )
    {
      // will keep track of the location in the current window
      size_t pane( 0);

      // to hold the moment of this window
      linal::Vector3D structural_moment( 0, 0 ,0);
      linal::Vector3D experiment_moment( 0, 0 ,0);

      storage::List< restraint::AccessibilityAAAssignment> assignments;

      // iterate over the current window of assignments
      for
      (
        storage::List< restraint::AccessibilityAAAssignment>::const_iterator assignment_itr( ITR);
        assignment_itr != ITR_END && pane < WINDOW_SIZE;
        ++assignment_itr, ++pane
      )
      {
        // try to find the accessibility in the data
        storage::Pair< bool, double> access
        (
          assignment_itr->GetAccessibilityByEnvironment( ENVIRONMENT_TYPE)
        );

        // true if the type of accessibility is not found in the current assignment or the residue is not defined
        if( !assignment_itr->GetAABase().IsDefined() || !access.First())
        {
          continue;
        }

        const linal::Vector3D exposure
        (
          CalculateSingleHydrophobicMoment( *assignment_itr->GetAABase(), assignment_itr->GetExposureValue())
        );
        const linal::Vector3D accessibility
        (
          CalculateSingleHydrophobicMoment( *assignment_itr->GetAABase(), access.Second())
        );

        structural_moment += exposure;
        experiment_moment += accessibility;
        assignments.PushBack( *assignment_itr);
      }

      const AccessibilityHydrophobicMoment::Window moment_window( structural_moment, experiment_moment, assignments);

      return moment_window;
    }

    linal::Vector3D AccessibilityHydrophobicMoment::CalculateSingleHydrophobicMoment
    (
      const biol::AABase &AA_BASE, const double HYDROPHOBICITY
    )
    {
      // get unit vector in direction of CA->(first side chain atom)
      const linal::Vector3D &ca_coords( AA_BASE.GetCA().GetCoordinates());
      const linal::Vector3D &cb_coords( AA_BASE.GetFirstSidechainAtom().GetCoordinates());
      const linal::Vector3D direction( ( cb_coords - ca_coords).Normalize());

      // calculate the moment
      const linal::Vector3D moment( direction * HYDROPHOBICITY);

      return moment;
    }

    linal::Vector3D AccessibilityHydrophobicMoment::CalculateMomentInXYPlane( const linal::Vector3D &MOMENT, const assemble::SSE &SSE)
    {
      // dihedral to x-axis
      const double x_dihedral
      (
        linal::Dihedral
        (
          MOMENT + SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_Z),
          SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_Z),
          SSE.GetCenter(),
          SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_X)
        )
      );

      BCL_MessageDbg( "dihedral angle between total moment and x-axis is " + util::Format()( x_dihedral));

      // dihedral to y-axis
      const double y_dihedral
      (
        linal::Dihedral
        (
          MOMENT + SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_Z),
          SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_Z),
          SSE.GetCenter(),
          SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_Y)
        )
      );

      BCL_MessageDbg( "dihedral angle between total moment and y-axis is " + util::Format()( y_dihedral));

      // moment in x-y plane
      const double distance( linal::Distance( SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_X), SSE.GetCenter()));
      BCL_MessageDbg( "distance is " + util::Format()( distance));

      const linal::Vector3D moment_in_xy_plane
      (
        linal::CoordinatesAngle
        (
          SSE.GetCenter(),
          SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_X),
          SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_Y),
          distance,
          x_dihedral,
          y_dihedral
        )
      );

      return moment_in_xy_plane;
    }

    std::ostream &AccessibilityHydrophobicMoment::ShowHydrophobicMomentWindows
    (
      const storage::List< AccessibilityHydrophobicMoment::Window> &WINDOW_LIST,
      std::ostream &OSTREAM,
      const assemble::SSE &SSE,
      const std::string &TAG,
      const util::Color &COLOR
    )
    {
      OSTREAM << "from pymol.cgo import *" << '\n';
      OSTREAM << "from pymol import cmd" << '\n';

      // iterate through the list
      for
      (
        storage::List< AccessibilityHydrophobicMoment::Window>::const_iterator
          window_itr( WINDOW_LIST.Begin()), window_itr_end( WINDOW_LIST.End());
        window_itr != window_itr_end;
        ++window_itr
      )
      {
        ShowHydrophobicMomentWindow( *window_itr, OSTREAM, SSE, TAG, COLOR);
      }
      const std::string radius( util::Format()( 0.1));
      {
        // x-axis
        const std::string x_start( util::Format()( SSE.GetCenter().X()));
        const std::string y_start( util::Format()( SSE.GetCenter().Y()));
        const std::string z_start( util::Format()( SSE.GetCenter().Z()));
        const linal::Vector3D xaxis( SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_X));
        const std::string x_end( util::Format()( xaxis.X()));
        const std::string y_end( util::Format()( xaxis.Y()));
        const std::string z_end( util::Format()( xaxis.Z()));
        const std::string x_axis
        (
          "xaxis" + TAG + "=[CYLINDER," + x_start + "," + y_start + "," + z_start + "," + x_end + "," + y_end + "," + z_end + ","
          + radius +
          "," + util::Format()( util::GetColors().e_Yellow->X()) +
          "," + util::Format()( util::GetColors().e_Yellow->Y()) +
          "," + util::Format()( util::GetColors().e_Yellow->Z()) +
          "," + util::Format()( util::GetColors().e_Yellow->X())   +
          "," + util::Format()( util::GetColors().e_Yellow->Y())   +
          "," + util::Format()( util::GetColors().e_Yellow->Z())   +
          "]\n"
        );
        OSTREAM << x_axis << '\n';
        OSTREAM << "cmd.load_cgo( xaxis" + TAG + ", 'x-axis" + TAG + "')" << '\n';
      }
      {
        // y-axis
        const std::string x_start( util::Format()( SSE.GetCenter().X()));
        const std::string y_start( util::Format()( SSE.GetCenter().Y()));
        const std::string z_start( util::Format()( SSE.GetCenter().Z()));
        const linal::Vector3D yaxis( SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_Y));
        const std::string x_end( util::Format()( yaxis.X()));
        const std::string y_end( util::Format()( yaxis.Y()));
        const std::string z_end( util::Format()( yaxis.Z()));
        const std::string y_axis
        (
          "yaxis" + TAG + "=[CYLINDER," + x_start + "," + y_start + "," + z_start + "," + x_end + "," + y_end + "," + z_end + ","
          + radius +
          "," + util::Format()( util::GetColors().e_Yellow->X()) +
          "," + util::Format()( util::GetColors().e_Yellow->Y()) +
          "," + util::Format()( util::GetColors().e_Yellow->Z()) +
          "," + util::Format()( util::GetColors().e_Yellow->X()) +
          "," + util::Format()( util::GetColors().e_Yellow->Y()) +
          "," + util::Format()( util::GetColors().e_Yellow->Z()) +
          "]\n"
        );
        OSTREAM << y_axis << '\n';
        OSTREAM << "cmd.load_cgo( yaxis" + TAG + ", 'y-axis" + TAG + "')" << '\n';
      }
      {
        // z-axis
        const std::string x_start( util::Format()( SSE.GetCenter().X()));
        const std::string y_start( util::Format()( SSE.GetCenter().Y()));
        const std::string z_start( util::Format()( SSE.GetCenter().Z()));
        const linal::Vector3D yaxis( SSE.GetCenter() + SSE.GetAxis( coord::GetAxes().e_Z));
        const std::string x_end( util::Format()( yaxis.X()));
        const std::string y_end( util::Format()( yaxis.Y()));
        const std::string z_end( util::Format()( yaxis.Z()));
        const std::string y_axis
        (
          "zaxis" + TAG + "=[CYLINDER," + x_start + "," + y_start + "," + z_start + "," + x_end + "," + y_end + "," + z_end + ","
          + radius +
          "," + util::Format()( util::GetColors().e_Yellow->X()) +
          "," + util::Format()( util::GetColors().e_Yellow->Y()) +
          "," + util::Format()( util::GetColors().e_Yellow->Z()) +
          "," + util::Format()( util::GetColors().e_Yellow->X()) +
          "," + util::Format()( util::GetColors().e_Yellow->Y()) +
          "," + util::Format()( util::GetColors().e_Yellow->Z()) +
          "]\n"
        );
        OSTREAM << y_axis << '\n';
        OSTREAM << "cmd.load_cgo( zaxis" + TAG + ", 'z-axis" + TAG + "')" << '\n';
      }

      return OSTREAM;
    }

    std::ostream &AccessibilityHydrophobicMoment::ShowHydrophobicMomentWindow
    (
      const AccessibilityHydrophobicMoment::Window &WINDOW,
      std::ostream &OSTREAM,
      const assemble::SSE &SSE,
      const std::string &TAG,
      const util::Color &COLOR
    )
    {
      BCL_MessageDbg( "num resis is " + util::Format()( WINDOW.GetAccessibilities().GetSize()));
      OSTREAM << "moments" + TAG + "=[]" << '\n';

      linal::Vector3D moment( 0, 0, 0);
      const std::string color_r_start( util::Format()( COLOR->X()));
      const std::string color_g_start( util::Format()( COLOR->Y()));
      const std::string color_b_start( util::Format()( COLOR->Z()));
      const std::string color_r_end(   util::Format()( COLOR->X()));
      const std::string color_g_end(   util::Format()( COLOR->Y()));
      const std::string color_b_end(   util::Format()( COLOR->Z()));

      // iterate through the assignments to show the total window moment on each residue
      for
      (
        storage::List< restraint::AccessibilityAAAssignment>::const_iterator
          itr( WINDOW.GetAccessibilities().Begin()), itr_end( WINDOW.GetAccessibilities().End());
        itr != itr_end; ++itr
      )
      {
        BCL_Assert( itr->GetAABase().IsDefined(), "aa ptr is not defined");
        const biol::AABase &aa( *itr->GetAABase());

        // moment line for calculated exposure
        {
          const std::string radius( util::Format()( 0.1));
          const std::string alpha_a( "ALPHA"), alpha_b( util::Format()( 1));
          const linal::Vector3D current_moment( WINDOW.GetCalculatedMoment());
          const linal::Vector3D end_point( current_moment + aa.GetCA().GetCoordinates());

          // draw line indicating moment
          const std::string x_start( util::Format()( aa.GetCA().GetCoordinates().X()));
          const std::string y_start( util::Format()( aa.GetCA().GetCoordinates().Y()));
          const std::string z_start( util::Format()( aa.GetCA().GetCoordinates().Z()));
          const std::string x_end(   util::Format()( end_point.X()));
          const std::string y_end(   util::Format()( end_point.Y()));
          const std::string z_end(   util::Format()( end_point.Z()));
          const std::string moment_string
          (
            "cur_mom=[ " + alpha_a + "," + alpha_b + ",CYLINDER," + x_start + "," + y_start + "," + z_start + ","
            + x_end + "," + y_end + "," + z_end + "," + radius +
            "," + color_r_start +
            "," + color_g_start +
            "," + color_b_start +
            "," + color_r_end +
            "," + color_g_end +
            "," + color_b_end +
            "]\n" +
            "moments" + TAG + ".extend(cur_mom)\n"
          );

          OSTREAM << moment_string << '\n';
        }
        // moment line for experimental accessibility
        {
          const std::string radius( util::Format()( 0.2));
          const std::string alpha_a( "ALPHA"), alpha_b( util::Format()( 0.5));
          const linal::Vector3D current_moment( WINDOW.GetExperimentMoment());
          const linal::Vector3D end_point( current_moment + aa.GetCA().GetCoordinates());

          // draw line indicating moment
          const std::string x_start( util::Format()( aa.GetCA().GetCoordinates().X()));
          const std::string y_start( util::Format()( aa.GetCA().GetCoordinates().Y()));
          const std::string z_start( util::Format()( aa.GetCA().GetCoordinates().Z()));
          const std::string x_end(   util::Format()( end_point.X()));
          const std::string y_end(   util::Format()( end_point.Y()));
          const std::string z_end(   util::Format()( end_point.Z()));
          const std::string moment_string
          (
            "cur_mom=[ " + alpha_a + "," + alpha_b + ",CYLINDER," + x_start + "," + y_start + "," + z_start + ","
            + x_end + "," + y_end + "," + z_end + "," + radius +
            "," + color_r_start +
            "," + color_g_start +
            "," + color_b_start +
            "," + color_r_end +
            "," + color_g_end +
            "," + color_b_end +
            "]\n" +
            "moments" + TAG + ".extend(cur_mom)\n"
          );

          OSTREAM << moment_string << '\n';
        }
      } //< end iteration
      OSTREAM << "cmd.load_cgo( moments" + TAG + ", 'moments" + TAG + "')" << '\n';

      // calculated exposure moment at sse center
      {
        OSTREAM << "calculated_exposure_centered" + TAG + "=[]" << '\n';
        const std::string radius( util::Format()( 0.1));
        const std::string alpha_a( "ALPHA"), alpha_b( util::Format()( 1));
        const linal::Vector3D current_moment( WINDOW.GetCalculatedMoment());
        const linal::Vector3D end_point( current_moment + SSE.GetCenter());

        // draw line indicating moment
        const std::string x_start( util::Format()( SSE.GetCenter().X()));
        const std::string y_start( util::Format()( SSE.GetCenter().Y()));
        const std::string z_start( util::Format()( SSE.GetCenter().Z()));
        const std::string x_end(   util::Format()( end_point.X()));
        const std::string y_end(   util::Format()( end_point.Y()));
        const std::string z_end(   util::Format()( end_point.Z()));
        const std::string moment_string
        (
          "cur_mom=[ " + alpha_a + "," + alpha_b + ",CYLINDER," + x_start + "," + y_start + "," + z_start + ","
          + x_end + "," + y_end + "," + z_end + "," + radius +
          "," + color_r_start +
          "," + color_g_start +
          "," + color_b_start +
          "," + color_r_end +
          "," + color_g_end +
          "," + color_b_end +
          "]\n" +
          "calculated_exposure_centered" + TAG + ".extend(cur_mom)\n"
        );

        OSTREAM << moment_string << '\n';

        OSTREAM << "cmd.load_cgo( calculated_exposure_centered" + TAG + ", 'calculated_exposure_centered" + TAG + "')" << '\n';
      }

      // experimental accessibility moment at sse center
      {
        OSTREAM << "experimental_access_centered" + TAG + "=[]" << '\n';
        const std::string radius( util::Format()( 0.2));
        const std::string alpha_a( "ALPHA"), alpha_b( util::Format()( 0.5));
        const linal::Vector3D current_moment( WINDOW.GetExperimentMoment());
        const linal::Vector3D end_point( current_moment + SSE.GetCenter());

        // draw line indicating moment
        const std::string x_start( util::Format()( SSE.GetCenter().X()));
        const std::string y_start( util::Format()( SSE.GetCenter().Y()));
        const std::string z_start( util::Format()( SSE.GetCenter().Z()));
        const std::string x_end(   util::Format()( end_point.X()));
        const std::string y_end(   util::Format()( end_point.Y()));
        const std::string z_end(   util::Format()( end_point.Z()));
        const std::string moment_string
        (
          "cur_mom=[ " + alpha_a + "," + alpha_b + ",CYLINDER," + x_start + "," + y_start + "," + z_start + ","
          + x_end + "," + y_end + "," + z_end + "," + radius +
          "," + color_r_start +
          "," + color_g_start +
          "," + color_b_start +
          "," + color_r_end +
          "," + color_g_end +
          "," + color_b_end +
          "]\n" +
          "experimental_access_centered" + TAG + ".extend(cur_mom)\n"
        );

        OSTREAM << moment_string << '\n';

        OSTREAM << "cmd.load_cgo( experimental_access_centered" + TAG + ", 'experimental_access_centered" + TAG + "')" << '\n';
      }

      // calculated exposure moment at sse center in x-y plane
      {
        OSTREAM << "calculated_exposure_xy" + TAG + "=[]" << '\n';
        const std::string radius( util::Format()( 0.1));
        const std::string alpha_a( "ALPHA"), alpha_b( util::Format()( 1));
        const linal::Vector3D current_moment( CalculateMomentInXYPlane( WINDOW.GetCalculatedMoment(), SSE));
        const linal::Vector3D end_point( current_moment);

        // draw line indicating moment
        const std::string x_start( util::Format()( SSE.GetCenter().X()));
        const std::string y_start( util::Format()( SSE.GetCenter().Y()));
        const std::string z_start( util::Format()( SSE.GetCenter().Z()));
        const std::string x_end(   util::Format()( end_point.X()));
        const std::string y_end(   util::Format()( end_point.Y()));
        const std::string z_end(   util::Format()( end_point.Z()));
        const std::string moment_string
        (
          "cur_mom=[ " + alpha_a + "," + alpha_b + ",CYLINDER," + x_start + "," + y_start + "," + z_start + ","
          + x_end + "," + y_end + "," + z_end + "," + radius +
          "," + color_r_start +
          "," + color_g_start +
          "," + color_b_start +
          "," + color_r_end +
          "," + color_g_end +
          "," + color_b_end +
          "]\n" +
          "calculated_exposure_xy" + TAG + ".extend(cur_mom)\n"
        );

        OSTREAM << moment_string << '\n';

        OSTREAM << "cmd.load_cgo( calculated_exposure_xy" + TAG + ", 'calculated_exposure_xy" + TAG + "')" << '\n';
      }

      // experimental accessibility moment at sse center in x-y plane
      {
        OSTREAM << "experimental_access_xy" + TAG + "=[]" << '\n';
        const std::string radius( util::Format()( 0.2));
        const std::string alpha_a( "ALPHA"), alpha_b( util::Format()( 0.5));
        const linal::Vector3D current_moment( CalculateMomentInXYPlane( WINDOW.GetExperimentMoment(), SSE));
        const linal::Vector3D end_point( current_moment);

        // draw line indicating moment
        const std::string x_start( util::Format()( SSE.GetCenter().X()));
        const std::string y_start( util::Format()( SSE.GetCenter().Y()));
        const std::string z_start( util::Format()( SSE.GetCenter().Z()));
        const std::string x_end(   util::Format()( end_point.X()));
        const std::string y_end(   util::Format()( end_point.Y()));
        const std::string z_end(   util::Format()( end_point.Z()));
        const std::string moment_string
        (
          "cur_mom=[ " + alpha_a + "," + alpha_b + ",CYLINDER," + x_start + "," + y_start + "," + z_start + ","
          + x_end + "," + y_end + "," + z_end + "," + radius +
          "," + color_r_start +
          "," + color_g_start +
          "," + color_b_start +
          "," + color_r_end +
          "," + color_g_end +
          "," + color_b_end +
          "]\n" +
          "experimental_access_xy" + TAG + ".extend(cur_mom)\n"
        );

        OSTREAM << moment_string << '\n';

        OSTREAM << "cmd.load_cgo( experimental_access_xy" + TAG + ", 'experimental_access_xy" + TAG + "')" << '\n';
      }

      return OSTREAM;
    }

    double AccessibilityHydrophobicMoment::CalculateMomentMagnitudeAgreement
    (
      const restraint::AccessibilityProfileAssignment &ASSIGNMENT
    ) const
    {
      // map for each chain a map that has for the moment/residue for a given sse
      // the order of the sses as ordered by the moment/residue should agree between the accessibility and espozure
      storage::Map< char, storage::Map< double, util::SiPtr< const assemble::SSE> > > exposure_chain_moment_sse;
      storage::Map< char, storage::Map< double, util::SiPtr< const assemble::SSE> > > accessibility_chain_moment_sse;

      for
      (
        storage::Map< util::SiPtr< const assemble::SSE>, storage::List< restraint::AccessibilityAAAssignment> >::const_iterator
          itr( ASSIGNMENT.GetSSEAssignments().Begin()), itr_end( ASSIGNMENT.GetSSEAssignments().End());
        itr != itr_end; ++itr
      )
      {
//        // make sure the siptr is defined
//        BCL_Assert( itr->first.IsDefined(), "SiPtr is not defined");
//        BCL_MessageDbg( m_AccessibilityType.GetName() + " scoring sse " + util::Format()( itr->first->GetIdentification()));
//
//        const util::SiPtr< const assemble::SSE> &current_sse( itr->first);
//        storage::VectorND< 2, storage::Map< util::SiPtr< const biol::AABase>, double> > exposure_accessibility( CalculateExposureAccessibility( itr->second, itr->first->GetType()));
//        BCL_Assert( exposure_accessibility.First().GetSize() == exposure_accessibility.Second().GetSize(), "sizes differ");
//        if( exposure_accessibility.First().IsEmpty())
//        {
//          BCL_MessageDbg( "no exposures of type " + m_AccessibilityType.GetName() + " to score");
//          continue;
//        }
//
//        const double num_data_aas( itr->second.GetSize());
//        linal::Vector3D exposure_moment( CalculateHydrophobicMoment( exposure_accessibility.First()).First());
//        BCL_Assert( exposure_chain_moment_sse[ current_sse->GetChainID()].Insert( std::pair< double, util::SiPtr< const assemble::SSE> >( exposure_moment.Norm() / num_data_aas, current_sse)).second, "could not insert moment, must already be same value");
//        linal::Vector3D accessibility_moment( CalculateHydrophobicMoment( exposure_accessibility.Second()).First());
//        BCL_Assert( accessibility_chain_moment_sse[ current_sse->GetChainID()].Insert( std::pair< double, util::SiPtr< const assemble::SSE> >( accessibility_moment.Norm() / num_data_aas, current_sse)).second, "could not insert moment, must already be same value");
      }

      double score( 0);
      // iterate through the chains
      for
      (
        storage::Map< char, storage::Map< double, util::SiPtr< const assemble::SSE> > >::const_iterator
          exp_itr( exposure_chain_moment_sse.Begin()),      exp_itr_end( exposure_chain_moment_sse.End()),
          epr_itr( accessibility_chain_moment_sse.Begin()), epr_itr_end( accessibility_chain_moment_sse.End());
        exp_itr != exp_itr_end && epr_itr != epr_itr_end;
        ++exp_itr, ++epr_itr
      )
      {
        // get reference on maps
        const storage::Map< double, util::SiPtr< const assemble::SSE> > &exp_map( exp_itr->second);
        const storage::Map< double, util::SiPtr< const assemble::SSE> > &epr_map( epr_itr->second);

        // iterate through the maps
        for
        (
          storage::Map< double, util::SiPtr< const assemble::SSE> >::const_iterator
            exp_map_itr( exp_map.Begin()), exp_map_itr_end( exp_map.End()),
            epr_map_itr( epr_map.Begin()), epr_map_itr_end( epr_map.End());
          exp_map_itr != exp_map_itr_end && epr_map_itr != epr_map_itr_end;
          ++exp_map_itr, ++epr_map_itr
        )
        {
          // sses are the same i.e. ordered the same
          if( exp_map_itr->second == epr_map_itr->second)
          {
            // give bonus
            --score;
          }
          else
          {
            // give penalty
            ++score;
          }
        }
      }

      return score;
    }

  ///////////////
  // operators //
  ///////////////

    double AccessibilityHydrophobicMoment::operator()( const restraint::AccessibilityProfileAssignment &ASSIGNMENT) const
    {
      // iterate through the profiles of the SSEs in ASSIGNMENT
      double score_sum( 0);

      // will hold the sum of the total number of restraint windows in the profile
      size_t window_sum( 0);

      // iterate over the sses and their associated restraints
      for
      (
        storage::Map
        <
          util::SiPtr< const assemble::SSE>, storage::List< restraint::AccessibilityAAAssignment>
        >::const_iterator
          itr( ASSIGNMENT.GetSSEAssignments().Begin()), itr_end( ASSIGNMENT.GetSSEAssignments().End());
        itr != itr_end;
        ++itr
      )
      {
        // make sure the sse siptr is defined and get a reference
        BCL_Assert( itr->first.IsDefined(), "SiPtr is not defined");
        const assemble::SSE &current_sse( *itr->first);
        BCL_MessageDbg( "scoring sse " + util::Format()( current_sse.GetIdentification()));

        // the list of assignments associated with the current sse
        const storage::List< restraint::AccessibilityAAAssignment> &current_sse_assignments( itr->second);

        // determine window size
        size_t window_size( util::GetUndefinedSize_t());
        storage::Map< biol::SSType, size_t>::const_iterator size_itr( m_WindowSizes.Find( current_sse.GetType()));
        if( size_itr != m_WindowSizes.End())
        {
          window_size = size_itr->second;
        }
        // true if the window size is not defined or is larger than the number of assignments
        if( !util::IsDefined( window_size) || window_size > current_sse_assignments.GetSize())
        {
          // set window size to the number of assignments for this sse
          window_size = current_sse_assignments.GetSize();
        }

        // get the list of moments for overlapping windows of the assignments
        const storage::List< AccessibilityHydrophobicMoment::Window> windows
        (
          CalculateHydrophobicMomentWindows( current_sse_assignments, window_size, m_AccessibilityType)
        );

        // no windows created
        if( windows.IsEmpty())
        {
          // go to next sse
          continue;
        }

        const size_t number_of_windows( windows.GetSize());
        BCL_MessageDbg( "number_of_windows " + util::Format()( number_of_windows));

        double current_score_sum( 0);

        // iterate through the list of overlapping windows of assignments
        for
        (
          storage::List< AccessibilityHydrophobicMoment::Window>::const_iterator
            window_itr( windows.Begin()), window_itr_end( windows.End());
          window_itr != window_itr_end;
          ++window_itr
        )
        {
          const double score( CalculateMomentAgreement( *window_itr, current_sse));
          BCL_MessageDbg
          (
            "current window " + util::Format()( window_itr->GetIdentification()) + " score is "
            + util::Format()( score)
          );
          current_score_sum += score;
        }

//        current_score_sum /= double( number_of_windows);

        score_sum += current_score_sum;
        window_sum += number_of_windows;
      }

      // normalize score_sum by the number of windows if there was at least 1 window
      if( window_sum != 0)
      {
        score_sum /= double( window_sum);
      }

      return score_sum;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AccessibilityHydrophobicMoment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AccessibilityType,       ISTREAM);
      io::Serialize::Read( m_WindowSizes      ,       ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AccessibilityHydrophobicMoment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AccessibilityType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WindowSizes,       OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AccessibilityHydrophobicMoment::Window::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_CalculatedMoment  ,       ISTREAM);
      io::Serialize::Read( m_ExperimentalMoment,       ISTREAM);
      io::Serialize::Read( m_AccessibilityAssignments, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AccessibilityHydrophobicMoment::Window::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // read members
      io::Serialize::Write( m_CalculatedMoment        , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExperimentalMoment      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AccessibilityAssignments, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &AccessibilityHydrophobicMoment::WriteDetailedSchemeAndValues
    (
      const restraint::AccessibilityProfileAssignment &ASSIGNMENT,
      std::ostream &OSTREAM
    ) const
    {
      OSTREAM << "AccessibilityHydrophobicMoment::WriteDetailedSchemeAndValues\n";
      OSTREAM << "nr sses is " << ASSIGNMENT.GetSSEAssignments().GetSize() << '\n';
      for
      (
        storage::Map
        <
          util::SiPtr< const assemble::SSE>,
          storage::List< restraint::AccessibilityAAAssignment>,
          assemble::SSELessThanNoOverlap
        >::const_iterator
          sse_itr( ASSIGNMENT.GetSSEAssignments().Begin()), sse_itr_end( ASSIGNMENT.GetSSEAssignments().End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // make sure the sse siptr is defined and get a reference
        BCL_Assert( sse_itr->first.IsDefined(), "SiPtr is not defined");
        const assemble::SSE &current_sse( *sse_itr->first);
        BCL_MessageDbg( "scoring sse " + util::Format()( current_sse.GetIdentification()));

        // the list of assignments associated with the current sse
        const storage::List< restraint::AccessibilityAAAssignment> &current_sse_assignments( sse_itr->second);

        // determine window size
        size_t window_size( util::GetUndefinedSize_t());
        storage::Map< biol::SSType, size_t>::const_iterator size_itr( m_WindowSizes.Find( current_sse.GetType()));
        if( size_itr != m_WindowSizes.End())
        {
          window_size = size_itr->second;
        }
        // true if the window size is not defined or is larger than the number of assignments
        if( !util::IsDefined( window_size) || window_size > current_sse_assignments.GetSize())
        {
          // set window size to the number of assignments for this sse
          window_size = current_sse_assignments.GetSize();
        }
        storage::List< AccessibilityHydrophobicMoment::Window> window_list
        (
          AccessibilityHydrophobicMoment::CalculateHydrophobicMomentWindows
          (
            current_sse_assignments, window_size, m_AccessibilityType
          )
        );

        std::string sse_id( sse_itr->first->GetIdentification());
        {
          util::StringReplacement space_replacer( util::StringReplacement::e_Any, " ", "");
          space_replacer.ReplaceEachIn( sse_id);
        }
        {
          util::StringReplacement space_replacer( util::StringReplacement::e_Any, "<==>", "");
          space_replacer.ReplaceEachIn( sse_id);
        }

        if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug)
        {
          io::OFStream write;
          io::File::MustOpenOFStream
          (
            write, "show_hydrophobic_moment_" + sse_id + "_" + m_AccessibilityType.GetString() + ".py"
          );

          AccessibilityHydrophobicMoment::ShowHydrophobicMomentWindows
          (
            window_list,
            write,
            *sse_itr->first,
            std::string( "sse_id"),
            util::GetColors().e_Cyan
          );
        }

        OSTREAM << "for sse " << current_sse.GetIdentification() << " nr assignments is "
                << ASSIGNMENT.GetSSEAssignments().GetSize() << " and nr windows is " << window_list.GetSize() << '\n';

        // iterate over the window list and print out the windows
        // iterate through the list of overlapping windows of assignments
        for
        (
          storage::List< AccessibilityHydrophobicMoment::Window>::const_iterator
            window_itr( window_list.Begin()), window_itr_end( window_list.End());
          window_itr != window_itr_end;
          ++window_itr
        )
        {
          const double score( CalculateMomentAgreement( *window_itr, current_sse));
          OSTREAM << window_itr->GetIdentification() << '\n';
          const linal::Vector3D exposure_moment_in_xyplane
          (
            CalculateMomentInXYPlane( window_itr->GetCalculatedMoment(), current_sse)
          );
          const linal::Vector3D accessibility_moment_in_xyplane
          (
            CalculateMomentInXYPlane( window_itr->GetExperimentMoment(), current_sse)
          );

          const double xaxis_moment_angle
          (
            linal::ProjAngle( current_sse.GetCenter(), exposure_moment_in_xyplane, accessibility_moment_in_xyplane)
          );
          OSTREAM << "exposure_moment_in_xyplane\n" << exposure_moment_in_xyplane
                  << "\naccessibility_moment_in_xyplane\n" << accessibility_moment_in_xyplane
                  << "\nxaxis_moment_angle " << xaxis_moment_angle
                  << " score: " << score << '\n';
        }
      }

      return OSTREAM;
    }

     double AccessibilityHydrophobicMoment::CalculateMomentAgreement
     (
       const AccessibilityHydrophobicMoment::Window &WINDOW, const assemble::SSE &SSE
     ) const
     {
       const linal::Vector3D exposure_moment_in_xyplane
       (
         CalculateMomentInXYPlane( WINDOW.GetCalculatedMoment(), SSE)
       );
       const linal::Vector3D accessibility_moment_in_xyplane
       (
         CalculateMomentInXYPlane( WINDOW.GetExperimentMoment(), SSE)
       );

       const double xaxis_moment_angle
       (
         linal::ProjAngle( SSE.GetCenter(), exposure_moment_in_xyplane, accessibility_moment_in_xyplane)
       ); //math::g_Pi - math::Absolute( exposure_xaxis_moment_angle - accessibility_xaxis_moment_angle));

       BCL_MessageDbg
       (
         "the angle between the exposure and accessibility moments is "
         + util::Format()( xaxis_moment_angle)
       );

       BCL_MessageDbg
       (
         "sse " + SSE.GetIdentification() + " moment angle : " +
         util::Format()( xaxis_moment_angle) + " current_score : " + util::Format()( ScoreAngle( xaxis_moment_angle))
       );
       return ScoreAngle( xaxis_moment_angle);
     }

     //! @brief scores the angle formed between the moments of the calculated and experimental exposures and a center
     //! @param CALC_CENTER_EXP_ANGLE the angle formed from the the calculated exposure moment in the x-y plane, the
     //!        center of the corresponding sse, and the experimental accessibility moment in the x-y plan.
     //!        calculated->center->experimental
     //! @return double which is the score for the angle between the calculated and experimental moments
     double AccessibilityHydrophobicMoment::ScoreAngle( const double CALC_CENTER_EXP_ANGLE)
     {
//       double current_score( util::GetUndefinedDouble());
       static const double min_score_cutoff( 0);
       static const double bonus_score_cutoff( math::g_Pi);

       static const math::TrigonometricTransition fnc( min_score_cutoff, bonus_score_cutoff, double( -1), double( 0));
       const double current_score( fnc( CALC_CENTER_EXP_ANGLE));

//       if( CALC_CENTER_EXP_ANGLE <= min_score_cutoff)
//       {
//         current_score = -1.0;
//       }
//       else if( CALC_CENTER_EXP_ANGLE > min_score_cutoff && CALC_CENTER_EXP_ANGLE < bonus_score_cutoff)
//       {
//         current_score = fnc( CALC_CENTER_EXP_ANGLE);
//       }
//       else
//       {
//         current_score = 0;
//       }

       return current_score;
     }

  } // namespace score
} // namespace bcl
