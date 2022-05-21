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

#ifndef BCL_SCORE_ACCESSIBILITY_HYDROPHOBIC_MOMENT_H_
#define BCL_SCORE_ACCESSIBILITY_HYDROPHOBIC_MOMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "io/bcl_io.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_accessibility_aa_assignment.h"
#include "restraint/bcl_restraint_accessibility_profile_assignment.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AccessibilityHydrophobicMoment
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_score_accessibility_hydrophobic_moment.cpp @endlink
    //! @author alexanns
    //! @date Apr 7, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AccessibilityHydrophobicMoment :
      public math::FunctionInterfaceSerializable< restraint::AccessibilityProfileAssignment, double>
    {

    protected:

    //////////
    // data //
    //////////

      //! the type of environment the accessibility was measured in that should be scored
      restraint::AccessibilityAA::EnvironmentEnum m_AccessibilityType;

      //! the number of restraints included in each window the moment will be calculated for
      storage::Map< biol::SSType, size_t> m_WindowSizes;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Window
      //! @brief TODO: add a brief comment
      //! @details TODO: add an detailed description to this class
      //!
      //! @remarks example unnecessary
      //! @author alexanns
      //! @date Jan 31, 2012
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class BCL_API Window :
        public util::ObjectInterface
      {

      private:

      //////////
      // data //
      //////////

        //! the value of the moment calculated from structure
        linal::Vector3D m_CalculatedMoment;

        //! the experimentally measured moment
        linal::Vector3D m_ExperimentalMoment;

        //! the list of accessibilities associated with this window
        storage::List< restraint::AccessibilityAAAssignment> m_AccessibilityAssignments;

      public:

        //! single instance of that class
        static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief default constructor
        Window();

        //! @brief constructor taking member variable parameters
        //! @param CALCULATED_MOMENT the value of the moment calculated from structure
        //! @param EXPERIMENTAL_MOMENT the experimentally measured moment
        //! @param ACCESSIBILITY_ASSIGNMENTS the list of accessibilities associated with this window
        Window
        (
          const linal::Vector3D &CALCULATED_MOMENT,
          const linal::Vector3D &EXPERIMENTAL_MOMENT,
          const storage::List< restraint::AccessibilityAAAssignment> &ACCESSIBILITY_ASSIGNMENTS
        );

        //! @brief Clone function
        //! @return pointer to new Window
        Window *Clone() const;

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! @return the class name as const ref std::string
        const std::string &GetClassIdentifier() const;

        //! @brief a string description of this class
        //! @return std::string which gives the data of this class
        std::string GetIdentification() const;

        //! @brief a string description of this class
        //! @return std::string which gives the data of this class
        std::string GetIdentificationInLine() const;

        //! @brief gives the moment calculated from structure
        //! @return vector 3d which is the moment calculated from structure
        const linal::Vector3D &GetCalculatedMoment() const;

        //! @brief gives the moment determined from experiment
        //! @return vector 3d which is the moment determined from experiment
        const linal::Vector3D &GetExperimentMoment() const;

        //! @brief gives the list of accessibilities associated with this window
        //! @return storage::List< restraint::AccessibilityAAAssignment> which are the window's accessibilities
        const storage::List< restraint::AccessibilityAAAssignment> &GetAccessibilities() const;

      ////////////////
      // operations //
      ////////////////

      ///////////////
      // operators //
      ///////////////

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

      }; // class Window

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AccessibilityHydrophobicMoment();

      //! @brief constructor taking member variable parameters
      //! @param ENVIRONMENT the type of environment the accessibility was measured in that should be scored
      //! @param WINDOW_SIZES for sstype, number of restraints included in each window the moment will be calculated for
      AccessibilityHydrophobicMoment
      (
        const restraint::AccessibilityAA::EnvironmentEnum &ENVIRONMENT,
        const storage::Map< biol::SSType, size_t> WINDOW_SIZES
      );

      //! @brief Clone function
      //! @return pointer to new AccessibilityHydrophobicMoment
      AccessibilityHydrophobicMoment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      static storage::List< AccessibilityHydrophobicMoment::Window> CalculateHydrophobicMomentWindows
      (
        const storage::List< restraint::AccessibilityAAAssignment> &ASSIGNMENTS, const size_t WINDOW_SIZE,
        const restraint::AccessibilityAA::EnvironmentEnum &ENVIRONMENT_TYPE
      );

      static AccessibilityHydrophobicMoment::Window CalculateHydrophobicMomentWindow
      (
        storage::List< restraint::AccessibilityAAAssignment>::const_iterator ITR,
        storage::List< restraint::AccessibilityAAAssignment>::const_iterator ITR_END,
        const size_t WINDOW_SIZE,
        const restraint::AccessibilityAA::EnvironmentEnum &ENVIRONMENT_TYPE
      );

      static linal::Vector3D CalculateSingleHydrophobicMoment
      (
        const biol::AABase &AA_BASE, const double HYDROPHOBICITY
      );

      static std::ostream &ShowHydrophobicMomentWindows
      (
        const storage::List< AccessibilityHydrophobicMoment::Window> &WINDOW_LIST,
        std::ostream &OSTREAM,
        const assemble::SSE &SSE,
        const std::string &TAG,
        const util::Color &COLOR
      );

      static std::ostream &ShowHydrophobicMomentWindow
      (
        const AccessibilityHydrophobicMoment::Window &WINDOW,
        std::ostream &OSTREAM,
        const assemble::SSE &SSE,
        const std::string &TAG,
        const util::Color &COLOR
      );

      static linal::Vector3D CalculateMomentInXYPlane( const linal::Vector3D &MOMENT, const assemble::SSE &SSE);

      double CalculateMomentMagnitudeAgreement
      (
        const restraint::AccessibilityProfileAssignment &ASSIGNMENT
      ) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an accessibility profile that has been assigned and scores it
      //! @param ASSIGNMENT the profile assignment that is going to be scored
      //! @return double that is the score of the accessibility profile assignment
      double operator()( const restraint::AccessibilityProfileAssignment &ASSIGNMENT) const;

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

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
       std::ostream &WriteDetailedSchemeAndValues
      (
        const restraint::AccessibilityProfileAssignment &ASSIGNMENT,
        std::ostream &OSTREAM
      ) const;

      double CalculateMomentAgreement
      (
        const AccessibilityHydrophobicMoment::Window &WINDOW, const assemble::SSE &SSE
      ) const;

      //! @brief scores the angle formed between the moments of the calculated and experimental exposures and a center
      //! @param CALC_CENTER_EXP_ANGLE the angle formed from the the calculated exposure moment in the x-y plane, the
      //!        center of the corresponding sse, and the experimental accessibility moment in the x-y plan.
      //!        calculated->center->experimental
      //! @return double which is the score for the angle between the calculated and experimental moments
      static double ScoreAngle( const double CALC_CENTER_EXP_ANGLE);

    }; // class AccessibilityHydrophobicMoment

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_ACCESSIBILITY_HYDROPHOBIC_MOMENT_H_
