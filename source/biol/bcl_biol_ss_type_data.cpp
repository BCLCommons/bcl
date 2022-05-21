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
#include "biol/bcl_biol_ss_type_data.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSTypeData::s_Instance( GetObjectInstances().AddInstance( new SSTypeData()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined secondary structure type
    SSTypeData::SSTypeData() :
      m_OneLetterCode( ' '),
      m_IsStructured( false),
      m_RadialExtent( util::GetUndefined< double>()),
      m_AnglePerTurn( util::GetUndefined< double>()),
      m_RiseInZPerResidue( util::GetUndefined< double>()),
      m_IdealPhi( util::GetUndefined< double>()),
      m_IdealPsi( util::GetUndefined< double>()),
      m_FragmentLength( util::GetUndefined< size_t>()),
      m_ContactWindowRadius( util::GetUndefined< size_t>()),
      m_TransformationMatrixForResidues(),
      m_ThreeStatePrediction( util::GetUndefined< double>()),
      m_PhiRange(),
      m_PsiRange()
    {
    }

    //! @brief construct secondary structure type from provided data
    //! @param ONE_LETTER_CODE one letter description
    //! @param IS_STRUCTURED Whether it has a defined geometry (only helix and strand has it)
    //! @param RADIAL_EXTENT radial extensions in angstroms
    //! @param ANGLE_PER_TURN angular turn per residue
    //! @param RISE_IN_Z_PER_RESIDUE rise in z axis per residue
    //! @param IDEAL_PHI ideal phi angle
    //! @param IDEAL_PSI ideal psi angle
    //! @param FRAGMENT_LENGTH Length of fragments created from this SSType that describe the curvature
    //! @param CONTACT_WINDOW_RADIUS number of neighbor residues in contact prediction descriptors
    //! @param THREE_STATE_PREDICTION Vector that contains three state prediction
    //! @param BACKBONE_PHI_RANGE permissible range of backbone phi angle for this ss type
    //! @param BACKBONE_PSI_RANGE permissible range of backbone psi angle for this ss type
    SSTypeData::SSTypeData
    (
      const char ONE_LETTER_CODE,
      const bool IS_STRUCTURED,
      const double RADIAL_EXTENT,
      const double ANGLE_PER_TURN,
      const double RISE_IN_Z_PER_RESIDUE,
      const double IDEAL_PHI,
      const double IDEAL_PSI,
      const size_t FRAGMENT_LENGTH,
      const size_t CONTACT_WINDOW_RADIUS,
      const linal::Vector3D &THREE_STATE_PREDICTION,
      const math::Range< double> &BACKBONE_PHI_RANGE,
      const math::Range< double> &BACKBONE_PSI_RANGE
    ) :
      m_OneLetterCode( ONE_LETTER_CODE),
      m_IsStructured( IS_STRUCTURED),
      m_RadialExtent( RADIAL_EXTENT),
      m_AnglePerTurn( ANGLE_PER_TURN),
      m_RiseInZPerResidue( RISE_IN_Z_PER_RESIDUE),
      m_IdealPhi( IDEAL_PHI),
      m_IdealPsi( IDEAL_PSI),
      m_FragmentLength( FRAGMENT_LENGTH),
      m_ContactWindowRadius( CONTACT_WINDOW_RADIUS),
      m_TransformationMatrixForResidues
      (
        math::TransformationMatrix3D
        (
          linal::Vector3D( 0.0, 0.0, RISE_IN_Z_PER_RESIDUE)
        )( math::RotationMatrix3D( coord::GetAxes().e_Z, ANGLE_PER_TURN))
      ),
      m_ThreeStatePrediction( THREE_STATE_PREDICTION),
      m_PhiRange( BACKBONE_PHI_RANGE),
      m_PsiRange( BACKBONE_PSI_RANGE)
    {
    }

    //! @brief virtual copy constructor
    SSTypeData *SSTypeData::Clone() const
    {
      return new SSTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return one letter code
    //! @return one letter code
    const char SSTypeData::GetOneLetterCode() const
    {
      return m_OneLetterCode;
    }

    //! @brief return if it is's structured
    //! @return whether it's structured
    bool SSTypeData::IsStructured() const
    {
      return m_IsStructured;
    }

    //! @brief return radial extensions in angstroms
    //! @return radial extensions in angstroms
    const double SSTypeData::GetRadialExtent() const
    {
      return m_RadialExtent;
    }

    //! @brief return angle per turn
    //! @return angle per turn
    const double SSTypeData::GetAnglePerTurn() const
    {
      return m_AnglePerTurn;
    }

    //! @brief return rise in z axis per residue
    //! @return rise in z axis per residue
    const double SSTypeData::GetRiseInZPerResidue() const
    {
      return m_RiseInZPerResidue;
    }

    //! @brief return ideal phi angle
    //! @return ideal phi angle
    const double SSTypeData::GetIdealPhi() const
    {
      return m_IdealPhi;
    }

    //! @brief return ideal psi angle
    //! @return ideal psi angle
    const double SSTypeData::GetIdealPsi() const
    {
      return m_IdealPsi;
    }

    //! @brief return Length of fragments created from this SSType that describe the curvature
    //! @return Length of fragments created from this SSType that describe the curvature
    const size_t SSTypeData::GetFragmentLength() const
    {
      return m_FragmentLength;
    }

    //! @brief return number of neighboring residues on each side to be used in contact prediction descriptors
    //! @return number of neighboring residues on each side to be used in contact prediction descriptors
    const size_t SSTypeData::GetContactWindowRadius() const
    {
      return m_ContactWindowRadius;
    }

    //! @brief return total number of neighbor residues to be used in contact prediction descriptors
    //! @return total number of neighbor residues to be used in contact prediction descriptors
    const size_t SSTypeData::GetContactWindowLength() const
    {
      return 2 * m_ContactWindowRadius + 1;
    }

    //! @brief return transformation matrix for residues for this SSType
    //! @return transformation matrix for residues for this SSType
    const math::TransformationMatrix3D &SSTypeData::GetTransformationMatrixForResidues() const
    {
      return m_TransformationMatrixForResidues;
    }

    //! @brief returns three state prediction vector with this SSType with a prediction of 1.0
    //! @return three state prediction vector with this SSType with a prediction of 1.0
    const linal::Vector3D &SSTypeData::GetThreeStatePrediction() const
    {
      return m_ThreeStatePrediction;
    }

    //! @brief returns backbone phi range for this SSType
    //! @return range of backbone phi angles that are permissible
    const math::Range< double> SSTypeData::GetBackbonePhiRange() const
    {
      return m_PhiRange;
    }

    //! @brief returns backbone psi range for this SSType
    //! @return range of backbone psi angles that are permissible
    const math::Range< double> SSTypeData::GetBackbonePsiRange() const
    {
      return m_PsiRange;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_OneLetterCode,        ISTREAM);
      io::Serialize::Read( m_IsStructured,         ISTREAM);
      io::Serialize::Read( m_RadialExtent,         ISTREAM);
      io::Serialize::Read( m_AnglePerTurn,         ISTREAM);
      io::Serialize::Read( m_RiseInZPerResidue,    ISTREAM);
      io::Serialize::Read( m_IdealPhi,             ISTREAM);
      io::Serialize::Read( m_IdealPsi,             ISTREAM);
      io::Serialize::Read( m_FragmentLength,       ISTREAM);
      io::Serialize::Read( m_ContactWindowRadius,  ISTREAM);
      io::Serialize::Read( m_ThreeStatePrediction, ISTREAM);
      io::Serialize::Read( m_PhiRange,             ISTREAM);
      io::Serialize::Read( m_PsiRange,             ISTREAM);

      // initialize transformation
      m_TransformationMatrixForResidues = math::TransformationMatrix3D( linal::Vector3D( 0.0, 0.0, m_RiseInZPerResidue));
      m_TransformationMatrixForResidues( math::RotationMatrix3D( coord::GetAxes().e_Z, m_AnglePerTurn));

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &SSTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_OneLetterCode,        OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IsStructured,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RadialExtent,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AnglePerTurn,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RiseInZPerResidue,    OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IdealPhi,             OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IdealPsi,             OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FragmentLength,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ContactWindowRadius,  OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ThreeStatePrediction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PhiRange,             OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PsiRange,             OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
