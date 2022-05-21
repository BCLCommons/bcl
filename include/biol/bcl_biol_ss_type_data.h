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

#ifndef BCL_BIOL_SS_TYPE_DATA_H_
#define BCL_BIOL_SS_TYPE_DATA_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSTypeData
    //! @brief This is a low level helper class to store secondary structure type properties
    //! @details Each SSType that is enumerated in class SSTypes class, is represented with a single SSTypeData
    //! instances which stores various information related to SSE representation.
    //!
    //! @see @link example_biol_ss_type_data.cpp @endlink
    //! @author karakam
    //! @date 01.29.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSTypeData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      char   m_OneLetterCode;       //!< One letter description
      bool   m_IsStructured;        //!< Whether it has a defined geometry (only helix and strand has it)
      double m_RadialExtent;        //!< Radial extensions in angstroms
      double m_AnglePerTurn;        //!< Angular turn per residue
      double m_RiseInZPerResidue;   //!< Rise In Z axis per residue
      double m_IdealPhi;            //!< Phi angle in a ideal SSE of this type
      double m_IdealPsi;            //!< Psi angle in a ideal SSE of this type
      size_t m_FragmentLength;      //!< Length of fragments created from this SSType that describe the curvature

      //! number of neighboring residues included for a given residue in contact prediction descriptors.
      size_t m_ContactWindowRadius;

      //! transformation matrices to be applied to residues of these sses
      math::TransformationMatrix3D m_TransformationMatrixForResidues;

      //! three state prediction vector
      linal::Vector3D m_ThreeStatePrediction;

      //! permissible backbone phi range for this ss type
      math::Range< double> m_PhiRange;

      //! permissible backbone psi range for this ss type
      math::Range< double> m_PsiRange;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined secondary structure type
      SSTypeData();

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
      SSTypeData
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
      );

      //! @brief virtual copy constructor
      SSTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return one letter code
      //! @return one letter code
      const char GetOneLetterCode() const;

      //! @brief return if it is's structured
      //! @return whether it's structured
      bool IsStructured() const;

      //! @brief return radial extensions in angstroms
      //! @return radial extensions in angstroms
      const double GetRadialExtent() const;

      //! @brief return angle per turn
      //! @return angle per turn
      const double GetAnglePerTurn() const;

      //! @brief return rise in z axis per residue
      //! @return rise in z axis per residue
      const double GetRiseInZPerResidue() const;

      //! @brief return ideal phi angle
      //! @return ideal phi angle
      const double GetIdealPhi() const;

      //! @brief return ideal psi angle
      //! @return ideal psi angle
      const double GetIdealPsi() const;

      //! @brief return Length of fragments created from this SSType that describe the curvature
      //! @return Length of fragments created from this SSType that describe the curvature
      const size_t GetFragmentLength() const;

      //! @brief return number of neighboring residues on each side to be used in contact prediction descriptors
      //! @return number of neighboring residues on each side to be used in contact prediction descriptors
      const size_t GetContactWindowRadius() const;

      //! @brief return total number of neighbor residues to be used in contact prediction descriptors
      //! @return total number of neighbor residues to be used in contact prediction descriptors
      const size_t GetContactWindowLength() const;

      //! @brief return transformation matrix for residues for this SSType
      //! @return transformation matrix for residues for this SSType
      const math::TransformationMatrix3D &GetTransformationMatrixForResidues() const;

      //! @brief returns three state prediction vector with this SSType with a prediction of 1.0
      //! @return three state prediction vector with this SSType with a prediction of 1.0
      const linal::Vector3D &GetThreeStatePrediction() const;

      //! @brief returns backbone phi range for this SSType
      //! @return range of backbone phi angles that are permissible
      const math::Range< double> GetBackbonePhiRange() const;

      //! @brief returns backbone psi range for this SSType
      //! @return range of backbone psi angles that are permissible
      const math::Range< double> GetBackbonePsiRange() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; //class SSTypeData

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_SS_TYPE_DATA_H_

