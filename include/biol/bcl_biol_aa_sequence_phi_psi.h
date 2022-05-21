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

#ifndef BCL_BIOL_AA_SEQUENCE_PHI_PSI_H_
#define BCL_BIOL_AA_SEQUENCE_PHI_PSI_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASequencePhiPsi
    //! @brief Stores center residue coordinates and phi/psi values for a sequence
    //! @details Alternate representation of an AASequence in space.  It is defined by the N, CA, and C coords of the
    //!          center residue, along with phi/psi values along the sequence going from N to C-term.  Superimposition
    //!          of the coordinates and subsequent setting of the phi/psi values allows for any sequence to adopt this
    //!          conformation.
    //!
    //! @see @link example_biol_aa_sequence_phi_psi.cpp @endlink
    //! @author weinerbe
    //! @date Jan 25, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASequencePhiPsi :
      public coord::MovableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! N coordinates of middle residue
      linal::Vector3D m_N;

      //! Ca coordinates of middle residue
      linal::Vector3D m_CA;

      //! C coordinates of middle residue
      linal::Vector3D m_C;

      //! vector of phi-psi angles
      storage::Vector< storage::VectorND< 2, double> > m_Angles;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AASequencePhiPsi();

      //! @brief construct from an AASequence
      //! @param SEQUENCE sequence to be used
      AASequencePhiPsi( const AASequence &SEQUENCE);

      //! @brief construct from the 3 coordinates and the phi/psi angles
      //! @param N_COORDS coordinates for N atom of middle residue
      //! @param CA_COORDS coordinates for CA atom of middle residue
      //! @param C_COORDS coordinates for C atom of middle residue
      //! @param ANGLES vector of phi-psi angles
      AASequencePhiPsi
      (
        const linal::Vector3D &N_COORDS,
        const linal::Vector3D &CA_COORDS,
        const linal::Vector3D &C_COORDS,
        const storage::Vector< storage::VectorND< 2, double> > &ANGLES
      );

      //! @brief Clone function
      //! @return pointer to new AASequencePhiPsi
      AASequencePhiPsi *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the N coordinates of middle residue
      //! @return the N coordinates of middle residue
      const linal::Vector3D &GetN() const
      {
        return m_N;
      }

      //! @brief get the Ca coordinates of middle residue
      //! @return the Ca coordinates of middle residue
      const linal::Vector3D &GetCA() const
      {
        return m_CA;
      }

      //! @brief get the C coordinates of middle residue
      //! @return the C coordinates of middle residue
      const linal::Vector3D &GetC() const
      {
        return m_C;
      }

      //! @brief get the phi/psi angles
      //! @return the phi/psi angles
      const storage::Vector< storage::VectorND< 2, double> > &GetAngles() const
      {
        return m_Angles;
      }

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const;

      //! @brief returns whether the object is defined (first phi and last psi can be nan but everything else cannot)
      //! @return whether the object is defined
      bool IsDefined() const;

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

      //! @brief initializes member variables from an AASequence
      //! @param SEQUENCE sequence to be used
      void InitializeAngles( const AASequence &SEQUENCE);

    }; // class AASequencePhiPsi

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_AA_SEQUENCE_PHI_PSI_H_ 
