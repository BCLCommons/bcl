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
#include "biol/bcl_biol_aa_sequence_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AASequencePhiPsi::s_Instance
    (
      GetObjectInstances().AddInstance( new AASequencePhiPsi())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASequencePhiPsi::AASequencePhiPsi() :
      m_N(),
      m_CA(),
      m_C(),
      m_Angles()
    {
    }

    //! @brief construct from an AASequence
    //! @param SEQUENCE sequence to be used
    AASequencePhiPsi::AASequencePhiPsi( const AASequence &SEQUENCE) :
      m_N(),
      m_CA(),
      m_C(),
      m_Angles()
    {
      InitializeAngles( SEQUENCE);
    }

    //! @brief construct from the 3 coordinates and the phi/psi angles
    //! @param N_COORDS coordinates for N atom of middle residue
    //! @param CA_COORDS coordinates for CA atom of middle residue
    //! @param C_COORDS coordinates for C atom of middle residue
    //! @param ANGLES vector of phi-psi angles
    AASequencePhiPsi::AASequencePhiPsi
    (
      const linal::Vector3D &N_COORDS,
      const linal::Vector3D &CA_COORDS,
      const linal::Vector3D &C_COORDS,
      const storage::Vector< storage::VectorND< 2, double> > &ANGLES
    ) :
      m_N( N_COORDS),
      m_CA( CA_COORDS),
      m_C( C_COORDS),
      m_Angles( ANGLES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AASequencePhiPsi
    AASequencePhiPsi *AASequencePhiPsi::Clone() const
    {
      return new AASequencePhiPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequencePhiPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the geometric center of the object
    //! @return the geometric center of the object
    linal::Vector3D AASequencePhiPsi::GetCenter() const
    {
      return m_CA;
    }

    //! @brief returns whether the object is defined (first phi and last psi can be nan but everything else cannot)
    //! @return whether the object is defined
    bool AASequencePhiPsi::IsDefined() const
    {
      // if any of the coords are undefined
      if( !m_N.IsDefined() || !m_CA.IsDefined() || !m_C.IsDefined())
      {
        // return false
        return false;
      }

      // check if the first element psi value is undefined
      if( !util::IsDefined( m_Angles.FirstElement().Second()))
      {
        return false;
      }

      // check if the last element phi value is undefined
      if( !util::IsDefined( m_Angles.LastElement().First()))
      {
        return false;
      }

      // if the # of angles is greater than 2
      if( m_Angles.GetSize() > 2)
      {
        // get iterators
        storage::Vector< storage::VectorND< 2, double> >::const_iterator angle_itr( m_Angles.Begin());
        storage::Vector< storage::VectorND< 2, double> >::const_iterator angle_itr_end( m_Angles.End());

        // already checked first and last elements so skip them
        ++angle_itr;
        --angle_itr_end;

        // iterate over the remaining angles
        for( ; angle_itr != angle_itr_end; ++angle_itr)
        {
          // if either value is undefined
          if( !util::IsDefined( angle_itr->First()) || !util::IsDefined( angle_itr->Second()))
          {
            // return false
            return false;
          }
        }
      }

      // if this point is reached all values are defined
      return true;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Translation to be applied
    void AASequencePhiPsi::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_N.Translate( TRANSLATION);
      m_CA.Translate( TRANSLATION);
      m_C.Translate( TRANSLATION);
    }

    //! @brief transform the object by a given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AASequencePhiPsi::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      m_N.Transform( TRANSFORMATION_MATRIX_3D);
      m_CA.Transform( TRANSFORMATION_MATRIX_3D);
      m_C.Transform( TRANSFORMATION_MATRIX_3D);
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void AASequencePhiPsi::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      m_N.Rotate( ROTATION_MATRIX_3D);
      m_CA.Rotate( ROTATION_MATRIX_3D);
      m_C.Rotate( ROTATION_MATRIX_3D);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASequencePhiPsi::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_N, ISTREAM);
      io::Serialize::Read( m_CA, ISTREAM);
      io::Serialize::Read( m_C, ISTREAM);
      io::Serialize::Read( m_Angles, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AASequencePhiPsi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_N, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_C, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Angles, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initializes member variables from an AASequence
    //! @param SEQUENCE sequence to be used
    void AASequencePhiPsi::InitializeAngles( const AASequence &SEQUENCE)
    {
      // if the sequence is empty or not continuous
      if( SEQUENCE.GetData().IsEmpty())
      {
        // warn the user and return
        BCL_MessageStd
        (
          "Given sequence is empty, so no phi-psi angles will be determined: " +
            SEQUENCE.GetSequenceIdentification()
        );
        return;
      }

      // get the middle index
      const size_t middle_index( SEQUENCE.GetSize() / 2);

      // get the coordinates of the first amino acid
      m_N = SEQUENCE.GetData()( middle_index)->GetAtom( GetAtomTypes().N).GetCoordinates();
      m_CA = SEQUENCE.GetData()( middle_index)->GetAtom( GetAtomTypes().CA).GetCoordinates();
      m_C = SEQUENCE.GetData()( middle_index)->GetAtom( GetAtomTypes().C).GetCoordinates();

      // iterate through the AAs
      for
      (
        AASequence::const_iterator aa_itr( SEQUENCE.GetData().Begin()),
          aa_itr_end( SEQUENCE.GetData().End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // get the phi/psi angles
        m_Angles.PushBack( SEQUENCE.CalculatePhiPsi( ( *aa_itr)->GetSeqID()));
      }
    }

  } // namespace biol

} // namespace bcl
