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
#include "biol/bcl_biol_aa_sequence_flexibility.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! @brief conversion to a string from a SequenceDirection
    //! @param SEQUENCE_DIRECTION the sequence direction to get a string for
    //! @return a string representing that sequence direction
    const std::string &AASequenceFlexibility::GetSequenceDirectionName
    (
      const AASequenceFlexibility::SequenceDirection &SEQUENCE_DIRECTION
    )
    {
      static const std::string s_descriptors[] =
      {
        "n_terminal",
        "c_terminal",
        "bidirectional",
        GetStaticClassName< SequenceDirection>()
      };
      return s_descriptors[ size_t( SEQUENCE_DIRECTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new AASequenceFlexibility
    AASequenceFlexibility *AASequenceFlexibility::Clone() const
    {
      return new AASequenceFlexibility( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequenceFlexibility::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief convenience function to calculate the number of different phi/psi between two SSEs
    //! @param SEQUENCE_A first sequence of interest
    //! @param SEQUENCE_B second sequence  of interest
    //! @return number of different phi/psi pairs between given sequences
    size_t AASequenceFlexibility::GetNumberDifferentPhiPsi
    (
      const AASequence &SEQUENCE_A,
      const AASequence &SEQUENCE_B
    )
    {
      // make sure they are the same SSEs
      BCL_Assert
      (
        SEQUENCE_A.GetSize() == SEQUENCE_B.GetSize(),
        "The sequences have to be same size for comparison!\nSEQUENCE_A " +
          util::Format()( SEQUENCE_A.GetSize()) + "\t" + SEQUENCE_A.Sequence() + "\nvs\nSEQUENCE_B " +
          util::Format()( SEQUENCE_B.GetSize())
      );

      // initialize count
      size_t count( 0);

      const int seqid_begin( SEQUENCE_A.GetFirstAA()->GetSeqID());
      const int seqid_end( SEQUENCE_B.GetLastAA()->GetSeqID());

      // iterate over seqids
      for( int seqid( seqid_begin); seqid <= seqid_end; ++seqid)
      {
        // calculate phi psi from both
        const storage::VectorND< 2, double> phi_psi_a( SEQUENCE_A.CalculatePhiPsi( seqid));
        const storage::VectorND< 2, double> phi_psi_b( SEQUENCE_B.CalculatePhiPsi( seqid));

        // if first residue do not check
        const bool phi_match
        (
          seqid == seqid_begin ? true : math::EqualWithinTolerance( phi_psi_a.First(), phi_psi_b.First())
        );

        const bool psi_match
        (
          seqid == seqid_end ? true : math::EqualWithinTolerance( phi_psi_a.Second(), phi_psi_b.Second())
        );

        // increment if any failed
        count += !( phi_match && psi_match);
      }

      // end
      return count;
    }

    //! @brief function to calculate a Transformation from begin atom, end atom and the rotation angle
    //! @param ATOM_BEGIN Atom that forms the beginning point of the rotation axis
    //! @param ATOM_END Atom that forms the end point of the rotation axis
    //! @param ROTATION degrees of rotation
    math::TransformationMatrix3D AASequenceFlexibility::CalculateTransformation
    (
      const linal::Vector3D &ATOM_BEGIN,
      const linal::Vector3D &ATOM_END,
      const double ROTATION
    )
    {
      // initialize the transformation so that you move the ATOM_BEGIN to origin
      math::TransformationMatrix3D transform( -ATOM_BEGIN);

      // rotate around the vector connecting ATOM_BEGIN to ATOM_END
      transform( math::RotationMatrix3D( linal::Vector3D( ATOM_END - ATOM_BEGIN), ROTATION));

      // move ATOM_BEGIN to its original position
      transform( ATOM_BEGIN);

      // end
      return transform;
    }

    //! @brief convenience function to rotate a specific atom of a residue with given transformation
    //! @param AMINO_ACID Amino acid of interest
    //! @param ATOM_TYPE AtomType of interest
    //! @param TRANSFORMATION Transformation to be applied
    void AASequenceFlexibility::TransformAtom
    (
      AABase &AMINO_ACID,
      const AtomType &ATOM_TYPE,
      const math::TransformationMatrix3D &TRANSFORMATION
    )
    {
      // get the correponding atom from this amino acid
      const Atom &atom( AMINO_ACID.GetAtom( ATOM_TYPE));

      // if the type is defined, thus the atom of the type ATOM_TYPE exists
      if( atom.GetType().IsDefined())
      {
        // make a copy of the atom and apply the psi transformation and put it back
        Atom new_atom( atom);
        new_atom.Transform( TRANSFORMATION);
        AMINO_ACID.SetAtom( new_atom);
      }
    }

    //! @brief calculate the difference in phi and psi to be applied
    //! @param SEQUENCE Sequence to be bent
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @param PHI_PSI phi and psi angles to be set
    //! @return change in phi and psi angles to be applied
    storage::VectorND< 2, double> AASequenceFlexibility::CalculatePhiPsiChange
    (
      AASequence &SEQUENCE,
      const int SEQ_ID,
      const storage::VectorND< 2, double> &PHI_PSI
    )
    {
      // calculate the current phi, psi values for the given amino acid
      const storage::VectorND< 2, double> current_phi_psi( SEQUENCE.CalculatePhiPsi( SEQ_ID));

      // calculate the differences
      return
        storage::VectorND< 2, double>
        (
          PHI_PSI.First() - current_phi_psi.First(), PHI_PSI.Second() - current_phi_psi.Second()
        );
    }

    //! @brief set the phi psi of the given amino acid to given values
    //! @param SEQUENCE Sequence to be bent
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @param PHI_PSI phi and psi angles to be set
    //! @param SEQUENCE_DIRECTION enumerator whether the bending should be applied towards N-terminal, C-terminal or both directions
    //! @return transformations applied towards N-terminal and C-Terminal
    storage::VectorND< 2, math::TransformationMatrix3D> AASequenceFlexibility::SetPhiPsi
    (
      AASequence &SEQUENCE,
      const int SEQ_ID,
      const storage::VectorND< 2, double> &PHI_PSI,
      const AASequenceFlexibility::SequenceDirection &SEQUENCE_DIRECTION
    )
    {
      // calculate the change and call the corresponding change phi psi function
      return ChangePhiPsi( SEQUENCE, SEQ_ID, CalculatePhiPsiChange( SEQUENCE, SEQ_ID, PHI_PSI), SEQUENCE_DIRECTION);
    }

    //! @brief change the phi psi of the given amino acid by the given values
    //! @param SEQUENCE Sequence to be bent
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @param PHI_PSI_CHANGE change in phi and psi angles to be applied
    //! @param SEQUENCE_DIRECTION enumerator whether the bending should be applied towards N-terminal, C-terminal or both directions
    //! @return transformations applied towards N-terminal and C-Terminal
    storage::VectorND< 2, math::TransformationMatrix3D> AASequenceFlexibility::ChangePhiPsi
    (
      AASequence &SEQUENCE,
      const int SEQ_ID,
      const storage::VectorND< 2, double> &PHI_PSI_CHANGE,
      const AASequenceFlexibility::SequenceDirection &SEQUENCE_DIRECTION
    )
    {
      // get the phi and psi change
      double phi_change( PHI_PSI_CHANGE.First());
      const double psi_change( PHI_PSI_CHANGE.Second());

      // if PHI or PSI is undefined return
      if( !util::IsDefined( phi_change) && !util::IsDefined( psi_change))
      {
        BCL_MessageDbg
        (
          "The provided PHI:" + util::Format()( phi_change) +
          " and PSI: " + util::Format()( psi_change) + " is undefined!"
        );
        return
          storage::VectorND< 2, math::TransformationMatrix3D>( util::GetUndefined< math::TransformationMatrix3D>());
      }

      // find the residue
      const AASequence::iterator aa_itr( SEQUENCE.FindAABySeqID( SEQ_ID));

      // make sure the AA is found
      BCL_Assert
      (
        aa_itr != SEQUENCE.End(),
        "No Amino acid with seqid " + util::Format()( SEQ_ID) + " in sequence " + SEQUENCE.GetSequenceIdentification()
      );

      // proline can not be rotated around phi
      if( ( *aa_itr)->GetType() == GetAATypes().PRO)
      {
        phi_change = 0.0;
      }

      // store reference to coordinates of nitrogen, ca and c for this amino acid
      const linal::Vector3D &coord_N( ( *aa_itr)->GetAtom( GetAtomTypes().N).GetCoordinates());
      const linal::Vector3D &coord_CA( ( *aa_itr)->GetCA().GetCoordinates());
      const linal::Vector3D &coord_C( ( *aa_itr)->GetAtom( GetAtomTypes().C).GetCoordinates());

      // initialize the transformations for phi and psi
      math::TransformationMatrix3D phi_transform, psi_transform;

      // if phi change is non-zero
      if( util::IsDefined( phi_change) && phi_change != double( 0.0))
      {
        // if towards N terminal or both directions
        if( SEQUENCE_DIRECTION == e_NTerminal || SEQUENCE_DIRECTION == e_Bidirectional)
        {
          // calculate phi transformation, the change has to be negated since we rotate residues beforehand
          phi_transform = CalculateTransformation( coord_N, coord_CA, phi_change);
          TransformAtom( **aa_itr, GetAtomTypes().H, phi_transform);
        }
        // else towards C terminal
        else
        {
          // calculate the phi transformation
          phi_transform = CalculateTransformation( coord_N, coord_CA, -phi_change);

          // transform this amino acid, keep hydrogen in place
          const Atom hydrogen( ( *aa_itr)->GetAtom( GetAtomTypes().H));
          ( *aa_itr)->Transform( phi_transform);
          if( hydrogen.GetType().IsDefined())
          {
            ( *aa_itr)->SetAtom( hydrogen);
          }
        }
      }

      // if psi change is non-zero
      if( util::IsDefined( psi_change) && psi_change != double( 0.0))
      {
        // if towards C terminal or both directions
        if( SEQUENCE_DIRECTION == e_CTerminal || SEQUENCE_DIRECTION == e_Bidirectional)
        {
          // calculate psi transformation and transform oxygen for this residue
          psi_transform = CalculateTransformation( coord_CA, coord_C, -psi_change);
          TransformAtom( **aa_itr, GetAtomTypes().O, psi_transform);
        }
        // else towards N terminal
        else
        {
          // calculate psi transformation the angle has to be negated since we rotate residues towards N-terminal
          psi_transform = CalculateTransformation( coord_CA, coord_C, psi_change);
          // transform the residue but keep oxygen in place
          const Atom oxygen( ( *aa_itr)->GetAtom( GetAtomTypes().O));
          ( *aa_itr)->Transform( psi_transform);
          // reset oxygen
          ( *aa_itr)->SetAtom( oxygen);
        }
      }

      // calculate the total transformation from phi and psi
      math::TransformationMatrix3D total_transformation( phi_transform);
      total_transformation( psi_transform);

      // construct the return transformation applied to N-terminal and C-terminal
      storage::VectorND< 2, math::TransformationMatrix3D> transformation_pair;

      // switch over bending types
      switch( SEQUENCE_DIRECTION)
      {
        case e_NTerminal:
        {
          // update the transformation pair
          transformation_pair.First() = total_transformation;

          // iterate over the N terminal part of the sequence
          for( AASequence::iterator itr( SEQUENCE.Begin()); itr != aa_itr; ++itr)
          {
            // apply the total transformation to each amino acid
            ( *itr)->Transform( total_transformation);
          }
          break;
        }
        case e_CTerminal:
        {
          // update the transformation pair
          transformation_pair.Second() = total_transformation;

          // iterate over the C terminal part of the sequence
          for( AASequence::iterator itr( aa_itr + 1), itr_end( SEQUENCE.End()); itr != itr_end; ++itr)
          {
            // apply the total transformation to each amino acid
            ( *itr)->Transform( total_transformation);
          }
          break;
        }
        case e_Bidirectional:
        {
          // update the transformation pair
          transformation_pair.First() = phi_transform;
          transformation_pair.Second() = psi_transform;

          // iterate over the N terminal part of the sequence
          for( AASequence::iterator itr( SEQUENCE.Begin()); itr != aa_itr; ++itr)
          {
            // apply the total transformation to each amino acid
            ( *itr)->Transform( phi_transform);
          }

          // iterate over the C terminal part of the sequence
          for( AASequence::iterator itr( aa_itr + 1), itr_end( SEQUENCE.End()); itr != itr_end; ++itr)
          {
            // apply the total transformation to each amino acid
            ( *itr)->Transform( psi_transform);
          }
          break;
        }
        case s_NumberSequenceDirections:
        {
          break;
        }
      }

      // return the pair of transformations applied
      return transformation_pair;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASequenceFlexibility::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AASequenceFlexibility::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
