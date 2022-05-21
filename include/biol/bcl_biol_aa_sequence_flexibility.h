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

#ifndef BCL_BIOL_AA_SEQUENCE_FLEXIBILITY_H_
#define BCL_BIOL_AA_SEQUENCE_FLEXIBILITY_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASequenceFlexibility
    //! @brief class for bending AASequence by changing or setting phi/psi angles on amino acids
    //! @details This class provides various functions for setting or changing phi/psi angles of amino acids in sequence
    //! while also making sure the rest of the sequence is also correctly transformed. The enumerator BendingTypes
    //! determines whether the changes should be applied toward N-terminal, C-terminal or towards both directions
    //!
    //! @see @link example_biol_aa_sequence_flexibility.cpp @endlink
    //! @author karakam
    //! @date Jan 22, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASequenceFlexibility :
      public util::ObjectInterface
    {
    public:

    //////////
    // data //
    //////////

      //! enumerator for sequence direction
      enum SequenceDirection
      {
        e_NTerminal,
        e_CTerminal,
        e_Bidirectional,
        s_NumberSequenceDirections
      };

      //! @brief conversion to a string from a SequenceDirection
      //! @param SEQUENCE_DIRECTION the sequence direction to get a string for
      //! @return a string representing that sequence direction
      static const std::string &GetSequenceDirectionName( const SequenceDirection &SEQUENCE_DIRECTION);

      //! @brief enum class wrapper for Unit
      typedef util::WrapperEnum< SequenceDirection, &GetSequenceDirectionName, s_NumberSequenceDirections> DirectionEnum;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AASequenceFlexibility
      AASequenceFlexibility *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief convenience function to calculate the number of different phi/psi between two SSEs
      //! @param SEQUENCE_A first sequence of interest
      //! @param SEQUENCE_B second sequence  of interest
      //! @return number of different phi/psi pairs between given sequences
      static size_t GetNumberDifferentPhiPsi
      (
        const AASequence &SEQUENCE_A,
        const AASequence &SEQUENCE_B
      );

      //! @brief function to calculate a Transformation from begin atom, end atom and the rotation angle
      //! @param ATOM_BEGIN Atom that forms the beginning point of the rotation axis
      //! @param ATOM_END Atom that forms the end point of the rotation axis
      //! @param ROTATION degrees of rotation
      static math::TransformationMatrix3D CalculateTransformation
      (
        const linal::Vector3D &ATOM_BEGIN,
        const linal::Vector3D &ATOM_END,
        const double ROTATION
      );

      //! @brief convenience function to rotate a specific atom of a residue with given transformation
      //! @param AMINO_ACID Amino acid of interest
      //! @param ATOM_TYPE AtomType of interest
      //! @param TRANSFORMATION Transformation to be applied
      static void TransformAtom
      (
        AABase &AMINO_ACID,
        const AtomType &ATOM_TYPE,
        const math::TransformationMatrix3D &TRANSFORMATION
      );

      //! @brief calculate the difference in phi and psi to be applied
      //! @param SEQUENCE Sequence to be bent
      //! @param SEQ_ID sequence id of the amino acid of interest
      //! @param PHI_PSI phi and psi angles to be set
      //! @return change in phi and psi angles to be applied
      static storage::VectorND< 2, double> CalculatePhiPsiChange
      (
        AASequence &SEQUENCE,
        const int SEQ_ID,
        const storage::VectorND< 2, double> &PHI_PSI
      );

      //! @brief set the phi psi of the given amino acid to given values
      //! @param SEQUENCE Sequence to be bent 
      //! @param SEQ_ID sequence id of the amino acid of interest
      //! @param PHI_PSI phi and psi angles to be set 
      //! @param SEQUENCE_DIRECTION enumerator whether the bending should be applied towards N-terminal, C-terminal or both directions
      //! @return transformations applied towards N-terminal and C-Terminal 
      static storage::VectorND< 2, math::TransformationMatrix3D> SetPhiPsi
      (
        AASequence &SEQUENCE,
        const int SEQ_ID,
        const storage::VectorND< 2, double> &PHI_PSI,
        const SequenceDirection &SEQUENCE_DIRECTION = e_Bidirectional
      );

      //! @brief change the phi psi of the given amino acid by the given values
      //! @param SEQUENCE Sequence to be bent
      //! @param SEQ_ID sequence id of the amino acid of interest
      //! @param PHI_PSI_CHANGE change in phi and psi angles to be applied
      //! @param SEQUENCE_DIRECTION enumerator whether the bending should be applied towards N-terminal, C-terminal or both directions
      //! @return transformations applied towards N-terminal and C-Terminal
      static storage::VectorND< 2, math::TransformationMatrix3D> ChangePhiPsi
      (
        AASequence &SEQUENCE,
        const int SEQ_ID,
        const storage::VectorND< 2, double> &PHI_PSI_CHANGE,
        const SequenceDirection &SEQUENCE_DIRECTION = e_Bidirectional
      );

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

    }; // class AASequenceFlexibility

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_AA_SEQUENCE_FLEXIBILITY_H_ 
