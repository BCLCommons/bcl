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

#ifndef BCL_FOLD_LOOP_PARAMETERS_H_
#define BCL_FOLD_LOOP_PARAMETERS_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopParameters
    //! @brief Parameterizes loop conformations corresponding to their euclidean and sequence distance, rotation
    //! angles and dihedral angles.
    //!
    //! @see @link example_fold_loop_parameters.cpp @endlink
    //! @author fischea
    //! @date Dec 16, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LoopParameters :
      public util::SerializableInterface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! translational span of the loop
      linal::Vector3D m_Translation;

      //! rotation angles of the second anchor with respect to the first anchor
      linal::Vector3D m_RotationAngles;

      //! sequence ids of the anchor points of the loop
      storage::Vector< int> m_Anchors;

      //! chain id of the loop
      char m_ChainID;

      //! dihedral angles of the loop
      storage::Vector< double> m_DihedralAngles;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      LoopParameters();

      //! @brief construct from members
      //! @param TRANSLATION translational span of the loop
      //! @param ROTATION rotation angles of the second anchor with respect to the first anchor
      //! @param ANCHORS sequence ids of the anchor points of the loop
      //! @param CHAIN_ID chain id of the loop
      //! @param ANGLES dihedral angles of the loop
      LoopParameters
      (
        const linal::Vector3D &TRANSLATION,
        const linal::Vector3D &ROTATION,
        const storage::Vector< int> &ANCHORS,
        const char &CHAIN_ID,
        const storage::Vector< double> &ANGLES
      );

      //! @brief construct from sequence and chain IDs
      //! @param ANCHORS sequence IDs of the anchors
      //! @param CHAIN_ID chain ID of the loop
      LoopParameters( const storage::Vector< int> &ANCHORS, const char &CHAIN_ID) :
        m_Anchors( ANCHORS),
        m_ChainID( CHAIN_ID)
      {
      }

      //! @brief clone function
      //! @return pointer to a new LoopParameters
      LoopParameters *Clone() const;

      //! @brief creates a parameterization of a loop from the given anchor points and dihedral angles
      //! @param FIRST_ANCHOR first anchor point of the loop
      //! @param SECOND_ANCHOR second anchor point of the loop
      //! @param ANGLES dihedral angles of the loop
      //! @return parameterization of the loop defined by the given arguments
      static util::ShPtr< LoopParameters> Create
      (
        const biol::AABase &FIRST_ANCHOR,
        const biol::AABase &SECOND_ANCHOR,
        const storage::Vector< double> &ANGLES = storage::Vector< double>()
      );

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the local coordinate system of the given amino acid
      //! @param AMINO_ACID amino aid to return the local coordinate system for
      //! @return local coordinate system of the given amino acid
      static storage::VectorND< 3, linal::Vector3D> GetLocalCoordinateSystem( const biol::AABase &AMINO_ACID);

      //! @brief computes the Euler angles between the given frames according to z-x'-z'' convention
      //! @param FRAME_A the global coordinate system with axes being the row vectors
      //! @param FRAME_B the local coordinate system with axes being the row vectors
      //! @return Euler angles between the coordinate systems according to z-x'-z'' convention
      static linal::Vector3D GetEulerAngles
      (
        const storage::VectorND< 3, linal::Vector3D> &FRAME_A,
        const storage::VectorND< 3, linal::Vector3D> &FRAME_B
      );

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the translational span of the loop
      //! @return translational span of the loop
      const linal::Vector3D &GetTranslation() const;

      //! @brief returns the rotation angles of the second anchor with respect to the first anchor
      //! @return rotation angles of the second anchor with respect to the first anchor
      const linal::Vector3D &GetRotation() const;

      //! @brief returns the sequence ids of the anchor points of the loop
      //! @return sequence ids of the anchor points of the loop
      const storage::Vector< int> &GetAnchors() const;

      //! @brief returns the chain of the anchor points of the loop
      //! @return chain id of the anchor points of the loop
      const char &GetChainID() const;

      //! @brief returns the sequence distance spanned by the loop
      //! @return sequence distance spanned by the loop
      size_t GetSequenceDistance() const;

      //! @brief returns the dihedral angles of the loop
      //! @return dihedral angles of the loop
      const storage::Vector< double> &GetAngles() const;

      //! @brief returns the sequence direction of the loop
      //! @return sequence direction of the loop
      biol::AASequenceFlexibility::SequenceDirection GetSequenceDirection() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class LoopParameters

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOOP_PARAMETERS_H_
