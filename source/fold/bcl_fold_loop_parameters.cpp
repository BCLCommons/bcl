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
#include "fold/bcl_fold_loop_parameters.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> LoopParameters::s_Instance
    (
      GetObjectInstances().AddInstance( new LoopParameters())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopParameters::LoopParameters() :
      m_Translation(),
      m_RotationAngles(),
      m_Anchors(),
      m_ChainID(),
      m_DihedralAngles()
    {
    }

    //! @brief construct from members
    //! @param TRANSLATION translational span of the loop
    //! @param ROTATION rotation angles of the second anchor with respect to the first anchor
    //! @param ANCHORS sequence ids of the anchor points of the loop
    //! @param CHAIN_ID chain id of the loop
    //! @param ANGLES dihedral angles of the loop
    LoopParameters::LoopParameters
    (
      const linal::Vector3D &TRANSLATION,
      const linal::Vector3D &ROTATION,
      const storage::Vector< int> &ANCHORS,
      const char &CHAIN_ID,
      const storage::Vector< double> &ANGLES
    ) :
      m_Translation( TRANSLATION),
      m_RotationAngles( ROTATION),
      m_Anchors( ANCHORS),
      m_ChainID( CHAIN_ID),
      m_DihedralAngles( ANGLES)
    {
      // make sure that the correct number of dihedral angles was provided for this length of the loop
      BCL_Assert
      (
        m_DihedralAngles.IsEmpty() || m_DihedralAngles.GetSize() == 2 * GetSequenceDistance() + 2,
        "Number of dihedral angles should either be 0 or two times the sequence distance plus two angles."
      );
    }

    //! @brief clone function
    //! @return pointer to a new LoopParameters
    LoopParameters *LoopParameters::Clone() const
    {
      return new LoopParameters( *this);
    }

    //! @brief creates a parameterization of a loop from the given anchor points and dihedral angles
    //! @param FIRST_ANCHOR first anchor point of the loop
    //! @param SECOND_ANCHOR second anchor point of the loop
    //! @param ANGLES dihedral angles of the loop
    //! @return parameterization of the loop defined by the given arguments
    util::ShPtr< LoopParameters> LoopParameters::Create
    (
      const biol::AABase &FIRST_ANCHOR,
      const biol::AABase &SECOND_ANCHOR,
      const storage::Vector< double> &ANGLES
    )
    {
      // compute the local coordinate system of the first anchor
      const storage::VectorND< 3, linal::Vector3D> coord_system( GetLocalCoordinateSystem( FIRST_ANCHOR));

      // compute the translation
      const linal::Vector3D translation_global
      (
        SECOND_ANCHOR.GetCA().GetCoordinates() - FIRST_ANCHOR.GetCA().GetCoordinates()
      );
      const linal::Vector3D translation_local
      (
        translation_global * coord_system.First(),  // local x-coordinate
        translation_global * coord_system.Second(), // local y-coordinate
        translation_global * coord_system.Third()   // local z-coordinate
      );

      // compute the rotation angles for the local coordinate system
      const storage::VectorND< 3, linal::Vector3D> coord_system_2( GetLocalCoordinateSystem( SECOND_ANCHOR));
      const linal::Vector3D rotation_angles( GetEulerAngles( coord_system, coord_system_2));

      // get the sequence ids of the anchors
      storage::Vector< int> anchor_ids;
      anchor_ids.PushBack( FIRST_ANCHOR.GetSeqID());
      anchor_ids.PushBack( SECOND_ANCHOR.GetSeqID());

      // get the chain id of the anchors
      const char &chain_id( FIRST_ANCHOR.GetChainID());

      // create the parameterization wrapper
      util::ShPtr< LoopParameters> sp_param
      (
        new LoopParameters( translation_local, rotation_angles, anchor_ids, chain_id, ANGLES)
      );

      return sp_param;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &LoopParameters::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the translational span of the loop
    //! @return translational span of the loop
    const linal::Vector3D &LoopParameters::GetTranslation() const
    {
      return m_Translation;
    }

    //! @brief returns the rotation angles of the second anchor with respect to the first anchor
    //! @return rotation angles of the second anchor with respect to the first anchor
    const linal::Vector3D &LoopParameters::GetRotation() const
    {
      return m_RotationAngles;
    }

    //! @brief returns the sequence ids of the anchor points of the loop
    //! @return sequence ids of the anchor points of the loop
    const storage::Vector< int> &LoopParameters::GetAnchors() const
    {
      return m_Anchors;
    }

    //! @brief returns the chain of the anchor points of the loop
    //! @return chain id of the anchor points of the loop
    const char &LoopParameters::GetChainID() const
    {
      return m_ChainID;
    }

    //! @brief returns the sequence distance spanned by the loop
    //! @return sequence distance spanned by the loop
    size_t LoopParameters::GetSequenceDistance() const
    {
      return std::abs( m_Anchors( 1) - m_Anchors( 0)) - 1;
    }

    //! @brief returns the dihedral angles of the loop
    //! @return dihedral angles of the loop
    const storage::Vector< double> &LoopParameters::GetAngles() const
    {
      return m_DihedralAngles;
    }

    //! @brief returns the sequence direction of the loop
    //! @return sequence direction of the loop
    biol::AASequenceFlexibility::SequenceDirection LoopParameters::GetSequenceDirection() const
    {
      // determine the sequence direction based of the sequence ids of the anchors
      if( m_Anchors( 0) < m_Anchors( 1))
      {
        return biol::AASequenceFlexibility::e_CTerminal;
      }
      else
      {
        return biol::AASequenceFlexibility::e_NTerminal;
      }
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &LoopParameters::GetAlias() const
    {
      static const std::string s_name( "LoopParameters");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopParameters::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Stores the conformation and orientation of a loop.");
      serializer.AddInitializer
      (
        "translation",
        "translation vector between the anchor residues of the loop",
        io::Serialization::GetAgent( &m_Translation)
      );
      serializer.AddInitializer
      (
        "rotation",
        "rotation angles in extrinsic x-y-z convention between the anchor residues of the loop",
        io::Serialization::GetAgent( &m_RotationAngles)
      );
      serializer.AddInitializer
      (
        "anchor ids",
        "sequence id of the loop's anchor residues",
        io::Serialization::GetAgentWithSizeLimits( &m_Anchors, 2, 2)
      );
      serializer.AddInitializer
      (
        "chain id",
        "chain id of the loop",
        io::Serialization::GetAgent( &m_ChainID)
      );
      serializer.AddInitializer
      (
        "dihedral angles",
        "dihedral angles of the loop including n-terminal psi-angle and c-terminal phi-angle",
        io::Serialization::GetAgent( &m_DihedralAngles)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns the local coordinate system of the given amino acid
    //! @param AMINO_ACID amino aid to return the local coordinate system for
    //! @return local coordinate system of the given amino acid
    storage::VectorND< 3, linal::Vector3D> LoopParameters::GetLocalCoordinateSystem
    (
      const biol::AABase &AMINO_ACID
    )
    {
      // get the coordinates to compute the local coordinate system
      const linal::Vector3D &ca_coord( AMINO_ACID.GetAtom( biol::GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D &c_coord( AMINO_ACID.GetAtom( biol::GetAtomTypes().C).GetCoordinates());
      const linal::Vector3D &o_coord( AMINO_ACID.GetAtom( biol::GetAtomTypes().O).GetCoordinates());

      // compute the orthonormal base of the local coordinate system
      const linal::Vector3D z_axis( ( c_coord - ca_coord).Normalize());
      const linal::Vector3D ca_o( o_coord - ca_coord);
      const linal::Vector3D y_axis( ( ca_o - ( ca_o * z_axis) * z_axis).Normalize());
      const linal::Vector3D x_axis( linal::CrossProduct( z_axis, y_axis));
      const storage::VectorND< 3, linal::Vector3D> base( x_axis, y_axis, z_axis);

      return base;
    }

    //! @brief computes the Euler angles between the given frames according to z-x'-z'' convention
    //! @param FRAME_A the global coordinate system with axes being the row vectors
    //! @param FRAME_B the local coordinate system with axes being the row vectors
    //! @return Euler angles between the coordinate systems according to z-x'-z'' convention
    linal::Vector3D LoopParameters::GetEulerAngles
    (
      const storage::VectorND< 3, linal::Vector3D> &FRAME_A,
      const storage::VectorND< 3, linal::Vector3D> &FRAME_B
    )
    {
      // convert the coordinate systems into matrices
      linal::Matrix3x3< double> frame_a( 0.0);
      frame_a.ReplaceRow( 0, FRAME_A.First());
      frame_a.ReplaceRow( 1, FRAME_A.Second());
      frame_a.ReplaceRow( 2, FRAME_A.Third());
      linal::Matrix3x3< double> frame_b( 0.0);
      frame_b.ReplaceRow( 0, FRAME_B.First());
      frame_b.ReplaceRow( 1, FRAME_B.Second());
      frame_b.ReplaceRow( 2, FRAME_B.Third());

      // compute the Euler angles between the two frames
      return linal::Matrix3x3< double>::ComputeEulerAnglesXYZ( frame_a, frame_b);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
