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

// include header of this class
#include "fold/bcl_fold_mutate_protein_model_chain_move.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_vector.h"

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< util::ObjectInterface> MutateProteinModelChainMove::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelChainMove())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelChainMove::MutateProteinModelChainMove() :
      m_MaxTranslationXY( 5.0),
      m_MaxTranslationZ( 5.0),
      m_MaxRotationAngle( 0.0),
      m_MaxInternalRotationAngle( 0.0)
    {
      // nothing else to do
    }

    //!
    //!
    //! @param MAX_TRANSLATION_XY
    //! @param MAX_TRANSLATION_Z
    //! @param MAX_ROTATION_ANGLE
    //!
    MutateProteinModelChainMove::MutateProteinModelChainMove
    (
      const double MAX_TRANSLATION_XY,
      const double MAX_TRANSLATION_Z,
      const double MAX_ROTATION_ANGLE,
      const double MAX_INTERNAL_ROTATION_ANGLE
    ) :
      m_MaxTranslationXY( MAX_TRANSLATION_XY),
      m_MaxTranslationZ( MAX_TRANSLATION_Z),
      m_MaxRotationAngle( MAX_ROTATION_ANGLE),
      m_MaxInternalRotationAngle( MAX_INTERNAL_ROTATION_ANGLE)
    {
      // nothing else to do
    }

    //!
    //! @brief Clone function
    //! @return a pointer to a copy of this object
    //!
    MutateProteinModelChainMove *MutateProteinModelChainMove::Clone() const
    {
      return new MutateProteinModelChainMove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //!
    //! @brief Get the name of this class
    //! @return the name of this class
    //!
    const std::string &MutateProteinModelChainMove::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateProteinModelChainMove::GetAlias() const
    {
      static const std::string s_name( "MutateProteinModelChainMove");
      return s_name;
    }

    //!
    //!
    //! @return
    //!
    io::Serializer MutateProteinModelChainMove::GetSerializer() const
    {
      io::Serializer serializer;

      serializer.AddInitializer
      (
        "max translation xy",
        "allowed maximum translation in x, y direction for each step of mutate",
        io::Serialization::GetAgent( &m_MaxTranslationXY)
      );

      serializer.AddInitializer
      (
        "max translation z",
        "allowed maximum translation in z direction for each step of mutate",
        io::Serialization::GetAgent( &m_MaxTranslationZ)
      );

      serializer.AddInitializer
      (
        "max rotation angle",
        "allowed maximum rotation angle for each step of mutate",
        io::Serialization::GetAgent( &m_MaxRotationAngle)
      );

      serializer.AddInitializer
      (
        "max internal rotation angle",
        "allowed maximum angle for internal rotation of the ligand",
        io::Serialization::GetAgent( &m_MaxInternalRotationAngle)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //!
    //!
    //! @param PROTEIN_MODEL
    //! @return
    //!
    char MutateProteinModelChainMove::GetLigandChainID( const assemble::ProteinModel &COMPLEX)
    {
      // char that holds the ID of the movable chain
      char movable_chain_id( COMPLEX.GetChainIDs()[ 0]);

      // get all chains of this protein model
      util::ShPtrVector< assemble::Chain> all_chains( COMPLEX.GetChains());

      // determine the movable chain, which is the smallest chain
      size_t min_chain_size( COMPLEX.GetNumberAAs());
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          itr( all_chains.Begin()), itr_end( all_chains.End());
        itr != itr_end;
        ++itr
      )
      {
        size_t current_chain_size( ( *itr)->GetNumberAAs());
        if( current_chain_size < min_chain_size)
        {
          min_chain_size = current_chain_size;
          movable_chain_id = ( *itr)->GetChainID();
        }
      }

      // return the ID of the movable chain
      return movable_chain_id;
    }

    //!
    //! @brief Generate placement for the given Chain at a random orientation with respect to a located Chain in the model
    //! @param MOVABLE_CHAIN chain to be moved
    //! @param PROTEIN_MODEL the given protein model to which the chain is to be added
    //! @return
    //!
    storage::Pair< math::TransformationMatrix3D, bool> MutateProteinModelChainMove::PlaceLigand
    (
      const assemble::Chain &LIGAND,
      const assemble::ProteinModel &COMPLEX
    ) const
    {

      // vector for holding fixed chains
      util::ShPtrVector< assemble::Chain> receptor_chains;

      // add fixed chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          itr( COMPLEX.GetChains().Begin()), itr_end( COMPLEX.GetChains().End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetChainID() != LIGAND.GetChainID())
        {
          receptor_chains.PushBack( *itr);
        }
      }

      //
      util::ShPtr< assemble::ProteinModel> receptor( new assemble::ProteinModel( receptor_chains));

      //! @TODO
      BCL_Debug( receptor->GetCenterOfMass());

      //
//        math::TransformationMatrix3D transformation
//        (
//          coord::OrientationInterface::GenerateRandomTransformationAroundCenter
//          (
//            m_MaxTranslationXY,
//            m_MaxRotationAngle,
//            receptor->GetCenterOfMass()
//          )
//        );

      math::RotationMatrix3D rotation
      (
        math::RotationMatrix3D().SetRand( m_MaxRotationAngle)
      );

      math::TransformationMatrix3D transformation( rotation);

      BCL_Debug( transformation);

      //
      double translation_xy( static_cast< double>( std::rand()) / static_cast< double>( RAND_MAX) * m_MaxTranslationXY);

      linal::Vector3D translation( translation_xy, translation_xy, 0.0);
      transformation.SetTranslation( translation);

      BCL_Debug( transformation);

      //

      //
      // transformation( math::RotationMatrix3D().SetRand( m_MaxRotationAngle));
      // transformation( linal::Vector3D().SetRandomTranslation( linal::Vector3D( m_MaxTranslationXY, m_MaxTranslationXY, m_MaxTranslationZ)));

      //
      // transformation( sp_fixed_part->GetOrientation());

      //
      return storage::Pair< math::TransformationMatrix3D, bool>( transformation, true);
    }

    //! @brief
    //! @param COMPLEX
    void MutateProteinModelChainMove::PreDock( assemble::ProteinModel &COMPLEX, bool MEMBRANE)
    {
      // get ligand
      assemble::ProteinModel ligand( GetLigand( COMPLEX));

      // get receptor
      assemble::ProteinModel receptor( GetReceptor( COMPLEX));

      // if is membrane protein complex
      if( MEMBRANE)
      {
        // move the receptor and the ligand into the membrane
        util::ShPtr< biol::Membrane> sp_membrane
        (
          COMPLEX.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
        );

        // if the membrane is defined, move the ligand center to the membrane center
        if( sp_membrane.IsDefined())
        {
          // calculate ligand geometric center
          linal::Vector3D ligand_center( coord::CenterOfMass( ligand.GetAtomCoordinates(), true));

          // move the ligand to the origin
          ligand.Translate( sp_membrane->GetOrientation().GetOrigin() - ligand_center + linal::Vector3D( 100, 100, 0));
        }
      }

      // translate the ligand onto the receptor so that they become in contact
      ligand.Translate( ComputeReceptorLigandShortestDistance( receptor, ligand).second);

//      util::ShPtr< assemble::ProteinModel> sp_model( COMPLEX.HardCopy());
//      sp_model->Replace( sp_ligand);
//      COMPLEX = *sp_model;
    }

    //!
    //!
    //! @param PROTEIN_MODEL
    //! @return
    //!
    math::MutateResult< assemble::ProteinModel> MutateProteinModelChainMove::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {

      // get the ID of the movable chain
      const char ligand_chain_id( GetLigandChainID( PROTEIN_MODEL));

      // get the movable chain
      const util::ShPtr< assemble::Chain> sp_ligand( PROTEIN_MODEL.GetChain( ligand_chain_id));

      // get the placement for the movable chain
      const storage::Pair< math::TransformationMatrix3D, bool> placement( PlaceLigand( *sp_ligand, PROTEIN_MODEL));

      // return if there is no placement defined
      if( !placement.Second())
      {
        return math::MutateResult< assemble::ProteinModel>( util::ShPtr< assemble::ProteinModel>(), *this);
      }

      math::TransformationMatrix3D transform( placement.First());

      // a copy of PROTEIN_MODEL
      util::ShPtr< assemble::ProteinModel> sp_model( PROTEIN_MODEL.HardCopy());

      // create a copy of the movable chain
      util::ShPtr< assemble::Chain> sp_ligand_copy( sp_ligand->HardCopy());

      // move the ligand to a new location
      sp_ligand_copy->Transform( transform);

      // rotate the ligand internally
      const linal::Vector3D ligand_com( coord::CenterOfMass( sp_ligand_copy->GetAtomCoordinates(), true));
      math::TransformationMatrix3D internal_rotation( -ligand_com);
      internal_rotation( math::RotationMatrix3D().SetRand( m_MaxInternalRotationAngle));
      internal_rotation( ligand_com);
      sp_ligand_copy->Transform( internal_rotation);

      // replace movable chain
      sp_model->Replace( sp_ligand_copy);

      // return the mutated model
      return math::MutateResult< assemble::ProteinModel>( sp_model, *this);
    }

    //!
    //!
    //! @param COMPLEX
    //! @return
    //!
    const util::ShPtr< assemble::Chain> &MutateProteinModelChainMove::GetLigand
    (
      const assemble::ProteinModel &COMPLEX
    )
    {
      return COMPLEX.GetChain( GetLigandChainID( COMPLEX));
    }

    //!
    //!
    //! @param COMPLEX
    //! @return
    //!
    util::ShPtrVector< assemble::Chain> MutateProteinModelChainMove::GetReceptor
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    )
    {
      // get the receptor chain id(s)
      std::string chain_ids( PROTEIN_MODEL.GetChainIDs());
      size_t ligand_chain_id_pos( chain_ids.find( GetLigandChainID( PROTEIN_MODEL)));
      chain_ids.erase( ligand_chain_id_pos, 1);

      // return the receptor chains
      return PROTEIN_MODEL.GetChains( chain_ids);
    }

    //! @brief computes the shortest distance between the ligand chain and all receptor chain(s)
    //! @param LIGAND the ligand chain
    //! @param RECEPTOR the receptor chain(s)
    //! @return
    std::pair< double, linal::Vector3D> MutateProteinModelChainMove::ComputeReceptorLigandShortestDistance
    (
      const assemble::ProteinModel &RECEPTOR,
      const assemble::ProteinModel &LIGAND
    )
    {
      // get all ligand atoms and all receptor atoms
      util::SiPtrVector< const biol::Atom> lig_atoms( LIGAND.GetAtoms());
      util::SiPtrVector< const biol::Atom> rec_atoms( RECEPTOR.GetAtoms());

      // maximum possible double values
      double shortest_distance
      (
        linal::Distance( ( *lig_atoms.Begin())->GetCoordinates(), ( *rec_atoms.Begin())->GetCoordinates())
      );
      linal::Vector3D shortest_distance_vector3d;

      // iterate over ligand atoms
      for
      (
        auto lig_atom_itr( lig_atoms.Begin()), lig_atom_itr_end( lig_atoms.End());
          lig_atom_itr != lig_atom_itr_end;
        ++lig_atom_itr
      )
      {
        // iterate over atoms in the receptor
        for
        (
          auto rec_atom_itr( rec_atoms.Begin()), rec_atom_itr_end( rec_atoms.End());
            rec_atom_itr != rec_atom_itr_end;
          ++rec_atom_itr
        )
        {
          // current distance
          double current_distance
          (
            linal::Distance( ( *lig_atom_itr)->GetCoordinates(), ( *rec_atom_itr)->GetCoordinates())
          );

          // update shortest distance
          if( current_distance < shortest_distance)
          {
            shortest_distance = current_distance;
            shortest_distance_vector3d = ( *rec_atom_itr)->GetCoordinates() - ( *lig_atom_itr)->GetCoordinates();
          }
        }
      }

      return std::pair< double, linal::Vector3D>( shortest_distance, shortest_distance_vector3d);
    }
  } // namespace fold
} // namespace bcl
