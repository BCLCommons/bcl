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
#include "fold/bcl_fold_mutate_membrane_chain_move.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// include other forward headers - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    //! single instance of this class
    const util::SiPtr< util::ObjectInterface> MutateMembraneChainMove::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateMembraneChainMove())
    );

    //! @brief default constructor
    MutateMembraneChainMove::MutateMembraneChainMove() :
      m_ReceptorChainIDs(),
      m_LigandChainIDs(),
      m_MaxTranslationX( 1.0),
      m_MaxTranslationY( 1.0),
      m_MaxTranslationZ( 1.0),
      m_MaxPhi( 1.0),
      m_MaxTheta( 1.0),
      m_MaxPsi( 1.0)
    {
    }

    //! @brief explicit constructor from given arguments
    //! @param RECEPTOR_CHAIN_IDS chain IDs of the receptor
    //! @param LIGAND_CHAIN_IDS chain IDs of the ligand
    //! @param MAX_TRANSLATION_X maximally allowed single-step translation along the x axis
    //! @param MAX_TRANSLATION_Y maximally allowed single-step translation along the y axis
    //! @param MAX_TRANSLATION_Z maximally allowed single-step translation along the y axis
    //! @param MAX_PHI maximally allowed single-step phi rotation angle
    //! @param MAX_THETA maximally allowed single-step theta rotation angle
    //! @param MAX_PSI maximally allowed single-step psi rotation angle
    MutateMembraneChainMove::MutateMembraneChainMove
    (
      const std::string &RECEPTOR_CHAIN_IDS,
      const std::string &LIGAND_CHAIN_IDS,
      double MAX_TRANSLATION_X,
      double MAX_TRANSLATION_Y,
      double MAX_TRANSLATION_Z,
      double MAX_PHI,
      double MAX_THETA,
      double MAX_PSI
    ) :
      m_ReceptorChainIDs( RECEPTOR_CHAIN_IDS),
      m_LigandChainIDs( LIGAND_CHAIN_IDS),
      m_MaxTranslationX( MAX_TRANSLATION_X),
      m_MaxTranslationY( MAX_TRANSLATION_Y),
      m_MaxTranslationZ( MAX_TRANSLATION_Z),
      m_MaxPhi( MAX_PHI),
      m_MaxTheta( MAX_THETA),
      m_MaxPsi( MAX_PSI)
    {
    }

    //! @brief Clone function
    //! @return a pointer to a copy of this object
    MutateMembraneChainMove *MutateMembraneChainMove::Clone() const
    {
      return new MutateMembraneChainMove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the name of this class
    //! @return the name of this class
    const std::string &MutateMembraneChainMove::GetClassIdentifier() const
    {
      return GetStaticClassName( this);
    }

    //! @brief get the chain IDs of the receptor
    //! @return the chain IDs of the receptor
    const std::string &MutateMembraneChainMove::GetReceptorChainIDs() const
    {
      return m_ReceptorChainIDs;
    }

    //! @brief get the chain IDs of the ligand
    //! @return the chain IDs of the ligand
    const std::string &MutateMembraneChainMove::GetLigandChainIDs() const
    {
      return m_LigandChainIDs;
    }

    //! @brief get the maximum translation along the x axis
    //! @return maximum translation along the x axis
    double MutateMembraneChainMove::GetMaxTranslationX() const
    {
      return m_MaxTranslationX;
    }

    //! @brief get the maximum translation along the y axis
    //! @return maximum translation along the y axis
    double MutateMembraneChainMove::GetMaxTranslationY() const
    {
      return m_MaxTranslationY;
    }

    //! @brief get the maximum translation along the z axis
    //! @return maximum translation along the z axis
    double MutateMembraneChainMove::GetMaxTranslationZ() const
    {
      return m_MaxTranslationZ;
    }

    //! @brief get the maximum of the phi Euler angle
    //! @return maximum of the phi Euler angle
    double MutateMembraneChainMove::GetMaxPhi() const
    {
      return m_MaxPhi;
    }

    //! @brief get the maximum of the theta Euler angle
    //! @return maximum of the theta Euler angle
    double MutateMembraneChainMove::GetMaxTheta() const
    {
      return m_MaxTheta;
    }

    //! @brief get the maximum of the psi Euler angle
    //! @return maximum of the psi Euler angle
    double MutateMembraneChainMove::GetMaxPsi() const
    {
      return m_MaxPsi;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateMembraneChainMove::GetAlias() const
    {
      static const std::string s_alias( "MutateMembraneChainMove");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateMembraneChainMove::GetSerializer() const
    {
      io::Serializer serializer;

      // add initializers
      serializer.AddInitializer
      (
        "receptor chain ids",
        "chain IDs for the receptor",
        io::Serialization::GetAgent( &m_ReceptorChainIDs)
      );
      serializer.AddInitializer
      (
        "ligand chain ids",
        "chain IDs for the ligand",
        io::Serialization::GetAgent( &m_LigandChainIDs)
      );
      serializer.AddInitializer
      (
        "max translation x",
        "maximum angstroms of translation along the x axis",
        io::Serialization::GetAgent( &m_MaxTranslationX),
        "1.0"
      );
      serializer.AddInitializer
      (
        "max translation y",
        "maximum angstroms of translation along the y axis",
        io::Serialization::GetAgent( &m_MaxTranslationY),
        "1.0"
      );
      serializer.AddInitializer
      (
        "max translation z",
        "maximum angstroms of translation along the z axis",
        io::Serialization::GetAgent( &m_MaxTranslationZ),
        "1.0"
      );
      serializer.AddInitializer
      (
        "max phi",
        "maximum randians of the phi Euler angle",
        io::Serialization::GetAgent( &m_MaxPhi),
        "1.0"
      );
      serializer.AddInitializer
      (
        "max theta",
        "maximum randians of the theta Euler angle",
        io::Serialization::GetAgent( &m_MaxTheta),
        "1.0"
      );
      serializer.AddInitializer
      (
        "max psi",
        "maximum randians of the psi Euler angle",
        io::Serialization::GetAgent( &m_MaxPsi),
        "1.0"
      );

      // return the serializer
      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @param LIGAND
    //! @return
    math::MutateResult< assemble::ProteinModel>
      MutateMembraneChainMove::operator()( const assemble::ProteinModel &MODEL) const
    {
      // make sure that receptor chain IDs and ligand chain IDs are set
      if( m_LigandChainIDs.empty() || m_ReceptorChainIDs.empty())
      {
        SetReceptorLigandChainIDs( MODEL);
      }

      // create a copy of the ligand
      const assemble::ProteinModel ligand
      (
        MODEL.GetChains( m_LigandChainIDs)
      );

      // get receptor
      const assemble::ProteinModel receptor
      (
        MODEL.GetChains( m_ReceptorChainIDs)
      );

      // try to get the membrane
      util::ShPtr< assemble::ProteinModelData> sp_data( MODEL.GetProteinModelData());
      util::ShPtr< biol::Membrane> sp_membrane( sp_data->GetData( assemble::ProteinModelData::e_Membrane));

      // make sure that the membrane is defined
      if( !sp_membrane.IsDefined())
      {
        BCL_MessageCrt( "No defined membrane is found!");
        return math::MutateResult< assemble::ProteinModel>( util::ShPtr< assemble::ProteinModel>(), *this);
      }

      //
      double global_rotation_angle
      (
        random::GetGlobalRandom().Double( math::Range< double>( -math::g_Pi, math::g_Pi))
      );

      // random rotaion around z axis
      math::RotationMatrix3D rotaion_z( linal::Vector3D( 0.0, 0.0, 1.0), global_rotation_angle);

      // random translations
      double translation_x
      (
        random::GetGlobalRandom().Double( math::Range< double>( -m_MaxTranslationX, m_MaxTranslationX))
      );
      double translation_y
      (
        random::GetGlobalRandom().Double( math::Range< double>( -m_MaxTranslationY, m_MaxTranslationY))
      );
      double translation_z
      (
        random::GetGlobalRandom().Double( math::Range< double>( -m_MaxTranslationZ, m_MaxTranslationZ))
      );
      linal::Vector3D translation( translation_x, translation_y, translation_z);

      // transformation matrix
      math::TransformationMatrix3D transformation3d( rotaion_z);
      transformation3d.SetTranslation( translation);

      // do global transformation on the ligand
      // perform deep copy of ligand without affecting other objects that might have pointer to this ligand
      util::ShPtr< assemble::ProteinModel> sp_ligand_hard_copy( ligand.HardCopy());
      sp_ligand_hard_copy->Transform( transformation3d);

      linal::Vector3D ligand_center( sp_ligand_hard_copy->GetCenter());
      sp_ligand_hard_copy->Translate( -ligand_center);

      // local random rotations
      double phi
      (
        random::GetGlobalRandom().Double( math::Range< double>( -m_MaxPhi, m_MaxPhi))
      );
      double theta
      (
        random::GetGlobalRandom().Double( math::Range< double>( -m_MaxTheta, m_MaxTheta))
      );
      double psi
      (
        random::GetGlobalRandom().Double( math::Range< double>( -m_MaxPsi, m_MaxPsi))
      );
      math::RotationMatrix3D rotation3D( phi, theta, psi);

      // transformation matrix
      transformation3d.SetUnit();
      transformation3d.SetRotation( rotation3D);

      // transform the ligand
      sp_ligand_hard_copy->Transform( transformation3d);
      sp_ligand_hard_copy->Translate( ligand_center);

      // create a docked model
      util::ShPtrVector< assemble::Chain> all_chains( receptor.GetChains());
      all_chains.Append( sp_ligand_hard_copy->GetChains());
      util::ShPtr< assemble::ProteinModel> sp_docked_model( new assemble::ProteinModel( all_chains));
      sp_docked_model->SetProteinModelData( sp_data);

      // return the ligand
      return math::MutateResult< assemble::ProteinModel>( sp_docked_model, *this);

    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the chain IDs for the ligand
    void MutateMembraneChainMove::SetReceptorLigandChainIDs( const assemble::ProteinModel &MODEL) const
    {
      // ligand chain IDs
      std::string ligand_chain_ids;

      // get all chains of this protein model
      util::ShPtrVector< assemble::Chain> all_chains( MODEL.GetChains());

      // number of amino acids in the entire model
      size_t min_chain_size( MODEL.GetNumberAAs());

      // determine the movable chain, which is the smallest chain
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          itr( all_chains.Begin()), itr_end( all_chains.End());
        itr != itr_end;
        ++itr
      )
      {
        // number of amino acids in current chain
        size_t current_chain_size( ( *itr)->GetNumberAAs());
        if( current_chain_size < min_chain_size)
        {
          min_chain_size = current_chain_size;
          ligand_chain_ids = ( *itr)->GetChainID();
        }
      }

      // set ligand chain IDs
      m_LigandChainIDs = ligand_chain_ids;

      // get all chain IDs
      std::string all_chain_ids( MODEL.GetChainIDs());

      // remove ligand chain IDs from all chain IDs
      for
      (
        std::string::const_iterator itr( ligand_chain_ids.begin()), itr_end( ligand_chain_ids.end());
          itr != itr_end;
        ++itr
      )
      {
        all_chain_ids.erase
        (
          std::remove( all_chain_ids.begin(), all_chain_ids.end(), *itr),
          all_chain_ids.end()
        );
      }

      // set receptor chain IDs
      m_ReceptorChainIDs = all_chain_ids;
    }

  } // namespace fold
} // namespace bcl

