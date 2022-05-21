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
#include "fold/bcl_fold_mutate_protein_model.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "coord/bcl_coord_move_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModel::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModel())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModel::MutateProteinModel() :
      m_Move(),
      m_Scheme( GetStaticClassName< MutateProteinModel>())
    {
    }

    //! @brief constructor from a Move and a scheme
    //! @param MOVE function that performs the move on the protein
    //! @param SCHEME scheme to be used
    MutateProteinModel::MutateProteinModel
    (
      const coord::MoveInterface &MOVE,
      const std::string &SCHEME
    ) :
      m_Move( MOVE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief clone
    MutateProteinModel *MutateProteinModel::Clone() const
    {
      return new MutateProteinModel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateProteinModel::GetAlias() const
    {
      static const std::string s_alias( "MutateProteinModel");
      return s_alias;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateProteinModel::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Changes a protein model using the provided mover.");
      serializer.AddInitializer
      (
        "mover",
        "mutates the protein model",
        io::Serialization::GetAgent( &m_Move)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModel::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      BCL_MessageVrb( "using Move: " + util::Format()( m_Move));

      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // skip if move not defined
      if( !m_Move.IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // hard copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.HardCopy());

      // make a result to be returned that stores this ShPtr
      math::MutateResult< assemble::ProteinModel> result( new_model, *this);

      // move the model
      m_Move->Move( *new_model);

      // end
      return result;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
