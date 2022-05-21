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
#include "fold/bcl_fold_mutate_membrane.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateMembrane::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateMembrane())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateMembrane::MutateMembrane() :
      m_Move(),
      m_Scheme( GetStaticClassName< MutateMembrane>())
    {
    }

    //! @brief constructor from a Move and a scheme
    //! @param MOVE function that performs the move on the membrane
    //! @param SCHEME scheme to be used
    MutateMembrane::MutateMembrane
    (
      const coord::MoveInterface &MOVE,
      const std::string &SCHEME
    ) :
      m_Move( MOVE.Clone()),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateMembrane
    MutateMembrane *MutateMembrane::Clone() const
    {
      return new MutateMembrane( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateMembrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief performs the mutate on the membrane associated with the passed protein model
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateMembrane::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // skip if move not defined
      if( !m_Move.IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // hard copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.HardCopy());

      // try to get the membrane
      util::ShPtr< assemble::ProteinModelData> sp_data( PROTEIN_MODEL.GetProteinModelData().HardCopy());
      util::ShPtr< biol::Membrane> sp_membrane( sp_data->GetData( assemble::ProteinModelData::e_Membrane).HardCopy());

      // skip if no membrane found
      if( !sp_membrane.IsDefined() || !sp_membrane->IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // apply the move
      m_Move->Move( *sp_membrane);

      // set the membrane data
      sp_data->Replace( assemble::ProteinModelData::e_Membrane, sp_membrane);
      new_model->SetProteinModelData( sp_data);

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateMembrane::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Move, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateMembrane::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Move, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
