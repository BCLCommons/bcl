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
#include "fold/bcl_fold_mutate_multimer.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "coord/bcl_coord_move_translate_random.h"
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

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateMultimer::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateMultimer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from members
    //! @param IS_DIHEDRAL is dihedral symmetry
    //! @param SCHEME scheme for the mutate
    MutateMultimer::MutateMultimer( const bool IS_DIHEDRAL, const std::string &SCHEME) :
      m_IsDihedral( IS_DIHEDRAL),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateMultimer
    MutateMultimer *MutateMultimer::Clone() const
    {
      return new MutateMultimer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateMultimer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateMultimer::GetAlias() const
    {
      static const std::string s_name( "MutateMultimer");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateMultimer::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Translates the protomer towards or away from the symmetry axis to change to multimer."
      );
      serializer.AddInitializer
      (
        "is dihedral",
        "whether symmetry is dihedral",
        io::Serialization::GetAgent( &m_IsDihedral)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief performs the mutate on the multiplier associated with the passed protein model
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateMultimer::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // get the multiplier
      util::ShPtr< assemble::ProteinModelData> sp_data( PROTEIN_MODEL.GetProteinModelData().HardCopy());
      util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        sp_data->GetData( assemble::ProteinModelData::e_Multiplier).HardCopy()
      );

      // return if no multiplier
      if( !sp_multiplier.IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // hard copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.HardCopy());

      // dihedral symmetry
      if( m_IsDihedral)
      {
        // translate in a random direction
        coord::MoveTranslateRandom random_translate( 0.0, 5.0, false);
        random_translate.Move( *new_model);
      }
      // cyclic symmetry
      else
      {
        // get the protein center
        const linal::Vector3D protein_center( PROTEIN_MODEL.GetCenter());

        // calculate footpoint on the axis
        const linal::Vector3D footpoint
        (
          linal::CalculateFootpoint
          (
            protein_center,
            sp_multiplier->GetCenter(),
            sp_multiplier->GetAxis( coord::GetAxes().e_Z)
          )
        );

        // get the vector from the footpoint to the protein
        const linal::Vector3D vector( footpoint - protein_center);

        // get a random number between -0.5 and 0.5
        const double distance_modifier( random::GetGlobalRandom().Double( math::Range< double>( -0.5, 0.5)));

        // apply a random translation toward or away from the axis
        new_model->Translate( distance_modifier * vector);
      }

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
