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
#include "assemble/bcl_assemble_biomolecule.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "math/bcl_math_mutate_result.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Biomolecule::s_Instance
    (
      GetObjectInstances().AddInstance( new Biomolecule())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Biomolecule::Biomolecule() :
      m_Threshold( util::GetUndefined< double>())
    {
    }

    //! @brief constructor from required members
    //! @param ATOM_TYPES atoms types to use for superimposition
    //! @param SUPERIMPOSE the superimposition measure to use
    //! @param THRESHOLD the threshold for which two chains are considered the same
    Biomolecule::Biomolecule
    (
      const storage::Set< biol::AtomType> &ATOM_TYPES,
      const util::ShPtr< quality::SuperimposeInterface> &SUPERIMPOSE,
      const double THRESHOLD
    ) :
      m_AtomTypes( ATOM_TYPES),
      m_Superimpose( SUPERIMPOSE),
      m_Threshold( THRESHOLD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Biomolecule
    Biomolecule *Biomolecule::Clone() const
    {
      return new Biomolecule( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Biomolecule::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &Biomolecule::GetScheme() const
    {
      static const std::string s_scheme( "biomolecule");
      return s_scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a PROTEIN_MODEL and returning a protein model with a protein model multiplier in the
    //! data, that defines a symmetric subunits
    //! @param PROTEIN_MODEL model with multiple chains
    //! @return MutateResult that results from finding the bio molecules
    math::MutateResult< ProteinModel> Biomolecule::operator()( const ProteinModel &PROTEIN_MODEL) const
    {
      util::ShPtr< ProteinModel> sp_new_model( new ProteinModel());
      util::ShPtr< ProteinModelData> sp_new_model_data( PROTEIN_MODEL.GetProteinModelData().HardCopy());

      // check if there is a multiplier present in the argument
      util::ShPtr< ProteinModelMultiplier> sp_multiplier( sp_new_model_data->GetData( ProteinModelData::e_Multiplier));
      if( sp_multiplier.IsDefined())
      {
        BCL_MessageCrt( "creating biomolecules from protein models with multiplie not implemented yet");
        return math::MutateResult< ProteinModel>( util::ShPtr< ProteinModel>(), *this);
      }
      storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > chain_transformations;

      storage::Set< char> inserted_chains;
      // iterate through chains of argument
      for
      (
        ProteinModel::const_iterator
          chain_itr_template( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr_template != chain_itr_end;
        ++chain_itr_template
      )
      {
        const char template_chain_id( ( *chain_itr_template)->GetChainID());
        if( inserted_chains.Contains( template_chain_id))
        {
          continue;
        }
        sp_new_model->Insert( *chain_itr_template);
        inserted_chains.Insert( template_chain_id);

        // insert identity transformation for template chain
        chain_transformations.PushBack
        (
          storage::Triplet< char, char, math::TransformationMatrix3D>
          (
            template_chain_id, template_chain_id, math::TransformationMatrix3D()
          )
        );

        const std::string template_sequence_string( ( *chain_itr_template)->GetSequence()->Sequence());
        const util::SiPtrVector< const linal::Vector3D> template_coords( ( *chain_itr_template)->GetAtomCoordinates( m_AtomTypes));

        // iterate through other chains to search for identical sequences that can be superimposed
        for
        (
          ProteinModel::const_iterator chain_itr( chain_itr_template + 1);
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          const std::string current_sequence_string( ( *chain_itr)->GetSequence()->Sequence());

          // skip chains of different sequence
          if( template_sequence_string != current_sequence_string)
          {
            BCL_MessageDbg( "found non-identical sequences:\n" + template_sequence_string + '\n' + current_sequence_string);
            continue;
          }

          const util::SiPtrVector< const linal::Vector3D> current_coords( ( *chain_itr)->GetAtomCoordinates( m_AtomTypes));

          // superimpose coordinates
          const double measure( m_Superimpose->CalculateMeasure( current_coords, template_coords));

          // is this superimposition within the threshold
          if( !m_Superimpose->GetComparisonFunction()( measure, m_Threshold))
          {
            continue;
          }
          const char current_chain_id( ( *chain_itr)->GetChainID());
          BCL_MessageDbg
          (
            "found superimposable chain pair: " + std::string( 1, template_chain_id)
            + '\n' + std::string( 1, current_chain_id)
          );

          // calculate the matrix
          const math::TransformationMatrix3D transformation( m_Superimpose->CalculateSuperimposition( current_coords, template_coords));
          chain_transformations.PushBack
          (
            storage::Triplet< char, char, math::TransformationMatrix3D>
            (
              template_chain_id, current_chain_id, transformation
            )
          );
          inserted_chains.Insert( current_chain_id);
        }
      }

      sp_multiplier = util::ShPtr< ProteinModelMultiplier>( new ProteinModelMultiplier( chain_transformations, *sp_new_model));
      sp_new_model_data->Insert( ProteinModelData::e_Multiplier, sp_multiplier);
      sp_new_model->SetProteinModelData( sp_new_model_data);

      return math::MutateResult< ProteinModel>( sp_new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Biomolecule::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AtomTypes  , ISTREAM);
      io::Serialize::Read( m_Superimpose, ISTREAM);
      io::Serialize::Read( m_Threshold  , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Biomolecule::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AtomTypes  , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Superimpose, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Threshold  , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
