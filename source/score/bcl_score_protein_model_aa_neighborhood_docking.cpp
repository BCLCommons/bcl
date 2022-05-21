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
#include "score/bcl_score_protein_model_aa_neighborhood_docking.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "fold/bcl_fold_default_flags.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_statistics.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelAANeighborhoodDocking::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelAANeighborhoodDocking())
    );

    //!
    //!
    ProteinModelAANeighborhoodDocking::ProteinModelAANeighborhoodDocking() :
      m_ScoreAANeighborhood(),
      m_Scheme( "aa_neighborhood_restraint")
    {
      // nothing else to do
    }

    //! @param SCORE_AA_NEIGHBORHOOD
    //! @param SCHEME
    ProteinModelAANeighborhoodDocking::ProteinModelAANeighborhoodDocking
    (
      const AANeighborhoodInterface &SCORE_AA_NEIGHBORHOOD,
      const std::string &SCHEME
    ) :
      m_ScoreAANeighborhood( SCORE_AA_NEIGHBORHOOD),
      m_Scheme( SCHEME)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new ProteinModelAANeighborhood copied from this one
    ProteinModelAANeighborhoodDocking *ProteinModelAANeighborhoodDocking::Clone() const
    {
      return new ProteinModelAANeighborhoodDocking( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelAANeighborhoodDocking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelAANeighborhoodDocking::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProteinModelAANeighborhoodDocking::GetAlias() const
    {
      static const std::string s_alias( "ProteinModelAANeighborhoodDocking");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelAANeighborhoodDocking::GetSerializer() const
    {
      // create a Serializer object
      io::Serializer serializer;

      // set description
      serializer.SetClassDescription
      (
        "AANeighborhoodDocking is a scoring term that penalizes docked models "
        "in which the neighborhood environment of interface residues deviates "
        "from the given environment restraint"
      );
      serializer.AddInitializer
      (
        "environment type",
        "the neighborhood environment to be used",
        io::Serialization::GetAgent( &m_ScoreAANeighborhood)
      );

      // return the serializer
      return serializer;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool ProteinModelAANeighborhoodDocking::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM
    )
    {
      return true;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the sum of exposures of all amino acids for the given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return the sum of exposures of all amino acids for the given ProteinModel
    double ProteinModelAANeighborhoodDocking::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // get the membrane
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      //
      const util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> >
        aa_neighbor_list_container_generater
        (
          assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
          (
            m_ScoreAANeighborhood->GetDistanceCutoff(),
            m_ScoreAANeighborhood->GetMinimalSequenceSeparation(),
            true,
            true
          )
        );

      // all neighbor list
      const assemble::AANeighborListContainer all_neighbor_lists
      (
        aa_neighbor_list_container_generater->operator ()( PROTEIN_MODEL)
      );

      // initialize score
      double score( 0.0);

      // initialize the number of scored residues
      size_t number_of_scored_residues( 0);

      // iterate all neighbor lists
      for
      (
        assemble::AANeighborListContainer::const_iterator
          list_itr( all_neighbor_lists.Begin()), list_itr_end( all_neighbor_lists.End());
        list_itr != list_itr_end;
        ++list_itr
      )
      {
        // get central amino acid
        const util::SiPtr< const biol::AABase> sp_current_aa( list_itr->second.GetCenterAminoAcid());

        // skip undefined aa types
        if( !sp_current_aa->GetType().IsDefined() || !sp_current_aa->GetType()->IsNaturalAminoAcid())
        {
          continue;
        }

        // check that the current amino acid has a defined coordinate
        if( !sp_current_aa->GetFirstSidechainAtom().GetCoordinates().IsDefined())
        {
          continue;
        }

        // skip residues that have undefined exposure prediction
        if( !util::IsDefined( sp_current_aa->GetExposurePrediction()))
        {
          continue;
        }

        // skip non-membrane residues
        if
        (
          sp_membrane->DetermineEnvironmentType( sp_current_aa->GetFirstSidechainAtom().GetCoordinates())
            != biol::GetEnvironmentTypes().e_MembraneCore
        )
        {
          continue;
        }

        // compute the deviation
        double absolute_deviation( m_ScoreAANeighborhood->operator ()( list_itr->second, sp_membrane));

        // skip undefined deviations
        if( !util::IsDefined( absolute_deviation))
        {
          continue;
        }

        // update score
        score += math::Sqr( absolute_deviation);

        // update the number of scored residues
        number_of_scored_residues += 1;
      }

      // compute the RMSD-type score
      return math::Sqrt( score / number_of_scored_residues);
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &ProteinModelAANeighborhoodDocking::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ScoreAANeighborhood, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream object
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @return returns the output stream
    std::ostream &ProteinModelAANeighborhoodDocking::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ScoreAANeighborhood, OSTREAM);
      io::Serialize::Write( m_Scheme, OSTREAM);

      // return the stream object
      return OSTREAM;
    }

    //! @brief reads the exposure predictions from a file
    //! @param ISTREAM stream to read from
    //! @param PROTEIN_MODEL protein model to contain predictions
    void ProteinModelAANeighborhoodDocking::ReadPredictions( std::istream &ISTREAM, assemble::ProteinModel &PROTEIN_MODEL)
    {
      // vector that stores seq_id exposure pairs
      storage::Map< char, storage::Vector< storage::Pair< int, std::string> > > chain_exposure;

      // read residue seq_id and its corresponding exposure
      char chain_id;
      int seq_id;
      std::string exposure;
      while( ISTREAM >> chain_id >> seq_id >> exposure)
      {
        chain_exposure[ chain_id].PushBack( storage::Pair< int, std::string>( seq_id, exposure));
      }

      //iterate over all chains in the sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // get exposures for the current chain
        const storage::Vector< storage::Pair< int, std::string> >
          current_chain_exposure( chain_exposure[ ( *chain_itr)->GetChainID()]);

        // read predictions for the sequence
        ReadPredictions( current_chain_exposure, *( *chain_itr)->GetSequence());
      }
    }

    //! @brief reads the exposure predictions from a file
    //! @param ISTREAM stream to read from
    //! @param SEQUENCE sequence to contain predictions
    void ProteinModelAANeighborhoodDocking::ReadPredictions
    (
      const storage::Vector< storage::Pair< int, std::string> > &CHAIN_EXPOSURE,
      biol::AASequence &SEQUENCE
    )
    {
      // assign exposure value to the right residue
      for
      (
        storage::Vector< storage::Pair< int, std::string> >::const_iterator
          exposure_itr( CHAIN_EXPOSURE.Begin()), exposure_itr_end( CHAIN_EXPOSURE.End());
        exposure_itr != exposure_itr_end;
        ++exposure_itr
      )
      {
        const size_t current_aa_index( exposure_itr->First() - 1);
        util::ShPtr< biol::AABase> sp_current_aa( SEQUENCE.GetAA( current_aa_index));
        sp_current_aa->SetExposurePrediction
        (
          util::ConvertStringToNumericalValue< double>( exposure_itr->Second())
        );
//          // iterate over the sequence
//          for
//          (
//            biol::AASequence::iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
//            aa_itr != aa_itr_end; ++aa_itr
//          )
//          {
//            // if no exposure value for current residue, set it to undefined
//            if( ( *aa_itr)->GetSeqID() != exposure_itr->First())
//            {
//              ( *aa_itr)->SetExposurePrediction( util::GetUndefinedDouble());
//            }
//            else
//            {
//              ( *aa_itr)->SetExposurePrediction( util::ConvertStringToNumericalValue< double>( exposure_itr->Second()));
//            }
//          }
      }
    }

  } // namespace score
} // namespace bcl

