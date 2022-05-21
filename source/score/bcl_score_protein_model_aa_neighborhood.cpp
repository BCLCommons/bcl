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
#include "score/bcl_score_protein_model_aa_neighborhood.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "fold/bcl_fold_default_flags.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    //! @brief conversion to a string from a Type
    //! @param TYPE the type to get a string for
    //! @return a string representing that type
    const std::string &ProteinModelAANeighborhood::GetTypeName( const NormalizationType &TYPE)
    {
      static const std::string s_descriptors[ s_NumberTypes + 1] =
      {
        "none",
        "normalize",
        "rmsd",
        "explained",
        "correlation_and_relative_error",
        GetStaticClassName< ProteinModelAANeighborhood>()
      };
      return s_descriptors[ size_t( TYPE)];
    }

    //! @brief conversion to a string from a Type
    //! @param TYPE the type to get a string for
    //! @return a string representing that type
    const std::string &ProteinModelAANeighborhood::GetTypeSuffix( const NormalizationType &TYPE)
    {
      static const std::string s_descriptors[ s_NumberTypes + 1] =
      {
        "",
        "_norm",
        "_rmsd",
        "_explained",
        "_corr",
        GetStaticClassName< ProteinModelAANeighborhood>()
      };
      return s_descriptors[ size_t( TYPE)];
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinModelAANeighborhood::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelAANeighborhood())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelAANeighborhood::ProteinModelAANeighborhood() :
      m_ScoreAANeighborhood(),
      m_AANeighborListContainerGenerator(),
      m_Normalization( e_None),
      m_Scheme()
    {
    }

    //! @brief constructor from ShPtr to aa exposure function
    //! @param SP_SCORE_AA_NEIGHBORHOOD ShPtr to AA neighborhood scoring function
    //! @param TYPE normalization type to use
    //! @param CONSIDER_DIFFERENT_CHAIN whether neighbors from different chains should be considered
    //! @param SCORE_TYPE score type
    //! @param READABLE_SCHEME scheme that is more human readable
    ProteinModelAANeighborhood::ProteinModelAANeighborhood
    (
      const util::ShPtr< AANeighborhoodInterface> &SP_SCORE_AA_NEIGHBORHOOD,
      const TypeEnum &TYPE,
      const bool CONSIDER_DIFFERENT_CHAIN,
      const ProteinModel::Type &SCORE_TYPE,
      const std::string &READABLE_SCHEME
    ) :
      m_ScoreAANeighborhood( SP_SCORE_AA_NEIGHBORHOOD),
      m_AANeighborListContainerGenerator
      (
        SP_SCORE_AA_NEIGHBORHOOD.IsDefined() ?
          assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
          (
            SP_SCORE_AA_NEIGHBORHOOD->GetDistanceCutoff(),
            SP_SCORE_AA_NEIGHBORHOOD->GetMinimalSequenceSeparation(),
            CONSIDER_DIFFERENT_CHAIN,
            true
          )
          :
          util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> >()
      ),
      m_Normalization( TYPE),
      m_Scheme( SP_SCORE_AA_NEIGHBORHOOD->GetScheme()),
      m_ScoreType( SCORE_TYPE),
      m_ReadableScheme( READABLE_SCHEME)
    {
      m_Scheme += GetTypeSuffix( TYPE);
      if( m_ReadableScheme.empty())
      {
        m_ReadableScheme = GetScheme();
      }
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new ProteinModelAANeighborhood copied from this one
    ProteinModelAANeighborhood *ProteinModelAANeighborhood::Clone() const
    {
      return new ProteinModelAANeighborhood( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelAANeighborhood::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelAANeighborhood::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the sum of exposures of all amino acids for the given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return the sum of exposures of all amino acids for the given ProteinModel
    double ProteinModelAANeighborhood::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // check if protein model has a membrane
      const util::SiPtr< const biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // collect AANeighborLists for all amino acids in the given protein model
      const assemble::AANeighborListContainer all_aa_neighbor_list
      (
        m_AANeighborListContainerGenerator->operator ()( PROTEIN_MODEL)
      );

      // count scored amino acids
      size_t number_scored_aas( 0);
      double score( 0);

      // vectors for storing contact numbers provided or computed from the folded model
      storage::Vector< double> predicted_contact_numbers;
      storage::Vector< double> model_contact_numbers;

      // iterate over all lists in the all_aa_neighbor_list
      assemble::AANeighborCount neighcount;
      for
      (
        assemble::AANeighborListContainer::const_iterator
          aa_itr( all_aa_neighbor_list.Begin()), aa_itr_end( all_aa_neighbor_list.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        const biol::AABase &current_aa( *aa_itr->second.GetCenterAminoAcid());

        // skip undefined aa types
        if( !current_aa.GetType().IsDefined() || !current_aa.GetType()->IsNaturalAminoAcid())
        {
          continue;
        }

        // check that the current amino acid has a defined coordinate
        if( !current_aa.GetFirstSidechainAtom().GetCoordinates().IsDefined())
        {
          continue;
        }

        // calculate exposure for the current amino acid
        const double current_score( m_ScoreAANeighborhood->operator()( aa_itr->second, sp_membrane));

        // average
        if( util::IsDefined( current_score) && current_score != double( 0))
        {
          if( m_Normalization == e_RMSD)
          {
            score += math::Sqr( current_score);
          }
          else if( m_Normalization == e_CorrelationPlusRelativeError)
          {
            // insert predicted contact number for the current amino acid
            predicted_contact_numbers.PushBack( current_aa.GetExposurePrediction());

            // insert model contact number
            model_contact_numbers.PushBack( neighcount( aa_itr->second));
          }
          else if( m_Normalization == e_FractionExplained)
          {
            score += ( current_score - current_aa.GetExposurePrediction()) / current_aa.GetExposurePrediction();
          }
          else
          {
            score += current_score;
          }
          ++number_scored_aas;
        }
      }

      // normalized final score over entire protein
      if( m_Normalization == e_Normalize && number_scored_aas > 0)
      {
        score /= double( number_scored_aas);
      }
      else if( m_Normalization == e_RMSD && number_scored_aas > 0)
      {
        score = math::Sqrt( score / double( number_scored_aas));
      }
      else if( m_Normalization == e_FractionExplained && number_scored_aas > 0)
      {
        score /= double( number_scored_aas);
      }
      else if( m_Normalization == e_CorrelationPlusRelativeError && number_scored_aas > 0)
      {
        // compute the Spearman (ranked) correlation coefficient between predicted contact numbers
        // and contact numbers derived from folded protein model
        const double spearman_correlation_coefficient( math::Statistics::CorrelationSpearman
          (
            predicted_contact_numbers.Begin(),
            predicted_contact_numbers.End(),
            model_contact_numbers.Begin(),
            model_contact_numbers.End()
          ));

        // get the number of potentially structured amino acids in the sse pool
        const util::ShPtr< assemble::SSEPool> sse_pool( PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool));
        const size_t number_potential_structured_aas( sse_pool->GetNumberPotentiallyStructuredAAs());

        // compute the completeness of the folded protein model
        const double completeness( number_scored_aas / double( number_potential_structured_aas));

        // total relative error of contact number computed based on the folded protein model
        double total_contact_number_relative_error( 0);
        for
        (
          storage::Vector< double>::const_iterator
            predicted_cn_itr( predicted_contact_numbers.Begin()), model_cn_itr( model_contact_numbers.Begin()),
            predicted_cn_itr_end( predicted_contact_numbers.End());
          predicted_cn_itr != predicted_cn_itr_end;
          ++predicted_cn_itr, ++model_cn_itr
        )
        {
          total_contact_number_relative_error +=
            math::Absolute( ( *predicted_cn_itr) - ( *model_cn_itr) / completeness) / ( *predicted_cn_itr);
        }

        // positive correlation is rewarding whereas deviation in averages is penalizing
        score = -spearman_correlation_coefficient * number_scored_aas + total_contact_number_relative_error * math::Sqr( completeness);
      }

      // end
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &ProteinModelAANeighborhood::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ScoreAANeighborhood, ISTREAM);
      io::Serialize::Read( m_AANeighborListContainerGenerator, ISTREAM);
      io::Serialize::Read( m_Normalization, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &ProteinModelAANeighborhood::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ScoreAANeighborhood, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AANeighborListContainerGenerator, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Normalization, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    ProteinModelAANeighborhood::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM
    ) const
    {
      // check if protein model has a membrane
      const util::SiPtr< const biol::Membrane> sp_membrane( PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane));

      // collect AANeighborLists for all amino acids in the given protein model
      const assemble::AANeighborListContainer all_aa_neighbor_list
      (
        m_AANeighborListContainerGenerator->operator ()( PROTEIN_MODEL)
      );

      // iterate over all lists in the all_aa_neighbor_list
      for
      (
        assemble::AANeighborListContainer::const_iterator
          aa_itr( all_aa_neighbor_list.Begin()), aa_itr_end( all_aa_neighbor_list.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        m_ScoreAANeighborhood->WriteDetailedSchemeAndValues( aa_itr->second, sp_membrane, OSTREAM) << '\n';
      }

      // end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
