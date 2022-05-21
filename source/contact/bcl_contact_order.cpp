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
#include "contact/bcl_contact_order.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_function_cached.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    math::FunctionCached< assemble::ProteinModel, double>
      Order::s_AbsoluteInstanceCached
      (
        util::ShPtr< Order>( new Order( e_Absolute, "raw")),
        &assemble::ProteinModel::GetDestructorSignal,
        &assemble::ProteinModel::GetChangeSignal
      );

    //! @brief NormalizationType as string
    //! @param NORMALIZATION_TYPE the NormalizationType
    //! @return the string for the NormalizationType
    const std::string &Order::GetNormalizationTypeDescriptor( const NormalizationType &NORMALIZATION_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "Absolute",
        "RelativeAAUsed",
        "RelativeSeqeunceLength",
        "RelativeSqrsAAUsed",
        "RelativeSqrSeqeunceLength",
        GetStaticClassName< NormalizationType>()
      };

      return s_descriptors[ NORMALIZATION_TYPE];
    }

  //////////
  // data //
  //////////

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &Order::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "co");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Order::Order() :
      m_Normalization( e_Absolute),
      m_Scheme( GetDefaultScheme()),
      m_CacheNeighborLists( false),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false, false
        )
      )
    {
    }

    //! @brief constructor from normalization type, scheme string and whether to use cache generating neighbor lists
    //! @brief NORMALIZE type of normalization
    //! @brief SCHEME scheme to report
    //! @brief CACHE whether to use cache generating neighbor lists
    Order::Order
    (
      const NormalizationType NORMALIZE,
      const std::string &SCHEME,
      const bool CACHE
    ) :
      m_Normalization( NORMALIZE),
      m_Scheme( SCHEME),
      m_CacheNeighborLists( CACHE),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false, CACHE
        )
      )
    {
    }

    //! @brief copy constructor
    Order *Order::Clone() const
    {
      return new Order( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Order::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief contact order for aa sequence
    //! @param AA_SEQUENCE sequence of interest
    //! @return contact order of sequence
    double Order::ContactOrder( const biol::AASequence &AA_SEQUENCE) const
    {
      // all amino acids
      const util::SiPtrVector< const biol::AABase> amino_acids( AA_SEQUENCE.GetMembers());

      // return normalized contact order
      return Normalize
      (
        ContactOrder
        (
          assemble::AANeighborListContainer
          (
            amino_acids, g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false
          )
        ),
        amino_acids.GetSize(),
        amino_acids.GetSize()
      );
    }

    //! @brief contact order for the given chain
    //! @param CHAIN chain of interest
    //! @return contact order for given chain
    double Order::ContactOrder( const assemble::Chain &CHAIN) const
    {
      // all amino acids
      const util::SiPtrVector< const biol::AABase> amino_acids( CHAIN.GetAminoAcids());

      // return normalized contact order
      return Normalize
      (
        ContactOrder
        (
          assemble::AANeighborListContainer
          (
            amino_acids, g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false
          )
        ),
        amino_acids.GetSize(),
        CHAIN.GetSequence()->GetSize()
      );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the contact order for a given protein model
    //! @param PROTEIN_MODEL protein model
    //! @return contact order for the given protein model
    double Order::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // currently, only one chain is defined
      if( PROTEIN_MODEL.GetChains().GetSize() != 1)
      {
        BCL_MessageVrb( "only single chain protein contact orders can be calculated!");
        return util::GetUndefined< double>();
      }
      // test whether this is the static instance. If so, it is necessary to compute the contact order
      if( this == &*s_AbsoluteInstanceCached.GetFunction())
      {
        return ContactOrder( m_NeighborGenerator->operator ()( PROTEIN_MODEL));
      }

      // return normalized contact order
      return Normalize
      (
        s_AbsoluteInstanceCached( PROTEIN_MODEL),
        PROTEIN_MODEL.GetChains().FirstElement()->GetNumberAAs(),
        PROTEIN_MODEL.GetChains().FirstElement()->GetSequence()->GetSize()
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Order::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Computes the contact order of a protein model allowing for various means of normalization"
      );
      serializer.AddInitializer
      (
        "normalization",
        "method of normalization",
        io::Serialization::GetAgent( &m_Normalization),
        GetNormalizationTypeDescriptor( e_Absolute)
      );
      serializer.AddInitializer
      (
        "cache neighbor lists",
        "whether to cache neighbor lists. This may speed up calculations if many methods make use of neighbor lists on "
        "the same protein model. If contact order is the only neighbor list being used, it may slow the computation "
        "down by ~10%",
        io::Serialization::GetAgent( &m_CacheNeighborLists),
        "False"
      );
      return serializer;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Order::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_NeighborGenerator =
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, s_ContactOrderSequencesDistanceCutoff, false, m_CacheNeighborLists
        );
      return true;
    }

  ////////////
  // helper //
  ////////////

    //! @brief normalize the contact order
    //! @param CONTACT_ORDER the unnormalized contact order
    //! @param NR_AAS_USED nr of amino acids used
    //! @param SEQUENCE_LENGTH the length of the amino acid sequence
    //! @return normalized contact order, depending of m_Normalization
    double Order::Normalize( const double CONTACT_ORDER, const size_t NR_AAS_USED, const size_t SEQUENCE_LENGTH) const
    {
      switch( m_Normalization)
      {
        case e_Absolute                 : return CONTACT_ORDER;
        case e_RelativeAAsUsed          : return NR_AAS_USED == 0 ? 0 : CONTACT_ORDER / double( NR_AAS_USED);
        case e_RelativeSequenceLength   : return SEQUENCE_LENGTH == 0 ? 0 : CONTACT_ORDER / double( SEQUENCE_LENGTH);
        case e_RelativeSqrAAsUsed       : return NR_AAS_USED == 0 ? 0 : math::Sqr( CONTACT_ORDER) / double( NR_AAS_USED);
        case e_RelativeSqrSequenceLength: return SEQUENCE_LENGTH == 0 ? 0 : math::Sqr( CONTACT_ORDER) / double( SEQUENCE_LENGTH);
        default                         : return util::GetUndefined< double>();
      }

      // should never reach this
      return util::GetUndefined< double>();
    }

    //! @brief calculate contact order form NeighborListContainer
    //! @param NEIGHBOR_LIST_CONTAINER container of neighbors
    //! @return the contact order
    double Order::ContactOrder( const assemble::AANeighborListContainer &NEIGHBOR_LIST_CONTAINER)
    {
      size_t seq_distance_sum( 0);
      size_t contacts( 0);

      // iterate over all aa neighbot list
      for
      (
        assemble::AANeighborListContainer::const_iterator
          itr( NEIGHBOR_LIST_CONTAINER.Begin()), itr_end( NEIGHBOR_LIST_CONTAINER.End());
        itr != itr_end;
        ++itr
      )
      {
        // sequence id of the center amino acid
        const int center_seq_id( itr->second.GetCenterAminoAcid()->GetSeqID());

        // iterate over neighbors
        for
        (
          assemble::AANeighborList::const_iterator neigh_itr( itr->second.Begin()), neigh_itr_end( itr->second.End());
          neigh_itr != neigh_itr_end;
          ++neigh_itr
        )
        {
          // calculate sequence distance
          const int seq_distance( neigh_itr->First()->GetSeqID() - center_seq_id);

          // assure that only one of the two pairs (stored redundant in AANeighborList) are used and test the distance
          if
          (
               seq_distance < int( s_ContactOrderSequencesDistanceCutoff)
            || neigh_itr->Second() >= g_ContactCbDistanceCutoff
          )
          {
            continue;
          }

          seq_distance_sum += size_t( seq_distance);
          contacts++;
        }
      }

      // if no contacts are found
      if( contacts == 0)
      {
        return double( 0.0);
      }

      // otherwise normalize
      return double( seq_distance_sum) / double( contacts);
    }

  } // namespace contact
} // namespace bcl
