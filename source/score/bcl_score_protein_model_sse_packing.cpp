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
#include "score/bcl_score_protein_model_sse_packing.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "assemble/bcl_assemble_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "linal/bcl_linal_vector_reference.h"
#include "signal/bcl_signal_signal.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_string_functions.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSEPacking::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelSSEPacking())
    );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CacheHelper
    //! @brief manages the cache to link protein models to their scores
    //!
    //! @remark example unnecessary
    //! @author mendenjl
    //! @date Mar 29, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class CacheHelper :
      public signal::Slots
    {
    private:
      //! @brief the cache that stores for each Argument the result
      //! the key is the pointer (address) of the argument, the data is a pointer to a result, NULL if the result does
      //! not exist (because the argument has changed)
      typedef std::map< assemble::ProteinModel *, const linal::Vector< double> *> CacheType;
      mutable CacheType m_Cache;

      //! Instance of the main class
      ProteinModelSSEPacking m_Packer;

    public:

      //! Default constructor
      CacheHelper() :
        m_Packer( ProteinModelSSEPacking::s_NumberTypes)
      {
      }

      //! Destructor
      ~CacheHelper()
      {
        // iterate over cache and delete all result pointer
        for( typename CacheType::iterator itr( m_Cache.begin()), itr_end( m_Cache.end()); itr != itr_end; ++itr)
        {
          // delete only if there is an actual object
          if( itr->second != NULL)
          {
            delete itr->second;
          }
        }
      }

      //! @brief return the only instance of this class
      static CacheHelper &GetInstance()
      {
        static CacheHelper s_instance;
        return s_instance;
      }

      //! @brief operator taking an ARGUMENT and returning a t_ResultType object
      //! If there is a result stored in the cache, return that
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @return function value of the given argument
      linal::Vector< double> operator()
      (
        const assemble::ProteinModel &MODEL,
        const assemble::SSEPool &MISSING_SSES,
        std::ostream &OSTREAM,
        const bool &DO_WRITE
      ) const
      {
        // search in cache
        typename CacheType::iterator itr( m_Cache.find( typename CacheType::key_type( &MODEL)));

        // no result found
        if( itr == m_Cache.end())
        {
          // try to insert
          const std::pair< typename CacheType::iterator, bool> insert_pair
          (
            m_Cache.insert
            (
              typename CacheType::value_type
              (
                typename CacheType::key_type( &MODEL),
                typename CacheType::mapped_type( m_Packer.Score( MODEL, MISSING_SSES, OSTREAM, DO_WRITE).Clone())
              )
            )
          );

          // check that insert was successful
          BCL_Assert( insert_pair.second, "could not insert the result!");

          // connect this to all signal handlers on ARGUMENT
          CacheHelper *non_const_this_ptr( const_cast< CacheHelper *>( this));

          // register this to the destructor signal with RemoveResultFromCache function
          MODEL.GetDestructorSignal().Connect( non_const_this_ptr, &CacheHelper::RemoveResultFromCache);
          // register this to the update signal with InvalidateResult function
          MODEL.GetChangeSignal().Connect( non_const_this_ptr, &CacheHelper::InvalidateResult);

          // return the result
          return *insert_pair.first->second;
        }
        else if( itr->second == NULL)
        {
          // result invalidated. Update it
          itr->second = m_Packer.Score( MODEL, MISSING_SSES, OSTREAM, DO_WRITE).Clone();
        }
        else if( DO_WRITE)
        {
          m_Packer.Score( MODEL, MISSING_SSES, OSTREAM, DO_WRITE);
        }
        // return the result
        return *itr->second;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief remove results for this object from the cache
      //! @param ARGUMENT argument that should be removed
      void RemoveResultFromCache( const assemble::ProteinModel &MODEL)
      {
        // find result for that Argument
        typename CacheType::iterator itr( m_Cache.find( typename CacheType::key_type( &MODEL)));

        // skip if there is no result cached
        if( itr == m_Cache.end())
        {
          return;
        }

        // copy of result - if result is deleted before the itr is removed from cache, if could be that the deletion of
        // the result triggers another call to this function
        const linal::Vector< double> *result( itr->second);

        // erase from cache
        m_Cache.erase( itr);

        // delete result
        delete result;
      }

      //! @brief delete a result if argument changes and function has to be reevaluated
      //! @param ARGUMENT the argument for which the result has to be forgotten
      void InvalidateResult( const assemble::ProteinModel &MODEL)
      {
        // find result for that Argument
        typename CacheType::iterator itr( m_Cache.find( typename CacheType::key_type( &MODEL)));

        // delete if there is a result cached
        if( itr != m_Cache.end() && itr->second != NULL)
        {
          // delete result
          const linal::Vector< double> *result( itr->second);
          itr->second = NULL;
          delete result;
        }
      }

    };

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    //! @param PSEUDOCOUNT pseudocount to use
    ProteinModelSSEPacking::ProteinModelSSEPacking( const Type &TYPE) :
      m_Type( TYPE)
    {
      if( TYPE == s_NumberTypes)
      {
        std::stringstream err_str;
        BCL_Assert
        (
          this->ReadInitializerSuccessHook( util::ObjectDataLabel(), err_str),
          "Couldn't read histogram file: " + err_str.str()
        );
      }
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelSSEPacking
    ProteinModelSSEPacking *ProteinModelSSEPacking::Clone() const
    {
      return new ProteinModelSSEPacking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelSSEPacking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelSSEPacking::GetScheme() const
    {
      static const std::string s_contact_names[ s_NumberTypes + 1] =
      {
         "sse_contact_type",
         "sse_contact_adjacency",
         "sse_orientation",
         "sse_interaction_angle_overlap",
         GetStaticClassName< Type>()
      };
      return s_contact_names[ size_t( m_Type)];
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &ProteinModelSSEPacking::GetAlias() const
    {
      return GetScheme();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelSSEPacking::GetSerializer() const
    {
      static const std::string s_desc[ s_NumberTypes + 1] =
      {
         "Propensity of each contact type, assuming that SSEs of the generic type are involved",
         "Propensity of the given SSE types to be in contact given that they are adjacent in sequence",
         "Orientational propensities given distance in # of SSEs apart the SSEs actually are",
         "Propensity of the SSE interaction to have the given angle and overlap with one another",
         GetStaticClassName< Type>()
      };
      io::Serializer serializer;
      serializer.SetClassDescription( std::string ( s_desc[ size_t( m_Type)]));
      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that scores the Protein model
    //! @param PROTEIN_MODEL the protein model for which all neighbor scores are calculated
    //! @return score
    double ProteinModelSSEPacking::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );

      // find all sse's not yet in the model. do not assume a linear loop between sses if another SSE could exist
      // between them, because that would invalidate the linear loop assumption
      util::SiPtrList< const assemble::SSE> non_overlapping
      (
        sp_pool.IsDefined()
        ? sp_pool->GetNonOverlappingSSEs( PROTEIN_MODEL)
        : util::SiPtrList< const assemble::SSE>()
      );

      // create an sse-pool with the missing SSEs
      assemble::SSEPool sse_pool( non_overlapping);

      // iterate through the chains of PROTEIN_MODEL
      std::stringstream oss;

      return CacheHelper::GetInstance().operator ()( PROTEIN_MODEL, sse_pool, oss, false)( size_t( m_Type));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ProteinModelSSEPacking::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM
    ) const
    {
      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );

      // find all sse's not yet in the model. do not assume a linear loop between sses if another SSE could exist
      // between them, because that would invalidate the linear loop assumption
      assemble::SSEPool non_overlapping
      (
        sp_pool.IsDefined()
        ? sp_pool->GetNonOverlappingSSEs( PROTEIN_MODEL)
        : util::SiPtrList< const assemble::SSE>()
      );

      CacheHelper::GetInstance().operator ()( PROTEIN_MODEL, non_overlapping, OSTREAM, true);
      // end
      return OSTREAM;
    }

    //! @brief helper function called by WriteDetailedSchemeAndValues and operator() so that the code remains in-sync
    //! @param CHAIN the chain of interest
    //! @param MISSING_SSES missing sses from the pool
    //! @param OSTREAM the output stream to write the detailed scheme to for this chain
    //! @param DO_WRITE set to true to actually write to the output stream; otherwise, nothing will be written
    //! @return the final score
    linal::Vector< double> ProteinModelSSEPacking::Score
    (
      const assemble::ProteinModel &MODEL,
      const assemble::SSEPool &MISSING_SSES,
      std::ostream &OSTREAM,
      const bool &DO_WRITE
    ) const
    {
      static const size_t n_categories( scorestat::ProteinModelPacking::s_NumberNaturalCategories);
      static const size_t n_non_weight_categories( scorestat::ProteinModelPacking::e_WeakInteraction);
      static const Type type_mapping[ n_categories] =
      {
        e_AdjacentSSEContactPropensity, // e_Adjacent
        e_AdjacentSSEContactPropensity, // e_NonAdjacent
        e_Orientation, // AdjP
        e_Orientation, // AdjAP
        e_Orientation, // SSE+1P
        e_Orientation, // SSE+1AP
        e_Orientation, // SSE+2P
        e_Orientation, // SSE+2AP
        e_Orientation, // SSE+3P
        e_Orientation, // SSE+3AP
        e_InteractionWeight,
        e_InteractionWeight,
        e_InteractionWeight
      };
      // sum of all scores
      linal::Vector< double> score( s_NumberTypes, 0.0);

      // need at least two sses
      if( MODEL.GetNumberSSEs() < 2)
      {
        return score;
      }
      BCL_Assert
      (
        m_ContactTypeEntropies.GetSize(),
        "No contact propensities loaded! The static instance of " + this->GetClassIdentifier() + " must be used. Or the "
        "histogram may not have been found"
      );

      // collect all non-coil sse's
      const util::SiPtrVector< const assemble::SSE> structured_sses
      (
        MODEL.GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
      );

      if( structured_sses.GetSize() < 2)
      {
        return score;
      }

      // initializes maps for storing histograms of distances between all types of amino acid pairs
      storage::Vector< storage::Vector< linal::Vector< size_t> > > packing_info_by_sse_contact_type
      (
        structured_sses.GetSize(),
        storage::Vector< linal::Vector< size_t> >
        (
          contact::GetTypes().GetEnumCount(),
          linal::Vector< size_t>( size_t( scorestat::ProteinModelPacking::s_NumberNaturalCategories), size_t( 0))
        )
      );
      storage::Vector< linal::Vector< size_t> > contact_type_counts
      (
        structured_sses.GetSize(),
        linal::Vector< size_t>( contact::GetTypes().GetEnumCount(), size_t( 0))
      );
      storage::Vector< linal::Vector< float> > sum_packing_info_by_sse_contact_type
      (
        contact::GetTypes().GetEnumCount(),
        linal::Vector< float>( size_t( scorestat::ProteinModelPacking::s_NumberNaturalCategories), float( 0))
      );
      linal::Vector< float> sum_ave_contact_type_counts( contact::GetTypes().GetEnumCount(), float( 0));
      storage::Vector< size_t> adj_noncontact_count( structured_sses.GetSize(), size_t( 0));
      storage::Vector< size_t> adj_contact_count( structured_sses.GetSize(), size_t( 0));
      storage::Vector< size_t> dist_contact_count( structured_sses.GetSize(), size_t( 0));

      // min probability to score two residues as being an almost-certain contact
      const double min_p( 0.04);

      const size_t min_contacting_atoms( m_PackingDefinition.GetMinAtomsInContact());

      assemble::VoxelGridAA interactions_detector( m_PackingDefinition.GetMinInteractionDistance());
      linal::Matrix< float> interactions_matrix
      (
        interactions_detector.GetSSEInteractionMatrix
        (
          structured_sses,
          MODEL.GetAminoAcids(),
          2,
          m_PackingDefinition.GetMinInteractionDistance(),
          false,
          min_p,
          false
        )
      );
      // handle SSE-pair interaction-weight and parallel/anti-parallel bias
      size_t a( 0);
      size_t aas_in_sses( 0);
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          itr_a( structured_sses.Begin()),
          itr_end( structured_sses.End());
        itr_a != itr_end;
        ++itr_a, ++a
      )
      {
        size_t b( a + 1);
        aas_in_sses += ( *itr_a)->GetSize();
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator itr_b( itr_a + 1);
          itr_b != itr_end;
          ++itr_b, ++b
        )
        {
          bool is_adjacent( itr_b == itr_a + 1 && ( *itr_a)->GetChainID() == ( *itr_b)->GetChainID());
          if( interactions_matrix( a, b) < min_contacting_atoms && !is_adjacent)
          {
            continue;
          }
          size_t effective_separation( b - a - 1);
          const int end_last_sse_seq_id( ( *itr_a)->GetLastAA()->GetPdbID());
          const int start_next_sse_seq_id( ( *itr_b)->GetFirstAA()->GetPdbID());

          // check for any overlapping sses in the pool
          for
          (
            auto itr_pool( MISSING_SSES.Begin()), itr_pool_end( MISSING_SSES.End());
            itr_pool != itr_pool_end;
            ++itr_pool
          )
          {
            const assemble::SSE &potential_sse( **itr_pool);
            if
            (
              potential_sse.GetFirstAA()->GetSeqID() > end_last_sse_seq_id
              && potential_sse.GetLastAA()->GetSeqID() < start_next_sse_seq_id
            )
            {
              // overlapping sse. inflate the number of residues between these sses to ensure that we don't check
              // the distance between the pair
              effective_separation += 1;
            }
          }
          if( effective_separation && interactions_matrix( a, b) < min_contacting_atoms)
          {
            continue;
          }
          const double max_loop_length
          (
            ( ( *itr_b)->GetFirstAA()->GetPdbID() - ( *itr_a)->GetLastAA()->GetPdbID()) * 2.56 + 2.11
          );
          const bool could_be_opposite_ori
          (
            ( *itr_b)->GetChainID() != ( *itr_a)->GetChainID()
            ||
            std::min
            (
              linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetLastAA()->GetCA().GetCoordinates()),
              std::min
              (
                linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_b)->GetFirstAA()->GetCA().GetCoordinates()),
                linal::Distance( ( *itr_a)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetFirstAA()->GetCA().GetCoordinates())
              ) + 2.34
            ) < max_loop_length
          );

          // for every pair of SSEs
          assemble::SSEGeometryPacking packing_ab( **itr_a, **itr_b, 0.0);
          if( !packing_ab.GetContactType().IsDefined())
          {
            continue;
          }
          const bool in_contact( interactions_matrix( a, b) >= min_contacting_atoms);
          if( DO_WRITE)
          {
            OSTREAM
              << ( *itr_a)->GetIdentification()
              << " " << ( *itr_b)->GetIdentification()
              << " " << packing_ab.GetContactType().GetName()
              << in_contact
              << ( is_adjacent ? "Adjacent" : "Non-Adjacent")
              << packing_ab.GetOrientationName( packing_ab.GetOrientation())
              << '\n';
          }
          const contact::Type rev_type( contact::GetTypes().Reverse( packing_ab.GetContactType()));
          if( in_contact)
          {
            ++contact_type_counts( a)( packing_ab.GetContactType());
            ++contact_type_counts( b)( rev_type);
          }
          m_PackingDefinition.AddPackingType
          (
            packing_ab,
            in_contact,
            effective_separation,
            false,
            could_be_opposite_ori,
            packing_info_by_sse_contact_type( a)( packing_ab.GetContactType())
          );
          m_PackingDefinition.AddPackingType
          (
            packing_ab,
            in_contact,
            effective_separation,
            false,
            could_be_opposite_ori,
            packing_info_by_sse_contact_type( b)( rev_type)
          );
        }
        const linal::Vector< size_t> &contact_counts( contact_type_counts( a));
        const float total_contact_counts( contact_counts.Sum());
        float norm[ n_categories] =
        {
          0.0, // e_Adjacent
          0.0, // e_NonAdjacent
          0.0, // AdjP
          0.0, // AdjAP
          0.0, // SSE+1P
          0.0, // SSE+1AP
          0.0, // SSE+2P
          0.0, // SSE+2AP
          0.0, // SSE+3P
          0.0, // SSE+3AP
          total_contact_counts,
          total_contact_counts,
          total_contact_counts
        };
        const storage::Vector< linal::Vector< size_t> > &p_info_this_sse( packing_info_by_sse_contact_type( a));
        for( size_t i( 0), n( contact_counts.GetSize()); i < n; ++i)
        {
          const linal::Vector< size_t> &contacts_counts_local( p_info_this_sse( i));
          for( size_t j( 0); j < n_non_weight_categories; ++j)
          {
            norm[ j] += contacts_counts_local( j);
          }
        }
        norm[ scorestat::ProteinModelPacking::e_AdjacentInContact]
          = norm[ scorestat::ProteinModelPacking::e_AdjacentNotInContact]
          += norm[ scorestat::ProteinModelPacking::e_AdjacentInContact];
        norm[ scorestat::ProteinModelPacking::e_AdjacentParallel]
          = norm[ scorestat::ProteinModelPacking::e_AdjacentAntiParallel]
          += norm[ scorestat::ProteinModelPacking::e_AdjacentParallel];
        norm[ scorestat::ProteinModelPacking::e_OneSSEApartParallel]
          = norm[ scorestat::ProteinModelPacking::e_OneSSEApartAntiParallel]
          += norm[ scorestat::ProteinModelPacking::e_OneSSEApartParallel];
        norm[ scorestat::ProteinModelPacking::e_TwoSSEApartParallel]
          = norm[ scorestat::ProteinModelPacking::e_TwoSSEApartAntiParallel]
          += norm[ scorestat::ProteinModelPacking::e_TwoSSEApartParallel];
        norm[ scorestat::ProteinModelPacking::e_ThreeOrMoreSSEApartParallel]
          = norm[ scorestat::ProteinModelPacking::e_ThreeOrMoreSSEApartAntiParallel]
          += norm[ scorestat::ProteinModelPacking::e_ThreeOrMoreSSEApartParallel];

        for( size_t i( 0), n( contact_counts.GetSize()); i < n; ++i)
        {
          if( total_contact_counts)
          {
            sum_ave_contact_type_counts( i) += float( contact_counts( i)) / total_contact_counts;
          }
          const linal::Vector< size_t> &contacts_counts_local( p_info_this_sse( i));
          for( size_t j( 0); j < n_categories; ++j)
          {
            if( contacts_counts_local( j))
            {
              sum_packing_info_by_sse_contact_type( i)( j) += float( contacts_counts_local( j)) / norm[ j];
            }
          }
        }
      }
      for( size_t i( 0), n_contact_types( contact::GetTypes().GetEnumCount()); i < n_contact_types; ++i)
      {
        contact::Type type( i);
        const linal::Vector< float> &spec_counts( sum_packing_info_by_sse_contact_type( i));
        if( sum_ave_contact_type_counts( i))
        {
          score( e_ContactType) += m_SSToContactTypeEntropy( i) * sum_ave_contact_type_counts( i);
          if( DO_WRITE)
          {
            OSTREAM << type.GetName() << " Count: " << sum_ave_contact_type_counts( i) << " sum-entropy: "
                    << sum_ave_contact_type_counts( i) * m_SSToContactTypeEntropy( i) << '\n';
          }
        }
        for( size_t j( 0); j < scorestat::ProteinModelPacking::s_NumberNaturalCategories; ++j)
        {
          if( spec_counts( j))
          {
            score( size_t( type_mapping[ j])) += m_ContactTypeEntropies( i)( j) * spec_counts( j);
            if( DO_WRITE)
            {
              OSTREAM << type.GetName() << " "
                      << scorestat::ProteinModelPacking::GetCategoryName( scorestat::ProteinModelPacking::Category( j))
                      << " Count: " << spec_counts( j) << " sum-entropy: "
                      << spec_counts( j) * m_ContactTypeEntropies( i)( j) << '\n';
            }
          }
        }
      }
      score *= double( aas_in_sses) / double( structured_sses.GetSize());

      // return the score sum
      return score;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
    bool ProteinModelSSEPacking::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( !command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, Score::AddHistogramPath( "ssepack_entropies.table"));
        if( !m_PackingDefinition.TryRead( util::ObjectDataLabel( input), ERROR_STREAM))
        {
          return false;
        }
        m_SSToContactTypeEntropy = m_PackingDefinition.ReadSSPairToContactEntropies( input);
        m_ContactTypeEntropies = m_PackingDefinition.ReadContactTypeEntropies( input);
        io::File::CloseClearFStream( input);
      }
      return true;
    }

  } // namespace score
} // namespace bcl
