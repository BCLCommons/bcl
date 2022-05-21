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
// include header for this class
#include "scorestat/bcl_scorestat_aa_distance.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AADistance::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AADistance())
    );

    //! @brief CategoryOptionName as string
    //! @param CATEGORY_OPTION the AADistanceCategoryOption
    //! @return the string for the CategoryOptionName
    const std::string &AADistance::GetCategoryOptionName( const AADistanceCategoryOption &CATEGORY_OPTION)
    {
      static const std::string s_names[] =
      {
          "OneChainNew",
          "OneChainOld",
          "AllChainNew",
          "AllChainOld",
          GetStaticClassName< AADistanceCategoryOption>()
      };
      return s_names[ CATEGORY_OPTION];
    }

    //! @brief CategoryFileName as string
    //! @param CATEGORY_OPTION the AADistanceCategoryOption
    //! @return the string for the CategoryFileName
    const std::string &AADistance::GetCategoryFileName( const AADistanceCategoryOption &CATEGORY_OPTION)
    {
      static const std::string s_names[] =
      {
          "aa_distances_one_chain_new.histograms",
          "aa_distances_one_chain_old.histograms",
          "aa_distances_all_chain_new.histograms",
          "aa_distances_all_chain_old.histograms",
          GetStaticClassName< AADistanceCategoryOption>()
      };
      return s_names[ CATEGORY_OPTION];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AADistance::AADistance() :
      m_CategoryOption( e_OneChainNew),
      m_BinSize( 1.0),
      m_ChainIds( ""),
      m_AADistSeqExcl( 6)
    {
    }

    //! @brief virtual copy constructor
    AADistance *AADistance::Clone() const
    {
      return new AADistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AADistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AADistance::GetOutFilePostfix() const
    {
      return GetCategoryFileName( m_CategoryOption);
    }

    //! @brief returns the binsize for the histogram
    //! @return the binsize for the histogram
    const double &AADistance::GetBinSize() const
    {
      return m_BinSize;
    }

    //! @brief returns chain id
    //! @return chain id
    const std::string &AADistance::GetChainId() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AADistance::GetAlias() const
    {
      static const std::string s_name( "AADistance");
      return s_name;
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string AADistance::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // initializes maps for storing histograms of distances between all types of amino acid pairs
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> > distance_map_one_old;
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> > distance_map_one_new;
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> > distance_map_all_old;
      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> > distance_map_all_new;

      // initialize maps for storing aa distances
      // iterate over all aa types to set the first aa
      for
      (
        biol::AATypes::const_iterator aa_1_itr( biol::GetAATypes().Begin()), aa_itr_end( biol::GetAATypes().End());
          aa_1_itr != aa_itr_end;
        ++aa_1_itr
      )
      {
        // iterate over all aa types to set the second aa
        for( biol::AATypes::const_iterator aa_2_itr( aa_1_itr); aa_2_itr != aa_itr_end; ++aa_2_itr)
        {
          // skip non-natural amino acids
          if( !( *aa_1_itr)->IsNaturalAminoAcid() || !( *aa_2_itr)->IsNaturalAminoAcid())
          {
            continue;
          }

          // initialize a histogram
          util::ShPtr< math::Histogram> histogram(
            new math::Histogram( double( 0), m_BinSize, size_t( 100 / m_BinSize)));

          // add aa pair-histogram to the map as a side effect of subscripting the map
          distance_map_one_old[ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] = histogram.HardCopy();
          distance_map_one_new[ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] = histogram.HardCopy();
          distance_map_all_old[ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] = histogram.HardCopy();
          distance_map_all_new[ std::pair< biol::AAType, biol::AAType>( *aa_1_itr, *aa_2_itr)] = histogram.HardCopy();
        }
      }

      // initialize histogram for storing disulfide bond statistics
      storage::VectorND< 2, math::Histogram> disulfide_bonds
      ( math::Histogram( double( 0), m_BinSize, size_t( 100 / m_BinSize)));

      // cutoff of disulfide bond length
      static const double disulfide_bond_cutoff( 2.5);

      // iterate through all models in the ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
          protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current protein model
        util::ShPtr< assemble::ProteinModel> sp_current_model( *protein_model_itr);
        const assemble::ProteinModel &protein_model( *sp_current_model);

        // get current pdb name before all the loops start
        const util::ShPtr< util::Wrapper< std::string> > &model_name_ptr(
          protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // instantiate sequences CaCb from pdb for each chain one
        BCL_MessageStd( "building sequences from pdb chains");
        util::SiPtrVector< const biol::AASequence> aa_sequences( protein_model.GetSequences());

        // number of chains
        BCL_MessageStd( "pdb has " + util::Format()( aa_sequences.GetSize()) + " chains.");

        // skip undefined pdbs
        if( aa_sequences.IsEmpty())
        {
          BCL_MessageCrt( "pdb does not contain chains.");
          continue;
        }

        // iterate through all chains
        for
        (
          util::SiPtrVector< const biol::AASequence>::const_iterator
            seq_itr( aa_sequences.Begin()), seq_itr_end( aa_sequences.End());
          seq_itr != seq_itr_end; ++seq_itr
        )
        {
          // skip undesired chains, if m_ChainIds is empty, then analyze all chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *seq_itr)->GetChainID()) == std::string::npos)
          {
            BCL_MessageDbg( std::string( "Skip chain: ") + ( *seq_itr)->GetChainID() + std::string( ", chain ids to use: ") + m_ChainIds);

            // continue to the next chain
            continue;
          }

          // iterate over the current chain to get the first amino acid residue
          for
          (
            biol::AASequence::const_iterator
              aa_1_itr( ( *seq_itr)->GetData().Begin()), aa_itr_end( ( *seq_itr)->GetData().End());
            aa_1_itr != aa_itr_end;
            ++aa_1_itr
          )
          {
            // proceed only if coordinates are given and it is a natural amino acid
            if
            (
              !( *aa_1_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined() ||
              !( *aa_1_itr)->GetType().IsDefined() || !( *aa_1_itr)->GetType()->IsNaturalAminoAcid()
            )
            {
              continue;
            }

            // iterate over the current chain and all other chains
            for( util::SiPtrVector< const biol::AASequence>::const_iterator seq_2_itr( seq_itr); seq_2_itr != seq_itr_end; ++seq_2_itr)
            {
              // iterate over the current chain to get the second amino acid residue to which the distance of the first amino acid residue is calculated
              for
              (
                biol::AASequence::const_iterator aa_2_itr( ( *seq_itr)->GetData().Begin() + 1);
                  aa_2_itr != ( *seq_itr)->GetData().End();
                ++aa_2_itr
              )
              {
                // proceed only if coordinates are given and it is a natural amino acid
                if( !( *aa_2_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined()
                    || !( *aa_2_itr)->GetType().IsDefined() || !( *aa_2_itr)->GetType()->IsNaturalAminoAcid())
                {
                  continue;
                }

                // check for disulfide bond
                if( ( *aa_1_itr)->GetType() == biol::GetAATypes().CYS && ( *aa_2_itr)->GetType() == biol::GetAATypes().CYS)
                {
                  // calculate distance between two cysteine residues
                  const double cys_cys_distance
                  (
                    biol::Distance( ( *aa_1_itr)->GetAtom( biol::GetAtomTypes().SG), ( *aa_2_itr)->GetAtom( biol::GetAtomTypes().SG))
                  );

                  // store the cys_cys_distance if it is defined and is within the range of disulfide bond length
                  if( util::IsDefined( cys_cys_distance) && cys_cys_distance < disulfide_bond_cutoff)
                  {
                    // aa pair distance
                    const double aa_aa_distance( biol::FirstSidechainAtomDistance( ( **aa_1_itr), ( **aa_2_itr)));

                    // insert into histogram
                    disulfide_bonds( 0).PushBack( aa_aa_distance);
                    disulfide_bonds( 1).PushBack( cys_cys_distance);

                    // skip the current amino acid residue since its distance has been calculated
                    continue;
                  }
                }

                // get sequence separation
                const size_t aa_aa_seq_separation( biol::SequenceSeparation( ( **aa_1_itr), ( **aa_2_itr)));

                // skip amino acid residues that are not from the current chain or are too close
                if( !util::IsDefined( aa_aa_seq_separation) || aa_aa_seq_separation < m_AADistSeqExcl)
                {
                  continue;
                }

                // distance between the first side chain atom of the first amino acid and that of the second amino acid
                const double aa_aa_distance( biol::FirstSidechainAtomDistance( ( **aa_1_itr), ( **aa_2_itr)));

                // distance weight
                const double aa_aa_distance_weight
                (
                  score::AAPairSidechainInteraction::WeightOfInteraction( ( **aa_1_itr), ( **aa_2_itr))
                );

                // add amino acid pair distance to the map
                biol::AAType first_type( std::min( ( *aa_1_itr)->GetType(), ( *aa_2_itr)->GetType()));
                biol::AAType second_type( std::max( ( *aa_1_itr)->GetType(), ( *aa_2_itr)->GetType()));
                if( seq_itr == seq_2_itr)
                {
                  distance_map_one_old[ std::pair< biol::AAType, biol::AAType>( first_type, second_type)]->PushBack( aa_aa_distance);
                  distance_map_one_new[ std::pair< biol::AAType, biol::AAType>( first_type, second_type)]->PushBack( aa_aa_distance, aa_aa_distance_weight);
                }
                else
                {
                  distance_map_all_old[ std::pair< biol::AAType, biol::AAType>( first_type, second_type)]->PushBack( aa_aa_distance);
                  distance_map_all_new[ std::pair< biol::AAType, biol::AAType>( first_type, second_type)]->PushBack( aa_aa_distance, aa_aa_distance_weight);
                }
              } // end of iterating the current chain to get the second amino acid residues
            } // end of iterating the current chain and all other chains
          } // end of iteration over the current chain to get the first amino acid residues
        } // end of iteration over all chains
      } // end of iteration over all models

      // write statistics
      std::ostringstream stream;

      // get the first residue
      for
      (
        biol::AATypes::const_iterator aa_1_itr( biol::GetAATypes().Begin()), aa_itr_end( biol::GetAATypes().End());
          aa_1_itr != aa_itr_end;
        ++aa_1_itr
      )
      {
        // get the second residue
        for
        (
          biol::AATypes::const_iterator aa_2_itr( aa_1_itr);
            aa_2_itr != aa_itr_end;
          ++aa_2_itr
        )
        {
          const std::pair< biol::AAType, biol::AAType> aa_pair( *aa_1_itr, *aa_2_itr);

          if( m_CategoryOption == e_AllChainNew)
          {
            // check to see if distance_map_all_new has the entry
            const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> >::const_iterator
              distance_map_all_chain_new_entry( distance_map_all_new.find( aa_pair));

            if( distance_map_all_chain_new_entry != distance_map_all_new.end() && !distance_map_all_chain_new_entry->second->IsEmpty())
            {
              stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
              stream << *distance_map_all_chain_new_entry->second;
            }
          }
          else if( m_CategoryOption == e_AllChainOld)
          {
            // check to see if distance_map_all_old has the entry
            const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> >::const_iterator
              distance_map_all_old_entry( distance_map_all_old.find( aa_pair));

            if( distance_map_all_old_entry != distance_map_all_old.end() && !distance_map_all_old_entry->second->IsEmpty())
            {
              stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
              stream << *distance_map_all_old_entry->second;
            }
          }
          else if( m_CategoryOption == e_OneChainNew)
          {
            // check to see if distance_map_one_new has the entry
            const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> >::const_iterator
              distance_map_one_new_entry( distance_map_one_new.find( aa_pair));

            if( distance_map_one_new_entry != distance_map_one_new.end() && !distance_map_one_new_entry->second->IsEmpty())
            {
              stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
              stream << *distance_map_one_new_entry->second;
            }
          }
          else if( m_CategoryOption == e_OneChainOld)
          {
            // check to see if distance_map_one_old has the entry
            const std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< math::Histogram> >::const_iterator
              distance_map_one_old_entry( distance_map_one_old.find( aa_pair));

            if( distance_map_one_old_entry != distance_map_one_old.end() && !distance_map_one_old_entry->second->IsEmpty())
            {
              stream << ( *aa_1_itr)->GetOneLetterCode() << " " << ( *aa_2_itr)->GetOneLetterCode() << "\n";
              stream << *distance_map_one_old_entry->second;
            }
          }
        }
      }

      return stream.str();

    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AADistance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes amino acid pair distance statistics.");

      parameters.AddInitializer
      (
        "category",
        "the distance category: OneChainNew, OneChainOld, AllChainNew, AllChainOld",
        io::Serialization::GetAgent( &m_CategoryOption),
        "OneChainNew"
      );

      parameters.AddInitializer
      (
        "bin_size",
        "the bin size for the histogram",
        io::Serialization::GetAgent( &m_BinSize),
        "1"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "seq_exlusion",
        "sequence exclusion: sequence separation below which aa distance is not considered",
        io::Serialization::GetAgent( &m_AADistSeqExcl),
        "6"
      );

      return parameters;
    }
  } // namespace scorestat
} // namespace bcl

