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
#include "descriptor/bcl_descriptor_aa_position_in_tm.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AAPositionInTM::s_AATMPositionInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AAPositionInTM()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAPositionInTM::AAPositionInTM()
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AAPositionInTM *AAPositionInTM::Clone() const
    {
      return new AAPositionInTM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAPositionInTM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAPositionInTM::GetAlias() const
    {
      static const std::string s_name( "AA_PositionInTM");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAPositionInTM::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A the element of interest
    //! @param STORAGE storage for the descriptor
    void AAPositionInTM::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if the sequence maps are empty, this indicates that this object is now operating over a new sequence,
      // so it is necessary to reload the files
      if( m_TMPosVec.IsEmpty())
      {
        LoadFiles();
      }

      const int seq_id( ELEMENT->GetData()->GetSeqID());
      storage::Map< int, size_t>::const_iterator map_itr( m_SeqIDToIndexMap.Find( seq_id));

      if( map_itr == m_SeqIDToIndexMap.End())
      {
        STORAGE( 0) = util::GetUndefined< float>();
        return;
      }

      const size_t vec_index( map_itr->second);
      STORAGE( 0) = m_TMPosVec( vec_index);

      return;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AAPositionInTM::SetObjectHook()
    {
      m_SeqIDToIndexMap.Reset();
      m_TMPosVec.Reset();
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference AAPositionInTM::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
    //! since the results are often in the cache
    void AAPositionInTM::LoadFiles()
    {
      // Create the extension
      const std::string extension( ".octo_topo");
      // Get protein model
      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

      // Check to make sure there is only a single chain in the given PDB file
      BCL_Assert
      (
        sp_protein_model->GetChains().GetSize() == 1,
        "More than one chain in the given protein model: "
         + util::Format()( sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile))
      );

      // Get the sequences
      util::ShPtr< biol::AASequence> sp_seq( sp_protein_model->GetChains().FirstElement()->GetSequence());

      // Get the filename
      util::ShPtr< util::Wrapper< std::string> > sp_filename_wrapper
      (
        sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );
      std::string pdb_filename( sp_filename_wrapper->GetData());

      // Remove the last extension
      std::string basename( io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( pdb_filename)));

      // Add my extension
      std::string octopus_filename( basename + extension);
      // Read in the file's data
      io::IFStream instream;
      io::File::MustOpenIFStream( instream, octopus_filename);

      // Read in prediction string
      std::string prediction_str( ReadPredictionsForAASequence( instream));

      size_t seq_size( sp_seq->GetSize());
      size_t pred_size( prediction_str.size());

      // Check to make sure that the lengths of both the read octopus predictions and the chain are the same
      BCL_Assert
      (
        pred_size == seq_size,
        "Prediction size and protein sequence size are not equal. Sequence size = " + util::Format()( seq_size)
        + " Prediction size = " + util::Format()( pred_size) + " for filename: " + pdb_filename
      );
      size_t def_seq_ids_index( 0);

      // Iterate through the entire protein sequence checking seq IDs and SS single character predictions
      for
      (
        biol::AASequence::const_iterator itr( sp_seq->Begin()), itr_end( sp_seq->End());
        itr < itr_end;
        ++itr, ++def_seq_ids_index
      )
      {
        const int temp_seq_id( ( *itr)->GetData()->GetSeqID());
        BCL_Assert( util::IsDefined( temp_seq_id), "Seq ID is undefined for " + pdb_filename);

        // Add entry to map
        m_SeqIDToIndexMap.Insert( storage::Pair< int, size_t>( temp_seq_id, def_seq_ids_index));
      }

      // Set vector sizes based on max seq ID
      m_TMPosVec.AllocateMemory( prediction_str.size());

      // Set temp size vector size
      storage::Vector< double> size_vector( prediction_str.size());
      size_t index( 0);

      // Set up sizes vector based on predictions string
      for
      (
        std::string::const_iterator itr( prediction_str.begin()), itr_end( prediction_str.end());
        itr < itr_end;
      )
      {
        size_t next_index( std::min( prediction_str.size(), prediction_str.find_first_not_of( *itr, index)));
        size_t temp_size( next_index - index);

        // the desired size is 1 less, due to the normalization
        temp_size = temp_size + 1;

        // Iterate internally until the current SSE is finished
        for( size_t pos( 0); index < next_index; ++itr, ++index, ++pos)
        {
          size_vector( index) = temp_size;
        }
      }

      BCL_MessageDbg( "AAPositionInTM sizes vector: " + util::Format()( size_vector));

      // Set up index tracker
      index = 0;
      // Set up internal TM position
      size_t tm_pos_index( 0);
      // Set up last outside char
      char last_outchar( 'i');
      // set up the direction
      int direction( 0);

      // Iterate through prediction string setting values for TM position index vector
      for
      (
        std::string::const_iterator itr( prediction_str.begin()), itr_end( prediction_str.end());
        itr < itr_end;
        ++itr, ++index
      )
      {
        // If "M" divide by size of TM portion + 1 which is set within the vector
        if( *itr == 'M')
        {
          if( last_outchar == 'o')
          {
            tm_pos_index = size_vector( index);
            direction = -1;
          }
          else if( last_outchar == 'i')
          {
            tm_pos_index = 0;
            direction = 1;
          }
          tm_pos_index += direction;
          m_TMPosVec.PushBack( tm_pos_index / size_vector( index));
        }
        // If o set to 1.0
        else if( *itr == 'o')
        {
          m_TMPosVec.PushBack( 1.0);
        }
        // Otherwise assume i and set to 0.0
        else if( *itr == 'i')
        {
          m_TMPosVec.PushBack( 0.0);
        }
        last_outchar = *itr;
      }

      // DEBUG print out the final vectors so I can check them manually!
      BCL_MessageDbg( "TM Position vector: " + util::Format()( m_TMPosVec));
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM and save to member
    //! string with all i, M, or o predictions from Octopus
    //! @param ISTREAM input stream
    std::string AAPositionInTM::ReadPredictionsForAASequence( std::istream &ISTREAM) const
    {
      // initialize necessary variables
      std::string line;
      bool read_flag( false);
      std::string pred_string( "");

      // while reading lines
      while( std::getline( ISTREAM, line))
      {
        const std::string trimmed_line( util::TrimString( line));

        // if the line contains identifier to start reading
        if( !read_flag)
        {
          if( trimmed_line == "SPOCTOPUS predicted topology:" || trimmed_line == "OCTOPUS predicted topology:")
          {
            read_flag = true;
          }
        }
        // line contains predictions
        else if( read_flag)
        {
          const std::string legitimate_octopus_chars( "MRior");
          const std::string octopus_char_mapping( "Moioi");
          // iterate over the string
          for
          (
            std::string::const_iterator char_itr( trimmed_line.begin()),
              char_itr_end( trimmed_line.end());
            char_itr != char_itr_end; ++char_itr
          )
          {
            // if this is a transmembrane helix
            const size_t char_pos( legitimate_octopus_chars.find( *char_itr));
            if( char_pos < legitimate_octopus_chars.size())
            {
              pred_string += octopus_char_mapping[ char_pos];
            }
          }
        }
      }

      // end
      return pred_string;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPositionInTM::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        GetAlias() + " based on octopus predictions. Eventually, this will be extended to support any sse method."
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
