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

// include header file
#include "bcl_descriptor_mapped_sequence.h"
// include from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_reference.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_ReturnType>
    MappedSequence< t_ReturnType>::MappedSequence() :
      m_IdKeyDelimiter( '\0'),
      m_NumberOutputs( 0)
    {
    }

    //! @brief virtual copy constructor
    template< typename t_ReturnType>
    MappedSequence< t_ReturnType> *MappedSequence< t_ReturnType>::Clone() const
    {
      return new MappedSequence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return constant reference to the class name
    template< typename t_ReturnType>
    const std::string &MappedSequence< t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns data label
    //! @return constant reference to data label
    template< typename t_ReturnType>
    const std::string &MappedSequence< t_ReturnType>::GetAlias() const
    {
      static const std::string s_name( "MappedSequence");
      return s_name;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_ReturnType>
    size_t MappedSequence< t_ReturnType>::GetNormalSizeOfFeatures() const
    {
      return m_NumberOutputs;
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    template< typename t_ReturnType>
    size_t MappedSequence< t_ReturnType>::GetNormalDimension() const
    {
      return m_KeyProperty->GetType().GetDimension();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_ReturnType>
    io::Serializer MappedSequence< t_ReturnType>::GetSerializer() const
    {

      io::Serializer parameters;
      parameters.SetClassDescription
      (
         "Retrieves data from a file, indexed by a key property."
         "\n    uses an id descriptor as a key to return the value(s) following that key on a row in the file file"
         "\n    Input files have rows with a key, followed by spaces or tabs, followed by a value string, e.g.:"
         "\n    2   5.0 6.0 12.0"
         "\n    1   3.0, 14.0, 18.0"
         "\n    If using this descriptor in multi-dimensional mode, be sure to specify a delimiter other than space"
        +
        type::Compare< t_ReturnType, float>::e_Same
        ? std::string( "\n    Note that commas are automatically stripped from the given inputs")
        : std::string()
      );

      parameters.AddInitializer
      (
        "file extension",
        "extension of a file with every line containing a key followed by a value",
        io::Serialization::GetAgent( &m_FileExtension)
      );

      parameters.AddInitializer
      (
        "key",
        "descriptor key (id type) to calculate the key for a given input",
        io::Serialization::GetAgent( &m_KeyProperty)
      );

      parameters.AddInitializer
      (
        "size",
        "size of the current descriptor",
        io::Serialization::GetAgent( &m_NumberOutputs)
      );

      parameters.AddOptionalInitializer
      (
        "delimiter",
        "Delimiter between the id and key. "
        "If not given, the ids must be fixed width, based on the # of characters as the given key descriptor"
        "If the delimiter is given, the read-in string will be tokenized, such that repeated spacing is ignored",
        io::Serialization::GetAgent( &m_IdKeyDelimiter)
      );

      return parameters;

    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
    //! since the results are often in the cache
    template< typename t_ResultType>
    void MappedSequence< t_ResultType>::LoadFile()
    {
      // reset the storage map
      m_Map.Reset();

      // get current object (protein chain)
      const util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model_with_cache( this->GetCurrentObject());

      // get the filename
      const util::ShPtr< util::Wrapper< std::string> > sp_filename_wrapper
      (
        sp_protein_model_with_cache->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );

      const std::string pdb_filename( sp_filename_wrapper->GetData());

      // get the basename
      const std::string basename( io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( pdb_filename)));

      // add user provided filename extension
      const std::string filename( basename + m_FileExtension);

      // open the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, filename);

      // create key and value to store data read from each line
      std::string key;
      linal::Vector< t_ResultType> value;

      // skip empty files
      if( input.eof()) {
        BCL_MessageDbg( util::Format()( filename + " is empty."));
        return;
      }

      // ignore comment lines
      util::ChopHeader( input);

      // scan in the file line by line
      while( !input.eof())
      {
        // get next line from input stream and store it in current line
        std::string current_line;
        std::getline( input, current_line);

        // skip empty lines
        if( current_line.find_first_of( m_IdKeyDelimiter) == std::string::npos)
        {
          BCL_MessageDbg( util::Format()( "Skip empty line."));
          continue;
        }

        // get key and value from current line
        if( !ReadKeyAndValue( current_line, key, value, util::GetLogger()))
        {
          BCL_Exit( "Could not read the current line: " + current_line, -1);
        }
        m_NumberOutputs = value.GetSize();
        m_Map[ key] = value;
      }
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_ReturnType>
    void MappedSequence< t_ReturnType>::SetObjectHook()
    {
      m_KeyProperty->SetObject( *this->GetCurrentObject());
      m_Map.Reset();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief clean the given key, based on the delimiter
    //! @param KEY the given key
    //! @return the key, duplicate spaces removed if present; all tabs converted to spaces
    template< typename t_ReturnType>
    std::string MappedSequence< t_ReturnType>::CleanKey( const std::string &KEY) const
    {
      if( m_IdKeyDelimiter == '\0')
      {
        return KEY;
      }

      std::string key_cleaned( util::TrimString( KEY));

      // replace tabs with spaces
      if( key_cleaned.find_first_of( '\t') != std::string::npos)
      {
        std::replace( key_cleaned.begin(), key_cleaned.end(), '\t', ' ');
      }

      // replace duplicate spaces with a single space
      if( key_cleaned.find( "  ") != std::string::npos)
      {
        std::string new_key;
        new_key.reserve( key_cleaned.size());
        new_key += key_cleaned[ 0];
        for( size_t i( 1), key_size( key_cleaned.size()); i < key_size; ++i)
        {
          if( key_cleaned[ i] != ' ' || new_key[ new_key.size() - 1] != ' ')
          {
            new_key += key_cleaned[ i];
          }
        }
        key_cleaned = new_key;
      }

      return key_cleaned;
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_ReturnType>
    void MappedSequence< t_ReturnType>::RecalculateImpl
    (
      const Iterator< biol::AABase> &ITR,
      linal::VectorReference< t_ReturnType> &STORAGE
    )
    {
      // compute the key property
      const linal::VectorConstReference< char> key_linal( m_KeyProperty->operator ()( ITR));

      // create a string with the key, pre-cleaned
      std::string key( CleanKey( std::string( key_linal.Begin(), key_linal.GetSize())));

      // if the sequence maps are empty, this indicates that this object is now operating over a new sequence,
      // so it is necessary to reload the files
      if( m_Map.IsEmpty())
      {
        LoadFile();
      }

      // check for the key in the map
      typename storage::Map< std::string, linal::Vector< t_ReturnType> >::const_iterator itr( m_Map.Find( key));

      // handle the case that the key is located in the map
      if( itr != m_Map.End())
      {
        STORAGE.CopyValues( itr->second);
      }
      else
      {
        STORAGE = util::GetUndefined< t_ReturnType>();
      }
    }

    namespace
    {

      //! @brief read a string into a linal vector of chars
      //! @param STRING the string to read
      //! @param STORAGE the vector of characters
      void StringToStorage( const std::string &STRING, linal::Vector< char> &STORAGE)
      {
        STORAGE = ( linal::Vector< char>( STRING.begin(), STRING.end()));
      }

      //! @brief read a string into a linal vector of floats
      //! @param STRING the string to read
      //! @param STORAGE the vector of floats
      void StringToStorage( const std::string &STRING, linal::Vector< float> &STORAGE)
      {
        // read in the value into a vector
        const storage::Vector< float> values_storage( util::SplitStringToNumerical< float>( STRING, ",\n\t\r "));

        // copy into the linal vector
        STORAGE = linal::Vector< float>( values_storage.Begin(), values_storage.End());
      }
    }

    //! @brief given a line, read the next key and value into passed-in objects
    //! @param LINE the next line
    //! @param KEY storage for the key
    //! @param VALUES storage for the values
    //! @param ERR_STREAM stream to write errors to
    //! @return true if a value and key were each read
    template< typename t_ReturnType>
    bool MappedSequence< t_ReturnType>::ReadKeyAndValue
    (
      const std::string &LINE,
      std::string &KEY,
      linal::Vector< t_ReturnType> &VALUES,
      std::ostream &ERR_STREAM
    ) const
    {
      // read in values as string
      std::string values;
      if( m_IdKeyDelimiter == '\0') // null delimiter
      {
        // fixed width
        if( LINE.size() < m_KeyProperty->GetSizeOfFeatures())
        {
          ERR_STREAM << "Key on line " << LINE << " is too short, should be exactly "
                     << m_KeyProperty->GetSizeOfFeatures() << " characters\n";
          return false;
        }

        // extract the key
        KEY = LINE.substr( 0, m_KeyProperty->GetSizeOfFeatures());

        // get everything else
        values = LINE.substr( m_KeyProperty->GetSizeOfFeatures());
      }
      else // other delimiters
      {
        const size_t delimiter_pos( LINE.find( m_IdKeyDelimiter));
        if( delimiter_pos == std::string::npos)
        {
          ERR_STREAM << "Delimiter " << m_IdKeyDelimiter << " was not found on line: " << LINE << '\n';
          return false;
        }

        // extract the key
        KEY = CleanKey( LINE.substr( 0, delimiter_pos));

        // extract everything else
        if( delimiter_pos + 1 == std::string::npos)
        {
          ERR_STREAM << "Found line with key but no values on line: " << LINE << '\n';
          return false;
        }

        values = LINE.substr( delimiter_pos + 1);
      }

        // convert string to storage
        StringToStorage( values, VALUES);
        if( VALUES.GetSize() != m_NumberOutputs && m_NumberOutputs)
        {
          if( type::Compare< t_ReturnType, float>::e_Same)
          {
            ERR_STREAM << "Found line with incorrect # of values, should have had ";
          }
          else // if( type::Compare< t_ReturnType, char>::e_Equal)
          {
            ERR_STREAM << "Found line with incorrect # of value characters, should have had ";
          }
          ERR_STREAM << m_NumberOutputs << " as the first line did, but instead had " << VALUES.GetSize() << '\n';
          ERR_STREAM << "Line was: " << LINE << '\n';
          return false;
        }
      return true;
    }
  } // namespace descriptor
} // namespace bcl

