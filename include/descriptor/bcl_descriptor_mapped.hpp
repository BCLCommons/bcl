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

// include header of this class
#include "bcl_descriptor_mapped.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector_reference.h"
#include "storage/bcl_storage_list.h"
#include "type/bcl_type_compare.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType, typename t_ReturnType>
    Mapped< t_DataType, t_ReturnType>::Mapped() :
      m_IdKeyDelimiter( '\0'),
      m_NumberOutputs( 0),
      m_Default()
    {
    }

    //! @brief virtual copy constructor
    template< typename t_DataType, typename t_ReturnType>
    Mapped< t_DataType, t_ReturnType> *Mapped< t_DataType, t_ReturnType>::Clone() const
    {
      return new Mapped( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Mapped< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Mapped< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_name( "Mapped");
      return s_name;
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    template< typename t_DataType, typename t_ReturnType>
    size_t Mapped< t_DataType, t_ReturnType>::GetNormalDimension() const
    {
      return m_KeyProperty->GetType().GetDimension();
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    template< typename t_DataType, typename t_ReturnType>
    Type::Symmetry Mapped< t_DataType, t_ReturnType>::GetSymmetry() const
    {
      return m_KeyProperty->GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    template< typename t_DataType, typename t_ReturnType>
    bool Mapped< t_DataType, t_ReturnType>::ConsiderRepeatedElements() const
    {
      return m_KeyProperty->ConsiderRepeatedElements();
    }

  //////////////////////
  // input and output //
  //////////////////////

    namespace
    {
      //! @brief read a string into a linal vector of chars
      //! @param STRING the string to read
      //! @param STORAGE the vector of characters
      void StringToStorage( const std::string &STRING, linal::Vector< char> &STORAGE)
      {
        STORAGE = linal::Vector< char>( STRING.begin(), STRING.end());
      }

      //! @brief read a string into a linal vector of floats
      //! @param STRING the string to read
      //! @param STORAGE the vector of floats
      void StringToStorage( const std::string &STRING, linal::Vector< float> &STORAGE)
      {
        // read in the value into a vector
        storage::Vector< float> values_storage( util::SplitStringToNumerical< float>( STRING, ",\n\t\r "));
        // copy into the linal vector
        STORAGE = linal::Vector< float>( values_storage.Begin(), values_storage.End());
      }

      //! @brief read a string into a linal vector of chars
      //! @param STRING the string to read
      //! @param STORAGE the vector of characters
      void Copy( const std::string &STRING, linal::VectorReference< char> &STORAGE)
      {
        if( STRING.size() == STORAGE.GetSize())
        {
          std::copy( STRING.begin(), STRING.end(), STORAGE.Begin());
        }
        else
        {
          STORAGE = char( 0);
        }
      }

      //! @brief read a string into a linal vector of floats
      //! @param STRING the string to read
      //! @param STORAGE the vector of floats
      void Copy( const linal::Vector< float> &VALUES, linal::VectorReference< float> &STORAGE)
      {
        if( VALUES.GetSize() == STORAGE.GetSize())
        {
          STORAGE.CopyValues( VALUES);
        }
        else
        {
          STORAGE = util::GetUndefined< float>();
        }
      }
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType, typename t_ReturnType>
    bool Mapped< t_DataType, t_ReturnType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // reset the map
      m_MapPtr.Reset();

      // open the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);

      // check that the file was not empty
      if( input.eof())
      {
        ERR_STREAM << "Mapped file was empty!\n";
        return false;
      }

      // create a key and value string
      std::string key;
      linal::Vector< t_ReturnType> value;

      // create std::string "current_line" to hold the current line of the mapped file
      std::string current_line;

      // read in till the first non-blank line of input
      // skip empty lines; particularly important if the file was generated on Windows and read on Linux/Mac,
      // in which case every other call to getline returns an empty line
      while( !input.eof() && current_line.find_first_not_of( " \t\r") == std::string::npos)
      {
        // get the next line in "ISTREAM" and put it into "current_line"
        std::getline( input, current_line);
      }

      // read the key and value
      if( !ReadKeyAndValue( current_line, key, value, ERR_STREAM))
      {
        return false;
      }
      m_NumberOutputs = value.GetSize();
      io::File::CloseClearFStream( input);

      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType, typename t_ReturnType>
    io::Serializer Mapped< t_DataType, t_ReturnType>::GetSerializer() const
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
        "file",
        "file with every line containing a key followed by a value",
        io::Serialization::GetAgentInputFilename( &m_Filename)
      );

      parameters.AddInitializer
      (
        "key",
        "descriptor key (id type) to calculate the key for a given input",
        io::Serialization::GetAgent( &m_KeyProperty)
      );

      parameters.AddOptionalInitializer
      (
        "delimiter",
        "Delimiter between the id and key. "
        "If not given, the ids must be fixed width, based on the # of characters as the given key descriptor"
        "If the delimiter is given, the read-in string will be tokenized, such that repeated spacing is ignored",
        io::Serialization::GetAgent( &m_IdKeyDelimiter)
      );

      parameters.AddOptionalInitializer
      (
        "default",
        "Value(s) returned if the key is not present in the file",
        io::Serialization::GetAgent( &m_Default)
      );
      return parameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_DataType, typename t_ReturnType>
    void Mapped< t_DataType, t_ReturnType>::SetObjectHook()
    {
      m_KeyProperty->SetObject( *this->GetCurrentObject());
    }

    //! @brief clean the given key, based on the delimiter
    //! @param KEY the given key
    //! @return the key, duplicate spaces removed if present; all tabs converted to spaces
    template< typename t_DataType, typename t_ReturnType>
    std::string Mapped< t_DataType, t_ReturnType>::CleanKey( const std::string &KEY) const
    {
      // handle null key, in which case the id property is fixed width
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
        for( size_t i( 0), key_size( key_cleaned.size()); i < key_size; ++i)
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
    template< typename t_DataType, typename t_ReturnType>
    void Mapped< t_DataType, t_ReturnType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< t_ReturnType> &STORAGE
    )
    {
      if( !m_MapPtr.IsDefined())
      {
        LoadFile();
      }
      // compute the key property
      const linal::VectorConstReference< char> key_linal( m_KeyProperty->operator()( ITR));

      // create a string with the key, pre-cleaned
      std::string key( CleanKey( std::string( key_linal.Begin(), key_linal.GetSize())));

      // check for the key in the map
      typename storage::Map< std::string, linal::Vector< t_ReturnType> >::const_iterator
        itr( m_MapPtr->Find( key));

      // handle the case that the key is located in the map
      if( itr != m_MapPtr->End())
      {
        STORAGE.CopyValues( itr->second);
      }
      else
      {
        BCL_MessageStd( key + " was not found in " + m_Filename + " available keys: " + util::Join( ", ", m_MapPtr->GetKeysAsVector()));
        Copy( m_Default, STORAGE);
      }
    }

    //! @brief given a line, read the next key and value into passed-in objects
    //! @param LINE the next line
    //! @param KEY storage for the key
    //! @param VALUES storage for the values
    //! @return true iff a value and key were each read
    template< typename t_DataType, typename t_ReturnType>
    bool Mapped< t_DataType, t_ReturnType>::ReadKeyAndValue
    (
      const std::string &LINE,
      std::string &KEY,
      linal::Vector< t_ReturnType> &VALUES,
      std::ostream &ERR_STREAM
    ) const
    {
      std::string values;
      if( m_IdKeyDelimiter == '\0')
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

        // extract everything else
        values = LINE.substr( m_KeyProperty->GetSizeOfFeatures());
      }
      else
      {
        // non-null delimiter, find it on the line
        const size_t delimiter_pos( LINE.find( m_IdKeyDelimiter));
        if( delimiter_pos == std::string::npos)
        {
          ERR_STREAM << "Delimiter " << m_IdKeyDelimiter << " was not found on line: " << LINE << '\n';
          return false;
        }

        // extract the key
        KEY = CleanKey( LINE.substr( 0, delimiter_pos));

        // extract everything else
        if( delimiter_pos + 1 == LINE.size())
        {
          ERR_STREAM << "Found line with key but no values: " << LINE << '\n';
          return false;
        }
        values = LINE.substr( delimiter_pos + 1);
      }

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

    //! @brief Load the file
    template< typename t_DataType, typename t_ReturnType>
    void Mapped< t_DataType, t_ReturnType>::LoadFile()
    {
      s_FilesMutex.Lock();
      storage::Map< std::string, linal::Vector< t_ReturnType> > &map( s_Files[ m_Filename]);
      m_MapPtr = util::SiPtr< const storage::Map< std::string, linal::Vector< t_ReturnType> > >( map);

      if( !map.IsEmpty())
      {
        // file was already loaded
        s_FilesMutex.Unlock();
        return;
      }

      // open the file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);

      // create a key and value string
      std::string key;
      linal::Vector< t_ReturnType> value;

      // check that the file was not empty
      BCL_Assert( !input.eof(), "Mapped file was empty!\n");

      std::ostringstream errs;

      // read in all the line of "ISTREAM"
      while( !input.eof())
      {
        // create std::string "current_line" to hold the current line of the accessibility file
        std::string current_line;

        // get the next line in "ISTREAM" and put it into "current_line"
        std::getline( input, current_line);

        // skip empty lines; particularly important if the file was generated on Windows and read on Linux/Mac,
        // in which case every other call to getline returns an empty line
        if( current_line.find_first_not_of( " \t") == std::string::npos)
        {
          continue;
        }

        // read the key and value
        BCL_Assert( ReadKeyAndValue( current_line, key, value, errs), errs.str());
        if( map.Has( key) && map[ key] != value)
        {
          BCL_MessageVrb
          (
            "Warning, Mapped file: " + m_Filename + " contains lines with different values for key: "
            + key + "; using last reported values (" + current_line + ")"
          );
        }
        map[ key] = value;
      }
      s_FilesMutex.Unlock();
      io::File::CloseClearFStream( input);
    }

  } // namespace descriptor
} // namespace bcl
