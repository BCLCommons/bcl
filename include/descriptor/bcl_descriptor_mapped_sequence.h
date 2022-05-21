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

#ifndef BCL_DESCRIPTOR_MAPPED_SEQUENCE_H_
#define BCL_DESCRIPTOR_MAPPED_SEQUENCE_H_

// include namespace header
#include "bcl_descriptor.h"

// include from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "biol/bcl_biol_aa_base.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MappedSequence
    //! @brief uses sequence id as a key to return neighbor count or neighbor vector for that key
    //! that is given in a file. Has a parameter the id - value delimiter, which separates the id
    //! property from the value.
    //!
    //! Input files must have at least one space/tab per row, with the key before the first space/tab.
    //! Everything after the first space/tab on each row will be the value.
    //!
    //!
    //! It is an error if the value property is missed for a particular key
    //! @see @link example_descriptor_mapped_sequence.cpp @endlink
    //! @author lib14
    //! @date March 5, 2015
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType>
    class MappedSequence :
      public Base< biol::AABase, t_ReturnType>
    {
    private:

    //////////
    // data //
    //////////

      // extension of the file to be read in
      std::string m_FileExtension;

      // ID/Key delimiter for the file
      char m_IdKeyDelimiter;

      // number of values per key
      size_t m_NumberOutputs;

      // implementation of the property that retrieves the key property for the protein
      util::Implementation< Base< biol::AABase, char> > m_KeyProperty;

      // map of protein 5-letter id to neighbor count/neighbor vector data file
      storage::Map< std::string, linal::Vector< t_ReturnType> > m_Map;

    public:

      // single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MappedSequence();

      //! @brief virtual copy constructor
      MappedSequence *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return constant reference to the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns data label
      //! @return constant reference to data label
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const;

      //! @brief return parameters for data members that are set up from the labels
      //! @return parameters for data members that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
      //! since the results are often in the cache
      void LoadFile();

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

    private:

    ////////////////
    // operations //
    ////////////////

      //! @brief clean the given key, based on the delimiter
      //! @param KEY the given key
      //! @return the key, duplicate spaces removed if present; all tabs converted to spaces
      std::string CleanKey( const std::string &KEY) const;

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< biol::AABase> &ITR,
        linal::VectorReference< t_ReturnType> &STORAGE
      );

      //! @brief given a line, read the next key and value into passed-in objects
      //! @param LINE the next line
      //! @param KEY storage for the key
      //! @param VALUES storage for the values
      //! @param ERR_STREAM stream to write errors to
      //! @return true if a value and key were each read
      bool ReadKeyAndValue
      (
        const std::string &LINE,
        std::string &KEY,
        linal::Vector< t_ReturnType> &Value,
        std::ostream &ERR_STREAM
      ) const;
    }; // end of class MappedSequence

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API MappedSequence< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API MappedSequence< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MAPPED_SEQUENCE_H_
