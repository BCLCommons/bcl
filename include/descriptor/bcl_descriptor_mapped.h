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

#ifndef BCL_DESCRIPTOR_MAPPED_H_
#define BCL_DESCRIPTOR_MAPPED_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "model/bcl_model_retrieve_interface.h"
#include "sched/bcl_sched_mutex.h"
#include "type/bcl_type_chooser.h"
#include "type/bcl_type_compare.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Mapped
    //! @brief uses one property value as a key to return a value for that key that is given in a file
    //! Has as a parameter the id - value delimiter, which separates the id property from the value
    //! If not given, the id's must be fixed width and have exactly the same # of characters as the given id property
    //! If the delimiter is given, however, the id columns are interpreted as a string of words.
    //! Commas are allowed in the value section for numeric types, and will be silently stripped from input.
    //!
    //! Input files must have at least one space/tab per row, with the key before the first space/tab
    //! Everything after the first space/tab on each line will be the value
    //! Note that if t_ReturnType is char, all outputs will be padded to the length of the longest row
    //! in the first row, e.g.
    //! ID  Value
    //! 2   5.0 6.0 12.0
    //! 1   3.0 14.0 18.0
    //!
    //! All properties can be strings too, e.g.:
    //! Letter  NextLetter
    //! Alpha   Beta
    //! Beta    Gamma
    //! Gamma   Phi
    //!
    //! It is an error if the value property is omitted for a particular key
    //! For multi-dimensional descriptors, the input file keys should be in identical order to those returned from
    //! GenerateDataSet with the same id descriptor set
    //!
    //! @see @link example_descriptor_mapped.cpp @endlink
    //! @author mendenjl
    //! @date Mar 14, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Mapped :
      public Base< t_DataType, t_ReturnType>
    {
    private:

      typedef typename type::Chooser< type::Compare< t_ReturnType, float>::e_Same, linal::Vector< float>, std::string>::Type
          t_Aggregation;

    //////////
    // data //
    //////////

      std::string m_Filename; //!< name of the file containing the table

      char        m_IdKeyDelimiter; //!< ID/Key delimiter for the file

      size_t      m_NumberOutputs; //!< # of t_ReturnTypes returned per t_DataType iterator

      t_Aggregation m_Default; //!< Default value, returned if the key is not found

      // implementation of the property that retrieves the key property for the molecule
      util::Implementation< Base< t_DataType, char> > m_KeyProperty;

      util::SiPtr< const storage::Map< std::string, linal::Vector< t_ReturnType> > > m_MapPtr; //!< Map from key to value

      //! Mutex to protect access to s_Files
      static sched::Mutex s_FilesMutex;

      //! Map from filename to key-value map
      static storage::Map< std::string, storage::Map< std::string, linal::Vector< t_ReturnType> > > s_Files;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Mapped();

      //! @brief virtual copy constructor
      Mapped *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_NumberOutputs;
      }

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const;

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const;

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

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
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< t_ReturnType> &STORAGE
      );

      //! @brief given a line, read the next key and value into passed-in objects
      //! @param LINE the next line
      //! @param KEY storage for the key
      //! @param VALUES storage for the values
      //! @param ERR_STREAM stream to write errors to
      //! @return true iff a value and key were each read
      bool ReadKeyAndValue
      (
        const std::string &LINE,
        std::string &KEY,
        linal::Vector< t_ReturnType> &VALUE,
        std::ostream &ERR_STREAM
      ) const;

      //! @brief Load the file
      void LoadFile();

    }; // class Mapped

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Mapped< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Mapped< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Mapped< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Mapped< char, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Mapped< chemistry::AtomConformationalInterface, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Mapped< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Mapped< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Mapped< char, char>;
  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MAPPED_H_
