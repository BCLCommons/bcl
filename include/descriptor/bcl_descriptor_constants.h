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

#ifndef BCL_DESCRIPTOR_CONSTANTS_H_
#define BCL_DESCRIPTOR_CONSTANTS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Constants
    //! @brief is an Interface derived class that returns a constant value independent of the argument passed
    //! @details This class can be used in global descriptors which are independent of are not solely determined by the
    //! argument It can also be used to provide an external information such as the oligomeric state of a protein.
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_constants.cpp @endlink
    //! @author mendenjl
    //! @date Jan 30, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Constants :
      public BaseSequence< t_DataType, t_ReturnType>
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector< t_ReturnType> m_Description;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Constants();

      //! @brief constructor from a description
      //! @param DESCRIPTION storage::Vector< t_ReturnType> that holds the constant description
      Constants( const t_ReturnType &DESCRIPTION);

      //! @brief constructor from a description
      //! @param DESCRIPTION storage::Vector< t_ReturnType> that holds the constant description
      Constants( const linal::Vector< t_ReturnType> &DESCRIPTION);

      //! @brief Clone function
      //! @return pointer to new Constants
      Constants< t_DataType, t_ReturnType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< t_ReturnType> &STORAGE);

    }; // class Constants

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Constants
    //! @brief Specialization for char return types
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_constants.cpp @endlink
    //! @author mendenjl
    //! @date Jan 30, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Constants< t_DataType, char> :
      public BaseSequence< t_DataType, char>
    {

    private:

    //////////
    // data //
    //////////

      std::string m_Description; //!< The constant values to return

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Constants()
      {
      }

      //! @brief Clone function
      //! @return pointer to new Constants
      Constants< t_DataType, char> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< char> &STORAGE);

    }; // class Constants

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Constants< chemistry::AtomConformationalInterface, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Constants< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Constants< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Constants< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Constants< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Constants< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Constants< char, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Constants< char, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_CONSTANTS_H_
