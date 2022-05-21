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

#ifndef BCL_DESCRIPTOR_NUMERIC_STRING_H_
#define BCL_DESCRIPTOR_NUMERIC_STRING_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "model/bcl_model_retrieve_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NumericString
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
    //! @see @link example_descriptor_numeric_string.cpp @endlink
    //! @author mendenjl
    //! @date Feb 17, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class NumericString :
      public Base< t_DataType, char>
    {
    private:

    //////////
    // data //
    //////////

      // property that returns the numeric value that will be converted to string
      util::Implementation< Base< t_DataType, float> > m_NumericProperty;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      NumericString *Clone() const;

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
      size_t GetNormalSizeOfFeatures() const;

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

      //! @brief hook that derived classes can override to add behavior after every time SetDimension is called
      void SetDimensionHook();

    private:

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< char> &STORAGE
      );

    }; // class NumericString

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API NumericString< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API NumericString< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API NumericString< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API NumericString< char>;
  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_NUMERIC_STRING_H_
