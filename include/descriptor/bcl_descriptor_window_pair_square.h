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

#ifndef BCL_DESCRIPTOR_WINDOW_PAIR_SQUARE_H_
#define BCL_DESCRIPTOR_WINDOW_PAIR_SQUARE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_pair.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WindowPairSquare
    //! @brief generates descriptions for a square window about two points
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //! @tparam t_ReturnType element type returned by the descriptor
    //!
    //! @see @link example_descriptor_window_pair_square.cpp @endlink
    //! @author mendenjl
    //! @date Feb 04, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class WindowPairSquare :
      public BasePair< t_DataType, t_ReturnType>
    {

    private:

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, t_ReturnType> > m_Descriptor;  //!< Actual descriptor implementation
      size_t                                                  m_Size;        //!< Expected size of the window

    public:

    //////////
    // data //
    //////////

      //! instances of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      WindowPairSquare();

      //! @brief constructor from a descriptor
      //! @param DESCRIPTOR descriptor to be used
      WindowPairSquare
      (
        const BasePair< t_DataType, t_ReturnType> &DESCRIPTOR,
        const size_t &SIZE
      );

      //! @brief constructor from a descriptor
      //! @param DESCRIPTOR descriptor to be used
      WindowPairSquare
      (
        const util::Implementation< Base< t_DataType, t_ReturnType> > &DESCRIPTOR,
        const size_t &SIZE
      );

      //! @brief Clone function
      //! @return pointer to new WindowPairSquare
      WindowPairSquare< t_DataType, t_ReturnType> *Clone() const;

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

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const;

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const;

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

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, t_ReturnType> > GetInternalDescriptors();

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT_A, ELEMENT_B: the element pair of interest
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        const iterate::Generic< const t_DataType> &ELEMENT_A,
        const iterate::Generic< const t_DataType> &ELEMENT_B,
        linal::VectorReference< t_ReturnType> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    }; // class WindowPairSquare

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API WindowPairSquare< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WindowPairSquare< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WindowPairSquare< char, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WindowPairSquare< char, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WINDOW_PAIR_SQUARE_H_
