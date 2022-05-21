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

#ifndef BCL_DESCRIPTOR_WINDOW_H_
#define BCL_DESCRIPTOR_WINDOW_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_window_alignment_type.h"
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
    //! @class Window
    //! @brief generates descriptions for a window, util::SiPtrVector, of arguments.
    //! @details This class allows generating descriptions that correspond to a provided Interface derived class for a
    //! util::SiPtrVector of arguments. It iterates over all arguments in the vector, gets the descriptions for each
    //! argument using m_Interface and returns the concatenated results.
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_window.cpp @endlink
    //! @author mendenjl
    //! @date Jan 31, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Window :
      public BaseElement< t_DataType, t_ReturnType>
    {

    private:

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, t_ReturnType> > m_Descriptor; //!< Actual descriptor implementation
      size_t m_Size;        //!< Window radius
      size_t m_SizeLeft;    //!< Size of the window to the left, not counting center element
      size_t m_SizeRight;   //!< Size of the window to the right, not counting center element

      //! True - window should reflect if it hits either end of the sequence
      //! False - just return NaNs if the boundary is exceeded
      bool m_Reflective;

      //! Alignment of the window relative to the current element
      WindowAlignmentEnum m_Alignment;

    public:

    //////////
    // data //
    //////////

      //! instances of that class
      static const util::SiPtr< const util::ObjectInterface> s_ReflectiveInstance;
      static const util::SiPtr< const util::ObjectInterface> s_OpaqueInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor from whether to be springy at the end and the internal descriptor, if desired
      Window
      (
        const bool &REFLECTIVE,
        const util::Implementation< Base< t_DataType, t_ReturnType> > &DESCRIPTOR
          = ( util::Implementation< Base< t_DataType, t_ReturnType> >())
      );

      //! @brief constructor from a descriptor
      //! @param DESCRIPTOR descriptor to be used
      Window
      (
        const BaseElement< t_DataType, t_ReturnType> &DESCRIPTOR,
        const size_t &SIZE
      );

      //! @brief constructor from a descriptor
      //! @param DESCRIPTOR descriptor to be used
      Window
      (
        const util::Implementation< Base< t_DataType, t_ReturnType> > &DESCRIPTOR,
        const size_t &SIZE
      );

      //! @brief Clone function
      //! @return pointer to new Window
      Window< t_DataType, t_ReturnType> *Clone() const;

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

      //! @brief get the indices of the returned values, relative to the central position
      //! @return the indices of the returned values, relative to the central position
      storage::Vector< int> GetRelativeIndices() const;

      //! @brief Set the window size
      //! @param RADIUS the radius of the window
      void SetRadius( const size_t &RADIUS);

      //! @brief Set the alignment of the window
      //! @param ALIGNMENT the new alignment of the window
      void SetAlignment( const WindowAlignmentType &ALIGNMENT);

      //! @brief get whether the window reflects at the end points
      //! @return true if the window reflects at the end points
      bool IsReflective() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, t_ReturnType> > GetInternalDescriptors();

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const t_DataType> &ELEMENT,
        linal::VectorReference< t_ReturnType> &STORAGE
      );

    }; // class Window

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Window< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Window< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Window< char, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Window< char, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WINDOW_H_
