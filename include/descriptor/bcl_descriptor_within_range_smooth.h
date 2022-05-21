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

#ifndef BCL_DESCRIPTOR_WITHIN_RANGE_SMOOTH_H_
#define BCL_DESCRIPTOR_WITHIN_RANGE_SMOOTH_H_

// include the namespace header
#include "bcl_descriptor.h"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "math/bcl_math_trigonometric_transition.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WithinRangeSmooth
    //! @details WithinRangeSmooth returns values between 0 and 1 of whether a descriptor's mean value
    //! @details is within a range, with a sinusoidal interpolation region on the ends (smooth cutoff)
    //! @details                ______
    //! @details looks like ___/      \____
    //!
    //! @see @link example_descriptor_within_range_smooth.cpp @endlink
    //! @author geanesar
    //! @date Feb 05, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class WithinRangeSmooth :
      public Base< t_DataType, float>
    {
    private:

    //////////
    // data //
    //////////

      //! the atom property encode in the 3D autocorrelation function
      util::Implementation< Base< t_DataType, float> > m_Descriptor;

      //! the beginning of the range
      float m_RangeBegin;

      //! the end of the range
      float m_RangeEnd;

      //! how wide the interpolation region is on the left
      float m_LeftWidth;

      //! how wide the interpolation region is on the right
      float m_RightWidth;

      //! whether the range is inclusive of the endpoints
      bool m_InclusiveRange;

      //! The trigonometric transition functions for the left side
      math::TrigonometricTransition m_LeftTransition;

      //! The trigonometric transition functions for the right side
      math::TrigonometricTransition m_RightTransition;

    public:

      //! single ascending instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      WithinRangeSmooth();

      //! @brief constructor
      //! @param RANGE_BEGIN where return values of 1 should start
      //! @param RANGE_END where return values of 1 should end
      //! @param LEFT_INTERP_BEGIN where sinusoidal interpolation should start to the left of the range
      //! @param RIGHT_INTERP_END where sinusoidal interpolation should end to the right of the range
      WithinRangeSmooth
      (
        const float &RANGE_BEGIN,
        const float &RANGE_END,
        const float &LEFT_WIDTH = 0.0,
        const float &RIGHT_WIDTH = 0.0,
        const bool &INCLUSIVE = true
      );

      //! @brief virtual copy constructor
      WithinRangeSmooth *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief set descriptor to be queried
      //! @param DESCRIPTOR the descriptor to use
      void SetDescriptor( const util::Implementation< Base< t_DataType, float> > &DESCRIPTOR);

      //! @brief get the descriptor that is queried
      //! @return an implementation of the descriptor
      const util::Implementation< Base< t_DataType, float> > &GetDescriptor() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void RecalculateImpl( const Iterator< t_DataType> &ITR, linal::VectorReference< float> &STORAGE);

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const;

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const;

    }; // class WithinRangeSmooth

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API WithinRangeSmooth< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WithinRangeSmooth< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WithinRangeSmooth< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WithinRangeSmooth< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WITHIN_RANGE_SMOOTH_H_
