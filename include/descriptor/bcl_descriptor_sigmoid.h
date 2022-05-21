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

#ifndef BCL_DESCRIPTOR_SIGMOID_H_
#define BCL_DESCRIPTOR_SIGMOID_H_

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
    //! @class Sigmoid
    //! @details Sigmoid a function of the form Amplitude / ( 1 + exp(-(x-x_offset)/slope)) + y_offset
    //!          For ease of user input, slope is computed by asking for a point on the curve other than that x_offset
    //! @see @link example_descriptor_sigmoid.cpp @endlink
    //! @author mendenjl
    //! @date Jun 21, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Sigmoid :
      public Base< t_DataType, float>
    {
    private:

    //////////
    // data //
    //////////

      //! the atom property encode in the 3D autocorrelation function
      util::Implementation< Base< t_DataType, float> > m_Descriptor;

      //! amplitude of the sigmoid
      float m_Amplitude;

      //! x-offset
      float m_XOffset;

      //! y-offset
      float m_YOffset;

      //! slope
      float m_Slope;

      //! any X on the curve other than at x-offset (this is used to compute the slope)
      float m_SlopeX;

      //! whether the range is inclusive of the endpoints
      float m_YAtSlopeX;

    public:

      //! single ascending instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Sigmoid();

      //! @brief constructor
      //! @param AMPLITUDE amplitude of the sigmoid
      //! @param X_OFFSET offset of the curve from 0 on x-axis; at this value, the sigmoid is 0.5
      //! @param Y_OFFSET offset of the curve from 0 on the y-axis
      //! @param SLOPE_X, SLOPE_Y any value on the curve
      Sigmoid
      (
        const float &AMPLITUDE,
        const float &X_OFFSET,
        const float &Y_OFFSET,
        const float &SLOPE_X,
        const float &SLOPE_Y
      );

      //! @brief virtual copy constructor
      Sigmoid *Clone() const;

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

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const;

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const;

    }; // class Sigmoid

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Sigmoid< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Sigmoid< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Sigmoid< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Sigmoid< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SIGMOID_H_
