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
#include "bcl_descriptor_sigmoid.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_combine.h"
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_trigonometric_transition.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, accepts bool of whether to auto-compute mean
    template< typename t_DataType>
    Sigmoid< t_DataType>::Sigmoid() :
      m_Amplitude( 1),
      m_XOffset( 0),
      m_YOffset( 0),
      m_Slope( 1),
      m_SlopeX( 1),
      m_YAtSlopeX( 1.0 / ( 1.0 + exp( -1.0)))
    {
    }

    //! @brief default constructor, accepts bool of whether to auto-compute mean
    template< typename t_DataType>
    Sigmoid< t_DataType>::Sigmoid
    (
      const float &AMPLITUDE,
      const float &X_OFFSET,
      const float &Y_OFFSET,
      const float &SLOPE_X,
      const float &SLOPE_Y
    ) :
      m_Amplitude( AMPLITUDE),
      m_XOffset( X_OFFSET),
      m_YOffset( Y_OFFSET),
      m_Slope( 0),
      m_SlopeX( SLOPE_X),
      m_YAtSlopeX( SLOPE_Y)
    {
      // Do not use trigonometric transitions if left/right widths are zero
      std::stringstream oss;
      BCL_Assert( ReadInitializerSuccessHook( util::ObjectDataLabel(), oss), oss.str());
    }

    //! @brief copy constructor
    template< typename t_DataType>
    Sigmoid< t_DataType> *Sigmoid< t_DataType>::Clone() const
    {
      return new Sigmoid( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &Sigmoid< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &Sigmoid< t_DataType>::GetAlias() const
    {
      static const std::string s_Name( "Sigmoid");

      return s_Name;
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    template< typename t_DataType>
    size_t Sigmoid< t_DataType>::GetNormalDimension() const
    {
      return m_Descriptor->GetType().GetDimension();
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t Sigmoid< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor->GetSizeOfFeatures();
    }

    //! @brief set descriptor to be queried
    //! @param DESCRIPTOR the descriptor to use
    template< typename t_DataType>
    void Sigmoid< t_DataType>::SetDescriptor( const util::Implementation< Base< t_DataType, float> > &DESCRIPTOR)
    {
      m_Descriptor = DESCRIPTOR;
    }

    //! @brief get the descriptor that is queried
    //! @return an implementation of the descriptor
    template< typename t_DataType>
    const util::Implementation< Base< t_DataType, float> > &Sigmoid< t_DataType>::GetDescriptor() const
    {
      return m_Descriptor;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool Sigmoid< t_DataType>::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      if( m_XOffset == m_SlopeX)
      {
        ERR_STREAM << " The given x-y pair cannot be at the X-offset, else there is no way to compute the slope!";
        return false;
      }
      if( m_YAtSlopeX == m_YOffset)
      {
        ERR_STREAM << " The given x-y pair cannot be at the Y-offset, else slope is infinite!";
        return false;
      }
      const double y_relative_dist( m_Amplitude / ( m_YAtSlopeX - m_YOffset) - 1.0);
      if( y_relative_dist <= 0.0)
      {
        ERR_STREAM << " the given y-point (" << m_YAtSlopeX << " is not on the sigmoid!";
        return false;
      }

      // Y' = A / ( 1 + exp(-(x_slope-x_offset)/slope)) + y_offset
      // -(x_slope-x_offset)/ln( A / (Y' - y_offset) - 1) = slope
      m_Slope = -( m_SlopeX - m_XOffset) / log( y_relative_dist);

      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer Sigmoid< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Amplitude / ( 1 + exp(-(x-x_offset)/slope)) + y_offset"
      );

      parameters.AddInitializer
      (
        "descriptor",
        "the property to consider",
        io::Serialization::GetAgent( &m_Descriptor)
      );

      parameters.AddInitializer
      (
        "amplitude",
        "amplitude of the sigmoid",
        io::Serialization::GetAgent( &m_Amplitude),
        "1.0"
      );

      parameters.AddInitializer
      (
        "x offset",
        "offset of the sigmoid midpoint along the x-axis",
        io::Serialization::GetAgent( &m_XOffset),
        "0"
      );

      parameters.AddInitializer
      (
        "y offset",
        "Offset of the sigmoid midpoint along the y-axis",
        io::Serialization::GetAgent( &m_YOffset),
        "0"
      );

      parameters.AddInitializer
      (
        "x",
        "Any x along the sigmoid other than x-offset, for which the y-value can be provided to compute the sigmoidal slope",
        io::Serialization::GetAgent( &m_SlopeX),
        "1"
      );

      parameters.AddInitializer
      (
        "y",
        "y on the desired sigmoidal curve",
        io::Serialization::GetAgent( &m_YAtSlopeX),
        "0.731"
      );

      return parameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > Sigmoid< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    template< typename t_DataType>
    Type::Symmetry Sigmoid< t_DataType>::GetSymmetry() const
    {
      return m_Descriptor->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    template< typename t_DataType>
    bool Sigmoid< t_DataType>::ConsiderRepeatedElements() const
    {
      return m_Descriptor->GetType().ConsiderRepeatedObjects();
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void Sigmoid< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // instantiate the property as a vector with indices that correspond to atoms
      const linal::VectorConstReference< float> &result( m_Descriptor->operator ()( ITR));
      size_t n_elements( result.GetSize());
      for( size_t i( 0); i < n_elements; ++i)
      {
        STORAGE( i) = m_YOffset + m_Amplitude / ( 1.0 + exp( -( result( i) - m_XOffset) / m_Slope));
      }
    }

  } // namespace descriptor
} // namespace bcl
