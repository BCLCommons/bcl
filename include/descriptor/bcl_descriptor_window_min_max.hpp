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
#include "bcl_descriptor_window_min_max.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_min_max.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor from whether to be springy at the end
    template< typename t_DataType>
    WindowMinMax< t_DataType>::WindowMinMax( const bool &REFLECTIVE) :
      m_Window( REFLECTIVE),
      m_InternalDescriptorSize( 1)
    {
    }

    //! @brief Clone function
    //! @return pointer to new WindowMinMax
    template< typename t_DataType>
    WindowMinMax< t_DataType> *WindowMinMax< t_DataType>::Clone() const
    {
      return new WindowMinMax( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &WindowMinMax< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &WindowMinMax< t_DataType>::GetAlias() const
    {
      static const std::string s_window_name( "WindowMinMax"), s_reflective_name( "ReflectingWindowMinMax");
      return m_Window.IsReflective() ? s_reflective_name : s_window_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t WindowMinMax< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return 2 * m_InternalDescriptorSize;
    }

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
    template< typename t_DataType>
    io::Serializer WindowMinMax< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters( m_Window.GetSerializer());
      parameters.SetClassDescription
      (
        "computes minimum (first value(s)) and maximum (second value(s)) of a descriptor over a window for " + GetStaticClassName< t_DataType>() +
        std::string
        (
          m_Window.IsReflective()
          ? " reflecting at boundaries"
          : " respecting the boundaries"
        ) + ".  For multi-value descriptors, the minimum and maximum are given consecutively for every column"
      );

      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void WindowMinMax< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      const linal::VectorConstReference< float> &window( m_Window( Iterator< t_DataType>( ELEMENT)));
      float *itr_results( STORAGE.Begin());
      for( size_t feature( 0); feature < m_InternalDescriptorSize; ++feature)
      {
        math::RunningMinMax< float> minmax;
        for
        (
          const float *itr_win( window.Begin() + feature), *itr_win_end( window.End());
          itr_win < itr_win_end;
          itr_win += m_InternalDescriptorSize
        )
        {
          if( util::IsDefined( *itr_win))
          {
            minmax += *itr_win;
          }
        }
        *itr_results++ = minmax.GetMin();
        *itr_results++ = minmax.GetMax();
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > WindowMinMax< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Window, &m_Window + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool WindowMinMax< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( !m_Window.ReadInitializerSuccessHook( LABEL, ERR_STREAM))
      {
        return false;
      }
      // get the window's relative indices
      m_InternalDescriptorSize = m_Window.GetSizeOfFeatures() / m_Window.GetRelativeIndices().GetSize();

      return true;
    }

  } // namespace descriptor
} // namespace bcl
