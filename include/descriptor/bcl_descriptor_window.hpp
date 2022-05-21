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
#include "bcl_descriptor_window.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType, typename t_ReturnType>
    Window< t_DataType, t_ReturnType>::Window
    (
      const bool &REFLECTIVE,
      const util::Implementation< Base< t_DataType, t_ReturnType> > &DESCRIPTOR
    ) :
      m_Descriptor( DESCRIPTOR),
      m_Size( 0),
      m_SizeLeft( 0),
      m_SizeRight( 0),
      m_Reflective( REFLECTIVE),
      m_Alignment( e_Center)
    {
    }

    //! @brief constructor from a descriptor
    //! @param DESCRIPTOR descriptor to be used
    template< typename t_DataType, typename t_ReturnType>
    Window< t_DataType, t_ReturnType>::Window
    (
      const BaseElement< t_DataType, t_ReturnType> &DESCRIPTOR,
      const size_t &SIZE
    ) :
      m_Descriptor( DESCRIPTOR),
      m_Size( SIZE),
      m_SizeLeft( m_Size),
      m_SizeRight( m_Size),
      m_Reflective( false),
      m_Alignment( e_Center)
    {
    }

    //! @brief constructor from a descriptor
    //! @param DESCRIPTOR descriptor to be used
    template< typename t_DataType, typename t_ReturnType>
    Window< t_DataType, t_ReturnType>::Window
    (
      const util::Implementation< Base< t_DataType, t_ReturnType> > &DESCRIPTOR,
      const size_t &SIZE
    ) :
      m_Descriptor( DESCRIPTOR),
      m_Size( SIZE),
      m_SizeLeft( m_Size),
      m_SizeRight( m_Size),
      m_Reflective( false),
      m_Alignment( e_Center)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Window
    template< typename t_DataType, typename t_ReturnType>
    Window< t_DataType, t_ReturnType> *Window< t_DataType, t_ReturnType>::Clone() const
    {
      return new Window( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Window< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Window< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_window_name( "Window"), s_reflective_name( "ReflectingWindow");
      return m_Reflective ? s_reflective_name : s_window_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType, typename t_ReturnType>
    size_t Window< t_DataType, t_ReturnType>::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor->GetSizeOfFeatures() * ( m_SizeLeft + m_SizeRight + 1);
    }

    //! @brief get the indices of the returned values, relative to the central position
    //! @return the indices of the returned values, relative to the central position
    template< typename t_DataType, typename t_ReturnType>
    storage::Vector< int> Window< t_DataType, t_ReturnType>::GetRelativeIndices() const
    {
      storage::Vector< int> indices;
      indices.AllocateMemory( GetNormalSizeOfFeatures());
      if( m_Alignment != e_JufoCenter)
      {
        for( size_t steps_left( 0); steps_left <= m_SizeLeft; ++steps_left)
        {
          indices.PushBack( -int( steps_left));
        }
        for( size_t steps_right( 1); steps_right <= m_SizeRight; ++steps_right)
        {
          indices.PushBack( int( steps_right));
        }
      }
      else
      {
        for( int steps( 0), total_steps( m_SizeRight + m_SizeLeft); steps <= total_steps; ++steps)
        {
          indices.PushBack( steps - int( m_SizeLeft));
        }
      }
      return indices;
    }

    //! @brief Set the window size
    //! @param RADIUS the radius of the window
    template< typename t_DataType, typename t_ReturnType>
    void Window< t_DataType, t_ReturnType>::SetRadius( const size_t &RADIUS)
    {
      m_Size = RADIUS;

      // call read initializer so that the internal data members get updated
      std::stringstream err_stream;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), err_stream);

      if( util::IsDefined( this->GetDimensionSetting()))
      {
        this->SetDimension( this->GetDimensionSetting());
      }
    }

    //! @brief Set the alignment of the window
    //! @param ALIGNMENT the new alignment of the window
    template< typename t_DataType, typename t_ReturnType>
    void Window< t_DataType, t_ReturnType>::SetAlignment( const WindowAlignmentType &ALIGNMENT)
    {
      m_Alignment = ALIGNMENT;

      // call read initializer so that the internal data members get updated
      std::stringstream err_stream;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), err_stream);

      if( util::IsDefined( this->GetDimensionSetting()))
      {
        this->SetDimension( this->GetDimensionSetting());
      }
    }

    //! @brief get whether the window reflects at the end points
    //! @return true if the window reflects at the end points
    template< typename t_DataType, typename t_ReturnType>
    bool Window< t_DataType, t_ReturnType>::IsReflective() const
    {
      return m_Reflective;
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
    template< typename t_DataType, typename t_ReturnType>
    io::Serializer Window< t_DataType, t_ReturnType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes descriptor windows for " + GetStaticClassName< t_DataType>() +
        std::string
        (
          m_Reflective
          ? " reflecting at boundaries"
          : " returning NaNs for positions beyond the end of the sequence"
        )
      );
      parameters.AddInitializer( "", "descriptor to use for the window", io::Serialization::GetAgent( &m_Descriptor));
      parameters.AddInitializer
      (
        "size",
        "desired size of the window (0 = just central value, 1 = central value and neighboring values, ...)",
        io::Serialization::GetAgent( &m_Size)
      );
      parameters.AddInitializer
      (
        "alignment",
        "Alignment of the window.  Use center to consider the window around this element; left to consider the window"
        " up until this element, and right to consider the window following this element. "
        "JufoCenter is strictly for compatiblity with the old JUFO; it is slower to calculate than normal Center, "
        "but yields the same set of values in a different order",
        io::Serialization::GetAgent( &m_Alignment),
        "Center"
      );

      return parameters;
    } // GetSerializer

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @class SpringIterator
    //! @brief an index iterator with a reflecting boundary condition
    //! Normal iterators have a terminal boundary condition; one must stop upon reaching the end
    //! This iterator class instead bounces off the ends of its container, e.g.
    //! if the values of the container are 1, 2, 3, 4 then iterating from the beginning would give this sequence:
    //! 1 2 3 4 3 2 1 2 3 4 ...
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API SpringIterator
    {
    private:

      size_t m_X;      //!< Current position
      long   m_Dx;     //!< Derivative of m_X; 1 if going up currently, -1 if going down
      size_t m_Length; //!< Overall length of the iteration

    public:

      //! @brief constructor from length, position, and whether or not to start by moving forward
      //! @param LENGTH the length of the container being iterated over
      //! @param POSITION current position in that container
      //! @param FORWARD true if we should start by moving forward, false otherwise
      explicit SpringIterator( const size_t &LENGTH, const size_t &POSITION, const bool &FORWARD) :
        m_X( POSITION),
        m_Dx( 1),
        m_Length( LENGTH)
      {
        if( FORWARD)
        {
          m_Dx = m_X ? 1 : -1;
        }
        else
        {
          m_Dx = m_X == m_Length - 1 ? 1 : -1;
        }
      }

      //! @brief increment the iterator, bounce off the ends if either of them is reached
      void operator++()
      {
        if( ( m_X += m_Dx) >= m_Length)
        {
          m_Dx = -m_Dx;
          m_X = m_Dx > 0 ? 1 : m_Length - 2;
        }
      }

      //! @brief Get the current position of the iterator
      const size_t &GetPosition() const
      {
        return m_X;
      }

    };

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType, typename t_ReturnType>
    void Window< t_DataType, t_ReturnType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< t_ReturnType> &STORAGE
    )
    {
      // current position in the storage vector reference
      size_t storage_pos( 0);

      // size of the sub feature
      size_t sub_feature_size( m_Descriptor->GetSizeOfFeatures());

      // copy the iterator
      Iterator< t_DataType> itr_seq( ELEMENT);

      // do the backward iteration
      // test whether boundary condition will be important for moving backwards
      const size_t steps_backward( std::min( m_SizeLeft, itr_seq.GetPosition()));
      // decrement to obtain previous values
      for
      (
        size_t window_position( 0);
        window_position <= steps_backward;
        ++window_position, --itr_seq, storage_pos += sub_feature_size
      )
      {
        STORAGE.CreateSubVectorReference( sub_feature_size, storage_pos).CopyValues( ( *m_Descriptor)( itr_seq));
      }

      if( steps_backward < m_SizeLeft)
      {
        if( !m_Reflective)
        {
          // just insert nans
          STORAGE.CreateSubVectorReference
          (
            sub_feature_size * ( m_SizeLeft - steps_backward),
            storage_pos
          ) = util::GetUndefined< t_ReturnType>();
          storage_pos += sub_feature_size * ( m_SizeLeft - steps_backward);
        }
        else
        {
          // create spring iterator for going reverse
          SpringIterator itr_spring_backward( itr_seq.GetSize(), 0, true);
          for
          (
            size_t window_position( steps_backward + 1);
            window_position <= m_SizeLeft;
            ++window_position, storage_pos += sub_feature_size
          )
          {
            const size_t old_pos( itr_spring_backward.GetPosition());
            ++itr_spring_backward;
            if( old_pos < itr_spring_backward.GetPosition())
            {
              ++itr_seq;
            }
            else if( old_pos > itr_spring_backward.GetPosition())
            {
              --itr_seq;
            }
            STORAGE.CreateSubVectorReference( sub_feature_size, storage_pos).CopyValues( ( *m_Descriptor)( itr_seq));
          }
        }
      }

      itr_seq = Iterator< t_DataType>( ELEMENT);
      // do the forward iteration
      // test whether boundary condition will be important for moving backwards
      const size_t steps_forward( std::min( m_SizeRight, itr_seq.GetSize() - itr_seq.GetPosition() - size_t( 1)));
      if( steps_forward == m_SizeRight)
      {
        ++itr_seq;
        // increment to obtain remaining values
        for
        (
          size_t window_position( 0);
          window_position < steps_forward;
          ++window_position, ++itr_seq, storage_pos += sub_feature_size
        )
        {
          STORAGE.CreateSubVectorReference( sub_feature_size, storage_pos).CopyValues( ( *m_Descriptor)( itr_seq));
        }
      }
      else if( !m_Reflective)
      {
        // increment to obtain remaining values
        for( ++itr_seq; itr_seq.NotAtEnd(); ++itr_seq, storage_pos += sub_feature_size)
        {
          STORAGE.CreateSubVectorReference( sub_feature_size, storage_pos).CopyValues( ( *m_Descriptor)( itr_seq));
        }
        STORAGE.CreateSubVectorReference( STORAGE.GetSize() - storage_pos, storage_pos) = util::GetUndefined< t_ReturnType>();
      }
      else
      {
        // create spring iterator for going reverse
        SpringIterator itr_spring_forward( itr_seq.GetSize(), itr_seq.GetPosition(), true);
        for
        (
          size_t window_position( 0);
          window_position < m_SizeRight;
          ++window_position, storage_pos += sub_feature_size
        )
        {
          const size_t old_pos( itr_spring_forward.GetPosition());
          ++itr_spring_forward;
          if( old_pos < itr_spring_forward.GetPosition())
          {
            ++itr_seq;
          }
          else if( old_pos > itr_spring_forward.GetPosition())
          {
            --itr_seq;
          }
          STORAGE.CreateSubVectorReference( sub_feature_size, storage_pos).CopyValues( ( *m_Descriptor)( itr_seq));
        }
      }

      // handle jufo-style alignment, in which the window indices are strictly ascending, starting from -window size
      if( m_Alignment == e_JufoCenter)
      {
        // reverse the values as needed to reorder the values as the old jufo would like them
        for( size_t forward( 0), center( ( m_SizeLeft + 1) * sub_feature_size); forward < center;)
        {
          size_t center_offset( center - sub_feature_size);
          if( center_offset == forward)
          {
            break;
          }
          for( size_t i( 0); i < sub_feature_size; ++i, ++forward)
          {
            std::swap( STORAGE( forward), STORAGE( center_offset + i));
          }
          center = center - sub_feature_size;
        }
      }

    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType, typename t_ReturnType>
    iterate::Generic< Base< t_DataType, t_ReturnType> > Window< t_DataType, t_ReturnType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, t_ReturnType> >( &m_Descriptor, &m_Descriptor + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType, typename t_ReturnType>
    bool Window< t_DataType, t_ReturnType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_Alignment == e_Center || m_Alignment == e_JufoCenter)
      {
        m_SizeLeft = m_SizeRight = m_Size;
      }
      else if( m_Alignment == e_Left)
      {
        m_SizeLeft = m_Size;
        m_SizeRight = 0;
      }
      else
      {
        m_SizeLeft = 0;
        m_SizeRight = m_Size;
      }
      if( m_Descriptor.IsDefined() && m_Descriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "Window descriptors only work with element-wise descriptors";
        return false;
      }
      return true;
    }

  } // namespace descriptor
} // namespace bcl
