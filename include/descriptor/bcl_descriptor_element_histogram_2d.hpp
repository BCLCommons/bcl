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
#include "bcl_descriptor_element_histogram_2d.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType>
    ElementHistogram2D< t_DataType>::ElementHistogram2D() :
      m_XAxis(),
      m_YAxis( 0),
      m_MinX( 0.0),
      m_MaxX( 1.0),
      m_MinY( 0.0),
      m_MaxY( 1.0),
      m_BinSizeX( 1.0),
      m_BinSizeY( 1.0),
      m_NumberBinsX( 1),
      m_NumberBinsY( 1),
      m_SmoothingDistance( 0.0),
      m_Catchall( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ElementHistogram2D
    template< typename t_DataType>
    ElementHistogram2D< t_DataType> *ElementHistogram2D< t_DataType>::Clone() const
    {
      return new ElementHistogram2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &ElementHistogram2D< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &ElementHistogram2D< t_DataType>::GetAlias() const
    {
      static const std::string s_alias( "ElementHistogram2D");
      return s_alias;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t ElementHistogram2D< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_NumberBinsX * m_NumberBinsY;
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
    io::Serializer ElementHistogram2D< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes a binary histogram using different descriptors for X and Y axis"
      );
      parameters.AddInitializer( "x", "descriptor to use for the x-axis", io::Serialization::GetAgent( &m_XAxis));
      parameters.AddInitializer( "y", "descriptor to use for the y-axis", io::Serialization::GetAgent( &m_YAxis));
      parameters.AddInitializer( "min x", "min value to plot for the x-axis", io::Serialization::GetAgent( &m_MinX));
      parameters.AddInitializer( "min y", "min value to plot for the y-axis", io::Serialization::GetAgent( &m_MinY));
      parameters.AddInitializer( "max x", "max value to plot for the x-axis", io::Serialization::GetAgent( &m_MaxX));
      parameters.AddInitializer( "max y", "max value to plot for the y-axis", io::Serialization::GetAgent( &m_MaxY));
      parameters.AddInitializer
      (
        "bin size x",
        "bin size for the x-axis",
        io::Serialization::GetAgentWithMin( &m_BinSizeX, 0.0)
      );
      parameters.AddInitializer
      (
        "bin size y",
        "bin size for the y-axis",
        io::Serialization::GetAgentWithMin( &m_BinSizeY, 0.0)
      );
      parameters.AddInitializer
      (
        "smoothing",
        "if non-zero, return a grid with points smoothed with a gaussian kernel, specifically "
        "Ae^(-BinDistance/SmoothingDistance), where BinDistance is the euclidean distance to center of the given bin, "
        "and A is chosen such that the histogram's sum still equals 1",
        io::Serialization::GetAgentWithMin( &m_SmoothingDistance, 0.0),
        "0"
      );
      parameters.AddInitializer
      (
        "catchall",
        "If set, values that fall outside the histogram boundaries will be placed into the nearest bin",
        io::Serialization::GetAgent( &m_Catchall),
        "False"
      );
      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void ElementHistogram2D< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      Iterator< t_DataType> itr_seq( ELEMENT);
      // get x and y-axis values
      const linal::VectorConstReference< float> x_values( m_XAxis->operator()( itr_seq));
      const linal::VectorConstReference< float> y_values( m_YAxis->operator()( itr_seq));

      if( x_values.GetSize() != y_values.GetSize())
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }
      STORAGE = 0;

      // insert each element into the histogram
      if( m_Catchall)
      {
        // no smoothing, just accumulate counts
        for( size_t col( 0), sub_feature_size( m_XAxis->GetSizeOfFeatures()); col < sub_feature_size; ++col)
        {
          const float &x( x_values( col)), &y( y_values( col));
          // compute x and y-position
          const size_t x_position( std::min( m_NumberBinsX - size_t( 1), size_t( ( std::max( x, m_MinX) - m_MinX) / m_BinSizeX)));
          const size_t y_position( std::min( m_NumberBinsY - size_t( 1), size_t( ( std::max( y, m_MinY) - m_MinY) / m_BinSizeY)));

          STORAGE( x_position * m_NumberBinsY + y_position) += 1.0;
          if( m_SmoothingDistance)
          {
            const float x_pos( ( std::min( std::max( x, m_MinX), m_MaxX) - m_MinX) / m_BinSizeX);
            const float y_pos( ( std::min( std::max( y, m_MinY), m_MaxY) - m_MinY) / m_BinSizeY);

            float x_bin_position( 0.5);
            float *itr_storage( STORAGE.Begin());
            for( size_t x_bin( 0); x_bin < m_NumberBinsX; ++x_bin, x_bin_position += 1.0)
            {
              const float x_distance_squared( math::Sqr( x_pos - x_bin_position));
              float y_bin_position( 0.5);
              for( size_t y_bin( 0); y_bin < m_NumberBinsY; ++y_bin, y_bin_position += 1.0, ++itr_storage)
              {
                const float y_distance_squared( math::Sqr( y_pos - y_bin_position));
                *itr_storage += 0.499999 * exp( -math::Sqrt( x_distance_squared + y_distance_squared) / m_SmoothingDistance);
              }
            }
          }
        }
      }
      else
      {
        for( size_t col( 0), sub_feature_size( m_XAxis->GetSizeOfFeatures()); col < sub_feature_size; ++col)
        {
          const float &x( x_values( col)), &y( y_values( col));

          // handle out-of-bounds values
          if( x < m_MinX || x > m_MaxX || y < m_MinY || y > m_MaxY)
          {
            continue;
          }

          // compute x and y-position
          const size_t x_position( std::min( m_NumberBinsX - size_t( 1), size_t( ( x - m_MinX) / m_BinSizeX)));
          const size_t y_position( std::min( m_NumberBinsY - size_t( 1), size_t( ( y - m_MinY) / m_BinSizeY)));

          STORAGE( x_position * m_NumberBinsY + y_position) += 1.0;

          if( m_SmoothingDistance)
          {
            // compute x and y-position
            const float x_position( ( x - m_MinX) / m_BinSizeX);
            const float y_position( ( y - m_MinY) / m_BinSizeY);

            float x_bin_position( 0.5);
            float *itr_storage( STORAGE.Begin());
            for( size_t x_bin( 0); x_bin < m_NumberBinsX; ++x_bin, x_bin_position += 1.0)
            {
              const float x_distance_squared( math::Sqr( x_position - x_bin_position));
              float y_bin_position( 0.5);
              for( size_t y_bin( 0); y_bin < m_NumberBinsY; ++y_bin, y_bin_position += 1.0, ++itr_storage)
              {
                const float y_distance_squared( math::Sqr( y_position - y_bin_position));
                *itr_storage += 0.499999 * exp( -math::Sqrt( x_distance_squared + y_distance_squared) / m_SmoothingDistance);
              }
            }
          }
        }
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > ElementHistogram2D< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_XAxis, &m_YAxis + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool ElementHistogram2D< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // check that both x and y axis are element-wise descriptors
      if
      (
        ( m_XAxis.IsDefined() && m_XAxis->GetType().GetDimension() != size_t( 1))
        || ( m_YAxis.IsDefined() && m_YAxis->GetType().GetDimension() != size_t( 1))
      )
      {
        ERR_STREAM << "ElementHistogram2D descriptors only work with element-wise descriptors";
        return false;
      }
      // check that both descriptors return an equal number of values
      if( m_XAxis.IsDefined() && m_YAxis.IsDefined() && m_XAxis->GetSizeOfFeatures() != m_YAxis->GetSizeOfFeatures())
      {
        ERR_STREAM << "ElementHistogram2D requires that the x-axis and y-axis descriptors return the same # of values";
        return false;
      }

      if( m_MinX > m_MaxX || m_MinY > m_MaxY)
      {
        ERR_STREAM << "min x must be < max x (and similarly for y-axis values)";
        return false;
      }
      if( !m_BinSizeX || !m_BinSizeY)
      {
        ERR_STREAM << "Non-zero bin sizes must be given for both the x and y axis";
        return false;
      }

      if( m_MinX == m_MaxX) // 1D histogram case
      {
        m_NumberBinsX = 1;
      }
      else
      {
        m_NumberBinsX = ( m_MaxX - m_MinX) / m_BinSizeX;
        if( m_MaxX - m_MinX - m_NumberBinsX * m_BinSizeX > m_BinSizeX / 2.0)
        {
          // rounding issue, add an extra bin
          ++m_NumberBinsX;
        }
      }
      if( m_MinY == m_MaxY) // 1D histogram case
      {
        m_NumberBinsY = 1;
      }
      else
      {
        m_NumberBinsY = ( m_MaxY - m_MinY) / m_BinSizeY;
        if( m_MaxY - m_MinY - m_NumberBinsY * m_BinSizeY > m_BinSizeY / 2.0)
        {
          // rounding issue, add an extra bin
          ++m_NumberBinsY;
        }
      }
      return true;
    }

  } // namespace descriptor
} // namespace bcl
