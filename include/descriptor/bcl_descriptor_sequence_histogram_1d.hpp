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
#include "bcl_descriptor_sequence_histogram_1d.h"
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
    SequenceHistogram1D< t_DataType>::SequenceHistogram1D() :
      m_Descriptor(),
      m_Min( 0.0),
      m_Max( 1.0),
      m_BinSize( 1.0),
      m_NumberBins( 1),
      m_SmoothingDistance( 0.0),
      m_Catchall( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SequenceHistogram1D
    template< typename t_DataType>
    SequenceHistogram1D< t_DataType> *SequenceHistogram1D< t_DataType>::Clone() const
    {
      return new SequenceHistogram1D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &SequenceHistogram1D< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &SequenceHistogram1D< t_DataType>::GetAlias() const
    {
      static const std::string s_alias( this->GetObjectName() + "Histogram1D");
      return s_alias;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t SequenceHistogram1D< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_NumberBins;
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
    io::Serializer SequenceHistogram1D< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "computes a histogram using a single descriptor");
      parameters.AddInitializer( "", "descriptor to use", io::Serialization::GetAgent( &m_Descriptor));
      parameters.AddInitializer( "min", "min value to consider", io::Serialization::GetAgent( &m_Min));
      parameters.AddInitializer( "max", "max value to consider", io::Serialization::GetAgent( &m_Max));
      parameters.AddInitializer( "bin size", "size of each bin", io::Serialization::GetAgentWithMin( &m_BinSize, 0.0));
      parameters.AddInitializer
      (
        "smoothing",
        "if non-zero, return a histogram with points smoothed with a gaussian kernel, specifically "
        "Ae^(-BinDistance/SmoothingDistance), where BinDistance is the euclidean distance to center of the given bin, "
        "and A is chosen such that the histogram's sum still equals 1",
        io::Serialization::GetAgentWithMin( &m_SmoothingDistance, 0.0),
        "0"
      );
      parameters.AddInitializer
      (
        "catchall",
        "If set, values that fall outside the boundaries will be placed into the nearest bin",
        io::Serialization::GetAgent( &m_Catchall),
        "False"
      );
      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void SequenceHistogram1D< t_DataType>::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      Iterator< t_DataType> itr_seq( this->GetEffectiveType());
      itr_seq.SetObject( *this->GetCurrentObject());
      // get x and y-axis values
      const linal::VectorConstReference< float> x_values( m_Descriptor->operator()( itr_seq));
      STORAGE = 0;

      // insert each element into the histogram
      if( m_Catchall)
      {
        for( size_t col( 0), sub_feature_size( m_Descriptor->GetSizeOfFeatures()); col < sub_feature_size; ++col)
        {
          // compute x position
          const size_t x_position
          (
            std::min
            (
              m_NumberBins - size_t( 1),
              size_t( ( std::max( x_values( col), m_Min) - m_Min) / m_BinSize)
            )
          );

          // handle basic histogram
          STORAGE( x_position) += 1.0;

          // handle smoothing
          // It is useful to limit the smoothing to a maximum effect of just under 0.5 for the
          // nearest bins to allow use of either accuracy or categorical max as an objective function
          if( m_SmoothingDistance)
          {
            const float &x( x_values( col));
            // compute x and y-position
            const float x_pos( ( x - m_Min) / m_BinSize);

            float x_bin_position( 0.5);
            float *itr_storage( STORAGE.Begin());
            for( size_t x_bin( 0); x_bin < m_NumberBins; ++x_bin, x_bin_position += 1.0, ++itr_storage)
            {
              *itr_storage += 0.4999 * exp( -math::Absolute( x_pos - x_bin_position) / m_SmoothingDistance);
            }
          }
        }
      }
      else
      {
        for( size_t col( 0), sub_feature_size( m_Descriptor->GetSizeOfFeatures()); col < sub_feature_size; ++col)
        {
          const float &x( x_values( col));

          // handle out-of-bounds values
          if( x < m_Min || x > m_Max)
          {
            continue;
          }

          // set the desired bin to 1
          STORAGE( std::min( m_NumberBins - size_t( 1), size_t( ( x - m_Min) / m_BinSize))) += 1.0;

          if( m_SmoothingDistance)
          {
            // handle smoothing
            // It is useful to limit the smoothing to a maximum effect of just under 0.5 for the
            // nearest bins to allow use of either accuracy or categorical max as an objective function
            const float x_pos( ( x - m_Min) / m_BinSize);

            float x_bin_position( 0.5);
            float *itr_storage( STORAGE.Begin());
            for( size_t x_bin( 0); x_bin < m_NumberBins; ++x_bin, x_bin_position += 1.0, ++itr_storage)
            {
              *itr_storage += 0.4999 * exp( -math::Absolute( x_pos - x_bin_position) / m_SmoothingDistance);
            }
          }
        }
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > SequenceHistogram1D< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool SequenceHistogram1D< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // check that both x and y axis are element-wise descriptors
      if( m_Descriptor.IsDefined() && m_Descriptor->GetType().GetDimension() != size_t( 0))
      {
        ERR_STREAM << "SequenceHistogram1D descriptors only work with whole-" << this->GetObjectName() << " descriptors";
        return false;
      }

      if( m_Min > m_Max)
      {
        ERR_STREAM << "min must be < max";
        return false;
      }
      if( !m_BinSize)
      {
        ERR_STREAM << "Non-zero bin sizes must be given";
        return false;
      }

      if( m_Min == m_Max) // Trivial histogram case; just catches values that fall into this bin
      {
        m_NumberBins = 1;
      }
      else
      {
        m_NumberBins = ( m_Max - m_Min) / m_BinSize;
        if( m_Max - m_Min - m_NumberBins * m_BinSize > m_BinSize / 2.0)
        {
          // rounding issue, add an extra bin
          ++m_NumberBins;
        }
      }
      return true;
    }

  } // namespace descriptor
} // namespace bcl
