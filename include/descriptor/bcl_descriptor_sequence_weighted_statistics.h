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

#ifndef BCL_DESCRIPTOR_SEQUENCE_WEIGHTED_STATISTICS_H_
#define BCL_DESCRIPTOR_SEQUENCE_WEIGHTED_STATISTICS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SequenceWeightedStatistics
    //! @brief Computes sequence-wide statistics (column-wise)
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_sequence_weighted_statistics.cpp @endlink
    //! @author mendenjl
    //! @date Feb 27, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template
    <
      typename t_DataType,
      typename t_Accumulator,
      const linal::Vector< float> &( t_Accumulator::*GetAccumulatedValues)() const
    >
    class SequenceWeightedStatistics :
      public BaseSequence< t_DataType, float>
    {

    //////////
    // data //
    //////////

      t_Accumulator                                    m_Accumulator; //!< Object that accumulates the statistics
      util::Implementation< Base< t_DataType, float> > m_Descriptor;  //!< property to calculate internally
      util::Implementation< Base< t_DataType, float> > m_Weight;      //!< property to use for weighting
      Type                                             m_Type;        //!< Cached type

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SequenceWeightedStatistics()
      {
      }

      //! @brief constructor from implementation
      SequenceWeightedStatistics
      (
        const util::Implementation< Base< t_DataType, float> > &DESCRIPTOR,
        const util::Implementation< Base< t_DataType, float> > &WEIGHT
      ) :
        m_Descriptor( DESCRIPTOR),
        m_Weight( WEIGHT)
      {
        BCL_Assert
        (
          ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
          "invalid initialization of " + GetAlias()
        );
      }

      //! @brief Clone function
      //! @return pointer to new SequenceWeightedStatistics
      SequenceWeightedStatistics< t_DataType, t_Accumulator, GetAccumulatedValues> *Clone() const
      {
        return new SequenceWeightedStatistics( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const
      {
        static std::string s_alias( Base< t_DataType, float>::GetObjectName() + "Weighted" + m_Accumulator.GetAliasFromFunctionPtr( GetAccumulatedValues));
        return s_alias;
      }

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      //! @note only one value is needed if this is a numeric descriptor, for char descriptors, assume that 99999 is the
      //! @note max, so 5 characters is sufficient
      size_t GetNormalSizeOfFeatures() const
      {
        return m_Descriptor->GetSizeOfFeatures();
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

    protected:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
      {
        if( !m_Descriptor.IsDefined())
        {
          // nothing to be done if the descriptor itself is undefined
          return true;
        }
        if( m_Descriptor->GetType().GetDimension() == size_t( 0))
        {
          ERR_STREAM << GetAlias() << " only makes sense when used on a descriptor that considers "
                     << this->GetElementName() << " of the " << this->GetObjectName();
          return false;
        }
        if( !m_Weight.IsDefined())
        {
          ERR_STREAM << GetAlias() << " requires a weighting descriptor";
          return false;
        }
        if( m_Weight->GetType().GetDimension() == size_t( 0))
        {
          ERR_STREAM << GetAlias() << " only makes sense when used on a descriptor that considers "
                     << this->GetElementName() << " of the " << this->GetObjectName();
          return false;
        }
        m_Type = m_Descriptor->GetType();
        m_Type.GeneralizeToHandle( m_Weight->GetType());
        m_Descriptor->SetDimension( m_Type.GetDimension());
        m_Weight->SetDimension( m_Type.GetDimension());
        if( m_Weight->GetSizeOfFeatures() != size_t( 1))
        {
          ERR_STREAM << GetAlias() << " requires a descriptor that returns only a single value, but "
                     << " " << m_Weight.GetString() << " returns " << m_Weight->GetSizeOfFeatures();
          return false;
        }

        return true;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        static const std::string s_lower_func
        (
          util::ToLower( m_Accumulator.GetAliasFromFunctionPtr( GetAccumulatedValues))
        );
        io::Serializer parameters;
        parameters.SetClassDescription
        (
          s_lower_func + " of a descriptor across the " + this->GetObjectName() + ", weighted by any other descriptor"
        );
        parameters.AddInitializer
        (
          "",
          "The descriptor to compute the " + s_lower_func + " of across the " + this->GetObjectName(),
          io::Serialization::GetAgent( &m_Descriptor)
        );
        parameters.AddInitializer
        (
          "weight",
          "The descriptor used to weight the primary descriptor",
          io::Serialization::GetAgent( &m_Weight)
        );
        return parameters;
      }

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE)
      {
        m_Descriptor->SetObject( *Base< t_DataType, float>::GetCurrentObject());
        m_Weight->SetObject( *Base< t_DataType, float>::GetCurrentObject());
        m_Accumulator.Reset();
        for
        (
          Iterator< t_DataType> itr( m_Type, *Base< t_DataType, float>::GetCurrentObject());
          itr.NotAtEnd();
          ++itr
        )
        {
          linal::VectorConstReference< float> result( m_Descriptor->operator()( itr));
          if( result.IsDefined())
          {
            linal::VectorConstReference< float> weight_vec( m_Weight->operator()( itr));
            if( weight_vec.GetSize() != size_t( 1))
            {
              STORAGE = util::GetUndefined< float>();
            }
            const float weight( weight_vec( 0));
            if( util::IsDefined( weight))
            {
              m_Accumulator.AddWeightedObservation( result, weight);
            }
          }
        }
        if( ( m_Accumulator.*GetAccumulatedValues)().GetSize())
        {
          STORAGE.CopyValues( ( m_Accumulator.*GetAccumulatedValues)());
        }
        else
        {
          STORAGE = util::GetUndefined< float>();
        }
      }

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class SequenceWeightedStatistics

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SEQUENCE_WEIGHTED_STATISTICS_H_
