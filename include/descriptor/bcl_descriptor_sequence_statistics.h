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

#ifndef BCL_DESCRIPTOR_SEQUENCE_STATISTICS_H_
#define BCL_DESCRIPTOR_SEQUENCE_STATISTICS_H_

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
    //! @class SequenceStatistics
    //! @brief Computes sequence-wide statistics (column-wise)
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_sequence_statistics.cpp @endlink
    //! @author mendenjl
    //! @date Feb 06, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template
    <
      typename t_DataType,
      typename t_Accumulator,
      const linal::Vector< float> &( t_Accumulator::*GetAccumulatedValues)() const
    >
    class SequenceStatistics :
      public BaseSequence< t_DataType, float>
    {

    //////////
    // data //
    //////////

      t_Accumulator                                    m_Accumulator; //!< Object that accumulates the statistics
      util::Implementation< Base< t_DataType, float> > m_Descriptor; //!< property to calculate internally

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SequenceStatistics()
      {
      }

      //! @brief constructor from implementation
      SequenceStatistics( const util::Implementation< Base< t_DataType, float> > &DESCRIPTOR) :
        m_Descriptor( DESCRIPTOR)
      {
      }

      //! @brief Clone function
      //! @return pointer to new SequenceStatistics
      SequenceStatistics< t_DataType, t_Accumulator, GetAccumulatedValues> *Clone() const
      {
        return new SequenceStatistics( *this);
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
        static std::string s_alias( Base< t_DataType, float>::GetObjectName() + m_Accumulator.GetAliasFromFunctionPtr( GetAccumulatedValues));
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
        if( m_Descriptor.IsDefined() && m_Descriptor->GetType().GetDimension() == size_t( 0))
        {
          ERR_STREAM << GetAlias() << " only makes sense if it is used on a descriptor that considers "
                     << this->GetElementName() << "s in a " << this->GetObjectName();
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
          "Returns the " + s_lower_func + " of the given descriptor across a " + this->GetObjectName()
        );
        parameters.AddInitializer
        (
          "",
          "The descriptor to compute the " + s_lower_func + " of across the " + this->GetObjectName(),
          io::Serialization::GetAgent( &m_Descriptor)
        );

        return parameters;
      }

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE)
      {
        m_Descriptor->SetObject( *Base< t_DataType, float>::GetCurrentObject());
        m_Accumulator.Reset();
        for
        (
          Iterator< t_DataType> itr( m_Descriptor->GetType(), *Base< t_DataType, float>::GetCurrentObject());
          itr.NotAtEnd();
          ++itr
        )
        {
          linal::VectorConstReference< float> result( m_Descriptor->operator()( itr));
          if( result.IsDefined())
          {
            m_Accumulator += result;
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

    }; // class SequenceStatistics

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SEQUENCE_STATISTICS_H_
