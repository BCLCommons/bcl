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

#ifndef BCL_DESCRIPTOR_STATISTICS_H_
#define BCL_DESCRIPTOR_STATISTICS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
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
    //! @class Statistics
    //! @brief Computes statistics on individual descriptor rows
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_statistics.cpp @endlink
    //! @author mendenjl
    //! @date Feb 06, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template
    <
      typename t_DataType,
      typename t_Accumulator,
      const float &( t_Accumulator::*GetAccumulatedValues)() const
    >
    class Statistics :
      public Base< t_DataType, float>
    {

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, float> > m_Descriptor; //!< property to calculate internally

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Statistics
      Statistics< t_DataType, t_Accumulator, GetAccumulatedValues> *Clone() const
      {
        return new Statistics( *this);
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
        static std::string s_alias( "Descriptor" + t_Accumulator().GetAliasFromFunctionPtr( GetAccumulatedValues));
        return s_alias;
      }

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      //! @note only one value is needed if this is a numeric descriptor, for char descriptors, assume that 99999 is the
      //! @note max, so 5 characters is sufficient
      size_t GetNormalSizeOfFeatures() const
      {
        return size_t( 1);
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

      //! @brief True if a derived class expects all calls to SetDimension to also call SetDimension on internally-held classes
      //! False because the internal descriptors for models cannot be changed in size
      bool InjectDimensions() const
      {
        return false;
      }

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      bool DimensionIsWellDefined() const
      {
        return m_Descriptor.IsDefined() && m_Descriptor->DimensionIsWellDefined();
      }

    protected:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
      {
        if( m_Descriptor.IsDefined() && m_Descriptor->GetSizeOfFeatures() <= size_t( 1))
        {
          ERR_STREAM << GetAlias() << " only makes sense if it is used on a descriptor that returns multiple values, but was only given: "
                     << ( m_Descriptor.IsDefined() ? m_Descriptor.GetLabel().ToString() : " \"\"");
          return false;
        }
        return true;
      }

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< float> &STORAGE
      )
      {
        t_Accumulator accumulator; // Object that accumulates the statistics
        linal::VectorConstReference< float> result( m_Descriptor->operator()( ITR));
        for
        (
          linal::VectorConstReference< float>::const_iterator
            itr_result( result.Begin()), itr_result_end( result.End());
          itr_result != itr_result_end;
          ++itr_result
        )
        {
          accumulator += *itr_result;
        }
        STORAGE( 0) = ( accumulator.*GetAccumulatedValues)();
      }

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook()
      {
        m_Descriptor->SetObject( *Base< t_DataType, float>::GetCurrentObject());
      }

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const
      {
        return m_Descriptor->GetType().GetDimension();
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        static const std::string s_lower_func
        (
          util::ToLower( t_Accumulator().GetAliasFromFunctionPtr( GetAccumulatedValues))
        );
        io::Serializer parameters;
        parameters.SetClassDescription
        (
          "Returns the " + s_lower_func + " of a descriptor"
        );
        parameters.AddInitializer
        (
          "",
          "The descriptor on which to compute the " + s_lower_func,
          io::Serialization::GetAgent( &m_Descriptor)
        );

        return parameters;
      }

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class Statistics

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_STATISTICS_H_
