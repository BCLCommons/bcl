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

#ifndef BCL_DESCRIPTOR_REPLACE_UNDEFINED_VALUES_H_
#define BCL_DESCRIPTOR_REPLACE_UNDEFINED_VALUES_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_average.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ReplaceUndefinedValues
    //! @brief replaces the undefined values returned by another atom property with a given value
    //!
    //! @see @link example_descriptor_replace_undefined_values.cpp @endlink
    //! @author mendenjl
    //! @date Nov 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class ReplaceUndefinedValues :
      public Base< t_DataType, float>
    {
    private:

    //////////
    // data //
    //////////

      util::Implementation< Base< t_DataType, float> > m_Descriptor; //!< property to calculate internally

      //! Property to use to calculate the replacement value : must return a single value per sequence
      util::Implementation< Base< t_DataType, float> > m_Replacement;

      //! Flag; if set, the replacement value is automatically the descriptor mean of all non-NaN values
      //! or 0 if all values are NaN
      bool m_ReplaceWithMeanNonNaNValue;

    public:

    //////////
    // data //
    //////////

      //! instances of the class
      static const util::SiPtr< const util::ObjectInterface> s_DynamicReplacementInstance;
      static const util::SiPtr< const util::ObjectInterface> s_ReplaceWithMeanInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param REPLACE_MEAN_WITH_NON_NAN_VALUE true to ignore the replacement property and just use the mean non-nan
      //!        value
      ReplaceUndefinedValues( const bool &REPLACE_MEAN_WITH_NON_NAN_VALUE = false) :
        m_Descriptor(),
        m_Replacement(),
        m_ReplaceWithMeanNonNaNValue( REPLACE_MEAN_WITH_NON_NAN_VALUE)
      {
      }

      //! @brief constructor from SOURCE descriptor (that could return NaNs) and replacement descriptor
      ReplaceUndefinedValues
      (
        const util::Implementation< Base< t_DataType, float> > &SOURCE,
        const util::Implementation< Base< t_DataType, float> > &REPLACEMENT
      ) :
        m_Descriptor( SOURCE),
        m_Replacement( REPLACEMENT),
        m_ReplaceWithMeanNonNaNValue( !REPLACEMENT.IsDefined())
      {
      }

      //! @brief copy constructor
      ReplaceUndefinedValues *Clone() const
      {
        return new ReplaceUndefinedValues( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const
      {
        static std::string s_define_alias( "DefineNaN");
        static std::string s_mean_alias( "SetNaNToDefinedDescriptorMean");
        return m_ReplaceWithMeanNonNaNValue ? s_mean_alias : s_define_alias;
      }

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const
      {
        return m_Descriptor->GetSymmetry();
      }

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const
      {
        return m_Descriptor->ConsiderRepeatedElements();
      }

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to forward calls from SetObject on to all internal implementations
      //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
      //! that function is called
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors()
      {
        // if auto-replacing with the mean of the non-nan values, then do not return the replacement property
        return
          iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Replacement + !m_ReplaceWithMeanNonNaNValue);
      }

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_Descriptor->GetSizeOfFeatures();
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
        io::Serializer parameters;

        parameters.AddInitializer
        (
          "",
          "descriptor whose undefined/NaN values should be replaced",
          io::Serialization::GetAgent( &m_Descriptor)
        );

        if( !m_ReplaceWithMeanNonNaNValue)
        {
          // dynamically replace values with another property
          parameters.SetClassDescription
          (
            "replaces undefined/NaN values in a descriptor with another value or descriptor"
          );

          parameters.AddInitializer
          (
            "replacement",
            "undefined values will be replaced with this property value (usually constants)",
            io::Serialization::GetAgent( &m_Replacement),
            "Constant(0)"
          );
        }
        else
        {
          // always replace with the mean non-nan value
          parameters.SetClassDescription
          (
            "replaces undefined/NaN values in a descriptor with the DescriptorMean of the defined values"
          );
        }

        return parameters;
      }

    private:

      //! @brief True if a derived class expects all calls to SetDimension to also call SetDimension on internally-held classes
      //! @note users do not normally need to override this function
      bool InjectDimensions() const
      {
        return false;
      }

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      bool DimensionIsWellDefined() const
      {
        return false;
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
        linal::VectorConstReference< float> ref( m_Descriptor->operator()( ITR));
        if( ref.GetSize() != this->GetSizeOfFeatures())
        {
          STORAGE = util::GetUndefined< float>();
        }
        else
        {
          STORAGE.CopyValues( ref);
        }
        if( !STORAGE.IsDefined())
        {
          if( !m_ReplaceWithMeanNonNaNValue)
          {
            // replacing with another property
            linal::VectorConstReference< float> reference( ( *m_Replacement)( ITR));

            // check dimensions
            BCL_Assert( ( STORAGE.GetSize() % reference.GetSize()) == 0, "Wrong sizes");

            // replacing undefined values only makes sense if the replacement property returns a single property
            for
            (
              linal::VectorReference< float>::iterator itr_store( STORAGE.Begin()), itr_store_end( STORAGE.End());
              itr_store != itr_store_end;
            )
            {
              for
              (
                linal::VectorConstReference< float>::const_iterator itr( reference.Begin()), itr_end( reference.End());
                itr != itr_end;
                ++itr, ++itr_store
              )
              {
                if( !util::IsDefined( *itr_store))
                {
                  *itr_store = *itr;
                }
              }
            }
          }
          else
          {
            // replacing with the mean defined value
            math::RunningAverage< float> ave_non_nan( float( 0.0));

            // get the average defined value
            for
            (
              linal::VectorReference< float>::iterator itr_store( STORAGE.Begin()), itr_store_end( STORAGE.End());
              itr_store != itr_store_end;
              ++itr_store
            )
            {
              if( util::IsDefined( *itr_store))
              {
                ave_non_nan += *itr_store;
              }
            }

            // set undefineds to the average defined value
            for
            (
              linal::VectorReference< float>::iterator itr_store( STORAGE.Begin()), itr_store_end( STORAGE.End());
              itr_store != itr_store_end;
              ++itr_store
            )
            {
              if( !util::IsDefined( *itr_store))
              {
                *itr_store = ave_non_nan.GetAverage();
              }
            }
          }
        }
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
      {
        if( !m_ReplaceWithMeanNonNaNValue && m_Descriptor->GetSizeOfFeatures() % m_Replacement->GetSizeOfFeatures())
        {
          ERR_STREAM
            << "Primary descriptor must return a number of values that is an integral multiple of the "
            << "# returned by the replacement!";
          return false;
        }
        return true;
      }

    }; // class ReplaceUndefinedValues

    BCL_EXPIMP_TEMPLATE template class BCL_API ReplaceUndefinedValues< char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ReplaceUndefinedValues< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ReplaceUndefinedValues< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ReplaceUndefinedValues< chemistry::AtomConformationalInterface>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_REPLACE_UNDEFINED_VALUES_H_
