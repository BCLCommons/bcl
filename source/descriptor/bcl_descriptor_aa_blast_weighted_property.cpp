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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_aa_blast_weighted_property.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"
#include "sched/bcl_sched_mutex.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////
  // data //
  //////////

    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t method( 0); method < AABlastWeightedProperty::s_NumberMethods; ++method)
        {
          last_instance =
            util::Enumerated< Base< biol::AABase, float> >::AddInstance
            (
              new AABlastWeightedProperty( static_cast< AABlastWeightedProperty::Method>( method))
            );
        }
        return last_instance;
      }
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AABlastWeightedProperty::s_Instance( AddInstances());

    //! @brief get the string for the method
    //! @param METHOD the method to retrieve the name for
    const std::string &AABlastWeightedProperty::GetMethodName( const Method &METHOD)
    {
      static const std::string s_names[ s_NumberMethods + 1] =
      {
        "AA_BlastProbabilityWeighted",
        "AA_BlastLogProbabilityWeighted",
        "AA_BlastLogPSign",
        "AA_BlastLogPSignWeighted",
        "AA_BlastLogPTTest",
        "AA_BlastLogPPearson",
        "AA_BlastLogPSpearman",
        GetStaticClassName< Method>()
      };
      return s_names[ METHOD];
    }

    namespace
    {
      //! @brief create a matrix with one row for each amino acid property, one value for each natural aa
      template< typename t_DataType>
      linal::Matrix< t_DataType> CreateAAPropertyMatrix()
      {
        linal::Matrix< t_DataType> values
        (
          size_t( biol::AATypeData::s_NumberPropertyTypes),
          size_t( biol::AATypes::s_NumberStandardAATypes),
          float( 0.0)
        );
        typename linal::Matrix< t_DataType>::iterator itr_values( values.Begin());
        for
        (
          biol::AATypeData::PropertyTypeEnum local_enum;
          local_enum != biol::AATypeData::s_NumberPropertyTypes;
          ++local_enum
        )
        {
          for
          (
            biol::AATypes::const_iterator
              itr( biol::GetAATypes().Begin()),
              itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
            itr != itr_end;
            ++itr, ++itr_values
          )
          {
            // store the value in the matrix
            *itr_values = ( *itr)->GetAAProperty( local_enum);
          }
        }
        return values;
      }

      //! @brief get a constant vector reference to the values for a given aa property
      //! @param AA_PROPERTY the aa-property of interest
      //! @return a constant vector reference to the values of the property for each natural aa type
      template< typename t_DataType>
      linal::VectorConstReference< t_DataType> GetAAPropertyVector( const biol::AATypeData::PropertyTypeEnum &PROPERTY)
      {
        static const linal::Matrix< t_DataType> s_aa_property( CreateAAPropertyMatrix< t_DataType>());
        return s_aa_property.GetRow( PROPERTY);
      }

      //! @brief calculate the weighted average for every amino acid property
      storage::Vector< float> CalculateWeightedAverages()
      {
        storage::Vector< float> averages;
        averages.AllocateMemory( biol::AATypeData::s_NumberPropertyTypes);
        for
        (
          biol::AATypeData::PropertyTypeEnum local_enum;
          local_enum != biol::AATypeData::s_NumberPropertyTypes;
          ++local_enum
        )
        {
          math::RunningAverage< float> ave( float( 0.0));
          for
          (
            biol::AATypes::const_iterator
              itr( biol::GetAATypes().Begin()),
              itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
            itr != itr_end;
            ++itr
          )
          {
            // compute the average, weighted by natural prevalence, for this property
            ave.AddWeightedObservation
            (
              ( *itr)->GetAAProperty( local_enum),
              ( *itr)->GetAAProperty( biol::AATypeData::e_NaturalPrevalence)
            );
          }
          averages.PushBack( ave.GetAverage());
        }
        return averages;
      }

      //! @brief calculate the weighted average for every amino acid property
      std::pair< storage::Vector< float>, storage::Vector< float> > CalculateAveragesSDs()
      {
        std::pair< storage::Vector< float>, storage::Vector< float> > averages_sds;
        averages_sds.first.AllocateMemory( biol::AATypeData::s_NumberPropertyTypes);
        averages_sds.second.AllocateMemory( biol::AATypeData::s_NumberPropertyTypes);
        for
        (
          biol::AATypeData::PropertyTypeEnum local_enum;
          local_enum != biol::AATypeData::s_NumberPropertyTypes;
          ++local_enum
        )
        {
          math::RunningAverageSD< float> ave;
          linal::VectorConstReference< float> property( GetAAPropertyVector< float>( local_enum));
          for( const float *itr( property.Begin()), *itr_end( property.End()); itr != itr_end; ++itr)
          {
            // compute the average, weighted by natural prevalence, for this property
            ave += *itr;
          }
          averages_sds.first.PushBack( ave.GetAverage());
          averages_sds.second.PushBack( ave.GetStandardDeviation());
        }
        return averages_sds;
      }

      //! @brief compute the weighted average of a particular amino acid property
      //! @param PROPERTY the amino acid property of interest
      //! @return the property average weighted by natural prevalence
      float GetWeightedAverage( const biol::AATypeData::PropertyTypeEnum &PROPERTY)
      {
        static const storage::Vector< float> s_averages( CalculateWeightedAverages());
        return s_averages( PROPERTY);
      }

      //! @brief compute the average of a particular amino acid property
      //! @param PROPERTY the amino acid property of interest
      //! @return the property average
      float GetAverage( const biol::AATypeData::PropertyTypeEnum &PROPERTY)
      {
        static const storage::Vector< float> s_averages( CalculateAveragesSDs().first);
        return s_averages( PROPERTY);
      }

      //! @brief compute the standard deviation of a particular amino acid property
      //! @param PROPERTY the amino acid property of interest
      //! @return the property standard deviation
      float GetStandardDeviation( const biol::AATypeData::PropertyTypeEnum &PROPERTY)
      {
        static const storage::Vector< float> s_averages( CalculateAveragesSDs().second);
        return s_averages( PROPERTY);
      }

      //! @brief simple function to round a float to an integer, respecting sign, e.g. 0.5 -> 1, -0.5 -> -1
      //! @param F the number to round
      int RoundInt( const float &F)
      {
        int tmp( static_cast< int>( F));
        return tmp + ( F - tmp >= 0.5) - ( F - tmp <= -0.5);
      }
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param METHOD weighting method to use
    AABlastWeightedProperty::AABlastWeightedProperty( const Method &METHOD) :
      m_Property(),
      m_Method( METHOD)
    {
    }

    //! @brief constructor from a list of properties
    //! @param PROPERTY property to be used
    //! @param METHOD weighting method to use
    AABlastWeightedProperty::AABlastWeightedProperty
    (
      const biol::AATypeData::PropertyType &PROPERTY,
      const Method &METHOD
    ) :
      m_Property( PROPERTY),
      m_Method( METHOD)
    {
      ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief Clone function
    //! @return pointer to new AABlastProfile
    AABlastWeightedProperty *AABlastWeightedProperty::Clone() const
    {
      return new AABlastWeightedProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AABlastWeightedProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AABlastWeightedProperty::GetAlias() const
    {
      return GetMethodName( m_Method);
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AABlastWeightedProperty::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AABlastWeightedProperty::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      Iterator< biol::AABase> itr_descriptor( ELEMENT);
      linal::VectorConstReference< float> blast_desc( m_TypeCalculator->operator()( itr_descriptor));

      // iterate over the blast profile and the reference values
      const float *itr_blast( blast_desc.Begin());

      // get a reference to the aa-properties
      linal::VectorConstReference< float> aa_property( GetAAPropertyVector< float>( m_Property));
      if( m_Method == e_LogPOffset)
      {
        math::RunningAverage< float> ave( float( 0.0));
        // compute an offset to apply to the blast profile so as to make the weighted sum weight == 1
        const double offset
        (
          ( double( 1.0) - blast_desc.Sum()) / double( biol::AATypes::s_NumberStandardAATypes)
        );
        for
        (
          const float *itr_aa_values( aa_property.Begin()), *itr_aa_values_end( aa_property.End());
          itr_aa_values != itr_aa_values_end;
          ++itr_aa_values, ++itr_blast
        )
        {
          // add the weighted observation
          ave.AddWeightedObservation( *itr_aa_values, *itr_blast + offset);
        }
        STORAGE( 0) = ave.GetAverage() / double( biol::AATypes::s_NumberStandardAATypes);
      }
      else if( m_Method == e_LogPSign)
      {
        math::RunningAverage< float> ave( float( 0.0));
        for
        (
          const float *itr_aa_values( aa_property.Begin()), *itr_aa_values_end( aa_property.End());
          itr_aa_values != itr_aa_values_end;
          ++itr_aa_values, ++itr_blast
        )
        {
          if( *itr_blast >= double( 0.0))
          {
            // add the observation
            ave += *itr_aa_values;
          }
        }
        STORAGE( 0) = ave.GetAverage();
      }
      else if( m_Method == e_LogPSignDiff || m_Method == e_LogPSignTValueSign)
      {
        // perform the weighted average seperately for blast values that are above/below 0.  Ignore values at 0
        math::RunningAverageSD< float> ave_above_zero;
        math::RunningAverageSD< float> ave_below_zero;

        // count the # of values >= 0
        size_t n_ge_0( 0);
        for
        (
          const float *itr_aa_values( aa_property.Begin()), *itr_aa_values_end( aa_property.End());
          itr_aa_values != itr_aa_values_end;
          ++itr_aa_values, ++itr_blast
        )
        {
          if( *itr_blast >= double( 0.0))
          {
            // add the weighted observation
            ave_above_zero.AddWeightedObservation( *itr_aa_values, *itr_blast + 1.0);
            ++n_ge_0;
          }
          else
          {
            // add the weighted observation
            ave_below_zero.AddWeightedObservation( *itr_aa_values, -*itr_blast);
          }
        }
        if( m_Method == e_LogPSignDiff)
        {
          // handle the case that all blast values were at or below zero.  This should typically happen only with
          // unnatural AA types
          if( n_ge_0 == 0 || n_ge_0 == biol::AATypes::s_NumberStandardAATypes)
          {
            STORAGE( 0) = GetAverage( m_Property);
          }
          else
          {
            // values in the blast profile exist both above and below zero, take ave above + (ave_above - ave below)
            // This will elucidate evolutionary selection of values due to a particular property
            STORAGE( 0) = ave_above_zero.GetAverage() + 4.0 * ELEMENT->GetBlastProfile().GetAlignmentWeight() * ( ave_above_zero.GetAverage() - ave_below_zero.GetAverage());
          }
        }
        else
        {
          // handle the case that all blast values were at or below zero.  This should typically happen only with
          // unnatural AA types
          if( n_ge_0 == 0 || n_ge_0 == biol::AATypes::s_NumberStandardAATypes)
          {
            STORAGE( 0) = 0.0;
          }
          else
          {
            // compute t-test values on the values above and below zero
            float t_test_above
            (
              ( ave_above_zero.GetAverage() - GetAverage( m_Property))
              * math::Sqrt( double( n_ge_0))
              / GetStandardDeviation( m_Property)
            );
            if( t_test_above > 5.0)
            {
              t_test_above = 5.0;
            }
            else if( t_test_above < -5.0)
            {
              t_test_above = -5.0;
            }
            float t_test_below
            (
              ( ave_below_zero.GetAverage() - GetAverage( m_Property))
              * math::Sqrt( double( biol::AATypes::s_NumberStandardAATypes - n_ge_0))
              / GetStandardDeviation( m_Property)
            );
            if( t_test_below > 5.0)
            {
              t_test_below = 5.0;
            }
            else if( t_test_below < -5.0)
            {
              t_test_below = -5.0;
            }

            if( ( t_test_above > 0.0) != ( t_test_below > 0.0))
            {
              STORAGE( 0) = t_test_above - t_test_below;
            }
            else
            {
              STORAGE( 0) = t_test_above + t_test_below;
            }
          }
        }
      }
      else if( m_Method == e_Probability)
      {
        math::RunningAverage< float> ave( float( 0.0));
        // no log scale
        for
        (
          const float *itr_aa_values( aa_property.Begin()), *itr_aa_values_end( aa_property.End());
          itr_aa_values != itr_aa_values_end;
          ++itr_aa_values, ++itr_blast
        )
        {
          // add the weighted observation
          ave.AddWeightedObservation( *itr_aa_values, *itr_blast);
        }
        STORAGE( 0) = ave.GetAverage();
      }
      else if( m_Method == e_PearsonCorrelation)
      {
        math::RunningAverageSD< float> blast_stats;
        math::RunningAverage< float> product_stats;

        for
        (
          const float *itr_aa_values( aa_property.Begin()), *itr_aa_values_end( aa_property.End());
          itr_aa_values != itr_aa_values_end;
          ++itr_aa_values, ++itr_blast
        )
        {
          const int blastcol( RoundInt( *itr_blast));
          const float weight
          (
            blastcol < -15 || blastcol > 15
            ? 20.0
            : m_BlastProfileWeighting( blastcol + s_BlastProfileWeightingOffset)
          );
          blast_stats.AddWeightedObservation( *itr_blast, weight);
          product_stats.AddWeightedObservation( *itr_blast * *itr_aa_values, weight);
        }
        // pearson correlation = (AveProduct(A,B) - Mean(A) * Mean(B))/(Std(A)*Std(B))
        const float ave_property( GetAverage( m_Property)), std_property( GetStandardDeviation( m_Property));
        const float ave_blast( blast_stats.GetAverage()), std_blast( blast_stats.GetStandardDeviation());
        const float std_product( std_property * std_blast);
        if( math::EqualWithinTolerance( std_product, float( 0.0)))
        {
          STORAGE( 0) = 0.0;
        }
        else
        {
          STORAGE( 0) = ( product_stats.GetAverage() - ave_blast * ave_property) / std_product;
        }
      }
      else if( m_Method == e_SpearmanCorrelation)
      {
        STORAGE( 0) =
          math::Statistics::CorrelationSpearman
          (
            itr_blast,
            blast_desc.End(),
            GetAAPropertyVector< float>( m_Property).Begin(),
            GetAAPropertyVector< float>( m_Property).End()
          );
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to forward calls from SetObject on to all internal implementations
    //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
    //! that function is called
    iterate::Generic< Base< biol::AABase, float> > AABlastWeightedProperty::GetInternalDescriptors()
    {
      if( m_TypeCalculator.IsDefined())
      {
        return iterate::Generic< Base< biol::AABase, float> >( &m_TypeCalculator, &m_TypeCalculator + 1);
      }
      return iterate::Generic< Base< biol::AABase, float> >();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AABlastWeightedProperty::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        m_Method == e_LogPOffset
        ? "Uses the blast profile (log) values to weight an AA property, specifically "
          "Mean( amino acid property value * ( blast profile log probability for given AA type for this AA + offset)) over all AA types,"
          "where offset is = ( 1.0 - blast profile log probability sum ) / # AA types"
        : "Uses the blast profile probabilities to weight an AA property, specifically "
          "Mean( amino acid property value * blast profile probability for given AA type for this AA) over all AA types"
      );

      serializer.AddInitializer( "property", "the AA property of interest", io::Serialization::GetAgent( &m_Property));
      serializer.AddInitializer
      (
        "blast type",
        "Blast or AA type descriptor to consider (must return 20 values per AA). Defaults to AABlastProfile",
        io::Serialization::GetAgent( &m_TypeCalculator),
        ( m_Method == e_Probability ? "AABlastProbability" : "AABlastProfile")
      );
      return serializer;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool AABlastWeightedProperty::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_TypeCalculator.IsDefined() && m_TypeCalculator->GetSizeOfFeatures() != size_t( 20))
      {
        ERR_STREAM << "A descriptor that returns 20 values is required for the AA type descriptor";
        return false;
      }
      if( !m_TypeCalculator.IsDefined())
      {
        m_TypeCalculator = ( m_Method == e_Probability ? "AABlastProbability" : "AABlastProfile");
      }

      // just call the two core methods for initialization; under a mutex.  This will ensure that other operations do
      // not require a mutex
      static sched::Mutex s_mutex;
      s_mutex.Lock();
      BCL_Assert
      (
        util::IsDefined( GetWeightedAverage( m_Property)),
        "Propensity weighted average for " + m_Property.GetString() + " was undefined"
      );
      BCL_Assert
      (
        util::IsDefined( GetAverage( m_Property)),
        "Average for " + m_Property.GetString() + " was undefined"
      );
      BCL_Assert
      (
        util::IsDefined( GetStandardDeviation( m_Property)),
        "Standard deviation for " + m_Property.GetString() + " was undefined"
      );
      s_mutex.Unlock();
      return true;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AABlastWeightedProperty::SetObjectHook()
    {
      if( m_Method != e_PearsonCorrelation)
      {
        return;
      }

      static const size_t s_blastvector_size = 2 * s_BlastProfileWeightingOffset + 1;
      // These values are appropriate for uniref50, sans fragments, with 5 iterations of psiblast, e-value of 0.01.
      static const float s_defaultblastweight[] =
      {
        20.00, 20.0, 20.0, 20.0, 20.0, 20.0,
        17.64, 8.96, 5.22, 3.25, 2.48, 1.87,
         1.53, 1.26, 1.05, 1.00, 1.53, 2.56,
         4.74, 8.41, 16.3, 20.0, 20.0, 20.0,
        20.00, 20.0, 20.0, 20.0, 20.0, 20.0,
        20.00
      };
      static linal::Vector< float> s_default_blast_weighting_vec( s_blastvector_size, s_defaultblastweight);
      m_BlastProfileWeighting = s_default_blast_weighting_vec;

      // initialize the weights to be the correct size
      linal::Vector< size_t> histogram( s_blastvector_size, size_t( 0));
      for
      (
        Iterator< biol::AABase> itr_seq( this->GetCurrentObject()->GetIterator());
        itr_seq.NotAtEnd();
        ++itr_seq
      )
      {
        linal::VectorConstReference< float> blast_profile( m_TypeCalculator->operator()( itr_seq));
        for
        (
          linal::VectorConstReference< float>::const_iterator
            itr_bp( blast_profile.Begin()), itr_bp_end( blast_profile.End());
          itr_bp != itr_bp_end;
          ++itr_bp
        )
        {
          const int bp( RoundInt( *itr_bp));
          if( bp > s_BlastProfileWeightingOffset || bp < -s_BlastProfileWeightingOffset)
          {
            continue;
          }
          ++histogram( bp + s_BlastProfileWeightingOffset);
        }
      }
      // test whether at least 400 AAs were seen; if not, just use the default vector for better accuracy
      if( histogram.Sum() < 400)
      {
        m_BlastProfileWeighting = s_default_blast_weighting_vec;
        return;
      }

      // add pseudocount of 1
      histogram += size_t( 1);

      // find the max and max index
      const size_t histogram_max_index
      (
        std::distance( histogram.Begin(), std::max_element( histogram.Begin(), histogram.End()))
      );

      // smooth out the histogram to make it bitonic
      const float histogram_max_value( histogram( histogram_max_index));
      for( size_t prev( 0), i( 1); i < histogram_max_index; ++prev, ++i)
      {
        if( histogram( prev) > histogram( i))
        {
          histogram( i) = histogram( prev);
        }
      }
      for( size_t prev( s_blastvector_size - 1), i( s_blastvector_size - 2); i > histogram_max_index; --prev, --i)
      {
        if( histogram( prev) > histogram( i))
        {
          histogram( i) = histogram( prev);
        }
      }

      // now normalize such that the smallest weight is 1 and none are larger than 20.
      // weight(i) = min(20,max_hist_value/histogram(i))
      float sum_positive( 0), sum_negative( 0);
      for( size_t i( 0); i < s_blastvector_size; ++i)
      {
        m_BlastProfileWeighting( i) = std::min( float( 20.0), histogram_max_value / float( histogram( i)));
        ( i < 16 ? sum_negative : sum_positive) += m_BlastProfileWeighting( i) * histogram( i);
      }
      static const double desired_weight_ratio( 1);
      const double weight_ratio( desired_weight_ratio * sum_positive / ( sum_negative + 1));
      for( size_t i( 0); i < 16; ++i)
      {
        m_BlastProfileWeighting( i) *= weight_ratio;
      }
    }

  } // namespace descriptor
} // namespace bcl
