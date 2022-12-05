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
#include "chemistry/bcl_chemistry_score_function_generic.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  ///////////
  // Enums //
  ///////////

    //! @brief CalculationType as string
    //! @param CALCULATION_TYPE the calculation type whose name is desired
    //! @return the name as string
    const std::string &ScoreFunctionGeneric::GetCalculationTypeName( const CalculationType &CALCULATION_TYPE)
    {
      static const std::string s_Names[ size_t( s_NumberCalculationTypes) + 1] =
      {
        "Index",
        "Sum",
        "Mean",
        "Min",
        "Max",
        "NormMax",
        "SoftMax",
        "Entropy",
        GetStaticClassName< CalculationType>()
      };
      return s_Names[ CALCULATION_TYPE];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ScoreFunctionGeneric::s_Instance
    (
      GetObjectInstances().AddInstance( new ScoreFunctionGeneric())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ScoreFunctionGeneric::ScoreFunctionGeneric() :
        m_Descriptor( descriptor::CheminfoProperty()),
        m_CalculationType( ScoreFunctionGeneric::e_Index),
        m_PropertyIndex( 0),
        m_Invert( false),
        m_Normalize( false),
        m_Noise( 0.00000001)
    {
    }

    //! @brief constructor with parameters
    //! @param DESCRIPTOR the descriptor to use
    ScoreFunctionGeneric::ScoreFunctionGeneric
    (
      const descriptor::CheminfoProperty &DESCRIPTOR
    ) :
      m_Descriptor( DESCRIPTOR),
      m_CalculationType( ScoreFunctionGeneric::e_Index),
      m_PropertyIndex( 0),
      m_Invert( false),
      m_Normalize( false),
      m_Noise( 0.00000001)
    {
    }

    //! @brief constructor with all parameters
    //! @param DESCRIPTOR the descriptor to use
    //! @param INDEX the reference index for certain scores
    //! @param INVERT invert each value in the property array
    //! @param NORMALIZE normalize property array to sum of values
    //! @param NOISE add some small value to bins to avoid ln(0) = nan
    ScoreFunctionGeneric::ScoreFunctionGeneric
    (
      const descriptor::CheminfoProperty &DESCRIPTOR,
      const size_t INDEX,
      const bool INVERT,
      const bool NORMALIZE,
      const double NOISE
    ) :
      m_Descriptor( DESCRIPTOR),
      m_CalculationType( ScoreFunctionGeneric::e_Index),
      m_PropertyIndex( INDEX),
      m_Invert( INVERT),
      m_Normalize( NORMALIZE),
      m_Noise( NOISE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ScoreFunctionGeneric
    ScoreFunctionGeneric *ScoreFunctionGeneric::Clone() const
    {
      return new ScoreFunctionGeneric( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ScoreFunctionGeneric::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the class name when used in a dynamic context
    //! @return the class name when used in a dynamic context
    const std::string &ScoreFunctionGeneric::GetAlias() const
    {
      static const std::string s_alias( "ScoreFunctionGeneric");
      return s_alias;
    }

    //! @brief return the value at a single index
    const double ScoreFunctionGeneric::CalcIndexValue( const linal::Vector< float> &PROPERTIES) const
    {
      return PROPERTIES( m_PropertyIndex);
    }

    //! @brief return the sum of property values
    const double ScoreFunctionGeneric::CalcSum( const linal::Vector< float> &PROPERTIES) const
    {
      return PROPERTIES.Sum();
    }

    //! @brief return the mean descriptor value
    const double ScoreFunctionGeneric::CalcMean( const linal::Vector< float> &PROPERTIES) const
    {
      return PROPERTIES.Sum() / PROPERTIES.GetSize();
    }

    //! @brief return the min descriptor value
    const double ScoreFunctionGeneric::CalcMin( const linal::Vector< float> &PROPERTIES) const
    {
      return PROPERTIES.Min();
    }

    //! @brief return the max descriptor value
    const double ScoreFunctionGeneric::CalcMax( const linal::Vector< float> &PROPERTIES) const
    {
      return PROPERTIES.Max();
    }

    //! @brief return the maximum value after prop_i/sum_0-->N(prop)
    const double ScoreFunctionGeneric::CalcNormMax( const linal::Vector< float> &PROPERTIES) const
    {
      // need a non-const vector
      linal::Vector< float> properties( PROPERTIES);

      const double &sum( CalcSum( properties));
      properties.Normalize();
      properties( m_PropertyIndex) /= sum;
      return properties.Max();
    }

    //! @brief return the maximum value after exp(prop_i)/sum_0-->N(exp(prop))
    const double ScoreFunctionGeneric::CalcSoftMax( const linal::Vector< float> &PROPERTIES) const
    {
      // need a non-const vector
      linal::Vector< float> properties( PROPERTIES);

      // compute exp(value)
      for
      (
          auto itr( properties.Begin()), itr_end( properties.End());
          itr != itr_end;
          ++itr
      )
      {
        *itr = std::exp( *itr);
      }

      // obtain normalized exponentials
      Normalize( properties);

      // return the maximum exponential normalized value
      return properties.Max();
    }

    //! @brief return the entropy of the dataset
    const double ScoreFunctionGeneric::CalcEntropy( const linal::Vector< float> &PROPERTIES) const
    {
      // need a non-const vector
      linal::Vector< float> properties( PROPERTIES);

      // denominator for relative likelihood
      Normalize( properties);

      // compute p*ln(p) for each value, where p is a normalized bin value
      for
      (
          auto itr( properties.Begin()), itr_end( properties.End());
          itr != itr_end;
          ++itr
      )
      {
        // add some noise inside the natural log to avoid nan
        *itr = ( *itr) * ( std::log( *itr + m_Noise));
      }

      return -1.0 * CalcSum( properties);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param MOLECULE molecule that needs to scored
    //! @return score of MOLECULE
    double ScoreFunctionGeneric::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      // setup score function options
      if( m_Descriptor.IsDefined())
      {
        // use passed property
        linal::Vector< double> properties( m_Descriptor->SumOverObject( MOLECULE));

        // invert each element of the property vector
        if( m_Invert)
        {
          Invert( properties);
        }

        // return raw value after potential inversion if the array is of size 1
        if( properties.GetSize() == size_t( 1))
        {
          return properties( 0);
        }

        // normalize property vector by sum
        if( m_Normalize)
        {
          Normalize( properties);
        }

      }
      // no no score defined; return 0.0
      BCL_MessageVrb( "No score defined; returning 0.0!");
      return 0.0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief invert each value of the property vector
    void ScoreFunctionGeneric::Invert( linal::Vector< float> PROPERTIES) const
    {
      for
      (
          auto itr( PROPERTIES.Begin()), itr_end( PROPERTIES.End());
          itr != itr_end;
          ++itr
      )
      {
        *itr = 1.0 / *itr;
      }
    }

    //! @brief normalize property vector by sum of all values
    void ScoreFunctionGeneric::Normalize( linal::Vector< float> PROPERTIES) const
    {
      const double &sum( CalcSum( PROPERTIES));
      for
      (
          auto itr( PROPERTIES.Begin()), itr_end( PROPERTIES.End());
          itr != itr_end;
          ++itr
      )
      {
        *itr = *itr / sum;
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief gets a serializer for constructing this class in a dynamic context
    //! @return the serializer containing member data
    io::Serializer ScoreFunctionGeneric::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "Compute a chemical property and transform the resultant (potentially multi-dimensional) array "
        "using one of several calculation types to return a final score"
      );
      member_data.AddInitializer
      (
        "descriptor",
        "the descriptor to calculate; "
        "if multi-valued, this will be transformed with the specified calculation type.",
        io::Serialization::GetAgent( &m_Descriptor)
      );
      member_data.AddInitializer
      (
        "calculation_type",
        "transform multi-dimensional array output with one of these allowed "
        "operations; final output will be a scalar value",
        io::Serialization::GetAgent( &m_CalculationType),
        "Index"
      );
      member_data.AddInitializer
      (
        "property_index",
        "the index of interest in a multi-dimensional property; "
        "no effect if the calculation type does not require a reference index",
        io::Serialization::GetAgent( &m_PropertyIndex),
        "0"
      );
      member_data.AddInitializer
      (
        "invert",
        "invert property values; occurs prior to any normalization",
        io::Serialization::GetAgent( &m_Invert),
        "false"
      );
      member_data.AddInitializer
      (
        "normalize",
        "normalize property values; occurs after any inversion; "
        "note that this is a raw normalization, therefore if the sum of "
        "all values is 0 then the resulting normalized vector will be undefined;"
        "this behavior is kept intentionally because if all values in a multi- "
        "dimensional array are 0 there is typically an issue.",
        io::Serialization::GetAgent( &m_Normalize),
        "false"
      );
      member_data.AddInitializer
      (
        "noise",
        "noise added to each value prior to taking the natural logarithm of said value; "
        "this allows calculation types such as 'SoftMax' to be estimated even when "
        "some of the property bins are 0; "
        "this is specific to calculation types that make use of logarithms.",
        io::Serialization::GetAgent( &m_Noise),
        "0.00000001"
      );

      return member_data;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ScoreFunctionGeneric::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ScoreFunctionGeneric::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
