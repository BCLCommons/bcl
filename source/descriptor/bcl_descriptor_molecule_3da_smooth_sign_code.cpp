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
#include "descriptor/bcl_descriptor_molecule_3da_smooth_sign_code.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_molecule_3da_smooth.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule3DASmoothSignCode::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule3DASmoothSignCode()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Molecule3DASmoothSignCode::Molecule3DASmoothSignCode() :
      m_NumberSteps( 12),
      m_StepSize( 1.0),
      m_Temperature( 100.0),
      m_Smooth( true),
      m_Interpolate( true)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief constructor from number of steps, and mapped atom property
    Molecule3DASmoothSignCode::Molecule3DASmoothSignCode
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const float &TEMPERATURE,
      const bool &SMOOTH,
      const bool &SQRT
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_Temperature( TEMPERATURE),
      m_Smooth( SMOOTH),
      m_Interpolate( true),
      m_Sqrt( SQRT)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule3DASmoothSignCode
    Molecule3DASmoothSignCode *Molecule3DASmoothSignCode::Clone() const
    {
      return new Molecule3DASmoothSignCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule3DASmoothSignCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule3DASmoothSignCode::GetAlias() const
    {
      static const std::string s_name( "3daSmoothSign");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule3DASmoothSignCode::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float Molecule3DASmoothSignCode::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get temperature of code
    //! @return const float  temperature of 3DA code
    const float &Molecule3DASmoothSignCode::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule3DASmoothSignCode::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule3DASmoothSignCode::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomProperty,
          &m_AtomProperty + 1
        );
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT_A, ELEMENT_B the elements of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DASmoothSignCode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_A,
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_B,
      linal::VectorReference< float> &STORAGE
    )
    {
      m_DiscreteCode = 0.0;

      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT_A), iterator_b( ELEMENT_B);
      float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

      // calculate exponential function for the RDF equation
      if( ELEMENT_A == ELEMENT_B)
      {
        m_DiscreteCode( int( property_atom_a > 0.0)) += property_atom_a * property_atom_a;
      }
      else
      {
        float property_atom_b( m_AtomProperty->operator ()( iterator_b)( 0));
        const float distance( linal::Distance( ELEMENT_A->GetPosition(), ELEMENT_B->GetPosition()));
        Accumulate( distance, property_atom_a, property_atom_b);
      }

      // normalization for backwards compatibility
      m_DiscreteCode *= 0.5;
      Smooth( STORAGE);
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DASmoothSignCode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      m_DiscreteCode = 0.0;

      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT);
      float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

      // iterate over distances and the remaining property values on the upper triangle of this matrix
      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_b( ELEMENT.Begin());
        iterator_b.NotAtEnd();
        ++iterator_b
      )
      {
        const chemistry::AtomConformationalInterface &atom_b( *iterator_b( 0));
        //store distance between both atoms
        const float distance( linal::Distance( ELEMENT->GetPosition(), atom_b.GetPosition()));

        const float property_atom_b( m_AtomProperty->operator ()( iterator_b)( 0));

        // calculate exponential function for the RDF equation
        if( iterator_b( 0) == ELEMENT)
        {
          m_DiscreteCode( int( property_atom_a > 0.0)) += property_atom_b * property_atom_b;
        }
        else
        {
          Accumulate( distance, property_atom_a, property_atom_b);
        }
      }
      // this is for backwards compatibility
      m_DiscreteCode *= 0.5;
      Smooth( STORAGE);
    } // Recalculate

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DASmoothSignCode::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      m_DiscreteCode = 0.0;

      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_a( this->GetCurrentObject()->GetIterator());
        iterator_a.NotAtEnd();
        ++iterator_a
      )
      {
        const chemistry::AtomConformationalInterface &atom_a( *iterator_a( 0));
        const float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

        // skip 0 properties since some properties are very sparse
        if( !property_atom_a)
        {
          continue;
        }

        // differentiate between positive and negative atoms
        m_DiscreteCode( int( property_atom_a > 0)) += 0.5 * math::Sqr( property_atom_a);

        // iterate over distances and the remaining property values on the upper triangle of this matrix
        Iterator< chemistry::AtomConformationalInterface> iterator_b( iterator_a);
        for( ++iterator_b; iterator_b.NotAtEnd(); ++iterator_b)
        {
          const chemistry::AtomConformationalInterface &atom_b( *iterator_b( 0));
          //store distance between both atoms
          const float distance( linal::Distance( atom_a.GetPosition(), atom_b.GetPosition()));

          float property_atom_b( m_AtomProperty->operator ()( iterator_b)( 0));
          Accumulate( distance, property_atom_a, property_atom_b);
        }
      }
      Smooth( STORAGE);
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule3DASmoothSignCode::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes the smooth radial distribution function using a given atom property"
      );

      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the smooth radial distribution function",
        io::Serialization::GetAgent( &m_AtomProperty)
      );
      parameters.AddInitializer
      (
        "step size",
        "size of each step in angstroms",
        io::Serialization::GetAgentWithRange( &m_StepSize, 0.01, 100.0),
        "0.25"
      );
      parameters.AddInitializer
      (
        "temperature",
        "increasing temperature spreads autocorrelation across more distant bins",
        io::Serialization::GetAgentWithRange( &m_Temperature, 0.0, 1000.0),
        "100"
      );
      parameters.AddInitializer
      (
        "steps",
        "# of steps/bins (each of size = step size) used in the radial distribution function",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "48"
      );
      parameters.AddInitializer
      (
        "gaussian",
        "whether to apply gaussian smoothing to the final curve. "
        "If set to false, temperature is ignored, interpolation is linear, and no gaussian smoothing is performed",
        io::Serialization::GetAgent( &m_Smooth),
        "True"
      );
      parameters.AddInitializer
      (
        "interpolate",
        "whether to interpolate values to the two nearest points; if false, all weight will be applied to the nearest bin",
        io::Serialization::GetAgent( &m_Interpolate),
        "True"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule3DASmoothSignCode::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_DiscreteCode = linal::Vector< float>( GetNormalSizeOfFeatures(), 0.0);
      m_SmoothingCoefficients =
        Molecule3DASmooth::GetSmoothingCoefficientVector( m_NumberSteps, m_Temperature, m_StepSize);
      if( m_NumberSteps == 0 || m_StepSize == 0.0)
      {
        ERR_STREAM << "m_NumberSteps equals zero - 3DA code will be empty!";
        return false;
      }
      if( m_AtomProperty.IsDefined() && m_AtomProperty->GetNormalSizeOfFeatures() != 1)
      {
        ERR_STREAM
           << "Expected a property that returned 1 properties per atom, but property returns "
           << m_AtomProperty->GetNormalSizeOfFeatures()
           << " values per atom ( property was "
           << m_AtomProperty->GetString()
           << ")";
        return false;
      }
      return true;
    }

    //! @brief create a gaussian-smoothed signal from m_DiscreteCode and store it in STORAGE
    //! @param STORAGE storage for the gaussian-smoothed signal
    void Molecule3DASmoothSignCode::Smooth( linal::VectorReference< float> &STORAGE) const
    {
      if( !m_Smooth)
      {
        STORAGE.CopyValues( m_DiscreteCode);
        return;
      }

      // generate the actual 3DA smooth code by applying the smoothing kernel over the 3DA
      int smoothing_vector_base_position( 0);
      for
      (
        linal::Vector< float>::iterator itr_3DA_smooth( STORAGE.Begin()),
        itr_3DA_smooth_end( STORAGE.End());
        itr_3DA_smooth != itr_3DA_smooth_end;
        ++itr_3DA_smooth, --smoothing_vector_base_position
      )
      {
        float neg_neg( 0), pos_pos( 0), pos_neg( 0);
        int smoothing_vector_position( smoothing_vector_base_position);
        for
        (
          linal::Vector< float>::const_iterator itr_3DA_rough( m_DiscreteCode.Begin()),
          itr_3DA_rough_end( m_DiscreteCode.End());
          itr_3DA_rough != itr_3DA_rough_end;
          ++itr_3DA_rough, ++smoothing_vector_position
        )
        {
          const float exponential_factor( m_SmoothingCoefficients( std::abs( smoothing_vector_position)));
          neg_neg += *itr_3DA_rough * exponential_factor;
          pos_pos += *++itr_3DA_rough * exponential_factor;
          pos_neg += *++itr_3DA_rough * exponential_factor;
        }
        *itr_3DA_smooth = neg_neg;
        *++itr_3DA_smooth = pos_pos;
        *++itr_3DA_smooth = pos_neg;
      }
    }

    //! @brief add an observed distance/property value to m_DiscreteCode
    //! @param DISTANCE actual distance of the two atoms
    //! @param PROP_A property from atom A
    //! @param PROP_B property from atom B
    void Molecule3DASmoothSignCode::Accumulate
    (
      const float &DISTANCE,
      const float &PROP_A,
      const float &PROP_B
    )
    {
      // get the product of the properties
      float product( PROP_A * PROP_B);

      int sign_id( PROP_A > 0.0);

      // differentiate between combinations of positive and negative properties
      if( product < 0.0)
      {
        product = -product;
        sign_id = 2;
      }
      else if( product == 0.0)
      {
        // product == 0 -> nothing to do
        return;
      }

      // What God intends:
      if( m_Sqrt)
      {
        product = math::Sqrt( product);
      }

      // calculate exponential function for the RDF equation
      if( !m_Interpolate)
      {
        const size_t closest_step( float( DISTANCE + 0.5 * m_StepSize) / m_StepSize);
        if( closest_step < m_NumberSteps)
        {
          m_DiscreteCode( 3 * closest_step + sign_id) += product;
        }
      }
      else
      {
        // assign atom distance to the closest 3DA bin
        size_t closest_step_low( std::min( size_t( m_NumberSteps - 1), size_t( DISTANCE / m_StepSize)));
        size_t closest_step_high( closest_step_low + 1);

        const float dist_lower_step( m_StepSize * closest_step_low - DISTANCE);
        const float dist_higher_step( dist_lower_step + m_StepSize);

        float low_exp_factor
        (
          m_Smooth
          ? exp( -m_Temperature * math::Sqr( m_StepSize * closest_step_low - DISTANCE))
          : std::max( dist_higher_step, float( 0.0)) / m_StepSize
        );
        float high_exp_factor
        (
          m_Smooth
          ? exp( -m_Temperature * math::Sqr( m_StepSize * closest_step_high - DISTANCE))
          : ( dist_higher_step >= float( 0.0) ? -dist_lower_step / m_StepSize : 0.0)
        );
        if( closest_step_high == m_NumberSteps)
        {
          // generate rough 3DA for boundry distances
          m_DiscreteCode( 3 * closest_step_low + sign_id) += product * low_exp_factor;
        }
        else
        {
          // normalize rough 3DA code
          product /= ( low_exp_factor + high_exp_factor);
          // generate rough 3DA code for all other distances
          m_DiscreteCode( 3 * closest_step_low + sign_id) += product * low_exp_factor;
          m_DiscreteCode( 3 * closest_step_high + sign_id) += product * high_exp_factor;
        }
      }
    }

  } // namespace descriptor
} // namespace bcl
