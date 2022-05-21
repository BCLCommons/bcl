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
#include "descriptor/bcl_descriptor_molecule_3da_soft_max.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_molecule_3da_smooth.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_min_max.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule3DASoftMax::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule3DASoftMax()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Molecule3DASoftMax::Molecule3DASoftMax() :
      m_NumberSteps( 12),
      m_StepSize( 1.0),
      m_Temperature( 100.0),
      m_Smooth( true)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief constructor from number of steps, and mapped atom property
    Molecule3DASoftMax::Molecule3DASoftMax
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const float &TEMPERATURE
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_Temperature( TEMPERATURE),
      m_Smooth( true)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule3DASoftMax
    Molecule3DASoftMax *Molecule3DASoftMax::Clone() const
    {
      return new Molecule3DASoftMax( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule3DASoftMax::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule3DASoftMax::GetAlias() const
    {
      static const std::string s_name( "3daSoftMax");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule3DASoftMax::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float Molecule3DASoftMax::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get temperature of code
    //! @return const float  temperature of 3DA code
    const float &Molecule3DASoftMax::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule3DASoftMax::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule3DASoftMax::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomProperty,
          &m_AtomProperty + 1
        );
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DASoftMax::Calculate( linal::VectorReference< float> &STORAGE)
    {
      m_DiscreteCode = math::GetLowestBoundedValue< float>();

      // instantiate the property as a vector with indices that correspond to atoms
      const linal::Vector< float> property( m_AtomProperty->CollectValuesOnEachElementOfObject( *this->GetCurrentObject()));

      // determine the minimum / maximum values for any bin that for which a value was obtained
      // this is needed in case no values are ever received for a particular bin
      math::RunningMinMax< float> min_max_real;

      // iterate over all possible pairs of atoms
      // iterate properties and surface areas simultaneously
      linal::Vector< float>::const_iterator itr_prop_a( property.Begin());
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_a( this->GetCurrentObject()->GetIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a, ++itr_prop_a
      )
      {
        // compute the property value for the atom onto itself
        const float prop( math::Sqr( *itr_prop_a));

        // generate rough 3DA code for each atom onto itself
        m_DiscreteCode( 0) = std::max( m_DiscreteCode( 0), prop);

        // add the value to the min/max tracker
        min_max_real += prop;

        // iterate over all possible pairs of atoms
        // iterate properties and surface areas simultaneously
        linal::Vector< float>::const_iterator itr_prop_b( itr_prop_a + 1);

        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_b( itr_atoms_a);
        for( ++itr_atoms_b; itr_atoms_b.NotAtEnd(); ++itr_atoms_b, ++itr_prop_b)
        {
          // store distance between both atoms
          const float distance( linal::Distance( itr_atoms_a->GetPosition(), itr_atoms_b->GetPosition()));

          // determine the nominal bin
          size_t nominal_bin( distance / m_StepSize);

          // continue if the atoms are beyond this distance
          if( nominal_bin >= m_NumberSteps)
          {
            continue;
          }

          // scale by the atom property
          float prop( *itr_prop_a * ( *itr_prop_b));

          // generate rough 3DA for boundry distances
          float &three_da_value( m_DiscreteCode( nominal_bin));
          three_da_value = std::max( prop, three_da_value);

          // add the value to the min/max tracker
          min_max_real += prop;
        }
      }

      // for all bins for which no max was ever obtained, set it to the minimum value that was seen, provided that
      // all values seen were not identical
      const float default_value( min_max_real.GetMax() > min_max_real.GetMin() ? min_max_real.GetMin() : float( 0.0));

      // generate the actual 3DA smooth code by applying the smoothing kernel over the 3DA
      int smoothing_vector_base_position( 0);
      for
      (
        linal::Vector< float>::iterator itr_3DA_smooth( STORAGE.Begin()),
        itr_3DA_smooth_end( STORAGE.End()), itr_3DA_rough_outer( m_DiscreteCode.Begin());
        itr_3DA_smooth != itr_3DA_smooth_end;
        ++itr_3DA_smooth, --smoothing_vector_base_position, ++itr_3DA_rough_outer
      )
      {
        // handle the case where no values were ever seen for that bin
        if( *itr_3DA_rough_outer == math::GetLowestBoundedValue< float>())
        {
          *itr_3DA_smooth = default_value;
          continue;
        }
        else if( !m_Smooth)
        {
          *itr_3DA_smooth = *itr_3DA_rough_outer;
          continue;
        }
        float bin_smooth( 0);
        int smoothing_vector_position( smoothing_vector_base_position);
        for
        (
          linal::Vector< float>::const_iterator itr_3DA_rough( m_DiscreteCode.Begin()), itr_3DA_rough_end( m_DiscreteCode.End());
          itr_3DA_rough != itr_3DA_rough_end;
          ++itr_3DA_rough, ++smoothing_vector_position
        )
        {
          if( *itr_3DA_rough != math::GetLowestBoundedValue< float>())
          {
            bin_smooth += *itr_3DA_rough * m_SmoothingCoefficients( std::abs( smoothing_vector_position));
          }
        }
        *itr_3DA_smooth = bin_smooth;
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule3DASoftMax::GetSerializer() const
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
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule3DASoftMax::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_DiscreteCode = linal::Vector< float>( GetNormalSizeOfFeatures(), math::GetLowestBoundedValue< float>());
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
  } // namespace descriptor
} // namespace bcl
