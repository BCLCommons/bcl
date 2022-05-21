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
#include "descriptor/bcl_descriptor_molecule_3da_soft_max_sign.h"

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
    const util::SiPtr< const util::ObjectInterface> Molecule3DASoftMaxSign::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule3DASoftMaxSign()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Molecule3DASoftMaxSign::Molecule3DASoftMaxSign() :
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
    Molecule3DASoftMaxSign::Molecule3DASoftMaxSign
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
    //! @return pointer to new Molecule3DASoftMaxSign
    Molecule3DASoftMaxSign *Molecule3DASoftMaxSign::Clone() const
    {
      return new Molecule3DASoftMaxSign( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule3DASoftMaxSign::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule3DASoftMaxSign::GetAlias() const
    {
      static const std::string s_name( "3daSoftMaxSign");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule3DASoftMaxSign::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float Molecule3DASoftMaxSign::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get temperature of code
    //! @return const float  temperature of 3DA code
    const float &Molecule3DASoftMaxSign::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule3DASoftMaxSign::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule3DASoftMaxSign::GetInternalDescriptors()
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
    void Molecule3DASoftMaxSign::Calculate( linal::VectorReference< float> &STORAGE)
    {
      m_DiscreteCode = 0.0;

      // catch NUMBER_STEPS == 0
      if( m_NumberSteps == 0)
      {
        BCL_MessageCrt( "NUMBER_STEPS equals zero - 3DA code will be empty!");
        return;
      }
      else if( m_AtomProperty->GetNormalSizeOfFeatures() != 1)
      {
        BCL_MessageCrt
        (
          "Expected an atomic property that returned 1 property per atom, but property returns "
          + util::Format()( m_AtomProperty->GetNormalSizeOfFeatures())
          + " values per atom (atom property was "
          + util::Format()( m_AtomProperty->GetString())
          + ")"
        );
        return;
      }

      // instantiate the property as a vector with indices that correspond to atoms
      const linal::Vector< float> property( m_AtomProperty->CollectValuesOnEachElementOfObject( *this->GetCurrentObject()));

      // iterate over all possible pairs of atoms
      // iterate properties and surface areas simultaneously
      linal::Vector< float>::const_iterator itr_prop_a( property.Begin());
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface>
            itr_atoms_a( this->GetCurrentObject()->GetIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a, ++itr_prop_a
      )
      {
        // compute the property value for the atom onto itself
        const float prop( math::Sqr( *itr_prop_a));

        // differentiate between positive and negative atoms
        int sign_id_a( *itr_prop_a > 0);

        // generate rough 3DA code for each atom onto itself
        m_DiscreteCode( sign_id_a) = std::max( m_DiscreteCode( sign_id_a), prop);

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

          int sign_id( sign_id_a);

          // differentiate between combinations of positive and negative properties
          if( prop < 0)
          {
            prop = -prop;
            sign_id = 2;
          }

          // generate rough 3DA for boundry distances
          float &three_da_value( m_DiscreteCode( 3 * nominal_bin + sign_id));
          three_da_value = std::max( prop, three_da_value);
        }
      }

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
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule3DASoftMaxSign::GetSerializer() const
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
    bool Molecule3DASoftMaxSign::ReadInitializerSuccessHook
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
  } // namespace descriptor
} // namespace bcl
