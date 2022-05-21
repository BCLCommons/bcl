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
#include "descriptor/bcl_descriptor_molecule_3d_sign_distribution.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule3DSignDistribution::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule3DSignDistribution()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Molecule3DSignDistribution::Molecule3DSignDistribution() :
      m_NumberSteps( 12),
      m_StepSize( 1.0)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief constructor from number of steps, and mapped atom property
    Molecule3DSignDistribution::Molecule3DSignDistribution
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const CheminfoProperty &ATOM_CENTER_PROPERTY,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_CenterProperty( ATOM_CENTER_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule3DSignDistribution
    Molecule3DSignDistribution *Molecule3DSignDistribution::Clone() const
    {
      return new Molecule3DSignDistribution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule3DSignDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule3DSignDistribution::GetAlias() const
    {
      static const std::string s_name( "3dDistributionSign");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule3DSignDistribution::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float Molecule3DSignDistribution::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule3DSignDistribution::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule3DSignDistribution::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomProperty,
          &m_CenterProperty + 1
        );
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void Molecule3DSignDistribution::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      m_DiscreteCode = 0.0;
      util::SiPtr< const coord::OrientationInterface> orientation
      (
        this->GetCurrentObject()
      );

      math::RunningAverage< linal::Vector3D> weighted_center;
      iterate::Generic< const chemistry::AtomConformationalInterface> itr_atom( this->GetCurrentObject()->GetIterator());
      size_t position( 0);
      for
      (
        Iterator< chemistry::AtomConformationalInterface> itr_desc( itr_atom);
        itr_atom.NotAtEnd();
        ++itr_atom, ++itr_desc, ++position
      )
      {
        const float weight( ( *m_CenterProperty)( itr_desc)( 0));
        weighted_center.AddWeightedObservation( itr_atom->GetPosition(), weight);
      }
      if( weighted_center.GetWeight() <= float( 0.0))
      {
        STORAGE = float( 0.0);
        return;
      }

      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_a( this->GetCurrentObject()->GetIterator());
        iterator_a.NotAtEnd();
        ++iterator_a
      )
      {
        const chemistry::AtomConformationalInterface &atom_a( *iterator_a( 0));
        const float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));
        const double distance( linal::Distance( atom_a.GetPosition(), weighted_center.GetAverage()));
        Accumulate( distance, property_atom_a);
      }
      STORAGE.CopyValues( m_DiscreteCode);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule3DSignDistribution::GetSerializer() const
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
        "center property",
        "property used to weight the atom distances. Should usually be something that is always >= 0",
        io::Serialization::GetAgent( &m_CenterProperty)
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
        "steps",
        "# of steps/bins (each of size = step size) used in the radial distribution function",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "48"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule3DSignDistribution::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_DiscreteCode = linal::Vector< float>( GetNormalSizeOfFeatures(), 0.0);
      if( m_NumberSteps == 0 || m_StepSize == 0.0)
      {
        ERR_STREAM << "m_NumberSteps equals zero - 3DA code will be empty!";
        return false;
      }
      bool atom_prop_bad( m_AtomProperty.IsDefined() && m_AtomProperty->GetNormalSizeOfFeatures() != 1);
      if
      (
        atom_prop_bad
        ||
        ( m_CenterProperty.IsDefined() && m_CenterProperty->GetNormalSizeOfFeatures() != 1)
      )
      {
        ERR_STREAM
           << "Expected a property that returned 1 properties per atom, but property returns "
           << ( atom_prop_bad ? m_AtomProperty->GetNormalSizeOfFeatures() : m_CenterProperty->GetNormalSizeOfFeatures())
           << " values per atom ( property was "
           << ( atom_prop_bad ? m_AtomProperty->GetString() : m_CenterProperty->GetString())
           << ")";
        return false;
      }
      return true;
    }

    //! @brief add an observed distance/property value to m_DiscreteCode
    //! @param DISTANCE actual distance of the two atoms
    //! @param PROP_A property from atom A
    //! @param PROP_B property from atom B
    void Molecule3DSignDistribution::Accumulate
    (
      const float &DISTANCE,
      const float &PROP_A
    )
    {
      // calculate exponential function for the RDF equation
      // assign atom distance to the closest 3DA bin
      size_t closest_step_low( std::min( size_t( m_NumberSteps - 1), size_t( DISTANCE / m_StepSize)));
      size_t closest_step_high( closest_step_low + 1);
      size_t offset( PROP_A < 0.0 ? m_NumberSteps : 0);
      double prop_a( math::Absolute( PROP_A));

      if( closest_step_high == m_NumberSteps)
      {
        m_DiscreteCode( offset + closest_step_low) += PROP_A;
      }
      else
      {
        const float dist_lower_step( m_StepSize * closest_step_low - DISTANCE);
        const float dist_higher_step( dist_lower_step + m_StepSize);

        float low_exp_factor( std::max( dist_higher_step, float( 0.0)) / m_StepSize);
        float high_exp_factor( dist_higher_step >= float( 0.0) ? -dist_lower_step / m_StepSize : 0.0);

        if( closest_step_high == m_NumberSteps)
        {
          // generate rough 3DA for boundry distances
          m_DiscreteCode( offset + closest_step_low) += prop_a * low_exp_factor;
        }
        else
        {
          // normalize rough 3DA code
          double product( prop_a / ( low_exp_factor + high_exp_factor));
          // generate rough 3DA code for all other distances
          m_DiscreteCode( offset + closest_step_low) += product * low_exp_factor;
          m_DiscreteCode( offset + closest_step_high) += product * high_exp_factor;
        }
      }
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Molecule3DSignDistribution::SetObjectHook()
    {
    }

  } // namespace descriptor
} // namespace bcl
