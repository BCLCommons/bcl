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
#include "descriptor/bcl_descriptor_molecule_rdf_max_sign_code.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeRDFMaxSignCode::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeRDFMaxSignCode()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeRDFMaxSignCode::MoleculeRDFMaxSignCode() :
      m_NumberSteps( 0),
      m_StepSize( 0.0),
      m_Temperature( 0.0)
    {
    }

    //! @brief constructor from number of steps, and mapped atom property
    MoleculeRDFMaxSignCode::MoleculeRDFMaxSignCode
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const float &TEMPERATURE
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_Temperature( TEMPERATURE)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeRDFMaxSignCode
    MoleculeRDFMaxSignCode *MoleculeRDFMaxSignCode::Clone() const
    {
      return new MoleculeRDFMaxSignCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeRDFMaxSignCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeRDFMaxSignCode::GetAlias() const
    {
      static const std::string s_name( "RDFMaxSign");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t MoleculeRDFMaxSignCode::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float MoleculeRDFMaxSignCode::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get temperature of code
    //! @return const float  temperature of 3DA code
    const float &MoleculeRDFMaxSignCode::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &MoleculeRDFMaxSignCode::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > MoleculeRDFMaxSignCode::GetInternalDescriptors()
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
    void MoleculeRDFMaxSignCode::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // catch NUMBER_STEPS == 0
      if( m_NumberSteps == 0)
      {
        BCL_MessageCrt( "NUMBER_STEPS equals zero - RDF code will be empty!");
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

      // Iterate over all possible pairs of atoms
      // Iterate properties and surface areas simultaneously
      linal::Vector< float>::const_iterator itr_prop_a( property.Begin());

      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface>
          itr_atoms_a( this->GetCurrentObject()->GetIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a, ++itr_prop_a
      )
      {

        // scale by the atom property for the atom onto itself
        const float prop( 0.5 * math::Sqr( *itr_prop_a));

        // Differentiate between positive and negative atoms
        const int a_sign_id( *itr_prop_a > 0);

        // Generate the actual RDF code for the case of the atom being paired with itself
        for( size_t temp_steps( 0), bin_steps( a_sign_id); temp_steps < m_NumberSteps; ++temp_steps, bin_steps += 3)
        {
          float temp_value( prop * exp( -m_Temperature * math::Sqr( m_StepSize * temp_steps)));
          STORAGE( bin_steps) = std::max( temp_value, STORAGE( bin_steps));
        }

        linal::Vector< float>::const_iterator itr_prop_b( itr_prop_a + 1);

        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_b( itr_atoms_a);
        for( ++itr_atoms_b; itr_atoms_b.NotAtEnd(); ++itr_atoms_b, ++itr_prop_b)
        {
          //store distance between both atoms
          const float distance
          (
            linal::Distance
            (
              itr_atoms_a->GetPosition(),
              itr_atoms_b->GetPosition()
            )
          );

          // scale by the atom property
          float prop( *itr_prop_a * ( *itr_prop_b));

          //Differentiate between combinations of positive and negative properties
          int sign_id( a_sign_id);

          if( prop < 0)
          {
            prop = -prop;
            sign_id = 2;
          }

          // Generate the actual RDF code for all other unique pairs
          for( size_t temp_steps( 0), bin_steps( sign_id); temp_steps < m_NumberSteps; ++temp_steps, bin_steps += 3)
          {
            float temp_value = prop * exp( -m_Temperature * math::Sqr( m_StepSize * temp_steps - distance));
            STORAGE( bin_steps) = std::max( STORAGE( bin_steps), temp_value);
          }
        }
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeRDFMaxSignCode::GetSerializer() const
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
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
