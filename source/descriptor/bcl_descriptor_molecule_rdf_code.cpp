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
#include "descriptor/bcl_descriptor_molecule_rdf_code.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeRDFCode::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeRDFCode()
      )
    );

    //! @brief default constructor
    MoleculeRDFCode::MoleculeRDFCode() :
      m_NumberSteps( 48),
      m_StepSize( 0.1),
      m_Temperature( 0.0)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief constructor from number of steps, and mapped atom property
    MoleculeRDFCode::MoleculeRDFCode
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
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeRDFCode
    MoleculeRDFCode *MoleculeRDFCode::Clone() const
    {
      return new MoleculeRDFCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeRDFCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeRDFCode::GetAlias() const
    {
      static const std::string s_name( "RDF");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t MoleculeRDFCode::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float MoleculeRDFCode::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float MoleculeRDFCode::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &MoleculeRDFCode::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > MoleculeRDFCode::GetInternalDescriptors()
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
    void MoleculeRDFCode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT);
      float prop_a( m_AtomProperty->operator ()( iterator_a)( 0));

      // iterate over distances and the remaining property values on the upper triangle of this matrix
      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_b( ELEMENT.Begin());
        iterator_b.NotAtEnd();
        ++iterator_b
      )
      {
        if( iterator_b( 0) == ELEMENT)
        {
          // skip atom onto itself for backwards compatibility
          continue;
        }

        const chemistry::AtomConformationalInterface &atom_b( *iterator_b( 0));

        //store distance between both atoms
        const float distance( linal::Distance( ELEMENT->GetPosition(), atom_b.GetPosition()));

        float prop_b( m_AtomProperty->operator ()( iterator_b)( 0));

        // Scale property by  surface areas.
        const float prop( ( prop_a) * ( prop_b));

        // Generate the actual RDF code
        for( size_t temp_steps( 0); temp_steps < m_NumberSteps; ++temp_steps)
        {
          STORAGE( temp_steps) += prop * exp( -m_Temperature * math::Sqr( m_StepSize * temp_steps - distance));
        }
      }

      // normalization for backwards compatibility
      STORAGE *= 0.5;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeRDFCode::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes the radial distribution function using a given atom property"
      );

      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the radial distribution function",
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
      static bool s_normalized( false);
      parameters.AddOptionalInitializer
      (
        "normalized",
        "deprecated; no longer used",
        io::Serialization::GetAgent( &s_normalized)
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeRDFCode::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // catch m_NumberSteps == 0
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
