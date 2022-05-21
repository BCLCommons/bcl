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
#include "descriptor/bcl_descriptor_molecule_rdf_grid_code.h"

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
    const util::SiPtr< const util::ObjectInterface> MoleculeRDFGridCode::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeRDFGridCode()
      )
    );

    //! @brief default constructor
    MoleculeRDFGridCode::MoleculeRDFGridCode() :
      m_AtomProperty(),
      m_WeightProperty(),
      m_DistanceNumberSteps( 24),
      m_PropertyNumberSteps( 12),
      m_DistanceStepSize   ( 0.5),
      m_PropertyStepSize   ( 0.5),
      m_DistanceTemperature( 100.0),
      m_PropertyTemperature( 100.0)
    {
      std::stringstream dummy_stream;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), dummy_stream);
    }

    //! @brief constructor from number of steps, and mapped atom property
    MoleculeRDFGridCode::MoleculeRDFGridCode
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const CheminfoProperty &WEIGHT_PROPERTY,
      const size_t       &DISTANCE_NUMBER_STEPS,
      const size_t       &PROPERTY_NUMBER_STEPS,
      const float        &DISTANCE_STEP_SIZE,
      const float        &PROPERTY_STEP_SIZE,
      const float        &DISTANCE_TEMPERATURE,
      const float        &PROPERTY_TEMPERATURE
    ) :
      m_AtomProperty       ( ATOM_PROPERTY        ),
      m_WeightProperty     ( WEIGHT_PROPERTY      ),
      m_DistanceNumberSteps( DISTANCE_NUMBER_STEPS),
      m_PropertyNumberSteps( PROPERTY_NUMBER_STEPS),
      m_DistanceStepSize   ( DISTANCE_STEP_SIZE   ),
      m_PropertyStepSize   ( PROPERTY_STEP_SIZE   ),
      m_DistanceTemperature( DISTANCE_TEMPERATURE ),
      m_PropertyTemperature( PROPERTY_TEMPERATURE )
    {
      std::stringstream dummy_stream;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), dummy_stream);
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeRDFGridCode
    MoleculeRDFGridCode *MoleculeRDFGridCode::Clone() const
    {
      return new MoleculeRDFGridCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeRDFGridCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeRDFGridCode::GetAlias() const
    {
      static const std::string s_name( "RDFGrid");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > MoleculeRDFGridCode::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomProperty,
          &m_WeightProperty + 1
        );
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeRDFGridCode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // get the overall grid size
      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT);
      float prop_a( m_AtomProperty->operator ()( iterator_a)( 0));
      float weight_a( m_WeightProperty->operator ()( iterator_a)( 0));

      // iterate over distances and the remaining property values on the upper triangle of this matrix
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atom_b( ELEMENT.Begin());
        itr_atom_b.NotAtEnd();
        ++itr_atom_b
      )
      {
        if( itr_atom_b == ELEMENT)
        {
          const float weight( math::Sqr( weight_a));
          const float *itr_auto_rdf( m_AutoRDFMatrix.Begin());
          const float *itr_auto_rdf_end( m_AutoRDFMatrix.End());
          for( float *itr_rdf( STORAGE.Begin()); itr_auto_rdf != itr_auto_rdf_end; ++itr_auto_rdf, ++itr_rdf)
          {
            *itr_rdf += weight * ( *itr_auto_rdf);
          }
          continue;
        }

        Iterator< chemistry::AtomConformationalInterface> iterator_b( itr_atom_b);
        //store distance between both atoms
        const float distance( linal::Distance( ELEMENT->GetPosition(), itr_atom_b->GetPosition()));

        float prop_b( m_AtomProperty->operator ()( iterator_b)( 0));
        float weight_b( m_WeightProperty->operator ()( iterator_b)( 0));

        // store distance between properties
        const float property_distance( math::Absolute( ( prop_a) - ( prop_b)));

        // Scale property by  weight
        float weight( ( weight_a) * ( weight_b));

        if
        (
          distance > m_DistanceStepSize * m_DistanceNumberSteps
          || property_distance > m_PropertyNumberSteps * m_PropertyNumberSteps
        )
        {
          continue;
        }

        // Generate the actual RDF code
        float *itr_rdf( STORAGE.Begin());
        float sum_influence( 0.0);
        for( size_t dist_steps( 0); dist_steps < m_DistanceNumberSteps; ++dist_steps)
        {
          const float three_d_distance
          (
            -m_DistanceTemperature * math::Sqr( m_DistanceStepSize * dist_steps - distance)
          );
          for( size_t prop_steps( 0); prop_steps < m_PropertyNumberSteps; ++prop_steps)
          {
            const float effective_property_distance
            (
              -m_PropertyTemperature * math::Sqr( m_PropertyStepSize * prop_steps - property_distance)
            );
            sum_influence += exp( three_d_distance + effective_property_distance);
          }
        }
        if( sum_influence < 1.0e-8)
        {
          continue;
        }
        weight *= m_AutoRDFMatrixSum / sum_influence;
        for( size_t dist_steps( 0); dist_steps < m_DistanceNumberSteps; ++dist_steps)
        {
          const float three_d_distance
          (
            -m_DistanceTemperature * math::Sqr( m_DistanceStepSize * dist_steps - distance)
          );
          for( size_t prop_steps( 0); prop_steps < m_PropertyNumberSteps; ++prop_steps, ++itr_rdf)
          {
            const float effective_property_distance
            (
              -m_PropertyTemperature * math::Sqr( m_PropertyStepSize * prop_steps - property_distance)
            );
            *itr_rdf += weight * exp( three_d_distance + effective_property_distance);
          }
        }
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeRDFGridCode::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes the radial distribution function using a given atom property, see\n"
        "http://www.opus.ub.uni-erlangen.de/opus/volltexte/2007/736/pdf/MarkusHemmerDissertation.pdf, p. 65 for details"
      );

      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the radial distribution function",
        io::Serialization::GetAgent( &m_AtomProperty)
      );
      parameters.AddInitializer
      (
        "weight property",
        "property whose product will weight the radial distribution function",
        io::Serialization::GetAgent( &m_WeightProperty)
      );
      parameters.AddInitializer
      (
        "distance steps",
        "number of distance bins (each of size distance size) used in the radial distribution function",
        io::Serialization::GetAgentWithRange( &m_DistanceNumberSteps, 1, 1000000),
        "24"
      );
      parameters.AddInitializer
      (
        "property steps",
        "number of property bins (each of size property step size) used in the radial distribution function",
        io::Serialization::GetAgentWithRange( &m_PropertyNumberSteps, 1, 1000000),
        "12"
      );
      parameters.AddInitializer
      (
        "distance step size",
        "size of each step for the distance axis in angstroms",
        io::Serialization::GetAgentWithRange( &m_DistanceStepSize, 0.01, 100.0),
        "0.5"
      );
      parameters.AddInitializer
      (
        "property step size",
        "size of each step for the property axis in angstroms",
        io::Serialization::GetAgentWithRange( &m_PropertyStepSize, 0.01, 10000.0),
        "0.5"
      );
      parameters.AddInitializer
      (
        "distance temperature",
        "increasing temperature spreads autocorrelation across more distant bins",
        io::Serialization::GetAgentWithRange( &m_DistanceTemperature, 0.0, 1000.0),
        "100"
      );
      parameters.AddInitializer
      (
        "property temperature",
        "same as distance temperature but used to distribute values of atom property over more distant bins",
        io::Serialization::GetAgentWithRange( &m_PropertyTemperature, 0.0, 1000.0),
        "100"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeRDFGridCode::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_AutoRDFMatrix = linal::Matrix< float>( m_DistanceNumberSteps, m_PropertyNumberSteps, float( 0));

      // instantiate matrix to store coefficients of the rdf code for atoms on themselves
      m_AutoRDFMatrixSum = 0.0;
      for( size_t dist_steps( 0); dist_steps < m_DistanceNumberSteps; ++dist_steps)
      {
        for( size_t prop_steps( 0); prop_steps < m_PropertyNumberSteps; ++prop_steps)
        {
          m_AutoRDFMatrixSum += m_AutoRDFMatrix( dist_steps, prop_steps) =
            0.5 * exp
            (
              -m_DistanceTemperature * math::Sqr( m_DistanceStepSize * dist_steps)
              - m_PropertyTemperature * math::Sqr( m_PropertyStepSize * prop_steps)
            );
        }
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
      else if( m_WeightProperty.IsDefined() && m_WeightProperty->GetNormalSizeOfFeatures() != 1)
      {
        ERR_STREAM
          << "Expected a weight property that returned 1 property per atom, but property returns "
          << m_WeightProperty->GetNormalSizeOfFeatures()
          << " values per atom ( property was "
          << m_WeightProperty->GetString()
          << ")";
        return false;
      }

      return true;
    }

  } // namespace descriptor
} // namespace bcl
