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
#include "descriptor/bcl_descriptor_molecule_2da_smooth_sign_code.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_connectivity.h"
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
    const util::SiPtr< const util::ObjectInterface> Molecule2DASmoothSignCode::s_SmoothInstance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule2DASmoothSignCode( e_2DASmoothSign)
      )
    );

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule2DASmoothSignCode::s_SignInstance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule2DASmoothSignCode( e_2DASign)
      )
    );

    //! @brief default constructor
    Molecule2DASmoothSignCode::Molecule2DASmoothSignCode( const Method SMOOTH) :
      m_NumberSteps( 12),
      m_Temperature( 100.0),
      m_Smooth( SMOOTH == e_2DASmoothSign)
    {
      m_DiscreteCode = linal::Vector< float>( GetNormalSizeOfFeatures(), 0.0);
      if( m_Smooth)
      {
        m_SmoothingCoefficients = GetSmoothingCoefficientVector( m_NumberSteps, m_Temperature);
      }
    }

    //! @brief constructor from number of steps, and mapped atom property
    Molecule2DASmoothSignCode::Molecule2DASmoothSignCode
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS,
      const float &TEMPERATURE,
      const Method SMOOTH
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_Temperature( TEMPERATURE),
      m_Smooth( SMOOTH == e_2DASmoothSign)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule2DASmoothSignCode
    Molecule2DASmoothSignCode *Molecule2DASmoothSignCode::Clone() const
    {
      return new Molecule2DASmoothSignCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule2DASmoothSignCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule2DASmoothSignCode::GetAlias() const
    {
      static const std::string s_name_smooth( "2DASmoothSign");
      static const std::string s_name_sign( "2DASign");
      return m_Smooth ? s_name_smooth : s_name_sign;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule2DASmoothSignCode::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get temperature of code
    //! @return const float temperature of 3DA code
    const float &Molecule2DASmoothSignCode::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule2DASmoothSignCode::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule2DASmoothSignCode::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomProperty,
          &m_AtomProperty + 1
        );
    }

    //! @brief Get a simple pointer to a all smoothing coefficient vectors
    //!        Because the same parameters are often reused across many instances of this class, it is helpful
    //!        to have only a single smoothing coefficients vector
    //! @param NUMBER_STEPS number of steps needed
    //! @param TEMPERATURE gaussian-smoothing temperature
    linal::VectorConstReference< float> Molecule2DASmoothSignCode::GetSmoothingCoefficientVector
    (
      const size_t &NUMBER_STEPS,
      const float &TEMPERATURE
    )
    {
      static std::map< std::pair< size_t, std::pair< float, float> >, linal::Vector< float> > s_map;
      static sched::Mutex s_mutex;
      s_mutex.Lock();
      std::pair< size_t, std::pair< float, float> > query( NUMBER_STEPS, std::make_pair( TEMPERATURE, 1));
      std::map< std::pair< size_t, std::pair< float, float> >, linal::Vector< float> >::const_iterator
        itr( s_map.find( query));
      if( itr != s_map.end())
      {
        // already have the smoothing coefficients in the map
        s_mutex.Unlock();
        return linal::VectorConstReference< float>( itr->second);
      }
      // create a new vector with the smoothing coefficients
      float current_bin_position( 0);
      linal::Vector< float> smoothing_coefficients( NUMBER_STEPS);
      for( size_t step_number( 0); step_number < NUMBER_STEPS; ++step_number, current_bin_position += 1)
      {
        smoothing_coefficients( step_number) = exp( -TEMPERATURE * math::Sqr( current_bin_position));
      }
      itr = s_map.insert( std::make_pair( query, smoothing_coefficients)).first;
      s_mutex.Unlock();
      return linal::VectorConstReference< float>( itr->second);
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule2DASmoothSignCode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_A,
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_B,
      linal::VectorReference< float> &STORAGE
    )
    {
      m_DiscreteCode = 0.0;

      // create the graph, if necessary
      if( m_Graph.GetSize() == size_t( 0))
      {
        chemistry::ConformationGraphConverter graph_maker;
        util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());
        m_Graph = graph_maker( *conformation);
      }

      // catch m_NumberSteps == 0
      if( m_NumberSteps == 0)
      {
        BCL_MessageCrt( "Because m_NumberSteps equals zero, 2DASmoothSign code will be empty!");
        return;
      }

      //iterate over all possible pairs of atom property values
      size_t index_prop_a( ELEMENT_A.GetPosition());
      size_t index_prop_b( ELEMENT_B.GetPosition());

      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT_A);
      Iterator< chemistry::AtomConformationalInterface> iterator_b( ELEMENT_B);

      float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

      // determine the distance from this vertex to all vertices
      util::ShPtr< storage::Vector< size_t> > distances
      (
        graph::Connectivity::DistancesToOtherVertices( m_Graph, index_prop_a)
      );

      if( index_prop_a == index_prop_b)
      {
        m_DiscreteCode( int( property_atom_a > 0.0)) += property_atom_a * property_atom_a;
        return;
      }
      else
      {
        float property_atom_b( m_AtomProperty->operator ()( iterator_b)( 0));
        const size_t distance( distances->operator()( index_prop_b));
        Accumulate( distance, property_atom_a, property_atom_b);
      }
    } // Recalculate

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule2DASmoothSignCode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // create the graph, if needed
      if( m_Graph.GetSize() == size_t( 0))
      {
        chemistry::ConformationGraphConverter graph_maker;
        util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());
        m_Graph = graph_maker( *conformation);
      }

      // catch m_NumberSteps == 0
      if( m_NumberSteps == 0)
      {
        BCL_MessageCrt( "Because m_NumberSteps equals zero, 2DASmoothSign code will be empty code will be empty!");
        return;
      }

      m_DiscreteCode = 0.0;

      //iterate over all possible pairs of atom property values
      size_t index_prop_a( ELEMENT.GetPosition());

      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT);
      // determine the distance from this vertex to all vertices
      util::ShPtr< storage::Vector< size_t> > distances
      (
        graph::Connectivity::DistancesToOtherVertices( m_Graph, index_prop_a)
      );

      float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

      // get an iterator on the distances
      storage::Vector< size_t>::const_iterator itr_distances( distances->Begin());

      // iterate over distances and the remaining property values on the upper triangle of this matrix
      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_b( ELEMENT.Begin());
        iterator_b.NotAtEnd();
        ++iterator_b, ++itr_distances
      )
      {
        if( *itr_distances < m_NumberSteps)
        {
          const float property_atom_b( m_AtomProperty->operator ()( iterator_b)( 0));
          Accumulate( *itr_distances, property_atom_a, property_atom_b);
        }
      }

      // this is for backwards compatibility
      m_DiscreteCode *= 0.5;
      Smooth( STORAGE);
    } // Recalculate

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void Molecule2DASmoothSignCode::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      // create the graph, if necessary
      if( m_Graph.GetSize() == size_t( 0))
      {
        chemistry::ConformationGraphConverter graph_maker;
        util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());
        m_Graph = graph_maker( *conformation);
      }

      // catch m_NumberSteps == 0
      if( m_NumberSteps == 0)
      {
        BCL_MessageCrt( "Because m_NumberSteps equals zero, 2DASmoothSign code will be empty!");
        return;
      }

      m_DiscreteCode = 0.0;

      //iterate over all possible pairs of atom property values
      size_t index_a( 0);
      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_a( this->GetCurrentObject()->GetIterator());
        iterator_a.NotAtEnd();
        ++iterator_a, ++index_a
      )
      {
        float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

        // skip 0 properties
        if( !property_atom_a)
        {
          continue;
        }

        // determine the distance from this vertex to all vertices
        util::ShPtr< storage::Vector< size_t> > distances
        (
          graph::Connectivity::DistancesToOtherVertices( m_Graph, index_a)
        );

        //differentiate pos and neg atoms
        m_DiscreteCode( int( property_atom_a > 0)) += 0.5 * math::Sqr( property_atom_a);

        // iterate over distances and the remaining property values on the upper triangle of this matrix
        Iterator< chemistry::AtomConformationalInterface> iterator_b( iterator_a);
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_distances( distances->Begin() + index_a), itr_distances_end( distances->End());
          itr_distances != itr_distances_end;
          ++itr_distances, ++iterator_b
        )
        {
          float property_atom_b( m_AtomProperty->operator ()( iterator_b)( 0));

          if( *itr_distances < m_NumberSteps)
          {
            Accumulate( *itr_distances, property_atom_a, property_atom_b);
          }
        }
      }
      Smooth( STORAGE);
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Molecule2DASmoothSignCode::SetObjectHook()
    {
      m_Graph = graph::ConstGraph< size_t, size_t>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule2DASmoothSignCode::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes 2D (bond distance) autocorrelation of a specified atom property"
      );

      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the smooth radial distribution function",
        io::Serialization::GetAgent( &m_AtomProperty)
      );
      parameters.AddInitializer
      (
        "steps",
        "# of steps/bins (each of size = step size) used in the radial distribution function",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "48"
      );
      if( m_Smooth)
      {
        parameters.AddInitializer
        (
          "temperature",
          "increasing temperature spreads autocorrelation across more distant bins",
          io::Serialization::GetAgentWithRange( &m_Temperature, 0.0, 1000.0),
          "100"
        );
      }
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule2DASmoothSignCode::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_DiscreteCode = linal::Vector< float>( GetNormalSizeOfFeatures(), 0.0);
      m_SmoothingCoefficients =
          Molecule2DASmoothSignCode::GetSmoothingCoefficientVector( m_NumberSteps, m_Temperature);

      if( m_NumberSteps == 0)
      {
        ERR_STREAM << "m_NumberSteps equals zero - 2DA code will be empty!";
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
    void Molecule2DASmoothSignCode::Smooth( linal::VectorReference< float> &STORAGE) const
    {
      //catch for not running Smoothing protocol
      if( !m_Smooth)
      {
        STORAGE.CopyValues( m_DiscreteCode);
        return;
      }

      // generate the actual 2DA smooth code by applying the smoothing kernel over the 2DA
      int smoothing_vector_base_position( 0);

      for
      (
        linal::Vector< float>::iterator itr_2DA_smooth( STORAGE.Begin()),
        itr_2DA_smooth_end( STORAGE.End());
        itr_2DA_smooth != itr_2DA_smooth_end;
        ++itr_2DA_smooth, --smoothing_vector_base_position
      )
      {
        float neg_neg( 0), pos_pos( 0), pos_neg( 0);
        int smoothing_vector_position( smoothing_vector_base_position);
        for
        (
          linal::Vector< float>::const_iterator itr_2DA_rough( m_DiscreteCode.Begin()),
          itr_2DA_rough_end( m_DiscreteCode.End());
          itr_2DA_rough != itr_2DA_rough_end;
          ++itr_2DA_rough, ++smoothing_vector_position
        )
        {
          const float exponential_factor( m_SmoothingCoefficients( std::abs( smoothing_vector_position)));
          neg_neg += *itr_2DA_rough * exponential_factor;
          pos_pos += *++itr_2DA_rough * exponential_factor;
          pos_neg += *++itr_2DA_rough * exponential_factor;
        }

        *itr_2DA_smooth = neg_neg;
        *++itr_2DA_smooth = pos_pos;
        *++itr_2DA_smooth = pos_neg;
      }
    }

    //! @brief add an observed distance/property value to m_DiscreteCode
    //! @param DISTANCE actual distance of the two atoms
    //! @param PROP_A property from atom A
    //! @param PROP_B property from atom B
    void Molecule2DASmoothSignCode::Accumulate
    (
      const size_t &DISTANCE,
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
      m_DiscreteCode( 3 * DISTANCE + sign_id) += product;
    }

  } // namespace descriptor
} // namespace bcl
