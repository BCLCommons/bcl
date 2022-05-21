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
#include "descriptor/bcl_descriptor_molecule_2da_code.h"

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
    const util::SiPtr< const util::ObjectInterface> Molecule2DACode::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule2DACode()
      )
    );

    //! @brief default constructor
    Molecule2DACode::Molecule2DACode() :
      m_NumberSteps( 12)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief constructor from number of steps, and mapped atom property
    Molecule2DACode::Molecule2DACode
    (
      const size_t NUMBER_STEPS,
      const CheminfoProperty &ATOM_PROPERTY
    ) :
      m_NumberSteps( NUMBER_STEPS),
      m_AtomProperty( ATOM_PROPERTY)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule2DACode
    Molecule2DACode *Molecule2DACode::Clone() const
    {
      return new Molecule2DACode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule2DACode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule2DACode::GetAlias() const
    {
      static const std::string s_name( "2DA");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule2DACode::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule2DACode::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule2DACode::GetInternalDescriptors()
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
    void Molecule2DACode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_A,
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_B,
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
        BCL_MessageCrt( "m_NumberSteps equals zero - 2da code will be empty!");
        return;
      }

      //iterate over all possible pairs of atom property values
      size_t index_prop_a( ELEMENT_A.GetPosition());
      size_t index_prop_b( ELEMENT_B.GetPosition());

      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT_A);
      const float prop_a( m_AtomProperty->operator ()( iterator_a)( 0));
      if( index_prop_a == index_prop_b)
      {
        STORAGE( 0) += prop_a * prop_a;
        return;
      }

      // determine the distance from this vertex to all vertices
      util::ShPtr< storage::Vector< size_t> > distances
      (
        graph::Connectivity::DistancesToOtherVertices( m_Graph, index_prop_a)
      );

      const size_t distance( distances->operator()( index_prop_b));
      if( distance < m_NumberSteps)
      {
        Iterator< chemistry::AtomConformationalInterface> iterator_b( ELEMENT_B);
        const float prop_b( m_AtomProperty->operator ()( iterator_b)( 0));
        STORAGE( distance) += 0.5 * prop_a * prop_b;
      }
    } // Recalculate

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule2DACode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
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
        BCL_MessageCrt( "m_NumberSteps equals zero - 2da code will be empty!");
        return;
      }

      //iterate over all possible pairs of atom property values
      size_t index_prop_a( ELEMENT.GetPosition());

      // determine the distance from this vertex to all vertices
      util::ShPtr< storage::Vector< size_t> > distances
      (
        graph::Connectivity::DistancesToOtherVertices( m_Graph, index_prop_a)
      );

      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT);
      float prop_a( m_AtomProperty->operator ()( iterator_a)( 0));

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
          const float prop_b( m_AtomProperty->operator ()( iterator_b)( 0));
          STORAGE( *itr_distances) += prop_a * prop_b;
        }
      }

      // normalization for backwards compatibility
      STORAGE( 0) *= 2.0;
      STORAGE *= 0.5;
    } // Recalculate

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void Molecule2DACode::Calculate( linal::VectorReference< float> &STORAGE)
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
        BCL_MessageCrt( "m_NumberSteps equals zero - 2da code will be empty!");
        return;
      }

      //iterate over all possible pairs of atom property values
      size_t index_a( 0);
      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_a( this->GetCurrentObject()->GetIterator());
        iterator_a.NotAtEnd();
        ++iterator_a, ++index_a
      )
      {
        float prop_a( m_AtomProperty->operator ()( iterator_a)( 0));

        // skip 0 properties
        if( !prop_a)
        {
          continue;
        }

        // determine the distance from this vertex to all vertices
        util::ShPtr< storage::Vector< size_t> > distances
        (
          graph::Connectivity::DistancesToOtherVertices( m_Graph, index_a)
        );

        Iterator< chemistry::AtomConformationalInterface> iterator_b( iterator_a);
        // iterate over distances and the remaining property values on the upper triangle of this matrix
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_distances( distances->Begin() + index_a), itr_distances_end( distances->End());
          itr_distances != itr_distances_end;
          ++itr_distances, ++iterator_b
        )
        {
          if( *itr_distances < m_NumberSteps)
          {
            STORAGE( *itr_distances) += prop_a * m_AtomProperty->operator ()( iterator_b)( 0);
          }
        }
      }
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Molecule2DACode::SetObjectHook()
    {
      m_Graph = graph::ConstGraph< size_t, size_t>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule2DACode::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes 2D (bond distance) autocorrelation of a specified atom property"
      );

      // Atom_Constant
      parameters.AddInitializer
      (
        "steps",
        "# of steps; corresponds to maximum # of bonds between atoms that will be considered for autocorrelation",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "11"
      );
      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the 2D-autocorrelation",
        io::Serialization::GetAgent( &m_AtomProperty),
        "Atom_Identity"
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
    bool Molecule2DACode::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
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

  } // namespace descriptor
} // namespace bcl
