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
#include "descriptor/bcl_descriptor_molecule_2da_max_sign.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule2DAMaxSign::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule2DAMaxSign()
      )
    );

    //! @brief default constructor
    Molecule2DAMaxSign::Molecule2DAMaxSign() :
      m_NumberSteps( 12)
    {
    }

    //! @brief constructor from number of steps, and mapped atom property
    Molecule2DAMaxSign::Molecule2DAMaxSign
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule2DAMaxSign
    Molecule2DAMaxSign *Molecule2DAMaxSign::Clone() const
    {
      return new Molecule2DAMaxSign( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule2DAMaxSign::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule2DAMaxSign::GetAlias() const
    {
      static const std::string s_name( "2DAMaxSign");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule2DAMaxSign::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule2DAMaxSign::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule2DAMaxSign::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomProperty,
          &m_AtomProperty + 1
        );
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void Molecule2DAMaxSign::Calculate
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
        BCL_MessageCrt( "Because m_NumberSteps equals zero, 2DAMaxSign code will be empty!");
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

        // determine the distance from this vertex to all vertices
        util::ShPtr< storage::Vector< size_t> > distances
        (
          graph::Connectivity::DistancesToOtherVertices( m_Graph, index_a)
        );

        //differentiate pos and neg atoms
        m_DiscreteCode( int( property_atom_a > 0)) =
            std::max( math::Sqr( property_atom_a), m_DiscreteCode( int( property_atom_a > 0)));

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
      STORAGE.CopyValues( m_DiscreteCode);
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Molecule2DAMaxSign::SetObjectHook()
    {
      m_Graph = graph::ConstGraph< size_t, size_t>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule2DAMaxSign::GetSerializer() const
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
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule2DAMaxSign::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_DiscreteCode = linal::Vector< float>( GetNormalSizeOfFeatures(), 0.0);

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

    //! @brief add an observed distance/property value to m_DiscreteCode
    //! @param DISTANCE actual distance of the two atoms
    //! @param PROP_A property from atom A
    //! @param PROP_B property from atom B
    void Molecule2DAMaxSign::Accumulate
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

      m_DiscreteCode( 3 * DISTANCE + sign_id) = std::max( product, m_DiscreteCode( 3 * DISTANCE + sign_id));
    }

  } // namespace descriptor
} // namespace bcl
