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
#include "descriptor/bcl_descriptor_molecule_2da_max_min.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_min_max.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule2DAMaxMin::s_MaxInstance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule2DAMaxMin( e_2DAMax)
      )
    );

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule2DAMaxMin::s_MinInstance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule2DAMaxMin( e_2DAMin)
      )
    );

    //! @brief constructor, max vs min choice only
    Molecule2DAMaxMin::Molecule2DAMaxMin( const Method MAXMIN) :
      m_NumberSteps( 12),
      m_MaxOrMin( MAXMIN == e_2DAMax),
      m_SubValue( util::GetUndefined< float>())
    {
    }

    //! @brief constructor from number of steps, mapped atom property, and min/max choice
    Molecule2DAMaxMin::Molecule2DAMaxMin
    (
      const size_t NUMBER_STEPS,
      const CheminfoProperty &ATOM_PROPERTY,
      const Method MAXMIN,
      const float SUBVALUE
    ) :
      m_NumberSteps( NUMBER_STEPS),
      m_AtomProperty( ATOM_PROPERTY),
      m_MaxOrMin( MAXMIN == e_2DAMax),
      m_SubValue( SUBVALUE)
    {

      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule2DAMaxMin
    Molecule2DAMaxMin *Molecule2DAMaxMin::Clone() const
    {
      return new Molecule2DAMaxMin( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule2DAMaxMin::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule2DAMaxMin::GetAlias() const
    {
      static const std::string s_name_max( "2DAMax");
      static const std::string s_name_min( "2DAMin");
      return m_MaxOrMin ? s_name_max : s_name_min;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule2DAMaxMin::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule2DAMaxMin::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

    //! @brief get sub value of code
    //! @return sub value given to 2da code
    const float &Molecule2DAMaxMin::GetSubValue() const
    {
      return m_SubValue;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule2DAMaxMin::GetInternalDescriptors()
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
    void Molecule2DAMaxMin::Calculate( linal::VectorReference< float> &STORAGE)
    {
      float extreme;
      STORAGE = extreme = m_MaxOrMin ? math::GetLowestBoundedValue< float>() : math::GetHighestBoundedValue< float>();

      math::RunningMinMax< float> min_max_real;

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
        BCL_MessageCrt( "m_NumberSteps equals zero - 2DAMaxMinSign code will be empty!");
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
            const float product( prop_a * m_AtomProperty->operator ()( iterator_b)( 0));
            STORAGE( *itr_distances) =
                m_MaxOrMin ?
                    std::max( product, STORAGE( *itr_distances)) :
                    std::min( product, STORAGE( *itr_distances));
            min_max_real += product;
          }
        }
      }

      //Fix for "empty bin" values for molecules smaller than the max bin size
      std::replace
      (
        STORAGE.Begin(),
        STORAGE.End(),
        extreme,
        util::IsDefined( m_SubValue)         //DECIDE if using a sub value
        ? m_SubValue                       //USE the sub value to replace
        : (
            min_max_real.GetRange()         //DECIDE if min and max are different
            ? (
                m_MaxOrMin                  //DECIDE if doing a 2DAMax or 2DAMin
                ?   min_max_real.GetMin()    //USE the min to replace for 2DAMax
                :   min_max_real.GetMax()    //USE the max to replace for 2DAMin
              )
            : 0                             //USE zero to replace if range=0
          )
      );

    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Molecule2DAMaxMin::SetObjectHook()
    {
      m_Graph = graph::ConstGraph< size_t, size_t>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule2DAMaxMin::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_MaxOrMin ?
        "computes max for each bin of the 2DA of a specified atom property" :
        "computes min for each bin of the 2DA of a specified atom property"
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
      parameters.AddInitializer
      (
        "substitution_value",
        "value to replace empty bins unless given a nan, then uses worst found value",
        io::Serialization::GetAgent( &m_SubValue),
        "nan"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule2DAMaxMin::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_NumberSteps == 0)
      {
        ERR_STREAM << "m_NumberSteps equals zero - 2DAMaxMin code will be empty!";
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
