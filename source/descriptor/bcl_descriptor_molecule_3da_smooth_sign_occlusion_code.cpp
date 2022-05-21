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
#include "descriptor/bcl_descriptor_molecule_3da_smooth_sign_occlusion_code.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "descriptor/bcl_descriptor_molecule_3da_smooth.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule3DASmoothSignOcclusionCode::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule3DASmoothSignOcclusionCode()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Molecule3DASmoothSignOcclusionCode::Molecule3DASmoothSignOcclusionCode() :
      m_AtomProperty(),
      m_NumberSteps( 12),
      m_StepSize( 1.0),
      m_MaxOcclusions( 2),
      m_ClashDistance( 0.25)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

//! @brief constructor from number of steps, and mapped atom property
    Molecule3DASmoothSignOcclusionCode::Molecule3DASmoothSignOcclusionCode
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const bool &SQRT
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_MaxOcclusions( 2),
      m_ClashDistance( 0.25),
      m_Sqrt( SQRT)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule3DASmoothSignOcclusionCode
    Molecule3DASmoothSignOcclusionCode *Molecule3DASmoothSignOcclusionCode::Clone() const
    {
      return new Molecule3DASmoothSignOcclusionCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule3DASmoothSignOcclusionCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule3DASmoothSignOcclusionCode::GetAlias() const
    {
      static const std::string s_name( "3daClashSensitiveSign");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule3DASmoothSignOcclusionCode::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float Molecule3DASmoothSignOcclusionCode::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule3DASmoothSignOcclusionCode::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule3DASmoothSignOcclusionCode::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomProperty,
          &m_AtomProperty + 1
        );
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT_A, ELEMENT_B the elements of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DASmoothSignOcclusionCode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_A,
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_B,
      linal::VectorReference< float> &STORAGE
    )
    {
      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT_A), iterator_b( ELEMENT_B);
      float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

      // calculate exponential function for the RDF equation. 0.5 normalization is for backwards compatibility
      if( ELEMENT_A == ELEMENT_B)
      {
        STORAGE( int( property_atom_a > 0.0)) += 0.5 * property_atom_a * property_atom_a;
      }
      else
      {
        float property_atom_b( 0.5 * m_AtomProperty->operator ()( iterator_b)( 0));
        const float distance( linal::Distance( ELEMENT_A->GetPosition(), ELEMENT_B->GetPosition()));
        Accumulate( STORAGE, distance, property_atom_a, property_atom_b, CountClashes( *ELEMENT_A, *ELEMENT_B));
      }
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DASmoothSignOcclusionCode::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      Iterator< chemistry::AtomConformationalInterface> iterator_a( ELEMENT);
      float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

      // iterate over distances and the remaining property values on the upper triangle of this matrix
      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_b( ELEMENT.Begin());
        iterator_b.NotAtEnd();
        ++iterator_b
      )
      {
        const chemistry::AtomConformationalInterface &atom_b( *iterator_b( 0));
        //store distance between both atoms
        const float distance( linal::Distance( ELEMENT->GetPosition(), atom_b.GetPosition()));

        const float property_atom_b( m_AtomProperty->operator ()( iterator_b)( 0));

        // calculate exponential function for the RDF equation
        if( iterator_b( 0) == ELEMENT)
        {
          STORAGE( int( property_atom_a > 0.0)) += 0.5 * property_atom_b * property_atom_b;
        }
        else
        {
          Accumulate( STORAGE, distance, 0.5 * property_atom_a, property_atom_b, CountClashes( *ELEMENT, atom_b));
        }
      }
    } // Recalculate

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DASmoothSignOcclusionCode::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      for
      (
        Iterator< chemistry::AtomConformationalInterface> iterator_a( this->GetCurrentObject()->GetIterator());
        iterator_a.NotAtEnd();
        ++iterator_a
      )
      {
        const chemistry::AtomConformationalInterface &atom_a( *iterator_a( 0));
        const float property_atom_a( m_AtomProperty->operator ()( iterator_a)( 0));

        // skip 0 properties since some properties are very sparse
        if( !property_atom_a)
        {
          continue;
        }

        // differentiate between positive and negative atoms
        STORAGE( int( property_atom_a > 0)) += 0.5 * math::Sqr( property_atom_a);

        // iterate over distances and the remaining property values on the upper triangle of this matrix
        Iterator< chemistry::AtomConformationalInterface> iterator_b( iterator_a);
        for( ++iterator_b; iterator_b.NotAtEnd(); ++iterator_b)
        {
          const chemistry::AtomConformationalInterface &atom_b( *iterator_b( 0));
          //store distance between both atoms
          const float distance( linal::Distance( atom_a.GetPosition(), atom_b.GetPosition()));

          float property_atom_b( m_AtomProperty->operator ()( iterator_b)( 0));
          Accumulate( STORAGE, distance, property_atom_a, property_atom_b, CountClashes( atom_a, atom_b));
        }
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule3DASmoothSignOcclusionCode::GetSerializer() const
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
        "steps",
        "# of steps/bins (each of size = step size) used in the radial distribution function",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "48"
      );
      parameters.AddInitializer
      (
        "max clashes",
        "max # of intervening atoms between considered distance pairs",
        io::Serialization::GetAgentWithRange( &m_MaxOcclusions, 1, 1000000),
        "2"
      );
      parameters.AddInitializer
      (
        "clash distance",
        "distance between atom pair vector and another atom that will be considered a clash",
        io::Serialization::GetAgentWithRange( &m_ClashDistance, 0, 1),
        "0.15"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule3DASmoothSignOcclusionCode::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
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

    //! @brief add an observed distance/property value to m_DiscreteCode
    //! @param DISTANCE actual distance of the two atoms
    //! @param PROP_A property from atom A
    //! @param PROP_B property from atom B
    //! @param N_CLASHES number of clashes
    void Molecule3DASmoothSignOcclusionCode::Accumulate
    (
      linal::VectorReference< float> &STORAGE,
      const float &DISTANCE,
      const float &PROP_A,
      const float &PROP_B,
      const size_t &N_CLASHES
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
      else if( product == 0.0)
      {
        // product == 0 -> nothing to do
        return;
      }

      // What God intends:
      if( m_Sqrt)
      {
        product = math::Sqrt( product);
      }

      // comput the clash-related offset
      const size_t clash_offset( std::min( N_CLASHES, m_MaxOcclusions) * m_NumberSteps * 3);

      // assign atom distance to the closest 3DA bin
      size_t closest_step_low( std::min( size_t( m_NumberSteps - 1), size_t( DISTANCE / m_StepSize)));
      size_t closest_step_high( closest_step_low + 1);

      const float dist_lower_step( m_StepSize * closest_step_low - DISTANCE);
      const float dist_higher_step( dist_lower_step + m_StepSize);

      float low_exp_factor
      (
        std::max( dist_higher_step, float( 0.0)) / m_StepSize
      );
      float high_exp_factor
      (
        ( dist_higher_step >= float( 0.0) ? -dist_lower_step / m_StepSize : 0.0)
      );
      if( closest_step_high == m_NumberSteps)
      {
        // generate rough 3DA for boundry distances
        STORAGE( clash_offset + 3 * closest_step_low + sign_id) += product * low_exp_factor;
      }
      else
      {
        // normalize rough 3DA code
        product /= ( low_exp_factor + high_exp_factor);
        // generate rough 3DA code for all other distances
        STORAGE( clash_offset + 3 * closest_step_low + sign_id) += product * low_exp_factor;
        STORAGE( clash_offset + 3 * closest_step_high + sign_id) += product * high_exp_factor;
      }
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    //! In this case, updates the line segments
    void Molecule3DASmoothSignOcclusionCode::SetObjectHook()
    {
      m_BondLines.Reset();
      util::SiPtr< const chemistry::ConformationInterface> si_ptr_conf( this->GetCurrentObject());
      const chemistry::ConformationInterface &conformation( *si_ptr_conf);
      storage::Vector< graph::UndirectedEdge< size_t> > bonds
      (
        conformation.GetAdjacencyList( chemistry::ConfigurationalBondTypeData::e_BondOrder)
      );
      util::SiPtrVector< const linal::Vector3D> coords( conformation.GetAtomCoordinates());
      for
      (
        storage::Vector< graph::UndirectedEdge< size_t> >::const_iterator
          itr_bonds( bonds.Begin()), itr_bonds_end( bonds.End());
        itr_bonds != itr_bonds_end;
        ++itr_bonds
      )
      {
        m_BondLines.PushBack
        (
          coord::LineSegment3D
          (
            *coords( itr_bonds->GetIndexLow()),
            *coords( itr_bonds->GetIndexHigh())
          )
        );
      }
    }

    //! @brief determine the # of clashes between two atom pairs
    size_t Molecule3DASmoothSignOcclusionCode::CountClashes
    (
      const chemistry::AtomConformationalInterface &A,
      const chemistry::AtomConformationalInterface &B
    ) const
    {
      const coord::LineSegment3D line( A.GetPosition(), B.GetPosition());
      size_t n_clashes( 0);
      for
      (
        storage::Vector< coord::LineSegment3D>::const_iterator itr( m_BondLines.Begin()), itr_end( m_BondLines.End());
        itr != itr_end;
        ++itr
      )
      {
        if
        (
          A.GetPosition() != itr->GetStartPoint() && B.GetPosition() != itr->GetEndPoint() &&
          A.GetPosition() != itr->GetEndPoint() && B.GetPosition() != itr->GetStartPoint() &&
          coord::ShortestConnectionBetweenLineSegments3D( *itr, line).First().GetLength() < m_ClashDistance
        )
        {
          ++n_clashes;
        }
      }
      return n_clashes;
    }

  } // namespace descriptor
} // namespace bcl
