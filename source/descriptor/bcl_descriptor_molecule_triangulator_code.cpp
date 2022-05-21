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
#include "descriptor/bcl_descriptor_molecule_triangulator_code.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeTriangulatorCode::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeTriangulatorCode()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeTriangulatorCode::MoleculeTriangulatorCode() :
      m_NumberSteps( 0),
      m_StepSize( 0.0),
      m_Temperature( 0.0),
      m_CutOff( 0.0)
    {
    }

    //! @brief constructor from number of steps, and mapped atom property
    MoleculeTriangulatorCode::MoleculeTriangulatorCode
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const float &TEMPERATURE,
      const float &CUTOFF
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_Temperature( TEMPERATURE),
      m_CutOff( CUTOFF)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeTriangulatorCode
    MoleculeTriangulatorCode *MoleculeTriangulatorCode::Clone() const
    {
      return new MoleculeTriangulatorCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeTriangulatorCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeTriangulatorCode::GetAlias() const
    {
      static const std::string s_name( "Triangulator");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t MoleculeTriangulatorCode::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float MoleculeTriangulatorCode::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get temperature of code
    //! @return const float  temperature of 3DA code
    const float &MoleculeTriangulatorCode::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &MoleculeTriangulatorCode::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > MoleculeTriangulatorCode::GetInternalDescriptors()
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
    void MoleculeTriangulatorCode::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // catch NUMBER_STEPS == 0
      if( m_NumberSteps == 0)
      {
        BCL_MessageCrt( "NUMBER_STEPS equals zero - triangulator code will be empty!");
        return;
      }
      else if( m_AtomProperty->GetNormalSizeOfFeatures() != 1)
      {
        BCL_MessageCrt
        (
          "Needed a property that returned 1 value per atom, instead got "
          + util::Format()( m_AtomProperty->GetNormalSizeOfFeatures())
          + " values per atom (atom property was "
          + util::Format()( m_AtomProperty->GetString())
          + ")"
        );
        return;
      }

      // instantiate the property as a vector with indices that correspond to atoms
      const linal::Vector< float> property( m_AtomProperty->CollectValuesOnEachElementOfObject( *this->GetCurrentObject()));

      // remember where we are in the property vector so that we can lookup the associated property value
      size_t count_a( 0);

      // LOOP 1 - itr_atoms_a
      // iterate over all atoms
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_a( this->GetCurrentObject()->GetIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a, ++count_a
      )
      {
        // remember where we are in the property vector so that we can lookup the associated property value
        size_t count_b( count_a);

        // LOOP 2 - itr_atoms_b
        // iterate over all atoms
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_b( itr_atoms_a);
          itr_atoms_b.NotAtEnd();
          ++itr_atoms_b, ++count_b
        )
        {
          // remember where we are in the property vector so that we can lookup the associated property value
          size_t count_c( count_b);

          // LOOP 3 - itr_atoms_c
          // iterate over all atoms
          for
          (
             iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_c( itr_atoms_b);
             itr_atoms_c.NotAtEnd();
             ++itr_atoms_c, ++count_c
          )
          {
            // generate each vector ab and ac
            linal::Vector3D a( itr_atoms_a->GetPosition());
            linal::Vector3D b( itr_atoms_b->GetPosition());
            linal::Vector3D c( itr_atoms_c->GetPosition());
            linal::Vector3D ab( a - b);
            linal::Vector3D ac( a - c);

            // store area of triangulation
            const float area( 0.5 * ( linal::CrossProduct( ab, ac).Norm()));

            //check if CUTOFF applies and the distance is larger than the CUTOFF
            if( m_CutOff && ( area > m_CutOff))
            {
              continue;
            }

            const float prop( property( count_a) * property( count_b) * property( count_c));

            // generate the actual Triangulator
            for( size_t step( 0); step < m_NumberSteps; ++step)
            {
              STORAGE( step) += prop * exp( -m_Temperature * math::Sqr( m_StepSize * step - area));
            }
          }
        }
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeTriangulatorCode::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes triangular autocorrelation of a specified atom property."
        "This is much like RDF, but considers all triplets of atoms."
      );

      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the triangular autocorrelation",
        io::Serialization::GetAgent( &m_AtomProperty)
      );
      parameters.AddInitializer
      (
        "step size",
        "size of each step in angstroms",
        io::Serialization::GetAgentWithRange( &m_StepSize, 0.01, 100.0),
        "1.0"
      );
      parameters.AddInitializer
      (
        "cutoff",
        "max area (A^2) to consider; 0.0 to consider all areas",
        io::Serialization::GetAgentWithRange( &m_CutOff, 0.0, 100.0),
        "0.0"
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
        "11"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
