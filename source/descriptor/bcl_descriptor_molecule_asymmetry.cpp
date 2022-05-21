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
#include "descriptor/bcl_descriptor_molecule_asymmetry.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeAsymmetry::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeAsymmetry()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeAsymmetry::MoleculeAsymmetry() :
      m_NumberSteps( 0),
      m_StepSize( 0.0),
      m_Temperature( 0.0),
      m_SumProps( false)
    {
    }

    //! @brief constructor from members.
    //! @param ATOM_PROPERTY Atom property for which stereochemistry rdf-like code will be scaled
    //! @param NUMBER_STEPS number of steps in the stereochemistry rdf-like code
    //! @param STEP_SIZE size of steps in the stereochemistry rdf-like code
    //! @param TEMPERATURE smoothing factor for the stereochemistry rdf-like code
    //! @param SUM_PROPS indicates whether to use normal property coefficient calculation (product of properties) or
    //                  alternate property coefficient calculation (sum of properties). These different methods may
    //                  be better or worse than the other depending on the dataset so both should be tried.
    MoleculeAsymmetry::MoleculeAsymmetry
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const float &TEMPERATURE,
      const bool &SUM_PROPS
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_Temperature( TEMPERATURE),
      m_SumProps( SUM_PROPS)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeAsymmetry
    MoleculeAsymmetry *MoleculeAsymmetry::Clone() const
    {
      return new MoleculeAsymmetry( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeAsymmetry::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeAsymmetry::GetAlias() const
    {
      static const std::string s_name( "MolecularAsymmetry");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t MoleculeAsymmetry::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float MoleculeAsymmetry::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get temperature of code
    //! @return const float  temperature of 3DA code
    const float &MoleculeAsymmetry::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &MoleculeAsymmetry::GetChemInfoProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > MoleculeAsymmetry::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomProperty,
          &m_AtomProperty + 1
        );
    }

    //! @brief Calculate the geometric center of the molecule without considering hydrogens.
    //! @param MOLECULE = molecule for which center will be calculated
    //! @return Vector3D of the coordinates of the geometric center of the molecule.
    linal::Vector3D MoleculeAsymmetry::GetCenterOfMoleculeWithoutHydrogens
    (
      const chemistry::ConformationInterface &MOLECULE
    ) const
    {
      math::RunningAverage< linal::Vector3D> averager;

      // iterate through the list of AtomsWithPosition and average their positions, skipping hydrogens
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_a( MOLECULE.GetAtomsIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a
      )
      {
        if( !( ( *itr_atoms_a).GetElementType() == chemistry::ElementTypes::GetEnums().e_Hydrogen))
        {
          averager += ( *itr_atoms_a).GetPosition();
        }
      }
      return averager.GetAverage();
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeAsymmetry::Calculate( linal::VectorReference< float> &STORAGE)
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);
      linal::Vector3D molecular_center( GetCenterOfMoleculeWithoutHydrogens( molecule));

      // Normalization factor for the maximum asymmetry of a triangle. For detailed explanation and derivation of this
      // const please see the confluence page.
      // https://structbio.vanderbilt.edu:8443/download/attachments/5113530/096225_proof.JPG?version=1&modificationDate=1336668531664
      const float normalization_factor( .096225);

      // This descriptor relies on defining planes using 3 atoms, therefore if the molecule has less than 3 non-H
      // atoms then it possesses no molecular asymmetry.
      if( molecule.GetNumberAtoms() - molecule.GetNumberHydrogens() < 3)
      {
        BCL_MessageDbg( "Not enough non-hydrogen atoms to calculate molecular asymmetry.");
        return;
      }
      m_AtomProperty->SetDimension( 1);

      // Get the atom properties of interest that will be used to weight the coordinates.
      linal::Vector< float> atom_properties( m_AtomProperty->CollectValuesOnEachElementOfObject( molecule));
      BCL_Assert( m_AtomProperty->GetNormalSizeOfFeatures() == 1, "Atom properties require 1 value per atom!");
      BCL_Assert( atom_properties.GetSize() == molecule.GetNumberAtoms(), "atom properties aint right");

      linal::Vector< float>::const_iterator itr_atom_a_prop( atom_properties.Begin());
      storage::Vector< linal::Vector3D> atom_coordinates( 3);
      linal::Vector3D current_side_lengths;
      // Iterate over all heavy-atom triplets.
      for
      (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_a( molecule.GetAtomsIterator());
          itr_atoms_a.NotAtEnd();
          ++itr_atoms_a, ++itr_atom_a_prop
      )
      {
        if( ( *itr_atoms_a).GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          continue;
        }
        atom_coordinates( 0) = ( *itr_atoms_a).GetPosition();
        linal::Vector< float>::const_iterator itr_atom_b_prop( atom_properties.Begin());
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_b( molecule.GetAtomsIterator());
          itr_atoms_b != itr_atoms_a;
          ++itr_atoms_b, ++itr_atom_b_prop
        )
        {
          if( ( *itr_atoms_b).GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
          {
            continue;
          }

          atom_coordinates( 1) = itr_atoms_b->GetPosition();
          linal::Vector< float>::const_iterator itr_atom_c_prop( atom_properties.Begin());

          for
          (
            iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_c( molecule.GetAtomsIterator());
            itr_atoms_c != itr_atoms_b;
            ++itr_atoms_c, ++itr_atom_c_prop
          )
          {
            if( itr_atoms_c->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
            {
              continue;
            }
            atom_coordinates( 2) = itr_atoms_c->GetPosition();

            // Calculate the normal vector of the plane generated by the atom triplet.
            linal::Vector3D plane_normal_vector
            (
              linal::CrossProduct
              (
                atom_coordinates( 1) - atom_coordinates( 0),
                atom_coordinates( 2) - atom_coordinates( 0)
              )
            );
            // Calculate the lengths of all sides of the triangle created by the atom triplet.
            // This is what provides a measure of the asymmetry of the triangle. Triangles with any sides equal
            // will have no asymmetry score and the order of the sides will reflect which direction the asymmetry
            // is pointing from the vantage point of the center of the molecule.
            current_side_lengths( 0) = linal::Distance( atom_coordinates( 0), atom_coordinates( 1));
            current_side_lengths( 1) = linal::Distance( atom_coordinates( 1), atom_coordinates( 2));
            current_side_lengths( 2) = linal::Distance( atom_coordinates( 2), atom_coordinates( 0));

            // Calculate the stereochemistry/volume score of the current atom triplet.
            // This score is based on the volume of the tetrahedron formed by the triplet and the molecular center.
            // Volume is calculated in such a way as it possesses a sign depending on whether the normal vector points
            // towards or away from the center of the molecule. This is important as it ensures that the numbered
            // order of the atoms does not affect the result. Atom triplets that lie further from the center of the
            // molecule will have a greater volume and in turn a greater score.
            const double current_volume
            (
              linal::ScalarProduct( plane_normal_vector, molecular_center - atom_coordinates( 0)) / 6.0
            );

            double current_stereochemistry( current_volume);

            // Side-lengths are factored into stereochemistry/volume score so that it can reflect mirror chirality and
            // asymmetry. Atom triplets that are more asymmetrically distributed will have a greater score. Stereochemistry
            // of the atom triplet is reflected in the fact that if their orders are reversed in space then the two
            // triangles "pointing" in opposite directions will have equal but opposite scores.
            current_stereochemistry *= ( current_side_lengths( 0) - current_side_lengths( 1)) *
              ( current_side_lengths( 1) - current_side_lengths( 2)) *
              ( current_side_lengths( 2) - current_side_lengths( 0));

            // Normalize by the maximum possible stereochemistry score given by the maximum side length. For derivation
            // of normalization_factor please see confluence page. In brief, the product of the difference between the
            // side lengths (a-b)(b-c)(c-a) will be between 0 and some maximum value that is achieved when the triangle
            // is most asymmetric. When the triangle is symmetric in that two sides are the same or very close to
            // equal, this product will approach zero. However, as the sides become more asymmetrical, this product
            // approaches its maximum value. This maximum value has been derived to be .096225. The mathematical
            // derivation of this value can be found on the confluence page and specifically in this image:
            // https://structbio.vanderbilt.edu:8443/download/attachments/5113530/096225_proof.JPG?version=1&modificationDate=1336668531664
            const double max_side
            (
              math::Statistics::MaximumValue( current_side_lengths.Begin(), current_side_lengths.End())
            );
            float current_stereo_cuberoot;
            if( max_side != 0.0)
            {
              current_stereochemistry /= ( normalization_factor * math::Pow( max_side, 3.0));
              if( current_stereochemistry == 0)
              {
                current_stereo_cuberoot = 0;
              }
              else
              {
                current_stereo_cuberoot = math::Pow( math::Absolute( current_stereochemistry), 1.0 / 3.0) *
                  current_stereochemistry / math::Absolute( current_stereochemistry);
              }
            }
            else
            {
              current_stereochemistry = 0;
              current_stereo_cuberoot = 0;
            }
//            BCL_Debug( current_stereo_cuberoot);

            // Calculate atom property weighting coefficient for rdf-like code.
            float property_coefficient;

            // Two methods considered for calculating the property coefficient are by product and by sum. Each has
            // benefits and drawbacks (product is sensitive to values of 0 for example while sum is sensitive to
            // equal and opposite property values). Depending on the dataset, one method may outperform the other and
            // both should be tested.
            if( m_SumProps)
            {
              property_coefficient = ( *itr_atom_a_prop) + ( *itr_atom_b_prop) + ( *itr_atom_c_prop);
            }
            else
            {
              property_coefficient = ( *itr_atom_a_prop) * ( *itr_atom_b_prop) * ( *itr_atom_c_prop);
            }

            // Calculate the actual rdf-like asymmetry code. This bins intensities of the cube-root of all scores
            // except 0. Binning method is similar to that of an rdf-code but steps cover the range of cube-root scores.
            // Codes weighted with the property calc_Identity will be coded entirely based on spatial distribution of
            // the molecule.
            // Final code represents an intensity "sum" of positive stereochemistry/volume scores minus negative
            // stereochemistry/volume scores. This way, simple enantiomers will represent mirror codes due to the fact
            // that equivalent scores will have opposite signs between the two enantiomers.
            if( current_stereo_cuberoot < 0)
            {
              property_coefficient *= -1;
            }
//            BCL_Debug( property_coefficient);

            // If the stereochemistry/asymmetry score is 0, then it should be skipped and not included in the rdf-like curve
            // because this curve should specifically reflect asymmetry and not be dictated by those parts of the
            // molecule that are symmetrical.
            if( math::EqualWithinTolerance( current_stereo_cuberoot, 0))
            {
              continue;
            }

            for( size_t temp_steps( 0); temp_steps < m_NumberSteps; ++temp_steps)
            {
              STORAGE( temp_steps) += property_coefficient * exp
                (
                  -m_Temperature * math::Sqr( ( m_StepSize * temp_steps) - math::Absolute( current_stereo_cuberoot))
                );
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
    io::Serializer MoleculeAsymmetry::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates an rdf-like asymmetry vector for the molecule based on overall distribution of atoms and/or properties. "
        "Use of this descriptor must be cited as: "
        "Sliwoski, Gregory, et al. \"BCL:: EMAS—Enantioselective Molecular Asymmetry Descriptor for 3D-QSAR.\" "
        "Molecules 17.8 (2012): 9971-9989.\n"
        "Link:  www.http://meilerlab.org/index.php/publications/show/2012\n"
      );

      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the molecular asymmetry score",
        io::Serialization::GetAgent( &m_AtomProperty)
      );
      parameters.AddInitializer
      (
        "step size",
        "size of each step in angstroms",
        io::Serialization::GetAgentWithRange( &m_StepSize, 0.01, 100.0),
        "0.1"
      );
      parameters.AddInitializer
      (
        "temperature",
        "increasing temperature spreads intensity across more distant bins",
        io::Serialization::GetAgentWithRange( &m_Temperature, 0.0, 1000.0),
        "100"
      );
      parameters.AddInitializer
      (
        "steps",
        "# of steps/bins (each of size = step size) used in the radial distribution function",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "24"
      );
      parameters.AddInitializer
      (
        "sum properties",
        "Use summation method to weight properties",
        io::Serialization::GetAgent( &m_SumProps)
      );
      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
