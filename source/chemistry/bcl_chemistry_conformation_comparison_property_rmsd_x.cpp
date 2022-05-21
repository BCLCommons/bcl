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
#include "chemistry/bcl_chemistry_conformation_comparison_property_rmsd_x.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonPropertyRMSDX::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonPropertyRMSDX())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ConformationComparisonPropertyRMSDX::ConformationComparisonPropertyRMSDX() :
      m_Properties
      (
        size_t( 1),
        descriptor::CheminfoProperty( descriptor::Constants< AtomConformationalInterface, float>( 1.0))
      ),
      m_PropertyWeight( 5.0),
      m_Base( 30.0)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        this->TryRead( util::ObjectDataLabel(), util::GetLogger());
      }
    }

    //! virtual copy constructor
    ConformationComparisonPropertyRMSDX *ConformationComparisonPropertyRMSDX::Clone() const
    {
      return new ConformationComparisonPropertyRMSDX( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonPropertyRMSDX::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonPropertyRMSDX::GetAlias() const
    {
      static const std::string s_Name( "PropertyRMSDX");
      return s_Name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return the RMSD between first and second molecule
    double ConformationComparisonPropertyRMSDX::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      if( MOLECULE_A.GetNumberAtoms() == size_t( 0) || MOLECULE_B.GetNumberAtoms() == size_t( 0))
      {
        return math::GetHighestBoundedValue< double>();
      }

      // use the property weighting
      const double sqrt_weight( math::Sqrt( m_PropertyWeight));

      const size_t n_atoms_a( MOLECULE_A.GetNumberAtoms());
      const size_t n_atoms_b( MOLECULE_B.GetNumberAtoms());

      // get a constant reference to the vector; necessary because m_Properties is mutable but we don't want it to
      // change in this function
      const storage::Vector< descriptor::CheminfoProperty> &properties( m_Properties);

      // compute the properties
      storage::Vector< linal::Vector< float> >
      properties_a( properties.GetSize()), properties_b( properties.GetSize());
      for( size_t i( 0), n_properties( m_Properties.GetSize()); i < n_properties; ++i)
      {
        properties_a( i) = properties( i)->CollectValuesOnEachElementOfObject( MOLECULE_A);
        properties_a( i) *= sqrt_weight;
        properties_b( i) = properties( i)->CollectValuesOnEachElementOfObject( MOLECULE_B);
        properties_b( i) *= sqrt_weight;

        // compute the distance matrix
        if( properties_a( i).GetSize() != n_atoms_a)
        {
          BCL_MessageCrt
          (
            "molecule with name " + MOLECULE_A.GetName() + " ignored as property could not be calculated"
          );
          return math::GetHighestBoundedValue< double>();
        }
        if( properties_b( i).GetSize() != n_atoms_b)
        {
          BCL_MessageCrt
          (
            "molecule with name " + MOLECULE_B.GetName() + " ignored as property could not be calculated"
          );
          return math::GetHighestBoundedValue< double>();
        }

      }

      // keep track of the smallest distances for each molecule
      linal::Vector< float> lowest_value_a( n_atoms_a, math::GetHighestBoundedValue< float>());
      linal::Vector< float> lowest_value_b( n_atoms_b, math::GetHighestBoundedValue< float>());

      iterate::Generic< const AtomConformationalInterface>
        itr_a( MOLECULE_A.GetAtomsIterator()), itr_b( MOLECULE_B.GetAtomsIterator());
      size_t atom_a_index( 0);
      const size_t n_properties( m_Properties.GetSize());
      for
      (
        float *itr_lowest_a( lowest_value_a.Begin());
        itr_a.NotAtEnd();
        ++itr_a, ++itr_lowest_a, ++atom_a_index
      )
      {
        float *itr_lowest_b( lowest_value_b.Begin());
        size_t atom_b_index( 0);
        for( itr_b.Restart(); itr_b.NotAtEnd(); ++itr_b, ++itr_lowest_b, ++atom_b_index)
        {
          // compute the distance between the atoms
          double sq_distance( linal::SquareDistance( itr_a->GetPosition(), itr_b->GetPosition()));

          // add in the property distance
          for( size_t i( 0); i < n_properties; ++i)
          {
            sq_distance += math::Sqr( properties_a( i)( atom_a_index) - properties_b( i)( atom_b_index));
          }

          // track the smallest distances
          if( *itr_lowest_a > sq_distance)
          {
            *itr_lowest_a = sq_distance;
          }
          if( *itr_lowest_b > sq_distance)
          {
            *itr_lowest_b = sq_distance;
          }
        }
      }

      // store all distances in a common vector and sort
      storage::Vector< float> all_smallest_distances_sq( lowest_value_a.Begin(), lowest_value_a.End());
      all_smallest_distances_sq.Append( storage::Vector< float>( lowest_value_b.Begin(), lowest_value_b.End()));
      all_smallest_distances_sq.Sort( std::less< float>());

      // compute the log multiplier such that if only one atom aligns, the denominator evalutes to 1/base
      const double log_weight( ( m_Base - 1) / ( m_Base * log( 2.0 * m_Base)));

      // now compute RMSD-Base, e.g. RMSD / 1 + ln(sqrt(N/base))
      // decomposing the denominator = 1 + 0.5 * ln(N/base) = 1 + ln(N) / 2 - ln(Base) / 2
      // precompute 1 - ln(Base) / 2
      const double base_offset( 1.0 - log( m_Base) * log_weight);

      double best_rmsd_x( math::GetHighestBoundedValue< double>());
      double sum_so_far( 0.0);
      size_t n_aligned( 0);
      for
      (
        size_t n_atoms_aligned( 1), total_atoms( all_smallest_distances_sq.GetSize());
        n_atoms_aligned <= total_atoms;
        ++n_atoms_aligned
      )
      {
        // update sum of squared deviations
        sum_so_far += all_smallest_distances_sq( n_atoms_aligned - 1);
        // compute RMSD squared so far
        const double rmsd_squared( sum_so_far / double( n_atoms_aligned));
        // compute RMSD-Xall_smallest_distances_sq
        const double rmsd_x( math::Sqrt( rmsd_squared) / ( base_offset + log( double( n_atoms_aligned) * 0.5) * log_weight));

        // update min value
        if( rmsd_x < best_rmsd_x && n_atoms_aligned >= 2 * std::min( lowest_value_a.GetSize(), lowest_value_b.GetSize()))
        {
          n_aligned = n_atoms_aligned;
          best_rmsd_x = rmsd_x;
        }
      }

      return best_rmsd_x;
    }

    //! @brief prepare the class for comparing a conformation
    //! @param MOLECULE the molecule to prepare to compare
    void ConformationComparisonPropertyRMSDX::Prepare( const ConformationInterface &MOLECULE) const
    {
      // calculate all properties on the molecule; this avoids potential clashes if multiple threads calculate properties
      // for the same molecule and try to update the cache simultaneously
      for
      (
        storage::Vector< descriptor::CheminfoProperty>::iterator
          itr( m_Properties.Begin()), itr_end( m_Properties.End());
        itr != itr_end;
        ++itr
      )
      {
        ( *itr)->CollectValuesOnEachElementOfObject( MOLECULE);
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonPropertyRMSDX::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the root-mean-squared-deviation between positions and properties between molecules using "
        "the formula sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 + (p1-p2)^2), where p refers to the property value. "
        "The RMSD is then divided by 1 + a * ln( n-atoms-aligned mol / base). "
        "a is set such that if only 1 atom aligns, the denominator evaluates to 1/base (e.g. a = (base-1)/(base*ln(base)))"
        "Alignment assignments are removed, in increasing order of RMSD, to optimize this function. "
        "The default property and property weight scheme corresponds to the ChargeRMSD metric reported in "
        "Gregory et al., ACS Chemical Neuroscience. PMID: 24528109"
      );

      parameters.AddInitializer
      (
        "base",
        "tunable parameter that should reflect the typical # of atoms that should be aligned",
        io::Serialization::GetAgentWithMin( &m_Base, 1.0),
        "30.0"
      );
      parameters.AddInitializer
      (
        "properties",
        "atom properties to consider, use multiply(Constant(X),property y) for weighting",
        io::Serialization::GetAgent( &m_Properties),
        util::ObjectDataLabel
        (
          "Combine(Multiply(Atom_SigmaCharge,Atom_VDWVolume))"
        )
      );
      parameters.AddInitializer
      (
        "property_weight",
        "Weighting to give the property relative to the euclidean distance (0-> no weighting)",
        io::Serialization::GetAgentWithMin( &m_PropertyWeight, 0.0),
        "5"
      );
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
