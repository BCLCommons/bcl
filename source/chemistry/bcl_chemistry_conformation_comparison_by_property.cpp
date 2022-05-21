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
#include "chemistry/bcl_chemistry_conformation_comparison_by_property.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_molecule_complete.h"
#include "descriptor/bcl_descriptor_named.h"
#include "descriptor/bcl_descriptor_sequence_size.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonByProperty::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonByProperty( false))
    );
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonByProperty::s_TanimotoInstance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonByProperty( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param TANIMOTO whether to use tanimoto distance rather than norm
    ConformationComparisonByProperty::ConformationComparisonByProperty( const bool &TANIMOTO) :
      m_Property
      (
        descriptor::Named< AtomConformationalInterface, float>
        (
          descriptor::SequenceSize< AtomConformationalInterface, float>(),
          "NAtoms",
          "Number of atoms"
        )
      ),
      m_Tanimoto( TANIMOTO),
      m_CacheLabel( m_Property.GetLabel())
    {
    }

    //! virtual copy constructor
    ConformationComparisonByProperty *ConformationComparisonByProperty::Clone() const
    {
      return new ConformationComparisonByProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonByProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonByProperty::GetAlias() const
    {
      static const std::string s_Name( "PropertyDistance"), s_TanimotoName( "PropertyTanimotoDistance");
      return m_Tanimoto ? s_TanimotoName : s_Name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return the RMSD between first and second molecule
    double ConformationComparisonByProperty::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      if( MOLECULE_A.GetNumberAtoms() == size_t( 0) || MOLECULE_B.GetNumberAtoms() == size_t( 0))
      {
        return math::GetHighestBoundedValue< double>();
      }
      util::SiPtr< const linal::Vector< float> > properties_a_si_ptr;
      util::SiPtr< const linal::Vector< float> > properties_b_si_ptr;

      // check that the properties existed, if not, cache them now
      properties_a_si_ptr = MOLECULE_A.FindInCache( m_CacheLabel);
      if( !properties_a_si_ptr.IsDefined())
      {
        Prepare( MOLECULE_A);
        properties_a_si_ptr = MOLECULE_A.FindInCache( m_CacheLabel);
      }
      properties_b_si_ptr = MOLECULE_B.FindInCache( m_CacheLabel);
      if( !properties_b_si_ptr.IsDefined())
      {
        Prepare( MOLECULE_B);
        properties_b_si_ptr = MOLECULE_B.FindInCache( m_CacheLabel);
      }

      // compute the properties
      const linal::Vector< float> &properties_a( *properties_a_si_ptr),
                                  &properties_b( *properties_b_si_ptr);
      // compute the distance matrix
      if( properties_a.GetSize() != m_Property->GetSizeOfFeatures() || !properties_a.IsDefined())
      {
        BCL_MessageCrt
        (
          "molecule with name " + MOLECULE_A.GetName() + " ignored as property could not be calculated"
        );
        return math::GetHighestBoundedValue< double>();
      }
      if( properties_b.GetSize() != m_Property->GetSizeOfFeatures() || !properties_b.IsDefined())
      {
        BCL_MessageCrt
        (
          "molecule with name " + MOLECULE_B.GetName() + " ignored as property could not be calculated"
        );
        return math::GetHighestBoundedValue< double>();
      }

      if( !m_Tanimoto)
      {
        // return euclidean distance
        return linal::Distance( properties_a, properties_b);
      }

      // return tanimoto distance.  For this to be a legitimate distance, no property value can be less than 0
      if( math::Statistics::MinimumValue( properties_a.Begin(), properties_a.End()) < float( 0.0))
      {
        return math::GetHighestBoundedValue< double>();
      }
      if( math::Statistics::MinimumValue( properties_b.Begin(), properties_b.End()) < float( 0.0))
      {
        return math::GetHighestBoundedValue< double>();
      }

      float sum( 0.0);
      for
      (
        const float *itr_a( properties_a.Begin()), *itr_a_end( properties_a.End()), *itr_b( properties_b.Begin());
        itr_a != itr_a_end;
        ++itr_a, ++itr_b
      )
      {
        const float ab_min( std::min( *itr_a, *itr_b));
        const float ab_max( std::max( std::max( *itr_a, *itr_b), float( 0.0000001)));
        sum += ab_min / ab_max;
      }
      sum /= float( properties_a.GetSize());
      return 1.0 - sum;
    }

    //! @brief prepare the class for comparing a conformation
    //! @param MOLECULE the molecule to prepare to compare
    void ConformationComparisonByProperty::Prepare( const ConformationInterface &MOLECULE) const
    {
      // calculate all properties on the molecule; this avoids potential clashes if multiple threads calculate properties
      // for the same molecule and try to update the cache simultaneously
      descriptor::CheminfoProperty property_copy( m_Property);
      property_copy->SetDimension( 0);
      linal::Vector< float> property;
      if( MOLECULE.GetNumberValences())
      {
        property = property_copy->SumOverObject( MoleculeComplete( MOLECULE));
      }
      else
      {
        property = property_copy->SumOverObject( MOLECULE);
      }
      MOLECULE.Cache( m_CacheLabel, property);
      if( m_Tanimoto && math::Statistics::MinimumValue( property.Begin(), property.End()) < float( 0.0))
      {
        BCL_MessageCrt
        (
          "Property returned value < 0, which is not allowed for tanimoto distance, on " + MOLECULE.GetName()
          + " values were: " + util::Format()( property)
        );
      }
    }

    //! @brief prepare the class for comparing conformations in the given ensemble
    //! @param ENSEMBLE the ensemble to prepare to compare
    void ConformationComparisonByProperty::PrepareEnsemble( const FragmentEnsemble &ENSEMBLE) const
    {
      // calculate all properties on the molecule; this avoids potential clashes if multiple threads calculate properties
      // for the same molecule and try to update the cache simultaneously
      descriptor::CheminfoProperty property_copy( m_Property);
      property_copy->SetDimension( 0);
      for( FragmentEnsemble::const_iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End()); itr != itr_end; ++itr)
      {
        linal::Vector< float> property( property_copy->SumOverObject( *itr));
        itr->Cache( m_CacheLabel, property);
        if( m_Tanimoto && math::Statistics::MinimumValue( property.Begin(), property.End()) < float( 0.0))
        {
          BCL_MessageCrt
          (
            "Property returned value < 0, which is not allowed for tanimoto distance, on " + itr->GetName()
            + " values were: " + util::Format()( property)
          );
        }
      }
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ConformationComparisonByProperty::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_Property->SetDimension( 0);
      //m_CacheLabel = m_Property->GetCacheLabel();
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonByProperty::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the "
        + std::string( m_Tanimoto ? "Tanimoto/Jaccard distance" : " root-mean-squared-deviation")
        + " between properties between molecules"
      );

      parameters.AddInitializer
      (
        "",
        "property to consider",
        io::Serialization::GetAgent( &m_Property),
        "NAtoms"
      );

      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
