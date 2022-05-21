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
#include "descriptor/bcl_descriptor_aa_atom_property.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////
  // data //
  //////////

    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAAtomProperty( true, true, true));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAAtomProperty( false, true, true));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAAtomProperty( true, false, true));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAAtomProperty( true, true, false));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAAtomProperty( false, true, false));
        return util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAAtomProperty( true, false, false));
      }
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAAtomProperty::s_Instance( AddInstances());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from exposure measure
    //! @param CONSIDER_SIDE_CHAIN whether to consider side chain atoms in the calculation
    //! @param CONSIDER_BACK_BONE whether to consider back bone atoms in the calculation
    //! @param PROPERTY property to calculate
    AAAtomProperty::AAAtomProperty
    (
      const bool &CONSIDER_SIDE_CHAIN,
      const bool &CONSIDER_BACK_BONE,
      const bool &COMPUTE_MEAN,
      const CheminfoProperty &PROPERTY
    ) :
      m_Property( PROPERTY),
      m_ConsiderSideChain( CONSIDER_SIDE_CHAIN),
      m_ConsiderBackBone( CONSIDER_BACK_BONE),
      m_ComputeMean( COMPUTE_MEAN)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AAAtomProperty
    AAAtomProperty *AAAtomProperty::Clone() const
    {
      return new AAAtomProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAAtomProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAAtomProperty::GetAlias() const
    {
      static const std::string s_bb_name( "SumBackBoneAtomProperty"),
                               s_sc_name( "SumSideChainAtomProperty"),
                               s_name( "SumAtomProperty"),
                               s_mean_bb_name( "MeanBackBoneAtomProperty"),
                               s_mean_sc_name( "MeanSideChainAtomProperty"),
                               s_mean_name( "MeanAtomProperty");
      return
        m_ConsiderBackBone
        ? (
            m_ConsiderSideChain
            ? ( m_ComputeMean ? s_mean_name : s_name)
            : ( m_ComputeMean ? s_mean_bb_name : s_bb_name)
          )
        : ( m_ComputeMean ? s_mean_sc_name : s_sc_name);
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAAtomProperty::GetNormalSizeOfFeatures() const
    {
      return m_Property->GetNormalSizeOfFeatures();
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AAAtomProperty::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      if( m_AAAtomPropertyStorage.IsEmpty())
      {
        CalculatePropertyOnProtein();
      }

      // find the entry for given amino acid
      const storage::Map
      <
        util::SiPtr< const biol::AABase>,
        linal::Vector< float>,
        biol::AALessThanSeqID
      >::const_iterator itr
      (
        m_AAAtomPropertyStorage.Find( *ELEMENT)
      );

      // check if the sasa for this aa exists
      if( itr == m_AAAtomPropertyStorage.End())
      {
        // set storage to undefined
        STORAGE = util::GetUndefined< float>();
      }
      else
      {
        // copy the descriptor value
        STORAGE.CopyValues( itr->second);
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAAtomProperty::GetSerializer() const
    {
      std::string atom_types( m_ConsiderBackBone ? ( m_ConsiderSideChain ? "all" : "back bone") : "side chain");
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Takes the " + std::string( m_ComputeMean ? "mean" : "sum")
        + " of a chemical descriptor on " + atom_types + " atoms of each AA"
      );

      serializer.AddInitializer
      (
        "",
        "atom property to use on the complete protein-molecule",
        io::Serialization::GetAgent( &m_Property)
      );

      return serializer;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AAAtomProperty::SetObjectHook()
    {
      // reset the map to store the AA's for the new protein model
      m_AAAtomPropertyStorage.Reset();
    }

    //! @brief populate the map with the given molecular property
    void AAAtomProperty::CalculatePropertyOnProtein()
    {
      // get a pointer to the protein
      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());
      const chemistry::AAFragmentComplete &protein_molecule( sp_protein_model->GetChemicalRepresentation());

      // compute the values over the full protein.  These would be different if calculated only over each residue
      // individually
      linal::Vector< float> atom_property_values( m_Property->CollectValuesOnEachElementOfObject( protein_molecule));

      // get # values per atom
      const size_t n_values_per_atom( m_Property->GetNormalSizeOfFeatures());

      // determine which residue each atom belongs to

      // property value of all atoms belonging to each residue
      storage::Vector< math::RunningAverage< linal::Vector< float> > > residue_atom_values( protein_molecule.GetResidueSequence().GetSize());

      // determine indices of atoms for each residue
      size_t atom_index( 0);

      // keep a pointer to the start of the atom property
      const float *itr_values( atom_property_values.Begin());
      for
      (
        storage::Vector< size_t>::const_iterator
          itr( protein_molecule.GetAtomsResidueIndices().Begin()),
          itr_end( protein_molecule.GetAtomsResidueIndices().End());
        itr != itr_end;
        ++itr, ++atom_index, itr_values += n_values_per_atom
      )
      {
        if( protein_molecule.GetAtomsType( atom_index)->IsSideChain() == m_ConsiderSideChain || m_ConsiderBackBone)
        {
          residue_atom_values( *itr) += linal::VectorConstReference< float>( n_values_per_atom, itr_values);
        }
      }

      // iterate through the ensemble, calculate the associated descriptor, and save it in the map
      util::SiPtrVector< const biol::AABase>::const_iterator itr_res( protein_molecule.GetResidueSequence().Begin());
      for
      (
        storage::Vector< math::RunningAverage< linal::Vector< float> > >::const_iterator
          itr( residue_atom_values.Begin()), itr_end( residue_atom_values.End());
        itr != itr_end;
        ++itr, ++itr_res
      )
      {
        if( itr->GetAverage().GetSize() > size_t( 0))
        {
          m_AAAtomPropertyStorage[ *itr_res] = itr->GetAverage() * float( m_ComputeMean ? 1.0 : itr->GetWeight());
        }
        else
        {
          m_AAAtomPropertyStorage[ *itr_res] = linal::Vector< float>( n_values_per_atom, float( 0));
        }
        BCL_Assert( m_AAAtomPropertyStorage[ *itr_res].GetSize(), "Should not be empty");
      }
    }
  } // namespace descriptor
} // namespace bcl
