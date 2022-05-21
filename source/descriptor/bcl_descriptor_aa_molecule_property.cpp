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
#include "descriptor/bcl_descriptor_aa_molecule_property.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"

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
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAMoleculeProperty( true, true));
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAMoleculeProperty( false, true));
        return util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAMoleculeProperty( true, false));
      }
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAMoleculeProperty::s_Instance( AddInstances());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from exposure measure
    //! @param CONSIDER_SIDE_CHAIN whether to consider side chain atoms in the calculation
    //! @param CONSIDER_BACK_BONE whether to consider back bone atoms in the calculation
    //! @param PROPERTY property to calculate
    AAMoleculeProperty::AAMoleculeProperty
    (
      const bool &CONSIDER_SIDE_CHAIN,
      const bool &CONSIDER_BACK_BONE,
      const CheminfoProperty &PROPERTY
    ) :
      m_Property( PROPERTY),
      m_ConsiderSideChain( CONSIDER_SIDE_CHAIN),
      m_ConsiderBackBone( CONSIDER_BACK_BONE),
      m_IsCachedOnAA( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AAMoleculeProperty
    AAMoleculeProperty *AAMoleculeProperty::Clone() const
    {
      return new AAMoleculeProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAMoleculeProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAMoleculeProperty::GetAlias() const
    {
      static const std::string s_bb_name( "BackBoneChemDescriptor"),
                               s_sc_name( "SideChainChemDescriptor"),
                               s_name( "ChemDescriptor");
      return m_ConsiderBackBone ? ( m_ConsiderSideChain ? s_name : s_bb_name) : s_sc_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAMoleculeProperty::GetNormalSizeOfFeatures() const
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
    void AAMoleculeProperty::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // handle the case where the property has already been cached
      if( m_IsCachedOnAA)
      {
        if( m_ConsiderBackBone)
        {
          // determine whether the element is the first or last in the chain
          bool is_n_terminus( false), is_c_terminus( false);
          // check for c terminus
          {
            iterate::Generic< const biol::AABase> element_next( ELEMENT);
            ++element_next;
            if( !element_next.NotAtEnd() || element_next->GetChainID() != ELEMENT->GetChainID())
            {
              is_c_terminus = true;
            }
          }
          // check for n-terminus
          if( ELEMENT.GetPosition() != size_t( 0))
          {
            iterate::Generic< const biol::AABase> element_prev( ELEMENT);
            --element_prev;
            if( element_prev->GetChainID() != ELEMENT->GetChainID())
            {
              is_n_terminus = true;
            }
          }
          else
          {
            is_n_terminus = true;
          }
          STORAGE.CopyValues
          (
            ELEMENT->GetType()->GetFragment( is_c_terminus, is_n_terminus).GetFromCache( m_PropertyString)
          );
        }
        else
        {
          // side chain only, no need to worry about the terminus
          STORAGE.CopyValues
          (
            ELEMENT->GetType()->GetFragment( false, false).GetFromCache( m_PropertyString)
          );
        }
        return;
      }

      if( m_AAMoleculePropertyStorage.IsEmpty())
      {
        // not yet cached; calculate on the protein
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
        m_AAMoleculePropertyStorage.Find( *ELEMENT)
      );

      // check if the sasa for this aa exists
      if( itr == m_AAMoleculePropertyStorage.End())
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
    io::Serializer AAMoleculeProperty::GetSerializer() const
    {
      std::string atom_types( m_ConsiderBackBone ? ( m_ConsiderSideChain ? "all" : "back bone") : "side chain");
      io::Serializer serializer;
      serializer.SetClassDescription( "Uses a chemical descriptor on " + atom_types + " atoms of each AA");

      serializer.AddInitializer
      (
        "",
        "molecular property to use on each residue individually",
        io::Serialization::GetAgent( &m_Property)
      );

      return serializer;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AAMoleculeProperty::SetObjectHook()
    {
      // reset the map to store the AA's for the new protein model
      m_AAMoleculePropertyStorage.Reset();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool AAMoleculeProperty::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_PropertyString = m_Property.GetLabel();

      // first, determine whether the property is already cached on the amino acids directly; in this case, the property
      // is conformation independent and can be retrieved directly, so no need for generating the chemical representation
      // which is a very slow process
      if( biol::GetAATypes().ALA->GetFragment( false, false).IsCached( m_PropertyString))
      {
        BCL_MessageVrb( "Property " + m_PropertyString.ToString() + " is cached");
        m_IsCachedOnAA = true;
      }
      else
      {
        BCL_MessageVrb( "Property " + m_PropertyString.ToString() + " is NOT cached");
        m_IsCachedOnAA = false;
      }

      // for thread-safety, ensure that all aa molecules have already been calculated
      for
      (
        biol::AATypes::const_iterator itr( biol::GetAATypes().Begin()), itr_end( biol::GetAATypes().End());
        itr != itr_end;
        ++itr
      )
      {
        ( *itr)->GetFragment( false, false).IsCached( m_PropertyString);
      }

      // add the appropriate prefix and/or suffix to the descriptor, depending on what part of the AA it applies to
      if( m_ConsiderBackBone && !m_ConsiderSideChain)
      {
        m_PropertyString = "BackBone" + m_PropertyString.ToString();
      }
      else if( m_ConsiderSideChain && !m_ConsiderBackBone)
      {
        m_PropertyString = "SideChain" + m_PropertyString.ToString();
      }
      return true;
    }

    //! @brief populate the map with the given molecular property
    void AAMoleculeProperty::CalculatePropertyOnProtein()
    {
      // get a pointer to the protein
      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

      // create a chemical representation of the protein model
      const chemistry::AAFragmentComplete &protein_molecule( sp_protein_model->GetChemicalRepresentation());

      // get fragments, one per residue
      chemistry::FragmentEnsemble fragment_ensemble
      (
        protein_molecule.GetResiduesAsFragments( m_ConsiderBackBone, m_ConsiderSideChain)
      );

      // iterate through the ensemble, calculate the associated descriptor, and save it in the map
      util::SiPtrVector< const biol::AABase>::const_iterator itr_res( protein_molecule.GetResidueSequence().Begin());
      for
      (
        chemistry::FragmentEnsemble::const_iterator itr( fragment_ensemble.Begin()), itr_end( fragment_ensemble.End());
        itr != itr_end;
        ++itr, ++itr_res
      )
      {
        m_AAMoleculePropertyStorage[ *itr_res] = m_Property->SumOverObject( *itr);
      }
    }

  } // namespace descriptor
} // namespace bcl
