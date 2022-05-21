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
#include "descriptor/bcl_descriptor_aa_pair_type.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAPairType::s_AsymmetricInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAPairType( false))
    );
    const util::SiPtr< const util::ObjectInterface> AAPairType::s_SymmetricInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAPairType( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAPairType::AAPairType( const bool &SYMMETRIC) :
      m_Stencil( linal::MakeVector< size_t>( 1)),
      m_Alignment( e_Center),
      m_Symmetric( SYMMETRIC)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AAPairType
    AAPairType *AAPairType::Clone() const
    {
      return new AAPairType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAPairType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAPairType::GetAlias() const
    {
      static const std::string s_name( "AAPairType"), s_symmetric_name( "SymmetricAAPairType");
      return m_Symmetric ? s_symmetric_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAPairType::GetNormalSizeOfFeatures() const
    {
      static const size_t s_asymmetric_size
      (
        size_t( biol::AATypes::s_NumberStandardAATypes)
        * size_t( biol::AATypes::s_NumberStandardAATypes)
      );
      static const size_t s_symmetric_size( ( s_asymmetric_size + size_t( biol::AATypes::s_NumberStandardAATypes)) / 2);
      return m_Symmetric ? s_symmetric_size : s_asymmetric_size;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    namespace
    {
      //! @brief add a type pair to the current vector
      //! @param STORAGE the vector to add the pair to
      //! @param TYPE_A the central aa type
      //! @param TYPE_B the distal aa type
      //! @param WEIGHT weight with which to add the types
      void AddTypePairToVector
      (
        linal::VectorReference< float> &STORAGE,
        const biol::AAType &TYPE_A,
        const biol::AAType &TYPE_B,
        const float &WEIGHT
      )
      {
        static const size_t n_aa_types( biol::AATypes::s_NumberStandardAATypes);

        size_t type_a_id
        (
          !TYPE_A.IsDefined() || TYPE_A->IsNaturalAminoAcid()
          ? TYPE_A.GetIndex()
          : TYPE_A->GetParentType() != TYPE_A ? TYPE_A->GetParentType()->GetIndex() : util::GetUndefined< size_t>()
        );
        size_t type_b_id
        (
          !TYPE_B.IsDefined() || TYPE_B->IsNaturalAminoAcid()
          ? TYPE_B.GetIndex()
          : TYPE_B->GetParentType() != TYPE_B ? TYPE_B->GetParentType()->GetIndex() : util::GetUndefined< size_t>()
        );

        // determine if the given storage vector is assuming an asymmetric distribution
        const bool is_asymmetric( n_aa_types * n_aa_types == STORAGE.GetSize());

        // test for symmetric and type a < type b
        if( !is_asymmetric && type_a_id < type_b_id)
        {
          std::swap( type_a_id, type_b_id);
        }

        // handle the common case that both types are defined
        if( util::IsDefined( type_a_id) && util::IsDefined( type_b_id))
        {
          BCL_Assert( type_a_id < n_aa_types, "Invalid type a id: " + util::Format()( type_a_id) + " for aa: " + TYPE_A->GetName());
          BCL_Assert( type_b_id < n_aa_types, "Invalid type b id: " + util::Format()( type_b_id) + " for aa: " + TYPE_B->GetName());
          if( is_asymmetric)
          {
            STORAGE( n_aa_types * type_a_id + type_b_id) += WEIGHT;
          }
          else
          {
            STORAGE( ( ( type_a_id + 1) * type_a_id) / 2 + type_b_id) += WEIGHT;
          }
        }
        else if( util::IsDefined( type_a_id))
        {
          BCL_Assert( type_a_id < n_aa_types, "Invalid type a id: " + util::Format()( type_a_id) + " for aa: " + TYPE_A->GetName());

          float *aa_type_begin( STORAGE.Begin() + n_aa_types * type_a_id);
          linal::VectorReference< float> aa_type_a( n_aa_types, aa_type_begin);
          aa_type_a += WEIGHT / float( n_aa_types);
        }
        else if( util::IsDefined( type_b_id))
        {
          BCL_Assert( type_b_id < n_aa_types, "Invalid type b id: " + util::Format()( type_b_id) + " for aa: " + TYPE_B->GetName());
          if( is_asymmetric)
          {
            for( size_t i( 0); i < n_aa_types; ++i)
            {
              STORAGE( i * n_aa_types + type_b_id) += WEIGHT / float( n_aa_types);
            }
          }
          else
          {
            for( size_t i( 0); i < n_aa_types; ++i)
            {
              STORAGE( i * ( i + 1) / 2 + type_b_id) += WEIGHT / float( n_aa_types);
            }
          }
        }
        else
        {
          // no defined aa types
          STORAGE += WEIGHT / float( STORAGE.GetSize());
        }
      }
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AAPairType::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // get the central type
      const biol::AAType &central_type( ELEMENT->GetType());
      if( m_Alignment != e_Right)
      {
        // if calculating left relationships
        // copy the iterator
        Iterator< biol::AABase> itr_element( ELEMENT);
        for
        (
          size_t stencil_id( 0), stencil_size( m_Stencil.GetSize());
          stencil_id < stencil_size;
          ++stencil_id
        )
        {
          const int target_pdb_id( ELEMENT->GetPdbID() - m_Stencil( stencil_id));
          while( itr_element.GetPosition() && itr_element( 0)->GetPdbID() > target_pdb_id)
          {
            --itr_element;
          }
          biol::AAType type_b;
          if( itr_element( 0)->GetPdbID() == target_pdb_id)
          {
            type_b = itr_element( 0)->GetType();
          }
          const float weight
          (
            m_StencilWeights.GetSize() == m_Stencil.GetSize() ? m_StencilWeights( stencil_id) : float( 1.0)
          );
          AddTypePairToVector( STORAGE, central_type, type_b, weight);
        }
      }
      if( m_Alignment != e_Left)
      {
        // if calculating right relationships
        // copy the iterator
        Iterator< biol::AABase> itr_element( ELEMENT);
        for
        (
          size_t stencil_id( 0), stencil_size( m_Stencil.GetSize());
          stencil_id < stencil_size;
          ++stencil_id
        )
        {
          const int target_pdb_id( ELEMENT->GetPdbID() + m_Stencil( stencil_id));
          while( itr_element.NotAtEnd() && itr_element( 0)->GetPdbID() < target_pdb_id)
          {
            ++itr_element;
          }
          biol::AAType type_b;
          if( itr_element.NotAtEnd() && itr_element( 0)->GetPdbID() == target_pdb_id)
          {
            type_b = itr_element( 0)->GetType();
          }
          const float weight
          (
            m_StencilWeights.GetSize() == m_Stencil.GetSize() ? m_StencilWeights( stencil_id) : float( 1.0)
          );
          AddTypePairToVector( STORAGE, central_type, type_b, weight);
        }
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPairType::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Returns a (normally binary) vector representing a 20x20 matrix (one value for each amino acid type pair)"
        "where M(x,y) = 1 iff x is the central amino acid in the stencil and y is the distal. "
        "If multiple amino acids are specified by the stencil, the average will sum to 1."
        "Undefined x would cause M(x,y) = 1/20 for all valid x"
      );

      parameters.AddInitializer
      (
        "alignment",
        "Alignment of the window.  Use center to consider the window around this element; left to consider the window"
        " up until this element, and right to consider the window following this element. "
        "JufoCenter is strictly for compatiblity with the old JUFO; it is slower to calculate than normal Center, "
        "but yields the same set of values in a different order",
        io::Serialization::GetAgent( &m_Alignment),
        "Center"
      );
      parameters.AddInitializer
      (
        "stencil",
        "AA relative distances to consider (not counting the central residue)",
        io::Serialization::GetAgentContainerWithCheck
        (
          &m_Stencil,
          io::Serialization::GetAgentWithRange( size_t( 1), size_t( 12))
        ),
        "(1,5)"
      );
      parameters.AddOptionalInitializer
      (
        "weight",
        "weight for each component in the stencil (1 by default).  If given, should be the same size as stencil",
        io::Serialization::GetAgent( &m_StencilWeights)
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
