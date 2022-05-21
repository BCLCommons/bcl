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
#include "descriptor/bcl_descriptor_aa_triplet_type.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
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
    const util::SiPtr< const util::ObjectInterface> AATripletType::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATripletType( false))
    );
    const util::SiPtr< const util::ObjectInterface> AATripletType::s_SymmetricInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AATripletType( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AATripletType::AATripletType( const bool &SYMMETRIC) :
      m_Stencils(),
      m_Symmetric( SYMMETRIC)
    {
      m_Stencils.PushBack( linal::MakeVectorND< int>( -1, 1));
    }

    //! @brief Clone function
    //! @return pointer to new AATripletType
    AATripletType *AATripletType::Clone() const
    {
      return new AATripletType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AATripletType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AATripletType::GetAlias() const
    {
      static const std::string s_name( "AATripletType"), s_symmetric_name( "AASymmetricTripletType");
      return m_Symmetric ? s_symmetric_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AATripletType::GetNormalSizeOfFeatures() const
    {
      static const size_t s_asymmetric_size
      (
        size_t( biol::AATypes::s_NumberStandardAATypes)
        * size_t( biol::AATypes::s_NumberStandardAATypes)
        * size_t( biol::AATypes::s_NumberStandardAATypes)
      );
      static const size_t s_symmetric_size
      (
        size_t( biol::AATypes::s_NumberStandardAATypes)
        * size_t( biol::AATypes::s_NumberStandardAATypes)
        * ( size_t( biol::AATypes::s_NumberStandardAATypes) + 1) / 2
      );
      return m_Symmetric ? s_symmetric_size : s_asymmetric_size;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AATripletType::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      static const size_t n_aa_types( biol::AATypes::s_NumberStandardAATypes);
      static const size_t n_aa_types_sq( n_aa_types * n_aa_types);

      // get the central type
      const biol::AAType central_type( ELEMENT->GetType());
      if( central_type.GetIndex() >= n_aa_types)
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }
      const size_t central_offset( central_type.GetIndex());

      for
      (
        storage::Vector< linal::VectorND< int, 2> >::const_iterator
          itr_stencil( m_Stencils.Begin()), itr_stencil_end( m_Stencils.End());
        itr_stencil != itr_stencil_end;
        ++itr_stencil
      )
      {
        biol::AAType lhs_type, rhs_type;
        // if calculating left relationships
        // copy the iterator
        Iterator< biol::AABase> itr_element_lhs( ELEMENT), itr_element_rhs( ELEMENT);
        const int lhs_offset( itr_stencil->operator()( 0));
        const int rhs_offset( itr_stencil->operator()( 1));
        const int target_pdb_id_lhs( ELEMENT->GetPdbID() + lhs_offset);
        const int target_pdb_id_rhs( ELEMENT->GetPdbID() + rhs_offset);
        if( lhs_offset < 0)
        {
          while( itr_element_lhs.GetPosition() && itr_element_lhs( 0)->GetPdbID() > target_pdb_id_lhs)
          {
            --itr_element_lhs;
          }
        }
        else
        {
          while( itr_element_lhs.NotAtEnd() && itr_element_lhs( 0)->GetPdbID() < target_pdb_id_lhs)
          {
            ++itr_element_lhs;
          }
          if( !itr_element_lhs.NotAtEnd())
          {
            --itr_element_lhs;
          }
        }
        if( itr_element_lhs( 0)->GetPdbID() == target_pdb_id_lhs)
        {
          lhs_type = itr_element_lhs( 0)->GetType();
        }
        else
        {
          STORAGE = util::GetUndefined< float>();
          return;
        }
        if( rhs_offset < 0)
        {
          while( itr_element_rhs.GetPosition() && itr_element_rhs( 0)->GetPdbID() > target_pdb_id_rhs)
          {
            --itr_element_rhs;
          }
        }
        else
        {
          while( itr_element_rhs.NotAtEnd() && itr_element_rhs( 0)->GetPdbID() < target_pdb_id_rhs)
          {
            ++itr_element_rhs;
          }
          if( !itr_element_rhs.NotAtEnd())
          {
            --itr_element_rhs;
          }
        }
        if( itr_element_rhs( 0)->GetPdbID() == target_pdb_id_rhs)
        {
          rhs_type = itr_element_rhs( 0)->GetType();
        }
        else
        {
          STORAGE = util::GetUndefined< float>();
          return;
        }

        if( lhs_type.GetIndex() >= n_aa_types || rhs_type.GetIndex() >= n_aa_types)
        {
          STORAGE = util::GetUndefined< float>();
          return;
        }

        if( !m_Symmetric)
        {
          STORAGE( central_offset + lhs_type.GetIndex() * n_aa_types + rhs_type.GetIndex() * n_aa_types_sq) += 1.0;
        }
        else
        {
          const size_t lhs_type_id( std::min( lhs_type.GetIndex(), rhs_type.GetIndex()));
          const size_t rhs_type_id( std::min( lhs_type.GetIndex(), rhs_type.GetIndex()));
          STORAGE( central_offset + n_aa_types * ( lhs_type_id * ( lhs_type_id + 1) / 2 + rhs_type_id)) += 1.0;
        }
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AATripletType::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Returns a (normally binary) vector representing a 20x20x20 tensor (one value for each amino acid type triplet)"
        "where M(x,y,z) = 1 iff the sequence surrounding the current amino acid is XYZ"
      );
      parameters.AddInitializer
      (
        "",
        "Stencils to use; relative positions used to return the AA type.  The central AA's position (0) is always included",
        io::Serialization::GetAgentWithSizeLimits( &m_Stencils, size_t( 0), size_t( 10)),
        "((-1,1))"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool AATripletType::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      if( m_Stencils.IsEmpty())
      {
        m_Stencils.PushBack( linal::MakeVectorND< int>( -1, 1));
      }
      return true;
    }

  } // namespace descriptor
} // namespace bcl
