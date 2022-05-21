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

// include header of this class
#include "bcl_descriptor_outer_product.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    OuterProduct< t_DataType>::OuterProduct() :
      m_Descriptors( 2)
    {
    }

    //! @brief constructor from an operation
    template< typename t_DataType>
    OuterProduct< t_DataType>::OuterProduct
    (
      const util::Implementation< Base< t_DataType, float> > &LHS,
      const util::Implementation< Base< t_DataType, float> > &RHS
    ) :
      m_Descriptors( 2)
    {
      m_Descriptors( 0) = LHS;
      m_Descriptors( 1) = RHS;
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    OuterProduct< t_DataType> *OuterProduct< t_DataType>::Clone() const
    {
      return new OuterProduct( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &OuterProduct< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    template< typename t_DataType>
    const std::string &OuterProduct< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "OuterProduct");
      return s_name;
    }

    //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
    //! @note users do not normally need to override this function
    template< typename t_DataType>
    bool OuterProduct< t_DataType>::DimensionIsWellDefined() const
    {
      return m_Descriptors.GetSize() == size_t( 2) && m_Descriptors( 0).IsDefined() && m_Descriptors( 1).IsDefined()
             && ( m_Descriptors( 0)->DimensionIsWellDefined() || m_Descriptors( 1)->DimensionIsWellDefined());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer OuterProduct< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Outer product of two descriptors");

      // outer product
      parameters.AddInitializer
      (
        "lhs",
        "argument for the left hand side of the operation",
        io::Serialization::GetAgent( &m_Descriptors( 0))
      );
      parameters.AddInitializer
      (
        "rhs",
        "argument for the right hand side of the operation",
        io::Serialization::GetAgent( &m_Descriptors( 1))
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool OuterProduct< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // first, update type
      m_Type = m_Descriptors( 0)->GetType();

      // generalize the type to support all of the arguments to the operation
      for( iterate::Generic< Base< t_DataType, float> > itr( GetInternalDescriptors()); itr.NotAtEnd(); ++itr)
      {
        m_Type.GeneralizeToHandle( itr->GetType());
      }

      // set the dimension for all internally-held descriptors
      for( iterate::Generic< Base< t_DataType, float> > itr( GetInternalDescriptors()); itr.NotAtEnd(); ++itr)
      {
        itr->SetDimension( m_Type.GetDimension());
      }

      // check that both descriptors return the same # of values, or they can return a single value
      m_FeatureColumnsPerObject = m_Descriptors( 0)->GetSizeOfFeatures() * m_Descriptors( 1)->GetSizeOfFeatures();
      return true;
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > OuterProduct< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( m_Descriptors.Begin(), m_Descriptors.End());
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void OuterProduct< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // instantiate vector to store the result of the next property
      if( m_Descriptors.GetSize() < size_t( 2))
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }

      // initialize the descriptor return vector with the first descriptor
      linal::VectorConstReference< float> result_a( m_Descriptors( 0)->operator ()( ITR));
      linal::VectorConstReference< float> result_b( m_Descriptors( 1)->operator ()( ITR));

      const float *b_row_begin( result_b.Begin()), *b_row_end( result_b.End());
      float *itr_storage( STORAGE.Begin());
      for( const float *a_row( result_a.Begin()), *a_row_end( result_a.End()); a_row != a_row_end; ++a_row)
      {
        const float a_val( *a_row);
        for( const float *b_row( b_row_begin); b_row != b_row_end; ++b_row, ++itr_storage)
        {
          *itr_storage = a_val * *b_row;
        }
      }
    }

  } // namespace descriptor
} // namespace bcl
