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
#include "bcl_descriptor_binary_operation.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from an operation
    template< typename t_DataType>
    BinaryOperation< t_DataType>::BinaryOperation
    (
      const util::Implementation< math::AssignmentOperationInterface< float> > &OPERATION
    ) :
      m_Descriptors( 2),
      m_Op( OPERATION)
    {
    }

    //! @brief constructor from an operation
    template< typename t_DataType>
    BinaryOperation< t_DataType>::BinaryOperation
    (
      const util::Implementation< Base< t_DataType, float> > &LHS,
      const util::Implementation< math::AssignmentOperationInterface< float> > &OPERATION,
      const util::Implementation< Base< t_DataType, float> > &RHS
    ) :
      m_Descriptors( 2),
      m_Op( OPERATION)
    {
      m_Descriptors( 0) = LHS;
      m_Descriptors( 1) = RHS;
      // call read initializer success hook to setup the type
      std::stringstream st;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), st);
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    BinaryOperation< t_DataType> *BinaryOperation< t_DataType>::Clone() const
    {
      return new BinaryOperation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &BinaryOperation< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    template< typename t_DataType>
    const std::string &BinaryOperation< t_DataType>::GetAlias() const
    {
      return m_Op->GetVerb();
    }

    //! @brief Return the assignment operation for this descriptor
    //! @return the assignment operation for this descriptor
    template< typename t_DataType>
    const util::Implementation< math::AssignmentOperationInterface< float> > &
      BinaryOperation< t_DataType>::GetAssignmentOperation() const
    {
      return m_Op;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer BinaryOperation< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( m_Op->GetVerb() + "s two descriptors (binary true/false)");

      // determine whether the lhs and rhs are interchangeable
      if( m_Op->IsSymmetric())
      {
        // comparison ops should only take two objects at a time
        const size_t max_size( m_Op->GetVerb() != "Equal" && m_Op->GetVerb() != "NotEqual" ? 1000 : 2);
        // Multiply, Add, And, Or
        parameters.AddInitializer
        (
          "",
          "Descriptors to " + m_Op->GetVerb(),
          io::Serialization::GetAgentWithSizeLimits( &m_Descriptors, size_t( 2), max_size) // Require at least two descriptors
        );
      }
      else
      {
        // Subtract, divide, exponentiate
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
      }

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool BinaryOperation< t_DataType>::ReadInitializerSuccessHook
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
      m_FeatureColumnsPerObject = m_Descriptors( 0)->GetSizeOfFeatures();
      for( iterate::Generic< Base< t_DataType, float> > itr( GetInternalDescriptors()); itr.NotAtEnd(); ++itr)
      {
        // get number of values returned for this descriptor
        const size_t itr_features_size( itr->GetSizeOfFeatures());
        if( itr_features_size == m_FeatureColumnsPerObject || itr_features_size == size_t( 1))
        {
          continue;
        }
        if( m_FeatureColumnsPerObject != size_t( 1))
        {
          ERR_STREAM
            << "Cannot perform operation on descriptors because descriptors return different #s of features: "
            << m_Descriptors( 0).GetString() << " returns " << m_FeatureColumnsPerObject
            << " values, " << itr->GetString() << " returns "
            << itr_features_size;
          return false;
        }
        m_FeatureColumnsPerObject = itr_features_size;
      }
      return true;
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > BinaryOperation< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( m_Descriptors.Begin(), m_Descriptors.End());
    }

    //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
    //! @note users do not normally need to override this function
    template< typename t_DataType>
    bool BinaryOperation< t_DataType>::DimensionIsWellDefined() const
    {
      if( m_Descriptors.GetSize() >= size_t( 2))
      {
        for
        (
          typename storage::Vector< util::Implementation< Base< t_DataType, float> > >::const_iterator
            itr( m_Descriptors.Begin()), itr_end( m_Descriptors.End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->IsDefined() && ( *itr)->DimensionIsWellDefined())
          {
            return true;
          }
        }
      }
      return false;
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void BinaryOperation< t_DataType>::RecalculateImpl
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
      linal::VectorConstReference< float> result( m_Descriptors( 0)->operator ()( ITR));

      // check for vectors of incorrect size -> indicates descriptors that were not calculated
      if( result.GetSize() != m_Descriptors( 0)->GetSizeOfFeatures())
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }

      // copy the values into storage.  Since any argument may have just one element, handle that case specially
      if( result.GetSize() == size_t( 1))
      {
        STORAGE = result( 0);
      }
      else
      {
        STORAGE.CopyValues( result);
      }

      if( m_Op->HasZeroProperty() && STORAGE == float( 0.0))
      {
        // avoid evaluating right hand side
        return;
      }

      // get a reference to the operation
      const math::AssignmentOperationInterface< float> &operation( *m_Op);
      iterate::Generic< Base< t_DataType, float> > itr( GetInternalDescriptors());
      for( ++itr; itr.NotAtEnd(); ++itr)
      {
        // calculate the next property
        linal::VectorConstReference< float> rhs( itr->operator()( ITR));

        // check for vectors of incorrect size -> indicates descriptors that were not calculated
        if( rhs.GetSize() != itr->GetSizeOfFeatures())
        {
          STORAGE = util::GetUndefined< float>();
          return;
        }

        const float *itr_rhs( rhs.Begin());
        if( STORAGE.GetSize() == rhs.GetSize())
        {
          for // perform the operation on each element of the vector
          (
            linal::Vector< float>::iterator itr_result( STORAGE.Begin()), itr_result_end( STORAGE.End());
            itr_result != itr_result_end;
            ++itr_rhs, ++itr_result
          )
          {
            operation( *itr_result, *itr_rhs);
          }
        }
        else // if( m_Descriptors( 1)->GetSizeOfFeatures() == size_t( 1))
        {
          BCL_Assert
          (
            rhs.GetSize() == size_t( 1),
            "expected that rhs was the same size as the lhs descriptor, or that one of them was of size == 1"
          );
          const float rhs_value( rhs( 0));
          for // constant rhs; perform the operation by walking over the LHS
          (
            linal::Vector< float>::iterator itr_result( STORAGE.Begin()), itr_result_end( STORAGE.End());
            itr_result != itr_result_end;
            ++itr_result
          )
          {
            operation( *itr_result, rhs_value);
          }
        }
      }
    }

  } // namespace descriptor
} // namespace bcl
