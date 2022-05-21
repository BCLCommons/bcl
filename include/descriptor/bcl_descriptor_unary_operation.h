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

#ifndef BCL_DESCRIPTOR_UNARY_OPERATION_H_
#define BCL_DESCRIPTOR_UNARY_OPERATION_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "math/bcl_math_assignment_unary_interface.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UnaryOperation
    //! @brief perform a math::Assignment on two descriptors
    //!
    //! @see @link example_descriptor_unary_operation.cpp @endlink
    //! @author mendenjl
    //! @date Mar 11, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class UnaryOperation :
      public Base< t_DataType, float>
    {
    private:

    //////////
    // data //
    //////////

      //! descriptor to perform the operation on
      util::Implementation< Base< t_DataType, float> > m_Descriptor;

      //! the operation that is performed to merge the descriptors (e.g. Log, Exp, Abs)
      util::Implementation< math::AssignmentUnaryInterface> m_Op;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from an operation
      explicit UnaryOperation
      (
        const util::Implementation< math::AssignmentUnaryInterface> &OPERATION,
        const util::Implementation< Base< t_DataType, float> > &DESCRIPTOR =
          ( util::Implementation< Base< t_DataType, float> >())
      );

      //! @brief virtual copy constructor
      UnaryOperation< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief Return the assignment operation for this descriptor
      //! @return the assignment operation for this descriptor
      const util::Implementation< math::AssignmentUnaryInterface> &GetAssignmentOperation() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_Descriptor->GetNormalSizeOfFeatures();
      }

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const
      {
        return m_Descriptor->GetType().GetDimension();
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors();

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

    private:

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< float> &STORAGE
      );

    }; // class UnaryOperation

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API UnaryOperation< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API UnaryOperation< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API UnaryOperation< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API UnaryOperation< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_UNARY_OPERATION_H_
