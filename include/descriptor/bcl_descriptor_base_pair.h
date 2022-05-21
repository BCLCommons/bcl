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

#ifndef BCL_DESCRIPTOR_BASE_PAIR_H_
#define BCL_DESCRIPTOR_BASE_PAIR_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "bcl_descriptor_iterator.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BasePair
    //! @brief A base class for overall descriptors of a element pairs, e.g. AA contacts, molecular RDF functions
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Jan 04, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class BasePair :
      public Base< t_DataType, t_ReturnType>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new BasePair
      virtual BasePair< t_DataType, t_ReturnType> *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

    private:

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      virtual void RecalculateImpl( const Iterator< t_DataType> &ITR, linal::VectorReference< t_ReturnType> &STORAGE)
      {
        Calculate( ITR( 0), ITR( 1), STORAGE);
      }

      //! @brief calculate the descriptors
      //! @param ELEMENT_A, ELEMENT_B: the element pair of interest
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        const iterate::Generic< const t_DataType> &ELEMENT_A,
        const iterate::Generic< const t_DataType> &ELEMENT_B,
        linal::VectorReference< t_ReturnType> &STORAGE
      ) = 0;

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const
      {
        return size_t( 2);
      }

      //! @brief True if a derived class expects all calls to SetDimension to also call SetDimension on internally-held classes
      //! @note users do not normally need to override this function
      bool InjectDimensions() const
      {
        return false;
      }

    }; // template class BasePair

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_BASE_PAIR_H_
