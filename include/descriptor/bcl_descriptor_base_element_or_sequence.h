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

#ifndef BCL_DESCRIPTOR_BASE_ELEMENT_OR_SEQUENCE_H_
#define BCL_DESCRIPTOR_BASE_ELEMENT_OR_SEQUENCE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "bcl_descriptor_iterator.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BaseElementOrSequence
    //! @brief A base class for descriptors of elements (atoms or aas) of a sequence, and a different implementation for
    //!        sequences themselves.  This class is useful when there are computations performed on the whole vector
    //!        for example, smoothing, by the descriptor, since performance is severely degraded if this must be
    //!        performed for every atom.
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Feb 24, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class BaseElementOrSequence :
      public Base< t_DataType, t_ReturnType>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new BaseElementOrSequence
      virtual BaseElementOrSequence< t_DataType, t_ReturnType> *Clone() const = 0;

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
        if( ITR.GetType().GetDimension() == size_t( 1))
        {
          Calculate( ITR( 0), STORAGE);
        }
        else if( ITR.GetType().GetDimension() == size_t( 2))
        {
          Calculate( ITR( 0), ITR( 1), STORAGE);
        }
        else
        {
          Calculate( STORAGE);
        }
      }

      //! @brief True if a derived class expects all calls to SetDimension to also call SetDimension on internally-held classes
      //! @note users do not normally need to override this function
      bool InjectDimensions() const
      {
        return false;
      }

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      bool DimensionIsWellDefined() const
      {
        return false;
      }

      //! @brief calculate the descriptors
      //! @param ELEMENT_A, ELEMENT_B the elements of the sequence of interest
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        const iterate::Generic< const t_DataType> &ELEMENT_A,
        const iterate::Generic< const t_DataType> &ELEMENT_B,
        linal::VectorReference< t_ReturnType> &STORAGE
      ) = 0;

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        const iterate::Generic< const t_DataType> &ELEMENT,
        linal::VectorReference< t_ReturnType> &STORAGE
      ) = 0;

      //! @brief calculate the descriptors across a sequence
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        linal::VectorReference< t_ReturnType> &STORAGE
      ) = 0;

      //! @brief hook that derived classes can override to add behavior after every time SetDimension is called
      void SetDimensionHook()
      {
        BCL_Assert
        (
          this->GetDimensionSetting() <= size_t( 2),
          this->GetClassIdentifier() + " does not work with more than 2 " + GetStaticClassName< t_DataType>()
        );
      }

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const
      {
        const size_t dim_setting( this->GetDimensionSetting());
        return util::IsDefined( dim_setting) ? dim_setting : size_t( 1);
      }

      //! @brief True if the intended dimension can only be determined if this descriptor is used with other descriptors
      //! @note users do not normally need to override this function
      bool DimensionIsContextual() const
      {
        return true;
      }

    }; // template class BaseElementOrSequence

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_BASE_ELEMENT_OR_SEQUENCE_H_
