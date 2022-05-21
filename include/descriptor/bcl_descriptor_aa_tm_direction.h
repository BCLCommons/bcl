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

#ifndef BCL_DESCRIPTOR_AA_TM_DIRECTION_H_
#define BCL_DESCRIPTOR_AA_TM_DIRECTION_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa_base.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AATMDirection
    //! @brief binary descriptor: returns 1 if the AA is in the membrane and is exposed to a pore, 0 otherwise
    //!
    //! @see @link example_descriptor_aa_tm_direction.cpp @endlink
    //! @author mendenjl
    //! @date Sep 19, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AATMDirection :
      public BaseElement< biol::AABase, float>
    {

    private:

    //////////
    // data //
    //////////

      std::string m_Alias; //!< Alias for this class

      sspred::CIPhiPsi::TMDirectionTypeEnum m_Type;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AATMDirection( const sspred::CIPhiPsi::TMDirectionTypeEnum &DIRECTION, const std::string &ALIAS) :
        m_Alias( ALIAS),
        m_Type( DIRECTION)
      {
      }

      //! @brief Clone function
      //! @return pointer to new AATMDirection
      AATMDirection *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

    }; // class AATMDirection

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_TM_DIRECTION_H_
