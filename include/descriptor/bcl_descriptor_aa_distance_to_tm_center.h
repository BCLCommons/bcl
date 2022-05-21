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

#ifndef BCL_DESCRIPTOR_AA_DISTANCE_TO_TM_CENTER_H_
#define BCL_DESCRIPTOR_AA_DISTANCE_TO_TM_CENTER_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_type.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "iterate/bcl_iterate_generic.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AADistanceToTMCenter
    //! @brief AADistanceToTMCenter gives the distance of an AA from the center of the membrane,
    //!        assuming that the protein has already been transformed using the PDBTM membrane, if one existed
    //!
    //! @see @link example_descriptor_aa_distance_to_tm_center.cpp @endlink
    //! @author mendenjl
    //! @date Apr 19, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AADistanceToTMCenter :
      public BaseElement< biol::AABase, float>
    {
    public:

    //////////
    // data //
    //////////

    private:

      //! whether the current protein is a membrane protein
      bool m_HaveMembrane;

      //! Whether we have already checked whether this protein is a membrane protein
      bool m_HaveCheckedForMembrane;

      //! Whether to normalize by PDBTM membrane distance to the range 0-1 (membrane core), 1-2 (transition), and 3 (solution)
      bool m_MembraneNorm;

      //! Default value, if no membrane is present
      double m_DefaultValue;

      //! membrane and transition width, needed for normalization
      double m_MembraneWidth;
      double m_TransitionWidth;

      //! single instance of the class
      static const util::SiPtr< const util::ObjectInterface> s_AATMPositionInstance;
      static const util::SiPtr< const util::ObjectInterface> s_AATMNormPositionInstance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param NORM whether to take the membrane norm
      AADistanceToTMCenter( const bool &MEMBRANE_NORM);

      //! @brief virtual copy constructor
      AADistanceToTMCenter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT_A the element of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class AADistanceToTMCenter

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_DISTANCE_TO_TM_CENTER_H_
