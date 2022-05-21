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

#ifndef BCL_DESCRIPTOR_MUTATION_AA_PROPERTY_H_
#define BCL_DESCRIPTOR_MUTATION_AA_PROPERTY_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_aa_data.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutationAAProperty
    //! @brief Computes desired AA / protein descriptorsfor a given mutation
    //!
    //! @see @link example_descriptor_mutation_aa_property.cpp @endlink
    //! @author mendenjl
    //! @date Jan 22, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutationAAProperty :
      public BaseElement< biol::Mutation, float>
    {

    private:

    //////////
    // data //
    //////////

      //! property to be calculated on the atoms of interest
      util::Implementation< Base< biol::AABase, float> > m_Property;

      //! bool; whether to use the native protein (true), or the mutant (false)
      bool m_UseNative;

      //! mutation set
      util::SiPtr< const biol::ProteinMutationSet> m_ModelWMutations;

      //! SiPtr to set of protein models encompassing all mutations
      util::SiPtr< const storage::Map< biol::Mutation, util::ShPtr< assemble::ProteinModelWithMutations> > > m_MutationsMap;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_NativeInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from exposure measure
      //! @param NATIVE true if the native protein for the mutant should be used
      //! @param PROPERTY property to calculate
      MutationAAProperty
      (
        const bool &NATIVE,
        const util::Implementation< Base< biol::AABase, float> > &PROPERTY
          = util::Implementation< Base< biol::AABase, float> >()
      );

      //! @brief Clone function
      //! @return pointer to new MutationAAProperty
      MutationAAProperty *Clone() const;

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

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const biol::Mutation> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief populate the map with the given molecular property
      void CalculatePropertyOnProtein();

    }; // class MutationAAProperty

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MUTATION_AA_PROPERTY_H_
