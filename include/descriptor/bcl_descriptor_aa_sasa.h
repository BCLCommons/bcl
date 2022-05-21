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

#ifndef BCL_DESCRIPTOR_AA_SASA_H_
#define BCL_DESCRIPTOR_AA_SASA_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "assemble/bcl_assemble_aa_exposure_interface.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_aa_data.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASasa
    //! @brief Generates the solvent accessible surface area of the given amino acid
    //!
    //! @see @link example_descriptor_aa_sasa.cpp @endlink
    //! @author woetzen, mendenjl
    //! @date Apr 24, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASasa :
      public BaseElement< biol::AABase, float>
    {

    private:

    //////////
    // data //
    //////////

      //! AAExposure function to be used for calculations
      util::ShPtr< assemble::AAExposureInterface> m_AAExposure;

      //! class alias (necessary because the AAExposureInterfaces are not under the Enumerated<> framework)
      std::string m_Alias;

      //! class description (necessary because the AAExposureInterfaces are not under the Enumerated<> framework)
      std::string m_Description;

      //! ShPtr to AANeighborListContainer generator
      util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > m_AANeighborListContainerGenerator;

      //! store on double for each pointer to a AABase
      storage::Map< util::SiPtr< const biol::AABase>, float, biol::AALessThanSeqID> m_AASasaStorage;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from exposure measure
      //! @param EXPOSURE_MEASURE exposure measure to use in order to calculate SASA
      //! @param ALIAS alias for this measure
      //! @param DESCRIPTION description of what this class does
      AASasa
      (
        const assemble::AAExposureInterface &EXPOSURE_MEASURE,
        const std::string &ALIAS,
        const std::string &DESCRIPTION
      );

      //! @brief Clone function
      //! @return pointer to new AASasa
      AASasa *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief write neighbor counts and neighbor vectors
      //! @return no return
      static void WriteNCNV
      (
        const size_t MINIMAL_SEQUENCE_SEPARATION,
        std::ostream &OSTREAM,
        const char &CHAIN_ID,
        const assemble::ProteinModel &MODEL
      );

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
        const iterate::Generic< const biol::AABase> &ELEMENT,
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
      CachePreference GetNormalCachePreference() const
      {
        // auto-cache this descriptor, as its calculation is slow
        return e_PreferCache;
      }

    }; // class AASasa

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_SASA_H_
