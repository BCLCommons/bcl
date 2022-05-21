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

#ifndef BCL_DESCRIPTOR_BUSER_METRIC_H_
#define BCL_DESCRIPTOR_BUSER_METRIC_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_atom_environment_bender.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_molecule_environment.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BuserMetric
    //! @brief calculates the Buser similarity scores
    //! @details MolPrint2D definition is from:
    //!   Veber, et. al. 2002, Molecular Properties That Influence the Oral Bioavailability of Drug Candidates
    //!
    //! @see @link example_descriptor_buser_metric.cpp @endlink
    //! @author vuot
    //! @date Jun 12, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BuserMetric :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {
    private:
      typedef chemistry::AtomEnvironmentBender::AtomTypeEnum t_AtomTypeEnum;
      //! Single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      t_AtomTypeEnum m_AtomType;
      size_t m_ActiveNum;
      //! List of molecule environments of active compounds
      storage::Vector< chemistry::MoleculeEnvironment> m_Actives;
      //! filename of the file containing active compounds
      std::string m_FileName;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BuserMetric();

      //! @brief Clone function
      //! @return pointer to new MoleculeRotatableBonds
      BuserMetric *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);
    };
  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_BUSER_METRIC_H_
