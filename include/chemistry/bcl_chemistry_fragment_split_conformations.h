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

#ifndef BCL_CHEMISTRY_FRAGMENT_SPLIT_CONFORMATIONS_H_
#define BCL_CHEMISTRY_FRAGMENT_SPLIT_CONFORMATIONS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_split_interface.h"
#include "bcl_chemistry_rotamer_library_interface.h"
#include "bcl_chemistry_sample_conformations.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentSplitConformations
    //! @brief This class get all components of a molecule complex that is provided.
    //!
    //! @see @link example_chemistry_fragment_split_conformations.cpp @endlink
    //! @author kothiwsk
    //! @date Feb 03, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentSplitConformations :
      public FragmentSplitInterface
    {

    private:

      //! object that samples molecular conformations
      SampleConformations                m_SampleConformations;

      // because this class is often called multiple times for the same molecule, the results from each molecule are
      // temporarily cached locally.
      typedef storage::Triplet
      <
        storage::Vector< AtomType>,
        storage::Vector< sdf::BondInfo>,
        FragmentEnsemble
      > t_SplitCacheTriplet;

      //! cache list for this object
      mutable util::SiPtr< storage::List< t_SplitCacheTriplet> > m_LocalCache;

      //! mutex for m_LocalCache access
      mutable util::SiPtr< sched::Mutex> m_Mutex;

      //! cache for recently generated molecules and mutex to access this list
      static storage::Map
      <
        std::string,
        storage::Pair< storage::List< t_SplitCacheTriplet>, sched::Mutex>
      > s_ConformationCache;

      bool m_UseCache;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief
      FragmentSplitConformations();

      //! virtual copy constructor
      FragmentSplitConformations *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! get the minimum size of a component of interest
      const size_t GetMinSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief operator returns an ensemble of conformations
      //! @param MOLECULE molecule of interest
      //! @return enseble of conformations
      FragmentEnsemble operator()( const ConformationInterface &MOLECULE) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_SPLIT_CONFORMATIONS_H_
