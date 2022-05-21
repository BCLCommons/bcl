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

#ifndef BCL_DESCRIPTOR_MOLECULE_FRAGMENT_RESCALE_H_
#define BCL_DESCRIPTOR_MOLECULE_FRAGMENT_RESCALE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeFragmentRescale
    //! @brief Computes sequence-wide statistics (column-wise)
    //!
    //! @see @link example_descriptor_molecule_fragment_rescale.cpp @endlink
    //! @author mendenjl
    //! @date Sep 08, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeFragmentRescale :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    public:
    //////////
    // data //
    //////////

      //! all possible neighbor adapt functions
      enum Statistic
      {
        e_MinMax,
        e_ZScore,
        e_SubtractMean,
        e_DivideMean,
        s_NumberStatistics
      };

      //! @brief Statistic as string
      //! @param STAT the statistic
      //! @return the string for the stat
      static const std::string &GetStatisticName( const Statistic &STAT);

      //! @brief NeighborKernelEnum enum I/O helper
      typedef util::WrapperEnum< Statistic, &GetStatisticName, s_NumberStatistics> StatisticEnum;

    private:

      storage::Vector< StatisticEnum> m_Statistics;
      CheminfoProperty m_Descriptor; //!< property to calculate internally
      bool m_CacheConformations; //!< Whether to cache conformations. Should be true when rescaling many conformations of a few configurations
      util::Implementation< chemistry::FragmentSplitInterface> m_Splitter; //!< Splitter to use on molecule to get the subfragments

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeFragmentRescale();

      //! @brief constructor from implementation and desired statistics
      MoleculeFragmentRescale
      (
        const CheminfoProperty &DESCRIPTOR,
        const storage::Vector< StatisticEnum> &STATS
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeFragmentRescale
      MoleculeFragmentRescale *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      //! @note only one value is needed if this is a numeric descriptor, for char descriptors, assume that 99999 is the
      //! @note max, so 5 characters is sufficient
      size_t GetNormalSizeOfFeatures() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

    }; // class MoleculeFragmentRescale

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_FRAGMENT_RESCALE_H_
