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

#ifndef BCL_DESCRIPTOR_AA_BLAST_PROFILE_ENTROPY_H_
#define BCL_DESCRIPTOR_AA_BLAST_PROFILE_ENTROPY_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_type.h"
#include "bcl_descriptor_window_alignment_type.h"
#include "linal/bcl_linal_matrix.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_symmetric_matrix.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AABlastProfileEntropy
    //! @brief Obtains probabilities of pairs of AAs from a file
    //!
    //! @see @link example_descriptor_aa_blast_profile_entropy.cpp @endlink
    //! @author mendenjl
    //! @date Apr 02, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AABlastProfileEntropy :
      public BaseElement< biol::AABase, float>
    {

    private:

    //////////
    // data //
    //////////

      //! used to store the entropy for each SS type and each ss-type pair
      typedef linal::Vector< float> t_Entropy;
      typedef storage::VectorND< 5, t_Entropy> t_BlastProfileBinEntropies;
      typedef storage::Vector< storage::Vector< t_BlastProfileBinEntropies> > t_AAEntropyStorage;

      //! entropies, indexed by dominant aa type (0), target aa type (1), blast profile bin (5), ss type (3)
      util::SiPtr< const t_AAEntropyStorage> m_EntropiesPtr;

      //! Actual file name for the statistics file
      std::string m_Filename;

      //! penalty multiplier for when the actual aa type is not one of the dominant blast types
      float m_MutantWeight;

      //! true to pretend that the only dominant AA type is the native type, when possible
      bool m_NativeDominance;

      //! number of entropies in file
      size_t m_NumberEntropies;

      //! descriptor used to calculate AA blast profile; must return 20 values
      util::Implementation< Base< biol::AABase, float> > m_TypeCalculator;

      //! Mutex for access to s_EntropyMapStorage
      static sched::Mutex s_EntropyMapMutex;

      //! static map holding all the probabilities from each file; key is filename
      static storage::Map< std::string, t_AAEntropyStorage> s_EntropyMapStorage;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param INITIALIZE true if it is desired that this class be fully initialized; should be set to true unless
      //!        using this class within the code
      AABlastProfileEntropy( const bool &INITIALIZE = true);

      //! @brief virtual copy constructor
      AABlastProfileEntropy *Clone() const;

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

      //! @brief core function, computes the ss-specific blast entropy for the blast profile of the given AABase
      //! @param ELEMENT iterator to the AA of interest
      //! @return vector containing the values, one for each AA
      t_Entropy ComputeBlastSSEntropy( const iterate::Generic< const biol::AABase> &ITR);

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT: the element pair of interest
      //! @param STORAGE storage for the descriptor
      //! @return true, if the calculation was successful
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to forward calls from SetObject on to all internal implementations
      //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
      //! that function is called
      iterate::Generic< Base< biol::AABase, float> > GetInternalDescriptors();

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

      //! @brief get the bin for a particular blast profile value
      //! @param BLAST_PROFILE the blast profile value
      //! @return the bin for the blast profile value
      int GetBin( const int &BP_VALUE) const;
    }; // class AABlastProfileEntropy

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_BLAST_PROFILE_ENTROPY_H_
