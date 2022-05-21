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

#ifndef BCL_DESCRIPTOR_AA_TRIPLET_HELIX_LOCATION_TYPE_H_
#define BCL_DESCRIPTOR_AA_TRIPLET_HELIX_LOCATION_TYPE_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_type.h"
#include "bcl_descriptor_window_alignment_type.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector_nd.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AATripletHelixLocationType
    //! @brief Returns helix initiation, core, and termination propensities
    //!
    //! @see @link example_descriptor_aa_triplet_helix_location_type.cpp @endlink
    //! @author mendenjl
    //! @date May 20, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AATripletHelixLocationType :
      public BaseElement< biol::AABase, float>
    {
    public:

    //////////
    // data //
    //////////

      //! Type used to store SS propensities in the order
      //! HelixInitiation-2 HelixInitiation-1 HelixInitiation HelixInitiation+1 HelixInitiation+2
      //! HelixCore
      //! HelixTermination-2 HelixTermination-1 HelixTermination HelixTermination+1 HelixTermination+2
      typedef linal::VectorND< float, 13> Propensity;

      //! Type used for storage of propensities.
      //! Outer index is the hash id
      typedef storage::Vector< Propensity> TripletSSPropensityType;

    private:

    //////////
    // data //
    //////////

      //! Filename containing the transition propensities for each type
      std::string m_SequencesFilename;

      //! weighting for the blast-computed values
      float m_BlastWeight;

      //! Triplet probabilities
      util::SiPtr< const TripletSSPropensityType> m_TripletProbabilitiesPtr;

      //! descriptor used to calculate AA blast profile; must return 20 values
      util::Implementation< Base< biol::AABase, float> > m_TypeCalculator;

      //! effective descriptor used to calculate AA blast profile; differs from m_TypeCalculator in that the former
      //! may be undefined, but if a blast threshold was given, then this descriptor becomes AABlastProbability
      util::Implementation< Base< biol::AABase, float> > m_EffTypeCalculator;

      //! Mutex for access to s_ProbabilityMapStorage
      static sched::Mutex s_TripletProbabilitiesMutex;

      //! static map holding all the probabilities from each file
      //! key is label containing m_SequencesFilename
      static storage::Map< std::string, TripletSSPropensityType> s_ProbabilitiesStorage;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AATripletHelixLocationType() :
        m_BlastWeight( 0.0)
      {
      }

      //! @brief virtual copy constructor
      AATripletHelixLocationType *Clone() const;

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to forward calls from SetObject on to all internal implementations
      //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
      //! that function is called
      iterate::Generic< Base< biol::AABase, float> > GetInternalDescriptors();

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

      //! @brief load the initial propensities vector
      //! @param TRIPLETS_PROB storage for all quintuplets loaded
      void LoadInitialPropensities( TripletSSPropensityType &TRIPLETS_PROB);
    }; // class AATripletHelixLocationType

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_TRIPLET_HELIX_LOCATION_TYPE_H_
