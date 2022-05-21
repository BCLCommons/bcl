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

#ifndef BCL_DESCRIPTOR_KOHONEN_MAP_INFO_H_
#define BCL_DESCRIPTOR_KOHONEN_MAP_INFO_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "model/bcl_model_kohonen_network_average.h"
#include "model/bcl_model_retrieve_interface.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KohonenMapInfo
    //! @brief combines the output of descriptors into one vector
    //!
    //! @see @link example_descriptor_kohonen_map_info.cpp @endlink
    //! @author mendenjl
    //! @date Sep 24, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class KohonenMapInfo :
      public Base< t_DataType, float>
    {
    private:

    //////////
    // data //
    //////////

      size_t m_NumberOutputs; //!< Number of outputs

      //! generic type of the descriptor (depends on the result descriptor the models were trained on)
      Type m_Type;

      //! Model storage; to retrieve models and descriptors
      util::Implementation< model::RetrieveInterface> m_ModelStorage;

      //! pointer to vector of models used by this property
      mutable util::SiPtr< const storage::Vector< util::OwnPtr< model::KohonenNetworkAverage> > > m_Models;

      //! unique descriptor sets needed by each model
      mutable storage::Vector< Combine< t_DataType, float> > m_Properties;

      //! ID of each Combine used, one per model
      storage::Vector< size_t> m_DescriptorIDs;

      //! true if positions will be returned
      bool m_Positions;

      //! true if distances will be returned
      bool m_Distances;

      //! number of models in the storage.  This is useful to cache since m_Models will only be defined once this
      //! descriptor has been used to calculate values for an object
      size_t m_NumberModels;

      //! Mutex to protect access to s_Models
      static sched::Mutex s_ModelsMutex;

      //! Mutex to protect access to descriptors
      static sched::Mutex s_DescriptorsMutex;

      //! Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
      //! gigantic models repetitively
      static storage::Map< std::string, storage::Vector< util::OwnPtr< model::KohonenNetworkAverage> > > s_Models;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, accepts bool of whether to auto-compute mean
      KohonenMapInfo();

      //! @brief virtual copy constructor
      KohonenMapInfo *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the type of this descriptor
      //! @return the type of this descriptor (should ignore dimension setting)
      Type GetType() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_NumberOutputs;
      }

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const
      {
        return m_Type.GetDimension();
      }

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const
      {
        return m_Type.GetSymmetry();
      }

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const
      {
        return m_Type.ConsiderRepeatedObjects();
      }

      //! @brief return the number of models in the storage used by this predictor
      //! @return the number of models in the storage used by this predictor
      size_t GetNumberOfModels() const
      {
        return m_NumberModels;
      }

      //! @brief return the number of outputs returned by each model in the storage
      //! @return the number of outputs returned by each model in the storage
      size_t GetModelsNumberOutputs() const
      {
        return m_NumberOutputs;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief load the models: called the 1st time that recalculate is called
      void LoadModels();

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief True if a derived class expects all calls to SetDimension to also call SetDimension on internally-held classes
      //! False because the internal descriptors for models cannot be changed in size
      bool InjectDimensions() const
      {
        return false;
      }

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      bool DimensionIsWellDefined() const
      {
        return m_ModelStorage.IsDefined();
      }

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< float> &STORAGE
      );

    }; // class KohonenMapInfo

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API KohonenMapInfo< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API KohonenMapInfo< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API KohonenMapInfo< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API KohonenMapInfo< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_KOHONEN_MAP_INFO_H_
