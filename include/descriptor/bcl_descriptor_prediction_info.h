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

#ifndef BCL_DESCRIPTOR_PREDICTION_INFO_H_
#define BCL_DESCRIPTOR_PREDICTION_INFO_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "bcl_descriptor_prediction.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_roc_curve.h"
#include "model/bcl_model_retrieve_interface.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PredictionInfo
    //! @brief combines the output of descriptors into one vector
    //!
    //! @see @link example_descriptor_prediction_info.cpp @endlink
    //! @author mendenjl
    //! @date Jan 19, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class PredictionInfo :
      public Base< t_DataType, float>
    {
    public:

      //! all possible neighbor adapt functions
      enum Statistic
      {
        e_Min,
        e_Max,
        e_Mean,
        e_StandardDeviation,
        s_NumberStatistics
      };

      //! @brief Statistic as string
      //! @param STAT the statistic
      //! @return the string for the stat
      static const std::string &GetStatisticName( const Statistic &STAT);

      //! @brief NeighborKernelEnum enum I/O helper
      typedef util::WrapperEnum< Statistic, ( &GetStatisticName), s_NumberStatistics> StatisticEnum;

    private:

    //////////
    // data //
    //////////

      // prediction descriptor; used for computing the actual prediction
      Prediction< t_DataType> m_Prediction;

      // statistics to compute
      storage::Vector< StatisticEnum> m_Statistics;

      //! Measures to compute
      storage::Vector< math::ContingencyMatrixMeasures::MeasureEnum> m_Measures;

      //! roc-curves for the model storage
      storage::Vector< storage::Map< double, math::ContingencyMatrix> > m_RocCurves;

      //! piecewise functions for LocalPPV, if that measure is requested
      storage::Vector< math::PiecewiseFunction> m_LocalPPVs;

      //! Bool : whether statistics or metrics go first in the output vector (dictated by user-ordering)
      bool m_StatisticsBeforeMetrics;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      PredictionInfo *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the type of this descriptor
      //! @return the type of this descriptor (should ignore dimension setting)
      Type GetType() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

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
        return m_Prediction.GetNormalDimension();
      }

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const
      {
        return m_Prediction.GetSymmetry();
      }

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const
      {
        return m_Prediction.ConsiderRepeatedElements();
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief True if a derived class has a well-defined dimension (vs. having a dimension determined by the inputs)
      //! @note users do not normally need to override this function
      bool DimensionIsWellDefined() const
      {
        return m_Prediction.DimensionIsWellDefined();
      }

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors();

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief load the models: called the 1st time that recalculate is called
      void LoadModels();

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

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

    }; // class PredictionInfo

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API PredictionInfo< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API PredictionInfo< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API PredictionInfo< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API PredictionInfo< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_PREDICTION_INFO_H_
