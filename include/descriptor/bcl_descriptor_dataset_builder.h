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

#ifndef BCL_DESCRIPTOR_DATASET_BUILDER_H_
#define BCL_DESCRIPTOR_DATASET_BUILDER_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "chemistry/bcl_chemistry.fwd.hh"
#include "iterate/bcl_iterate.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "bcl_descriptor_combine.h"
#include "bcl_descriptor_type.h"
#include "model/bcl_model_feature_label_set.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DatasetBuilder
    //! @brief Mutable dataset containing features, results, identification, and rescale
    //!
    //! @see @link example_descriptor_dataset_builder.cpp @endlink
    //! @author mendenjl
    //! @date Jan 04, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class DatasetBuilder :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Type of dataset to build (e.g. atom-based pairwise)
      Type m_Type;

      //! Dataset descriptor for the features
      Combine< t_DataType, float> m_Feature;

      //! Dataset descriptor for the results
      Combine< t_DataType, float> m_Result;

      //! Identification for each row in the feature dataset
      //! This too has to be a matrix to maintain fixed width
      Combine< t_DataType, char> m_Identification;

      //! Feature label set for features
      util::ShPtr< model::FeatureLabelSet> m_FeatureLabels;

      //! Feature label set for features
      util::ShPtr< model::FeatureLabelSet> m_ResultLabels;

      //! Feature label set for features
      util::ShPtr< model::FeatureLabelSet> m_IdLabels;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DatasetBuilder();

      //! @brief constructor from info
      //! @param FEATURE descriptor for the features
      //! @param RESULT descriptor for the results
      //! @param ID descriptor for the ID
      DatasetBuilder
      (
        const util::ObjectDataLabel &FEATURE,
        const util::ObjectDataLabel &RESULT,
        const util::ObjectDataLabel &ID
      );

      //! @brief Clone function
      //! @return pointer to new DatasetBuilder
      DatasetBuilder *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the size of each feature in the dataset
      //! @return return the size of each feature in the dataset
      size_t GetFeatureSize() const;

      //! @brief return the size of each id (# characters) in the dataset
      //! @return return the size of each id (# characters) in the dataset
      size_t GetIdSize() const;

      //! @brief return the size of each result in the dataset
      //! @return return the size of each result in the dataset
      size_t GetResultSize() const;

      //! @brief set the feature code
      //! @param FEATURES the desired feature code
      void SetFeatureCode( const util::ObjectDataLabel &FEATURES);

      //! @brief set the results code
      //! @param RESULTS the desired results code
      void SetResultsCode( const util::ObjectDataLabel &RESULTS);

      //! @brief set the id code
      //! @param IDS the desired id code
      void SetIdCode( const util::ObjectDataLabel &IDS);

      //! @brief return the feature code
      //! @return return the feature code
      const Combine< t_DataType, float> &GetFeatureCode() const;

      //! @brief return the result code
      //! @return return the result code
      const Combine< t_DataType, float> &GetResultCode() const;

      //! @brief return the id code
      //! @return return the id code
      const Combine< t_DataType, char> &GetIdCode() const;

      //! @brief Type of dataset to build (e.g. atom-based pairwise)
      //! @return the type of dataset that will be built
      Type GetType() const;

      //! @brief Override the computed type of descriptor
      //! @param TYPE the actual type of descriptors desired
      void SetType( const Type &TYPE);

    ////////////////
    // operations //
    ////////////////

      //! @brief create a dataset from the given sequence interface
      //! @param SEQUENCE sequence of interest
      Dataset operator()( const SequenceInterface< t_DataType> &SEQUENCE);

      //! @brief create a dataset from the given sequence interfaces
      //! @param ITERATOR iterator to the sequence interfaces
      //! @param START start feature inside the given iterator
      //! @param MAX_FEATURES maximum number of features to generate
      Dataset operator()
      (
        const iterate::Generic< const SequenceInterface< t_DataType> > &ITERATOR,
        const size_t &START = 0,
        const size_t &MAX_FEATURES = std::numeric_limits< size_t>::max()
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class DatasetBuilder

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API DatasetBuilder< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API DatasetBuilder< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API DatasetBuilder< biol::Mutation>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_DATASET_BUILDER_H_
