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

#ifndef BCL_MODEL_RESCALE_FEATURE_DATA_SET_H_
#define BCL_MODEL_RESCALE_FEATURE_DATA_SET_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "math/bcl_math_linear_function.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RescaleFeatureDataSet
    //! @brief This class will rescale all feature vectors within a FeatureDataSet to within RANGE
    //!
    //! @author loweew
    //!
    //! @see @link example_model_rescale_feature_data_set.cpp @endlink
    //!
    //! @date 10/01/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RescaleFeatureDataSet :
      public util::ObjectInterface
    {

    public:

    //////////
    // enum //
    //////////

      enum Type
      {
        e_None,   //!< No rescaling, works even if there are nan's in the data
        e_MinMax, //!< Rescaled such that min/max fall at the edges of the given range, works w/ nan's in the data
        e_AveStd, //!< Rescaled such that ave +/- 2 * std fall at the edges of the given range, fails w/ nan's in the data
        s_NumberTypes
      };

      //! @brief Type as string
      //! @param TYPE the type
      //! @return the string for the kernel
      static const std::string &GetTypeName( const Type &TYPE);

      //! @brief TypeEnum enum I/O helper
      typedef util::WrapperEnum< Type, &GetTypeName, s_NumberTypes> TypeEnum;

    private:

    //////////
    // data //
    //////////

      //! the rescale function (ranges)
      storage::Vector< math::Range< float> > m_RescaleRanges;

      //! the range to scale to
      math::Range< float> m_Range;

      //! Cached rescale functions
      storage::Vector< math::LinearFunction> m_RescaleFunction;

      //! Cached descale functions
      storage::Vector< math::LinearFunction> m_DescaleFunction;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! global static function that allows altering the default null value
      //! The standard default (undefined float) is interpreted to mean picking the
      //! middle of the output range
      static float &GetDefaultValueForEmptyRangedColumns();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RescaleFeatureDataSet();

      //! @brief construct from feature result data set
      //! @param DATA the feature data from which to obtain the ranges
      //! @param RANGE features will be rescaled to this range
      //! @param USE_VARIANCE true to set the range such that 1 corresponds to 1 standard deviation away from the mean
      RescaleFeatureDataSet
      (
        const FeatureDataSetInterface< float> &DATA,
        const math::Range< float> &RANGE,
        const Type &TYPE = e_MinMax
      );

      //! @brief construct from a matrix const interface
      //! @param DATA the matrix from which to obtain the ranges
      //! @param RANGE features will be rescaled to this range
      //! @param TYPE type of rescaling to perform
      RescaleFeatureDataSet
      (
        const linal::MatrixConstInterface< float> &DATA,
        const math::Range< float> &RANGE,
        const Type &TYPE = e_MinMax
      );

      //! @brief construct from a rescaler and a function to be applied to all rescalings
      //! @param RESCALING the original rescaling function
      //! @param MODIFIERS modification function be applied to all scaling
      RescaleFeatureDataSet
      (
        const RescaleFeatureDataSet &RESCALING,
        const util::ShPtrVector< math::FunctionInterfaceSerializable< double, double> > &MODIFIERS
      );

      //! @brief Clone function
      //! @return pointer to new RescaleFeatureDataSet
      RescaleFeatureDataSet *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to the scale ranges
      //! @return a rescale range objects
      const storage::Vector< math::Range< float> > &GetRescaleRanges() const
      {
        return m_RescaleRanges;
      }

      //! @brief get the number of features
      //! @return the number of columns that should be in any matrix that is passed in
      size_t GetSize() const
      {
        return m_RescaleRanges.GetSize();
      }

      //! @brief rescale a value for a particular column using the appropriate rescale function
      //! @param COLUMN column that the value logically belongs to
      //! @param VALUE value to rescale
      //! @return the rescaled value
      float RescaleValue( const size_t &COLUMN, const float &VALUE) const;

      //! @brief descale a value for a particular column using the appropriate descale function
      //! @param COLUMN column that the value logically belongs to
      //! @param VALUE value to descale
      //! @return the descaled value
      float DescaleValue( const size_t &COLUMN, const float &VALUE) const;

      //! @brief get rescale output range
      //! @return rescale output range
      const math::Range< float> &GetRange() const
      {
        return m_Range;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief scales back to pre-rescaled values
      //! @param FEATURE the feature data set interface to scale back to original values
      //! @return FeatureDataSet of scaled values
      FeatureDataSet< float> DeScale( const FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief scales back to pre-rescaled values
      //! @param MATRIX the matrix to descale
      void DeScaleMatrix( linal::MatrixInterface< float> &MATRIX) const;

      //! @brief scale a MatrixInterface
      //! @param MATRIX the matrix to rescale
      void RescaleMatrix( linal::MatrixInterface< float> &MATRIX) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief scale a FeatureDataSet
      //! @param FEATURE the FeatureDataSetInterface to rescale
      //! @return rescaled FeatureDataSet to one range
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief test for equality
      //! @param OTHER other rescaling function
      //! @return true if the rescaling functions are equal
      bool operator !=( const RescaleFeatureDataSet &OTHER) const;

      //! @brief test for equality
      //! @param OTHER other rescaling function
      //! @return true if the rescaling functions are equal
      bool operator ==( const RescaleFeatureDataSet &OTHER) const;

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
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief determine m_RescaleRanges from the given input matrix
      //! @param DATA the matrix from which to obtain the ranges
      //! This option does handle nan's in the input data
      void GetInputRangesFromMatrix( const linal::MatrixConstInterface< float> &DATA);

      //! @brief determine m_RescaleRanges from the given input matrix
      //! @param DATA the matrix from which to obtain the ranges
      //! @param FRACTION (0,1) the minimum fraction of each column's values that should be covered by the chosen range
      //! This option does not handle nan's in the input data
      void GetInputRangesFromMatrixRange
      (
        const linal::MatrixConstInterface< float> &DATA,
        const float &PERCENT
      );

      //! @brief determine m_RescaleRanges from the given input matrix
      //! @param DATA the matrix from which to obtain the ranges
      //! This option does not handle nan's in the input data
      void GetInputRangesFromMatrixAveStd( const linal::MatrixConstInterface< float> &DATA);

      //! @brief construct rescaling and descaling functions
      void ConstructScalingFunctions();

    }; // class RescaleFeatureDataSet

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_RESCALE_FEATURE_DATA_SET_H_
