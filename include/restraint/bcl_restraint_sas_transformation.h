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

#ifndef BCL_RESTRAINT_SAS_TRANSFORMATION_H_
#define BCL_RESTRAINT_SAS_TRANSFORMATION_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_scattering_data.h"
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "score/bcl_score_sas_type.h"
#include "util/bcl_util_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasTransformation
    //! @brief performs transformation of computed SAXS curves into different scales
    //! @details scales compuated curve, normalize curves, log10 scale, take derivative of curves
    //!
    //! @see @link example_restraint_saxs_transformation.cpp @endlink
    //! @author putnamdk
    //! @date Sep 11, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasTransformation :
      public util::FunctionInterfaceSerializable
      <
        SasExperimentalAndCalculatedData,
        SasExperimentalAndCalculatedData
      >
    {

    public:

      enum TransformationType
      {
        e_Absolute,
        e_ScaleCalculatedProfile,
        e_NormalizeData,
        e_SetYScale,
        e_Log10Profiles,
        e_DerivativeProfiles,
        e_None,
        s_NumberTransformationTypes
      };

      //! @brief TransformationType as string
      //! @param TRANSFORMATION_TYPE the TransformationType
      //! @return the string for the TransformationType
      static const std::string &GetFunctionDescriptor( const TransformationType &TRANSFORMATION_TYPE);

      //! @brief TransformationType is used for I/O of ScoreFunction
      typedef util::WrapperEnum
          < TransformationType, &GetFunctionDescriptor, s_NumberTransformationTypes> TransformationTypeEnum;

    private:

    //////////
    // data //
    //////////

      //! container to hold the enums representing the transformations and order of transformatations to perform
      storage::Vector< TransformationTypeEnum> m_Transformations;

      //! flag to print transformations performed
      bool m_OutputTransformations;

      //! flag to use errors in scaling operations
      bool m_UseErrors;

      //! parameter for y_max value
      double m_Ymax;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SasTransformation();

      //! @brief constructor taking member variables
      //! @param OUTPUT_TRANSFORMATIONS flag to print transformations to a file
      //! @param USE_ERRORS DATA_TYPE flag to use errors
      //! @param Y_MAX the max value for the y scale
      SasTransformation
      (
        const storage::Vector< TransformationTypeEnum> &TRANSFORMATION_TYPES,
        const bool &OUTPUT_TRANSFORMATIONS,
        const bool &USE_ERRORS,
        const double &Y_MAX
      );

      //! @brief Clone function
      //! @return pointer to new AnalyzeAtomDistanceHeatmap
      SasTransformation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief performs the desired operation.
      //! @param takes a changeable reference to SasExperimentalAndCalculatedData

      void ScaleCalculatedProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const;
      void NormalizeData( SasExperimentalAndCalculatedData &DATA_OBJECT) const;
      void SetYScale( SasExperimentalAndCalculatedData &DATA_OBJECT) const;
      void LogProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const;
      void DerivativeProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const;
      void LogtoAbsolute( SasExperimentalAndCalculatedData &DATA_OBJECT) const;

      const bool &GetUseErrors() const
      {
        return m_UseErrors;
      }

      const storage::Vector< TransformationTypeEnum> &GetTransformationTypes() const
      {
        return m_Transformations;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief Transforms the data object into desired form
      //! @param ORIGINAL_DATA the raw data to perform transformations on
      //! @return Object that has been transformed
      SasExperimentalAndCalculatedData operator()( const SasExperimentalAndCalculatedData &ORIGINAL_DATA) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

    }; // class SasTransformation

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_TRANSFORMATION_H_
