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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "restraint/bcl_restraint_sas_transformation.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_scattering_data.h"
#include "restraint/bcl_restraint_saxs_opti_result.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    //! @brief ScoreFunction as string
    //! @param SCORE_FUNCTION the ScoreFunction
    //! @return the string for the ScoreFunction
    const std::string &SasTransformation::GetFunctionDescriptor( const TransformationType &TRANSFORMATION_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "Absolute",
        "Scale",
        "Normalize",
        "SetYMax",
        "Log10",
        "Derivative",
        "None",
        GetStaticClassName< TransformationType>()
      };

      return s_descriptors[ TRANSFORMATION_TYPE];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasTransformation::s_Instance
    (
      util::Enumerated
      <
        util::FunctionInterfaceSerializable
        <
          SasExperimentalAndCalculatedData,
          SasExperimentalAndCalculatedData
        >
      >::AddInstance( new SasTransformation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasTransformation::SasTransformation() : m_OutputTransformations( false), m_UseErrors( false)
    {
      storage::Vector< TransformationTypeEnum> default_transforms;

      default_transforms.PushBack( SasTransformation::e_ScaleCalculatedProfile);
      default_transforms.PushBack( SasTransformation::e_SetYScale);
      default_transforms.PushBack( SasTransformation::e_Log10Profiles);
      default_transforms.PushBack( SasTransformation::e_DerivativeProfiles);

      m_Transformations = default_transforms;
      m_Ymax = 1.0;
    }

    //! @brief constructor taking member variables
    //! @param OUTPUT_TRANSFORMATIONS flag to print transformations to a file
    //! @param USE_ERRORS DATA_TYPE flag to use errors
    //! @param Y_MAX the max value for the y scale
    SasTransformation::SasTransformation
    (
      const storage::Vector< TransformationTypeEnum> &TRANSFORMATION_TYPES,
      const bool &PRINT_TRANSFORMATIONS,
      const bool &USE_ERRORS,
      const double &Y_MAX
    )
      :
      m_Transformations( TRANSFORMATION_TYPES),
      m_OutputTransformations( PRINT_TRANSFORMATIONS),
      m_UseErrors( USE_ERRORS),
      m_Ymax( Y_MAX)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceHeatmap
    SasTransformation *SasTransformation::Clone() const
    {
      return new SasTransformation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasTransformation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SasTransformation::GetAlias() const
    {
      static const std::string s_Name( "Transform");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    void SasTransformation::ScaleCalculatedProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      const double scaling_factor
      (
        SasAnalysis::CalculateScalingWeight( DATA_OBJECT, m_UseErrors)
      );

      DATA_OBJECT.ScaleCalculatedData( scaling_factor);
    }

    void SasTransformation::NormalizeData( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.NormalizeData();
    }

    void SasTransformation::SetYScale( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.SetYScale( m_Ymax);
    }

    void SasTransformation::LogProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.Log10();
    }

    void SasTransformation::DerivativeProfile( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.Derivative();
    }

    void SasTransformation::LogtoAbsolute( SasExperimentalAndCalculatedData &DATA_OBJECT) const
    {
      DATA_OBJECT.LogtoAbsolute();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    SasExperimentalAndCalculatedData SasTransformation::operator()
     ( const SasExperimentalAndCalculatedData &ORIGINAL_DATA) const
    {
      SasExperimentalAndCalculatedData data_object( ORIGINAL_DATA);

      if( m_OutputTransformations)
      {
        data_object.WriteToGnuplotFileName( "raw.data");
      }

      // Iterate over the storage vector of Enums
      for
      (
          storage::Vector< TransformationTypeEnum>::const_iterator
            enum_itr( m_Transformations.Begin()),
            enum_itr_end( m_Transformations.End());
            enum_itr != enum_itr_end;
            ++enum_itr
      )
      {
        switch( *enum_itr)
        {
          case e_ScaleCalculatedProfile:
            ScaleCalculatedProfile( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "scaled.data");
            }
            break;
          case e_SetYScale:
            SetYScale( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "scaled_y_axis.data");
            }
            break;
          case e_NormalizeData:
            NormalizeData( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "normalized.data");
            }
            break;
          case e_Log10Profiles:
            LogProfile( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "logarithmic.data");
            }
            break;
          case e_DerivativeProfiles:
            DerivativeProfile( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "derivative.data");
            }
            break;
          case e_Absolute:
            LogtoAbsolute( data_object);
            if( m_OutputTransformations)
            {
              data_object.WriteToGnuplotFileName( "absolute.data");
            }
            break;
          case e_None:
            // Do nothing
            break;
          default:
            BCL_Assert
            (
              false,
              "Unknown Transformation Procedure"
            );
            break;
        }
      }
      return data_object;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SasTransformation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Transforms Computed SAXS Profile and Experimental Saxs Profile into different forms based on user input."
      );

      SasTransformation default_object;

      parameters.AddInitializer
      (
        "transformations",
        "computed scaling factor between experimental and calculated profiles and scale calculated profile",
        io::Serialization::GetAgent( &m_Transformations),
        io::Serialization::GetAgent( &default_object.m_Transformations)->GetLabel().ToString()
      );

      parameters.AddInitializer
      (
        "print_transformations",
        "Print transformations performed",
        io::Serialization::GetAgent( &m_OutputTransformations),
        "false"
      );

      parameters.AddInitializer
      (
        "use_errors",
        "use experimental errors in scaling saxs profiles",
        io::Serialization::GetAgent( &m_UseErrors)
      );

      parameters.AddInitializer
      (
        "y_max",
        "set max value on y scale for graphing",
        io::Serialization::GetAgent( &m_Ymax),
        "-9999"
      );

      return parameters;
    }
  } // namespace restraint
} // namespace bcl
